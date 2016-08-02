///////////////////////////////////////////////////////////////////////
//
//  File: OutputTecplot.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Dat file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <set>
#include <string>
using namespace std;

#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/PtsIO.h>

#include "OutputTecplot.h"

namespace Nektar
{
namespace FieldUtils
{

std::string TecplotZoneTypeMap[] = {
    "ORDERED",
    "LINESEG",
    "TRIANGLE",
    "QUADRILATERAL",
    "TETRAHEDRON",
    "BRICK",
    "POLYGON",
    "POLYHEDRON"
};

ModuleKey OutputTecplot::m_className[] =
{
    GetModuleFactory().RegisterCreatorFunction(ModuleKey(eOutputModule, "dat"),
                                               OutputTecplot::create,
                                               "Writes a Tecplot file."),
    GetModuleFactory().RegisterCreatorFunction(ModuleKey(eOutputModule, "plt"),
                                               OutputTecplot::createBinary,
                                               "Writes a Tecplot file in binary"
                                               " plt format.")
};

OutputTecplot::OutputTecplot(FieldSharedPtr f, bool binary) : OutputModule(f),
                                                              m_binary(binary)
{
    if (!f->m_setUpEquiSpacedFields)
    {
        m_requireEquiSpaced = true;
    }

    m_config["double"] =
        ConfigOption(true, "0", "Write double-precision data: more "
                                "accurate but more disk space required");
}

OutputTecplot::~OutputTecplot()
{
}

/**
 * @brief Helper function to write binary data to stream.
 */
template<typename T> void WriteStream(std::ostream &outfile, T data)
{
    T tmp = data;
    outfile.write(reinterpret_cast<char *>(&tmp), sizeof(T));
}

/**
 * @brief Specialisation of WriteStream to support writing std::string.
 *
 * Tecplot binary formats represent all strings by writing out their characters
 * as 32-bit integers, followed by a 32-bit integer null (0) character to denote
 * the end of the string.
 */
template<> void WriteStream(std::ostream &outfile, std::string data)
{
    // Convert string to array of int32_t
    for (std::string::size_type i = 0; i < data.size(); ++i)
    {
        char strChar = data[i];
        int32_t strCharInt = strChar;
        WriteStream(outfile, strCharInt);
    }

    // Now dump out zero character to terminate
    WriteStream(outfile, 0);
}

/**
 * @brief Specialisation of WriteStream to support writing Nektar::Array
 * datatype.
 */
template<typename T> void WriteStream(std::ostream &outfile,
                                      Array<OneD, T> data)
{
    outfile.write(reinterpret_cast<char *>(&data[0]),
                  data.num_elements() * sizeof(T));
}

/**
 * @brief Specialisation of WriteStream to support writing std::vector datatype.
 */
template<typename T> void WriteStream(std::ostream  &outfile,
                                      std::vector<T> data)
{
    outfile.write(reinterpret_cast<char *>(&data[0]),
                  data.size() * sizeof(T));
}

void OutputTecplot::Process(po::variables_map &vm)
{
    LibUtilities::PtsFieldSharedPtr fPts = m_f->m_fieldPts;

    m_doError = (vm.count("error") == 1) ? true : false;
    m_numBlocks = 0;

    // Do nothing if no expansion defined
    if (fPts == LibUtilities::NullPtsField && !m_f->m_exp.size())
    {
        return;
    }

    if (m_f->m_verbose)
    {
        if (m_f->m_comm->TreatAsRankZero())
        {
            cout << "OutputTecplot: Writing file..." << endl;
        }
    }

    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();

    int nprocs = m_f->m_comm->GetSize();
    int rank   = m_f->m_comm->GetRank();
    // Amend for parallel output if required
    if (nprocs != 1)
    {
        int dot       = filename.find_last_of('.');
        string ext    = filename.substr(dot, filename.length() - dot);
        string procId = "_P" + boost::lexical_cast<std::string>(rank);
        string start  = filename.substr(0, dot);
        filename      = start + procId + ext;
    }

    // Open output file
    ofstream outfile(filename.c_str(), m_binary ? ios::binary : ios::out);

    std::vector<std::string> var;
    bool writeHeader = true;

    if (fPts == LibUtilities::NullPtsField)
    {
        // Standard tensor-product element setup.
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> fDef =
            m_f->m_fielddef;

        if (fDef.size())
        {
            var = fDef[0]->m_fields;
        }

        // Calculate number of FE blocks
        m_numBlocks = GetNumTecplotBlocks();

        // Calculate coordinate dimension
        int nBases = m_f->m_exp[0]->GetExp(0)->GetNumBases();
        MultiRegions::ExpansionType HomoExpType = m_f->m_exp[0]->GetExpType();

        m_coordim = m_f->m_exp[0]->GetExp(0)->GetCoordim();

        if (HomoExpType == MultiRegions::e3DH1D)
        {
            int nPlanes = m_f->m_exp[0]->GetZIDs().num_elements();
            if (nPlanes == 1) // halfMode case
            {
                // do nothing
            }
            else
            {
                nBases += 1;
                m_coordim += 1;
                NekDouble tmp = m_numBlocks * (nPlanes - 1);
                m_numBlocks   = (int)tmp;
            }
        }
        else if (HomoExpType == MultiRegions::e3DH2D)
        {
            nBases += 2;
            m_coordim += 2;
        }

        m_zoneType = (TecplotZoneType)(2*(nBases-1) + 1);

        // Calculate connectivity
        CalculateConnectivity();

        // Set up storage for output fields
        m_fields = Array<OneD, Array<OneD, NekDouble> >(
            var.size() + m_coordim);

        // Get coordinates
        int totpoints = m_f->m_exp[0]->GetTotPoints();
        cout << totpoints << endl;

        for (int i = 0; i < m_coordim; ++i)
        {
            m_fields[i] = Array<OneD, NekDouble>(totpoints);
        }

        if (m_coordim == 1)
        {
            m_f->m_exp[0]->GetCoords(m_fields[0]);
        }
        else if (m_coordim == 2)
        {
            m_f->m_exp[0]->GetCoords(m_fields[0], m_fields[1]);
        }
        else
        {
            m_f->m_exp[0]->GetCoords(m_fields[0], m_fields[1], m_fields[2]);
        }

        if (var.size())
        {
            // Backward transform all data
            for (int i = 0; i < m_f->m_exp.size(); ++i)
            {
                if (m_f->m_exp[i]->GetPhysState() == false)
                {
                    m_f->m_exp[i]->BwdTrans(m_f->m_exp[i]->GetCoeffs(),
                                            m_f->m_exp[i]->UpdatePhys());
                }
            }

            for (int i = 0; i < m_f->m_exp.size(); ++i)
            {
                m_fields[i + m_coordim] = m_f->m_exp[i]->UpdatePhys();
            }
        }
    }
    else
    {
        m_coordim = fPts->GetDim();

        if (fPts->GetNpoints() == 0)
        {
            return;
        }

        // Grab connectivity information.
        fPts->GetConnectivity(m_conn);

        // Get fields and coordinates
        fPts->GetPts(m_fields);

        switch (fPts->GetPtsType())
        {
            case LibUtilities::ePtsFile:
            case LibUtilities::ePtsLine:
                m_numPoints.resize(1);
                m_numPoints[0] = fPts->GetNpoints();
                m_zoneType = eOrdered;
                break;
            case LibUtilities::ePtsPlane:
                m_numPoints.resize(2);
                m_numPoints[0] = fPts->GetPointsPerEdge(0);
                m_numPoints[1] = fPts->GetPointsPerEdge(1);
            case LibUtilities::ePtsBox:
                m_numPoints.resize(3);
                m_numPoints[2] = fPts->GetPointsPerEdge(2);
                m_zoneType     = eOrdered;
                break;
            case LibUtilities::ePtsTriBlock:
            {
                m_zoneType = eFETriangle;
                for (int i = 0; i < m_conn.size(); ++i)
                {
                    m_numBlocks += m_conn[i].num_elements() / 3;
                }
                break;
            }
            case LibUtilities::ePtsTetBlock:
            {
                m_zoneType = eFETetrahedron;
                for (int i = 0; i < m_conn.size(); ++i)
                {
                    m_numBlocks += m_conn[i].num_elements() / 4;
                }
                break;
            }
            default:
                ASSERTL0(false, "This points type is not supported yet.");
        }

        // Only write header if we're root or FE blocks
        writeHeader = m_zoneType != eOrdered || rank == 0;
    }

    if (writeHeader)
    {
        WriteTecplotHeader(outfile, var);
    }
    WriteTecplotZone(outfile);
    WriteTecplotConnectivity(outfile);

    cout << "Written file: " << filename << endl;
}

/**
 * Write Tecplot Files Header
 * @param   outfile Output file name.
 * @param   var                 variables names
 */
void OutputTecplot::WriteTecplotHeader(std::ofstream &outfile,
                                       std::vector<std::string> &var)
{
    if (m_binary)
    {
        // Version number
        outfile << "#!TDV112";

        // Int value of 1 for endian check
        WriteStream(outfile, 1);

        // We'll probably write a full solution field
        WriteStream(outfile, 0);

        // Title
        std::string title = "";
        WriteStream(outfile, title);

        // Number of variables
        WriteStream(outfile, (int)(m_coordim + var.size()));

        std::string coords[3] = { "x", "y", "z" };

        // Write names of space coordinates and variables
        for (int i = 0; i < m_coordim; ++i)
        {
            WriteStream(outfile, coords[i]);
        }

        for (int i = 0; i < var.size(); ++i)
        {
            WriteStream(outfile, var[i]);
        }
    }
    else
    {
        outfile << "Variables = x";

        if (m_coordim == 2)
        {
            outfile << ", y";
        }
        else if (m_coordim == 3)
        {
            outfile << ", y, z";
        }

        if (var.size())
        {
            outfile << var[0];
            for (int i = 1; i < var.size(); ++i)
            {
                outfile << ", " << var[i];
            }
        }

        outfile << std::endl << std::endl;
    }
}

/**
 * Write Tecplot Files Zone
 * @param   outfile    Output file name.
 * @param   expansion  Expansion that is considered
 */
void OutputTecplot::WriteTecplotZone(std::ofstream &outfile)
{
    int i, j;

    if (m_binary)
    {
        // Number of points in zone
        int nPoints = m_fields[0].num_elements();

        // Don't bother naming zone
        WriteStream(outfile, 299.0f); // Zone marker

        // Write same name as preplot
        std::string zonename = "ZONE 001";
        WriteStream(outfile, zonename);

        WriteStream(outfile, -1); // No parent zone
        WriteStream(outfile, -1); // No strand ID
        WriteStream(outfile, 0.0); // Solution time
        WriteStream(outfile, -1); // Unused, set to -1

        // Zone type: 1 = lineseg, 3 = quad, 5 = brick
        WriteStream(outfile, (int)m_zoneType);

        WriteStream(outfile, 0); // Data at nodes
        WriteStream(outfile, 0); // No 1-1 face neighbours
        WriteStream(outfile, 0); // No user-defined connections

        if (m_zoneType == eOrdered)
        {
            for (i = 0; i < m_numPoints.size(); ++i)
            {
                WriteStream(outfile, m_numPoints[i]);
            }

            for (i = m_numPoints.size(); i < 3; ++i)
            {
                WriteStream(outfile, 0);
            }
        }
        else
        {
            WriteStream(outfile, nPoints); // Total number of points
            WriteStream(outfile, m_numBlocks); // Number of blocks
        }

        WriteStream(outfile, 0); // Unused
        WriteStream(outfile, 0); // Unused
        WriteStream(outfile, 0); // Unused
        WriteStream(outfile, 0); // No auxiliary data names

        // Finalise header
        WriteStream(outfile, 357.0f);

        // Now start to write data section so that we can dump geometry
        // information

        // Data marker
        WriteStream(outfile, 299.0f);

        // Data format: either double or single depending on user
        // options
        bool useDoubles = m_config["double"].m_beenSet;

        for (int j = 0; j < m_fields.num_elements(); ++j)
        {
            WriteStream(outfile, useDoubles ? 2 : 1);
        }

        // No passive variables or variable sharing, no zone
        // connectivity sharing (we only dump one zone)
        WriteStream(outfile, 0);
        WriteStream(outfile, 0);
        WriteStream(outfile, -1);

        // Write out min/max of field data
        for (int j = 0; j < m_fields.num_elements(); ++j)
        {
            WriteStream(outfile, Vmath::Vmin(nPoints, m_fields[j], 1));
            WriteStream(outfile, Vmath::Vmax(nPoints, m_fields[j], 1));
        }

        // Now dump out field data.
        if (useDoubles)
        {
            // For doubles, we can just write data.
            for (i = 0; i < m_fields.num_elements(); ++i)
            {
                WriteStream(outfile, m_fields[i]);
            }
        }
        else
        {
            // For single precision, needs typecast first.
            for (i = 0; i < m_fields.num_elements(); ++i)
            {
                int nPoints = m_fields[0].num_elements();
                vector<float> tmp(nPoints);
                std::copy(&m_fields[i][0], &m_fields[i][nPoints-1], &tmp[0]);
                WriteStream(outfile, tmp);
            }
        }
    }
    else
    {
        // Write either points or finite element block
        if (m_zoneType != eOrdered)
        {
            outfile << "Zone, N=" << m_fields[0].num_elements() << ", E="
                    << m_numBlocks << ", F=FEBlock, ET="
                    << TecplotZoneTypeMap[m_zoneType] << std::endl;
        }
        else
        {
            std::string dirs[] = { "I", "J", "K" };
            outfile << "Zone";
            for (i = 0; i < m_numPoints.size(); ++i)
            {
                outfile << ", " << dirs[i] << "=" << m_numPoints[i];
            }
            outfile << ", F=POINT" << std::endl;
        }

        // Write out coordinates and field data
        for (j = 0; j < m_fields.num_elements(); ++j)
        {
            for (i = 0; i < m_fields[j].num_elements(); ++i)
            {
                if ((!(i % 1000)) && i)
                {
                    outfile << std::endl;
                }
                outfile << m_fields[j][i] << " ";
            }
            outfile << std::endl;
        }
    }
}

int OutputTecplot::GetNumTecplotBlocks()
{
    int returnval = 0;

    if (m_f->m_exp[0]->GetExp(0)->GetNumBases() == 1)
    {
        for (int i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
        {
            returnval += (m_f->m_exp[0]->GetExp(i)->GetNumPoints(0) - 1);
        }
    }
    else if (m_f->m_exp[0]->GetExp(0)->GetNumBases() == 2)
    {
        for (int i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
        {
            returnval += (m_f->m_exp[0]->GetExp(i)->GetNumPoints(0) - 1) *
                         (m_f->m_exp[0]->GetExp(i)->GetNumPoints(1) - 1);
        }
    }
    else
    {
        for (int i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
        {
            returnval += (m_f->m_exp[0]->GetExp(i)->GetNumPoints(0) - 1) *
                         (m_f->m_exp[0]->GetExp(i)->GetNumPoints(1) - 1) *
                         (m_f->m_exp[0]->GetExp(i)->GetNumPoints(2) - 1);
        }
    }

    return returnval;
}

// /**
//  * Write Tecplot Files Field
//  * @param   outfile    Output file name.
//  * @param   expansion  Expansion that is considered
//  */
// void OutputTecplot::WriteTecplotField(const int field, std::ofstream &outfile)
// {
//     int totpoints = m_f->m_exp[0]->GetTotPoints();

//     if (m_doError)
//     {
//         NekDouble l2err =
//             m_f->m_exp[0]->L2(m_f->m_exp[field]->UpdatePhys());

//         if (m_f->m_comm->GetRank() == 0)
//         {
//             cout << "L 2 error (variable "
//                  << m_f->m_fielddef[0]->m_fields[field] << ") : " << l2err
//                  << endl;
//         }
//     }
//     else if (!m_binary)
//     {
//         for (int i = 0; i < totpoints; ++i)
//         {
//             outfile << m_f->m_exp[field]->GetPhys()[i] << " ";
//             if ((!(i % 1000)) && i)
//             {
//                 outfile << std::endl;
//             }
//         }
//         outfile << std::endl;
//     }
//     else if (m_config["double"].m_beenSet)
//     {
//         WriteStream(outfile, m_f->m_exp[field]->UpdatePhys());
//     }
//     else
//     {
//         vector<float> tmp(totpoints);
//         Array<OneD, NekDouble> data = m_f->m_exp[field]->UpdatePhys();
//         std::copy(&data[0], &data[totpoints], &tmp[0]);
//             WriteStream(outfile, tmp);
//     }
// }

void OutputTecplot::WriteTecplotConnectivity(std::ofstream &outfile)
{
    if (m_binary)
    {
        for (int i = 0; i < m_conn.size(); ++i)
        {
            WriteStream(outfile, m_conn[i]);
        }
    }
    else
    {
        for (int i = 0; i < m_conn.size(); ++i)
        {
            const int nConn = m_conn[i].num_elements();
            for (int j = 0; j < nConn - 1; ++j)
            {
                outfile << m_conn[i][j] + 1 << " ";
            }
            outfile << m_conn[i][nConn-1] + 1 << endl;
        }
    }
}

void OutputTecplot::CalculateConnectivity()
{
    int i, j, k, l;
    int nbase = m_f->m_exp[0]->GetExp(0)->GetNumBases();
    int cnt   = 0;

    m_conn.resize(m_f->m_exp[0]->GetNumElmts());

    for (i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
    {
        cnt = m_f->m_exp[0]->GetPhys_Offset(i);

        if (nbase == 1)
        {
            int cnt2 = 0;
            int np0  = m_f->m_exp[0]->GetExp(i)->GetNumPoints(0);

            Array<OneD, int> conn(2 * (np0 - 1));

            for (k = 1; k < np0; ++k)
            {
                conn[cnt2++] = cnt + k;
                conn[cnt2++] = cnt + k - 1;
            }

            m_conn[i] = conn;
        }
        else if (nbase == 2)
        {
            int np0       = m_f->m_exp[0]->GetExp(i)->GetNumPoints(0);
            int np1       = m_f->m_exp[0]->GetExp(i)->GetNumPoints(1);
            int totPoints = m_f->m_exp[0]->GetTotPoints();
            int nPlanes   = 1;
            int cnt2      = 0;

            if (m_f->m_exp[0]->GetExpType() == MultiRegions::e3DH1D)
            {
                nPlanes = m_f->m_exp[0]->GetZIDs().num_elements();

                // default to 2D case for HalfMode when nPlanes = 1
                if (nPlanes > 1)
                {
                    totPoints = m_f->m_exp[0]->GetPlane(0)->GetTotPoints();

                    Array<OneD, int> conn(8 * (np1 - 1) * (np0 - 1) * (nPlanes - 1));

                    for (int n = 1; n < nPlanes; ++n)
                    {
                        for (j = 1; j < np1; ++j)
                        {
                            for (k = 1; k < np0; ++k)
                            {
                                conn[cnt2++] = cnt + (n - 1) * totPoints +
                                    (j - 1) * np0 + k - 1;
                                conn[cnt2++] = cnt + (n - 1) * totPoints +
                                    (j - 1) * np0 + k;
                                conn[cnt2++] = cnt + (n - 1) * totPoints +
                                    j * np0 + k;
                                conn[cnt2++] = cnt + (n - 1) * totPoints +
                                    j * np0 + k - 1;
                                conn[cnt2++] = cnt + n * totPoints +
                                    (j - 1) * np0 + k - 1;
                                conn[cnt2++] = cnt + n * totPoints +
                                    (j - 1) * np0 + k;
                                conn[cnt2++] = cnt + n * totPoints +
                                    j * np0 + k;
                                conn[cnt2++] = cnt + n * totPoints +
                                    j * np0 + k - 1;
                            }
                        }
                    }
                    m_conn[i] = conn;
                }
            }

            if (nPlanes == 1)
            {
                Array<OneD, int> conn(4 * (np0 - 1) * (np1 - 1));
                for (j = 1; j < np1; ++j)
                {
                    for (k = 1; k < np0; ++k)
                    {
                        conn[cnt2++] = cnt + (j - 1) * np0 + k - 1;
                        conn[cnt2++] = cnt + (j - 1) * np0 + k;
                        conn[cnt2++] = cnt + j * np0 + k;
                        conn[cnt2++] = cnt + j * np0 + k - 1;
                    }
                }
                m_conn[i] = conn;
            }
        }
        else if (nbase == 3)
        {
            int np0  = m_f->m_exp[0]->GetExp(i)->GetNumPoints(0);
            int np1  = m_f->m_exp[0]->GetExp(i)->GetNumPoints(1);
            int np2  = m_f->m_exp[0]->GetExp(i)->GetNumPoints(2);
            int cnt2 = 0;

            Array<OneD, int> conn(8 * (np0 - 1) * (np1 - 1) * (np2 - 1));

            for (j = 1; j < np2; ++j)
            {
                for (k = 1; k < np1; ++k)
                {
                    for (l = 1; l < np0; ++l)
                    {
                        conn[cnt2++] =
                            cnt + (j - 1) * np0 * np1 + (k - 1) * np0 + l - 1;
                        conn[cnt2++] =
                            cnt + (j - 1) * np0 * np1 + (k - 1) * np0 + l;
                        conn[cnt2++] =
                            cnt + (j - 1) * np0 * np1 + k * np0 + l;
                        conn[cnt2++] =
                            cnt + (j - 1) * np0 * np1 + k * np0 + l - 1;
                        conn[cnt2++] =
                            cnt + j * np0 * np1 + (k - 1) * np0 + l - 1;
                        conn[cnt2++] =
                            cnt + j * np0 * np1 + (k - 1) * np0 + l;
                        conn[cnt2++] =
                            cnt + j * np0 * np1 + k * np0 + l;
                        conn[cnt2++] =
                            cnt + j * np0 * np1 + k * np0 + l - 1;
                    }
                }
            }

            m_conn[i] = conn;
        }
        else
        {
            ASSERTL0(false, "Not set up for this dimension");
        }

    }
}
}
}
