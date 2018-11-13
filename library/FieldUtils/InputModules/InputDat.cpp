////////////////////////////////////////////////////////////////////////////////
//
//  File: InputDat.cpp
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
//  Description: Read tecplot dat for isontours in 2D triangular FE block format
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/PtsIO.h>

#include <tinyxml.h>

#include "InputDat.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey InputDat::m_className[1] = {
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "dat"),
        InputDat::create,
        "Reads Tecplot dat file for FE block triangular format.")
};

/**
 * @brief Set up InputDat object.
 *
 */
InputDat::InputDat(FieldSharedPtr f) : InputModule(f)
{
    m_allowedFiles.insert("dat");
}

/**
 *
 */
InputDat::~InputDat()
{
}

/**
 *
 */
void InputDat::Process(po::variables_map &vm)
{
    string line;
    std::ifstream datFile;

    // Open the file stream.
    string fname = m_f->m_inputfiles["dat"][0];

    datFile.open(fname.c_str());
    if (!datFile.good())
    {
        cerr << "Error opening file: " << fname << endl;
        abort();
    }

    // read variables
    // currently assume there are x y and z coordinates
    int dim = 3;
    vector<string> fieldNames;
    while (!datFile.eof())
    {
        getline(datFile, line);
        string linetest = line;
        boost::to_upper(linetest);
        if (linetest.find("VARIABLES") != string::npos)
        {
            std::size_t pos = line.find('=');
            pos++;

            // note this expects a comma separated list but
            // does not work for white space separated lists!
            bool valid = ParseUtils::GenerateVector(
                line.substr(pos), fieldNames);
            ASSERTL0(valid, "Unable to process list of field variable in "
                            " VARIABLES list:  " +
                                line.substr(pos));

            // remove coordinates from fieldNames
            fieldNames.erase(fieldNames.begin(), fieldNames.begin() + dim);

            break;
        }
    }


    // set up basic parameters
    int nfields = fieldNames.size();
    int totvars = dim + nfields;
    Array<OneD, Array<OneD, NekDouble> > pts(totvars);
    vector<Array<OneD, int> > ptsConn;
    vector< int> ptsPerEdge;
    LibUtilities::PtsType ptype;
    
    // read zones
    while (!datFile.eof())
    {
        getline(datFile, line);
        string linetest = line;
        boost::to_upper(linetest);
        // Read data from all zones
        if ((linetest.find("ZONE") != string::npos))
        {
            if(line.find(""))
            {
                ReadTecplotFEBlockZone(datFile, line, pts, ptsConn);
                ptype = LibUtilities::ePtsTriBlock;
            }
            else if(line.find("POINT"))
            {
                ReadTecplotDatZone(datFile, line, ptype, ptsPerEdge,pts);
            }
            else
            {
                ASSERTL0(false,"Unrecognised zone format");
            }
        }
    }

    datFile.close();

    // These are set after reading in case of multiple files.
    m_f->m_fieldPts = MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
        dim, fieldNames, pts);
    m_f->m_fieldPts->SetPtsType(ptype);
    if(ptype == LibUtilities::ePtsTriBlock)
    {
        m_f->m_fieldPts->SetConnectivity(ptsConn);
    }
    else
    {
        m_f->m_fieldPts->SetPointsPerEdge(ptsPerEdge);
    }

    // save field names
    m_f->m_variables = fieldNames;
}

/**
 *
 */
void InputDat::ReadTecplotFEBlockZone(std::ifstream &datFile,
                                      string &line,
                                      Array<OneD, Array<OneD, NekDouble> > &pts,
                                      vector<Array<OneD, int> > &ptsConn)
{
    ASSERTL0(line.find("FEBlock") != string::npos,
             "Routine only set up for FEBLock format");
    ASSERTL0(line.find("ET") != string::npos, "Routine only set up TRIANLES");

    // read the number of nodes
    stringstream s;
    string tag;
    int start, end;

    s.clear();
    s.str(line);
    tag = s.str();

    // read the number of vertices
    start     = tag.find("N=");
    end       = tag.find_first_of(',', start);
    int nvert = atoi(tag.substr(start + 2, end).c_str());

    // read the number of elements
    start     = tag.find("E=");
    end       = tag.find_first_of(',', start);
    int nelmt = atoi(tag.substr(start + 2, end).c_str());

    // set-up or extend m_pts array;
    int norigpts  = pts[0].num_elements();
    int totfields = pts.num_elements();
    Array<OneD, Array<OneD, NekDouble> > origpts(totfields);
    for (int i = 0; i < totfields; ++i)
    {
        origpts[i] = pts[i];
        pts[i]     = Array<OneD, NekDouble>(norigpts + nvert);
    }

    NekDouble value;
    for (int n = 0; n < totfields; ++n)
    {

        for (int i = 0; i < norigpts; ++i)
        {
            pts[n][i] = origpts[n][i];
        }
        for (int i = 0; i < nvert; ++i)
        {
            datFile >> value;
            pts[n][norigpts + i] = value;
        }
    }

    // read connectivity and add to list
    int intvalue;
    Array<OneD, int> conn(3 * nelmt);
    for (int i = 0; i < 3 * nelmt; ++i)
    {
        datFile >> intvalue;
        intvalue -= 1; // decrement intvalue by 1 for c array convention
        conn[i] = norigpts + intvalue;
    }
    ptsConn.push_back(conn);

    getline(datFile, line);
}

    
void InputDat::ReadTecplotDatZone(std::ifstream &datFile,
                                  string &line,
                                  LibUtilities::PtsType  &ptype,
                                  vector<int> &ptsperedge,
                                  Array<OneD, Array<OneD, NekDouble> > &pts)
{
    // read the number of nodes
    stringstream s;
    string tag;
    int start, end;
    int nj = 1;
    int nk = 1;
    
    s.clear();
    s.str(line);
    tag = s.str();

    // read the number of I Points
    start  = tag.find("I=");
    end    = tag.find_first_of(',', start);
    int ni = atoi(tag.substr(start + 2, end).c_str());
    ptsperedge.push_back(ni);
    ptype = LibUtilities::ePtsLine; 
    
    // read the number of J points
    start     = tag.find("J=");
    if(start != std::string::npos)
    {
        end = tag.find_first_of(',', start);
        nj  = atoi(tag.substr(start + 2, end).c_str());
        ptsperedge.push_back(nj);
        ptype = LibUtilities::ePtsPlane; 
    }

    // read the number of K points
    start     = tag.find("K=");
    if(start != std::string::npos)
    {
        end = tag.find_first_of(',', start);
        nk  = atoi(tag.substr(start + 2, end).c_str());
        ptsperedge.push_back(nk);
        ptype = LibUtilities::ePtsBox; 
    }

    // set-up or extend m_pts array;
    int npts  = ni*nj*nk;
    int totfields = pts.num_elements();
    for (int i = 0; i < totfields; ++i)
    {
        pts[i] = Array<OneD, NekDouble>(npts);
    }
    
    NekDouble value;
    for (int i = 0; i < npts; ++i)
    {
        for (int n = 0; n < totfields; ++n)
        {
            datFile >> value;
            pts[n][i] = value;
        }
    }
}
}
}
