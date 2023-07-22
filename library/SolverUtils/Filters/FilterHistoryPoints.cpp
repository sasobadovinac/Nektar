///////////////////////////////////////////////////////////////////////////////
//
// File: FilterHistoryPoints.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Outputs values at specific points during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
using namespace std;

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <SolverUtils/Filters/FilterHistoryPoints.h>
#include <boost/core/ignore_unused.hpp>
#include <boost/format.hpp>

namespace Nektar
{
namespace SolverUtils
{
std::string FilterHistoryPoints::className =
    GetFilterFactory().RegisterCreatorFunction("HistoryPoints",
                                               FilterHistoryPoints::create);

/**
 *
 */
FilterHistoryPoints::FilterHistoryPoints(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem> &pEquation, const ParamMap &pParams)
    : Filter(pSession, pEquation)
{
    // OutputFile
    auto it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        m_outputFile = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        m_outputFile = it->second;
    }
    if (m_outputFile.length() >= 4 &&
        m_outputFile.substr(m_outputFile.length() - 4) == ".his")
    {
        m_outputFile = m_outputFile.substr(0, m_outputFile.length() - 4);
    }

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        m_outputFrequency = round(equ.Evaluate());
    }

    // output all data into one file
    it = pParams.find("OutputOneFile");
    if (it == pParams.end())
    {
        m_outputOneFile = true;
    }
    else
    {
        std::string sOption = it->second.c_str();
        m_outputOneFile     = (boost::iequals(sOption, "true")) ||
                          (boost::iequals(sOption, "yes"));
    }

    // OutputPlane
    m_session->MatchSolverInfo("Homogeneous", "1D", m_isHomogeneous1D, false);
    if (m_isHomogeneous1D)
    {
        it = pParams.find("OutputPlane");
        if (it == pParams.end())
        {
            m_outputPlane = -1;
        }
        else
        {
            LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
            m_outputPlane = round(equ.Evaluate());
        }

        it = pParams.find("WaveSpace");
        if (it == pParams.end())
        {
            m_waveSpace = false;
        }
        else
        {
            std::string sOption = it->second.c_str();
            m_waveSpace         = (boost::iequals(sOption, "true")) ||
                          (boost::iequals(sOption, "yes"));
        }
    }

    // Points
    if (pParams.end() != (it = pParams.find("Points")))
    {
        m_pointNumMap[LibUtilities::ePtsFile] = Array<OneD, int>(1, 0);
        m_historyPointStream.str(it->second);
    }
    else if (pParams.end() != (it = pParams.find("line")))
    {
        vector<NekDouble> values;
        ASSERTL0(ParseUtils::GenerateVector(it->second, values),
                 "Failed to interpret line string");

        ASSERTL0(values.size() > 2, "line string should contain 2*Dim+1 values "
                                    "N,x0,y0,z0,x1,y1,z1");

        double tmp;
        ASSERTL0(std::modf(values[0], &tmp) == 0.0, "N is not an integer");
        ASSERTL0(values[0] > 1, "N is not a valid number");

        int dim  = (values.size() - 1) / 2;
        int npts = values[0];

        Array<OneD, int> num(1, npts);
        Array<OneD, NekDouble> delta(6, 0.);
        for (int i = 0; i < dim; ++i)
        {
            delta[i + 0] = values[i + 1];
            delta[i + 3] = (values[dim + i + 1] - values[i + 1]) / (npts - 1);
        }
        m_pointDatMap[LibUtilities::ePtsLine] = delta;
        m_pointNumMap[LibUtilities::ePtsLine] = num;
    }
    else if (pParams.end() != (it = pParams.find("plane")))
    {
        vector<NekDouble> values;
        ASSERTL0(ParseUtils::GenerateVector(it->second, values),
                 "Failed to interpret plane string");

        ASSERTL0(values.size() > 9,
                 "plane string should contain 4 Dim+2 values "
                 "N1,N2,x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3");

        double tmp;
        ASSERTL0(std::modf(values[0], &tmp) == 0.0, "N1 is not an integer");
        ASSERTL0(std::modf(values[1], &tmp) == 0.0, "N2 is not an integer");

        ASSERTL0(values[0] > 1, "N1 is not a valid number");
        ASSERTL0(values[1] > 1, "N2 is not a valid number");

        int dim = (values.size() - 2) / 4;

        Array<OneD, int> npts(3);
        npts[0] = values[0];
        npts[1] = values[0] * values[1];
        npts[2] = values[1];

        Array<OneD, NekDouble> delta(12, 0.);
        for (int i = 0; i < dim; ++i)
        {
            delta[i + 0] = values[2 + i];
            delta[i + 3] = values[2 + 3 * dim + i];
            delta[i + 6] = (values[2 + 1 * dim + i] - values[2 + 0 * dim + i]) /
                           (values[0] - 1);
            delta[i + 9] = (values[2 + 2 * dim + i] - values[2 + 3 * dim + i]) /
                           (values[0] - 1);
        }
        m_pointDatMap[LibUtilities::ePtsPlane] = delta;
        m_pointNumMap[LibUtilities::ePtsPlane] = npts;
    }
    else if (pParams.end() != (it = pParams.find("box")))
    {
        vector<NekDouble> values;
        ASSERTL0(ParseUtils::GenerateVector(it->second, values),
                 "Failed to interpret box string");

        ASSERTL0(values.size() == 9, "box string should contain 9 values "
                                     "N1,N2,N3,xmin,xmax,ymin,ymax,zmin,zmax");

        int dim = 3;
        Array<OneD, int> npts(3);
        npts[0] = values[0];
        npts[1] = values[0] * values[1];
        npts[2] = values[0] * values[1] * values[2];

        Array<OneD, NekDouble> delta(6, 0.);
        for (int i = 0; i < dim; ++i)
        {
            delta[i + 0] = values[3 + 2 * i];
            delta[i + 3] =
                (values[4 + 2 * i] - values[3 + 2 * i]) / (values[i] - 1);
        }
        m_pointDatMap[LibUtilities::ePtsBox] = delta;
        m_pointNumMap[LibUtilities::ePtsBox] = npts;
    }
    else
    {
        ASSERTL0(false, "Missing parameter 'Points'.");
    }
}

/**
 *
 */
FilterHistoryPoints::~FilterHistoryPoints()
{
}

bool FilterHistoryPoints::GetPoint(Array<OneD, NekDouble> gloCoord, int I)
{
    if (m_pointNumMap.count(LibUtilities::ePtsFile))
    {
        m_historyPointStream >> gloCoord[0] >> gloCoord[1] >> gloCoord[2];
        return !m_historyPointStream.fail();
    }
    else if (m_pointNumMap.count(LibUtilities::ePtsLine))
    {
        if (I >= m_pointNumMap[LibUtilities::ePtsLine][0])
        {
            return false;
        }
        Array<OneD, NekDouble> values = m_pointDatMap[LibUtilities::ePtsLine];
        Array<OneD, NekDouble> delta =
            m_pointDatMap[LibUtilities::ePtsLine] + 3;
        for (int n = 0; n < 3; ++n)
        {
            gloCoord[n] = values[n] + I * delta[n];
        }
        return true;
    }
    else if (m_pointNumMap.count(LibUtilities::ePtsPlane))
    {
        if (I >= m_pointNumMap[LibUtilities::ePtsPlane][1])
        {
            return false;
        }
        Array<OneD, NekDouble> values = m_pointDatMap[LibUtilities::ePtsPlane];
        Array<OneD, NekDouble> delta =
            m_pointDatMap[LibUtilities::ePtsPlane] + 6;
        Array<OneD, int> i(2);
        i[1]      = I / m_pointNumMap[LibUtilities::ePtsPlane][0];
        i[0]      = I - i[1] * m_pointNumMap[LibUtilities::ePtsPlane][0];
        double n1 = -1. + (NekDouble)m_pointNumMap[LibUtilities::ePtsPlane][2];
        for (int n = 0; n < 3; ++n)
        {
            gloCoord[n] = (values[n] + i[0] * delta[n]) * (1.0 - i[1] / n1) +
                          (values[3 + n] + i[0] * delta[3 + n]) * (i[1] / n1);
        }
        return true;
    }
    else if (m_pointNumMap.count(LibUtilities::ePtsBox))
    {
        if (I >= m_pointNumMap[LibUtilities::ePtsBox][2])
        {
            return false;
        }
        Array<OneD, NekDouble> values = m_pointDatMap[LibUtilities::ePtsBox];
        Array<OneD, NekDouble> delta = m_pointDatMap[LibUtilities::ePtsBox] + 3;
        Array<OneD, int> i(3);
        i[2] = I / m_pointNumMap[LibUtilities::ePtsBox][1];
        I -= i[2] * m_pointNumMap[LibUtilities::ePtsBox][1];
        i[1] = I / m_pointNumMap[LibUtilities::ePtsBox][0];
        i[0] = I - i[1] * m_pointNumMap[LibUtilities::ePtsBox][0];
        for (int n = 0; n < 3; ++n)
        {
            gloCoord[n] = values[n] + i[n] * delta[n];
        }
        return true;
    }
    return false;
}

/**
 *
 */
void FilterHistoryPoints::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    ASSERTL0(!m_historyPointStream.fail(), "No history points in stream.");

    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

    // All processes read all the history points
    Array<OneD, NekDouble> gloCoord(3, 0.0);
    std::vector<Array<OneD, const NekDouble>> GloCoords;
    const int coordim = pFields[0]->GetGraph()->GetSpaceDimension();

    int i = -1;
    while (GetPoint(gloCoord, ++i))
    {
        // Overwrite gloCoord[2] for 3DH1D using m_outputPlane
        // if it is defined
        if (m_isHomogeneous1D && (m_outputPlane != -1 || m_waveSpace))
        {
            int nplanes    = pFields[0]->GetHomogeneousBasis()->GetZ().size();
            NekDouble lhom = pFields[0]->GetHomoLen();
            int plane;
            if (m_outputPlane != -1)
            {
                plane = m_outputPlane;
            }
            else
            {
                // Pick plane near the point
                plane = round((gloCoord[coordim] * nplanes) / lhom);
            }
            if (m_waveSpace)
            {
                m_planeIDs.push_back(plane);
            }
            NekDouble Z = (pFields[0]->GetHomogeneousBasis()->GetZ())[plane];
            Z           = (Z + 1) * lhom / 2;
            if (fabs(gloCoord[coordim] - Z) >
                    NekConstants::kVertexTheSameDouble &&
                vComm->GetRank() == 0)
            {
                cout << "Resetting History point from z = " << gloCoord[coordim]
                     << " to z = " << Z << endl;
            }
            gloCoord[coordim] = Z;
        }

        // Save a copy of the global coordinate
        GloCoords.push_back(
            Array<OneD, const NekDouble>(gloCoord.size(), gloCoord.data()));
    }

    // Check which history points this process is responsible for
    const int vRank       = vComm->GetRank();
    const int vColumnRank = vComm->GetColumnComm()->GetRank();
    const int vRowRank    = vComm->GetRowComm()->GetRank();
    m_historyPointsSize   = GloCoords.size();
    Array<OneD, int> procList(m_historyPointsSize, -1);
    Array<OneD, int> idList(m_historyPointsSize, -1);

    Array<OneD, NekDouble> locCoord(pFields[0]->GetShapeDimension(), 0.);
    std::vector<Array<OneD, const NekDouble>> LocCoords;

    // For each history point, find the element that contains it. If no element
    // is available on this process, idList[i] = -1 and procList[i] = -1
    for (i = 0; i < m_historyPointsSize; ++i)
    {
        // Create a copy of the global coordinate
        gloCoord = GloCoords[i];

        // Determine the expansion and local coordinates
        if (m_isHomogeneous1D)
        {
            idList[i] = pFields[0]->GetPlane(0)->GetExpIndex(
                gloCoord, locCoord, NekConstants::kGeomFactorsTol);
        }
        else
        {
            idList[i] = pFields[0]->GetExpIndex(gloCoord, locCoord,
                                                NekConstants::kGeomFactorsTol);
        }

        for (int j = 0; j < locCoord.size(); ++j)
        {
            locCoord[j] = std::max(locCoord[j], -1.0);
            locCoord[j] = std::min(locCoord[j], 1.0);
        }

        // If this process owns an element that contains the history point, we
        // flag this by updating the corresponding element in procList
        if (idList[i] != -1 && vColumnRank == 0)
        {
            procList[i] = vRank;
        }

        // Save a copy of the local coordinates
        LocCoords.push_back(
            Array<OneD, const NekDouble>(locCoord.size(), locCoord.data()));
    }

    // Clear list before populating it
    m_historyList.clear();

    // Share process ranks among all processes. If multiple processes found
    // that they own a history point, procList will only contain the rank of
    // the process with the largest rank
    vComm->AllReduce(procList, LibUtilities::ReduceMax);
    for (i = 0; i < m_historyPointsSize; ++i)
    {
        // If process owns a history point, it saves a reference to its global
        // and local coordinates. Note that, if multiple processes have found
        // that they own a history point, which can happen if the point is
        // located on a partition boundary, only the process with the highest
        // rank will store the point.
        if (procList[i] == vRowRank)
        {
            // Note that the reference counting in the Array class prevents the
            // data stored in GloCoords and LocCoords from being deallocated
            // when they go out of scope, since m_historyList will increase the
            // reference counter by 1 when it gets a reference to the data.
            m_historyList.push_back(
                std::make_tuple(GloCoords[i], LocCoords[i], i, idList[i]));
        }
    }

    // Root process must have access to all coordinates for writing the file
    if (vComm->TreatAsRankZero())
    {
        m_historyPoints =
            Array<OneD, Array<OneD, const NekDouble>>(m_historyPointsSize);

        for (i = 0; i < m_historyPointsSize; ++i)
        {
            // Create a reference to the underlying array. Note that the data
            // stored in GloCoords won't be deallocated because the Array class
            // has a reference counter.
            m_historyPoints[i] = GloCoords[i];
        }
    }

    // Check that each history point is allocated to a process
    if (vComm->TreatAsRankZero())
    {
        for (i = 0; i < m_historyPointsSize; ++i)
        {
            gloCoord = m_historyPoints[i];

            // Write an error if no process owns history point
            ASSERTL0(procList[i] != -1,
                     "History point " +
                         boost::lexical_cast<std::string>(gloCoord[0]) + ", " +
                         boost::lexical_cast<std::string>(gloCoord[1]) + ", " +
                         boost::lexical_cast<std::string>(gloCoord[2]) +
                         " cannot be found in the mesh.");
        }

        m_session->MatchSolverInfo("Driver", "Adaptive", m_adaptive, false);
    }
    v_Update(pFields, time);
}

/**
 *
 */
void FilterHistoryPoints::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Only output every m_outputFrequency.
    if ((m_index++) % m_outputFrequency)
    {
        return;
    }

    const int numFields = pFields.size();
    const int coordim   = pFields[0]->GetGraph()->GetSpaceDimension();
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
    Array<OneD, NekDouble> data(m_historyPointsSize * numFields, 0.0);
    Array<OneD, NekDouble> physvals;
    Array<OneD, NekDouble> locCoord;

    // Pull out data values field by field
    Array<OneD, NekDouble> gloCoord(3, 0.0);
    for (int j = 0; j < numFields; ++j)
    {
        if (m_isHomogeneous1D)
        {
            Array<OneD, const unsigned int> planes = pFields[j]->GetZIDs();
            int nPlanes    = pFields[j]->GetHomogeneousBasis()->GetZ().size();
            NekDouble lHom = pFields[j]->GetHomoLen();
            for (auto &x : m_historyList)
            {
                ASSERTL0(pFields[j]->GetWaveSpace(),
                         "HistoryPoints in Homogeneous1D require that solution "
                         "is in wavespace");
                gloCoord        = std::get<0>(x);
                locCoord        = std::get<1>(x);
                const int indx  = std::get<2>(x);
                const int expId = std::get<3>(x);

                NekDouble value = 0.0;
                NekDouble BetaT =
                    2. * M_PI * fmod(gloCoord[coordim], lHom) / lHom;
                for (size_t n = 0; n < planes.size(); ++n)
                {
                    if (m_waveSpace && planes[n] != m_planeIDs[indx])
                    {
                        continue;
                    }
                    physvals = pFields[j]->GetPlane(n)->UpdatePhys() +
                               pFields[j]->GetPhys_Offset(expId);

                    // transform elemental data if required.
                    if (pFields[j]->GetPhysState() == false)
                    {
                        pFields[j]->GetPlane(n)->GetExp(expId)->BwdTrans(
                            pFields[j]->GetPlane(n)->GetCoeffs() +
                                pFields[j]->GetCoeff_Offset(expId),
                            physvals);
                    }
                    // Interpolate data
                    NekDouble coeff =
                        pFields[j]->GetPlane(n)->GetExp(expId)->StdPhysEvaluate(
                            locCoord, physvals);

                    if (m_waveSpace)
                    {
                        value = coeff;
                    }
                    else
                    {
                        if (planes[n] == 0)
                        {
                            value += coeff;
                        }
                        else if (planes[n] == 1)
                        {
                            value += cos(0.5 * nPlanes * BetaT) * coeff;
                        }
                        else if (planes[n] % 2 == 0)
                        {
                            NekDouble phase = (planes[n] >> 1) * BetaT;
                            value += cos(phase) * coeff;
                        }
                        else
                        {
                            NekDouble phase = (planes[n] >> 1) * BetaT;
                            value += -sin(phase) * coeff;
                        }
                    }
                }
                // store data
                data[indx * numFields + j] = value;
            }
        }
        else
        {
            for (auto &x : m_historyList)
            {
                locCoord        = std::get<1>(x);
                const int indx  = std::get<2>(x);
                const int expId = std::get<3>(x);

                physvals = pFields[j]->UpdatePhys() +
                           pFields[j]->GetPhys_Offset(expId);

                // transform elemental data if required.
                if (pFields[j]->GetPhysState() == false)
                {
                    pFields[j]->GetExp(expId)->BwdTrans(
                        pFields[j]->GetCoeffs() +
                            pFields[j]->GetCoeff_Offset(expId),
                        physvals);
                }

                // interpolate point
                data[indx * numFields + j] =
                    pFields[j]->GetExp(expId)->StdPhysEvaluate(locCoord,
                                                               physvals);
            }
        }
    }

    // Exchange history data
    // This could be improved to reduce communication but works for now
    vComm->AllReduce(data, LibUtilities::ReduceSum);

    // TODO: Why not only call this routine if we are rank 0?
    WriteData(vComm->GetRank(), data, numFields, time);
}

void FilterHistoryPoints::WriteData(const int &rank,
                                    const Array<OneD, NekDouble> &data,
                                    const int &numFields, const NekDouble &time)
{
    // Only the root process writes out history data
    if (rank == 0)
    {
        Array<OneD, NekDouble> gloCoord(3, 0.0);
        if (!m_outputOneFile || m_index == 1)
        {
            std::stringstream vOutputFilename;
            if (m_outputOneFile)
            {
                vOutputFilename << m_outputFile << ".his";
            }
            else
            {
                vOutputFilename << m_outputFile << "_" << m_outputIndex
                                << ".his";
            }
            ++m_outputIndex;
            if (m_adaptive)
            {
                m_outputStream.open(vOutputFilename.str().c_str(),
                                    ofstream::app);
            }
            else
            {
                m_outputStream.open(vOutputFilename.str().c_str());
            }
            m_outputStream << "# History data for variables (:";

            for (int i = 0; i < numFields; ++i)
            {
                m_outputStream << m_session->GetVariable(i) << ",";
            }

            if (m_isHomogeneous1D)
            {
                m_outputStream << ") at points:" << endl;
            }
            else
            {
                m_outputStream << ") at points:" << endl;
            }

            for (int i = 0; i < m_historyPoints.size(); ++i)
            {
                gloCoord = m_historyPoints[i];

                m_outputStream << "# " << boost::format("%6.0f") % i;
                m_outputStream << " " << boost::format("%25.19e") % gloCoord[0];
                m_outputStream << " " << boost::format("%25.19e") % gloCoord[1];
                m_outputStream << " " << boost::format("%25.19e") % gloCoord[2];
                m_outputStream << endl;
            }

            if (m_isHomogeneous1D)
            {
                if (m_waveSpace)
                {
                    m_outputStream << "# (in Wavespace)" << endl;
                }
            }
        }

        // Write data values point by point
        for (int k = 0; k < m_historyPoints.size(); ++k)
        {
            m_outputStream << boost::format("%25.19e") % time;
            for (int j = 0; j < numFields; ++j)
            {
                m_outputStream
                    << " "
                    << boost::format("%25.19e") % data[k * numFields + j];
            }
            m_outputStream << endl;
        }

        if (!m_outputOneFile)
        {
            m_outputStream.close();
            ;
        }
    }
}

/**
 *
 */
void FilterHistoryPoints::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(time);

    if (pFields[0]->GetComm()->GetRank() == 0 && m_outputOneFile)
    {
        m_outputStream.close();
    }
}

/**
 *
 */
bool FilterHistoryPoints::v_IsTimeDependent()
{
    return true;
}
} // namespace SolverUtils
} // namespace Nektar
