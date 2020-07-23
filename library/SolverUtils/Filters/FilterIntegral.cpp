///////////////////////////////////////////////////////////////////////////////
//
// File FilterIntegral.cpp
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
// Description: Outputs integrals of fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <SolverUtils/Filters/FilterIntegral.h>

#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string.hpp>

namespace Nektar
{
namespace SolverUtils
{
std::string FilterIntegral::className =
    GetFilterFactory().RegisterCreatorFunction("Integral",
                                               FilterIntegral::create);

FilterIntegral::FilterIntegral(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem> &pEquation, const ParamMap &pParams)
    : Filter(pSession, pEquation)
{
    std::string outName;

    // OutputFile
    auto it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        outName = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Empty parameter 'OutputFile'.");
        outName = it->second;
    }
    outName += ".int";

    // Composites (to calculate integrals on)
    it = pParams.find("Composites");
    ASSERTL0(it != pParams.end(), "Missing parameter 'Composites'.");
    ASSERTL0(it->second.length() > 0, "Empty parameter 'Composites'.");

    std::vector<std::string> splitComposite;
    boost::split(m_splitCompString, it->second, boost::is_any_of(","),
                 boost::token_compress_on);

    for (auto &comp : m_splitCompString)
    {
        boost::trim(comp);
        size_t first   = comp.find_first_of('[') + 1;
        size_t last    = comp.find_last_of(']') - 1;
        auto tmpString = comp.substr(first, last - first + 1);

        std::vector<unsigned int> tmpVec;
        bool parseGood = ParseUtils::GenerateSeqVector(tmpString, tmpVec);

        ASSERTL0(parseGood && !tmpVec.empty(),
                 "Unable to read composite regions index range for "
                 "FilterIntegral: " +
                     comp);

        m_compVector.emplace_back(tmpVec);
    }

    // OutputPrecision
    size_t precision;
    it = pParams.find("OutputPrecision");
    if (it == pParams.end())
    {
        precision = 7;
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Empty parameter 'OutputPrecision'.");
        precision = std::stoi(it->second);
    }

    // Lock equation system pointer
    auto equationSys = m_equ.lock();
    ASSERTL0(equationSys, "Weak pointer expired");

    m_numVariables = equationSys->GetNvariables();

    m_comm = pSession->GetComm();
    if (m_comm->GetRank() == 0)
    {
        m_outFile.open(outName);
        ASSERTL0(m_outFile.good(), "Unable to open: '" + outName + "'");
        m_outFile.setf(std::ios::scientific, std::ios::floatfield);
        m_outFile.precision(precision);
        m_outFile << "#Time";

        for (auto &compName : m_splitCompString)
        {
            for (size_t j = 0; j < m_numVariables; ++j)
            {
                std::string varName = equationSys->GetVariable(j);
                m_outFile << " " + compName + "_" + varName + "_integral";
            }
        }
        m_outFile << std::endl;
    }

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = 1;
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Empty parameter 'OutputFrequency'.");
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        m_outputFrequency = round(equ.Evaluate());
    }
}

void FilterIntegral::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{

    // Create map from element ID -> expansion list ID
    auto expList   = pFields[0]->GetExp();
    auto meshGraph = pFields[0]->GetGraph();

    for (size_t i = 0; i < expList->size(); ++i)
    {
        auto exp                                           = (*expList)[i];
        m_geomElmtIdToExpId[exp->GetGeom()->GetGlobalID()] = i;
    }

    // Create a map from geom ID -> trace expansion list ID
    std::map<size_t, size_t> geomIdToTraceId;
    auto trace = pFields[0]->GetTrace()->GetExp();
    for (size_t i = 0; i < trace->size(); ++i)
    {
        auto exp                                       = (*trace)[i];
        geomIdToTraceId[exp->GetGeom()->GetGlobalID()] = i;
    }

    // Get comp list dimension from first composite & element
    auto composites = pFields[0]->GetGraph()->GetComposites();
    size_t meshDim  = pFields[0]->GetGraph()->GetMeshDimension();

    for (int i = 0; i < m_compVector.size(); ++i)
    {
        // Check composite is present in the rank
        if (composites.find(m_compVector[i][0]) == composites.end())
        {
            continue;
        }

        std::vector<std::shared_ptr<SpatialDomains::Geometry>> geomVec =
            composites[m_compVector[i][0]]->m_geomVec;
        size_t dim =
            composites[m_compVector[i][0]]->m_geomVec[0]->GetShapeDim();

        // Vector of all geometry IDs contained within the composite list
        std::vector<size_t> compGeomIds;
        for (auto compNum : m_compVector[i])
        {
            ASSERTL0(composites.find(compNum) != composites.end(),
                     "In FilterIntegral defined composite C[" +
                         std::to_string(compNum) +
                         "] does not exist in the mesh.")

            auto compGeom = composites[compNum]->m_geomVec;

            for (auto &geom : compGeom)
            {
                compGeomIds.emplace_back(geom->GetGlobalID());
            }

            // Only check first element in each comp for dimension
            ASSERTL0(
                dim == compGeom[0]->GetShapeDim(),
                "Differing geometry dimensions specified in FilterIntegral '" +
                    m_splitCompString[i] + "'.");
        }

        std::vector<std::pair<LocalRegions::ExpansionSharedPtr, int>>
            tmpCompExp(compGeomIds.size());

        // If dimension of composite == dimension of mesh then we only need the
        // expansion of the element
        if (dim == pFields[0]->GetShapeDimension())
        {
            for (size_t j = 0; j < compGeomIds.size(); ++j)
            {
                tmpCompExp[j] = std::make_pair(
                    (*expList)[m_geomElmtIdToExpId[compGeomIds[j]]], -1);
            }
        }
        // however if the dimension is less we need the expansion of the element
        // containing the global composite geometry and the face/edge local ID
        // within that. 3D mesh -> 2D, 2D -> 1D.
        // @TODO: Restructure with the new dimension independent functions
        //        and check all is correct with this filter.
        else if (meshDim == 3 && dim == 2)
        {
            for (size_t j = 0; j < compGeomIds.size(); ++j)
            {
                LocalRegions::ExpansionSharedPtr exp =
                    (*trace)[geomIdToTraceId[compGeomIds[j]]];
                LocalRegions::Expansion2DSharedPtr exp2D =
                    std::dynamic_pointer_cast<LocalRegions::Expansion2D>(exp);

                LocalRegions::ExpansionSharedPtr leftAdjElmtExp =
                    std::dynamic_pointer_cast<LocalRegions::Expansion>(
                        exp2D->GetLeftAdjacentElementExp());
                int leftAdjElmtFace = exp2D->GetLeftAdjacentElementTrace();

                tmpCompExp[j] = std::make_pair(leftAdjElmtExp, leftAdjElmtFace);
            }
        }
        else if (meshDim == 2 && dim == 1)
        {
            for (size_t j = 0; j < compGeomIds.size(); ++j)
            {
                LocalRegions::ExpansionSharedPtr exp =
                    (*trace)[geomIdToTraceId[compGeomIds[j]]];
                LocalRegions::Expansion1DSharedPtr exp1D =
                    std::dynamic_pointer_cast<LocalRegions::Expansion1D>(exp);

                LocalRegions::ExpansionSharedPtr leftAdjElmtExp =
                    std::dynamic_pointer_cast<LocalRegions::Expansion>(
                        exp1D->GetLeftAdjacentElementExp());
                int leftAdjElmtEdge = exp1D->GetLeftAdjacentElementTrace();

                tmpCompExp[j] = std::make_pair(leftAdjElmtExp, leftAdjElmtEdge);
            }
        }
        else
        {
            ASSERTL0(false,
                     "FilterIntegral: Only composite dimensions equal to or "
                     "one lower than the mesh dimension are supported.")
        }

        m_compExpMap[i] = tmpCompExp;
    }

    v_Update(pFields, time);
}

void FilterIntegral::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    if (m_index++ % m_outputFrequency > 0)
    {
        return;
    }

    if (m_comm->GetRank() == 0)
    {
        m_outFile << time;
    }

    for (size_t j = 0; j < m_compVector.size(); ++j)
    {
        for (size_t i = 0; i < m_numVariables; ++i)
        {
            Array<OneD, const NekDouble> phys;
            phys = pFields[i]->GetPhys();

            NekDouble sum = 0.0;
            NekDouble c   = 0.0;

            // Check if composite is on the rank
            if (m_compExpMap.find(j) != m_compExpMap.end())
            {
                auto compExp   = m_compExpMap[j];
                size_t dim     = compExp[0].first->GetGeom()->GetShapeDim();
                size_t meshDim = pFields[i]->GetGraph()->GetMeshDimension();

                // Evaluate integral using improved Kahan–Babuška summation
                // algorithm to reduce numerical error from adding floating
                // points
                for (auto &expPair : compExp)
                {
                    NekDouble input = 0;
                    auto exp        = expPair.first;

                    if (meshDim == dim)
                    {
                        size_t offset = pFields[i]->GetPhys_Offset(
                            m_geomElmtIdToExpId[exp->GetGeom()->GetGlobalID()]);
                        input = exp->Integral(phys + offset);
                    }
                    else if (meshDim == 3 && dim == 2)
                    {
                        Array<OneD, NekDouble> facePhys;
                        exp->GetTracePhysVals(expPair.second, exp, phys,
                                             facePhys);
                        input =
                            pFields[i]
                                ->GetTrace()
                                ->GetExp(exp->GetGeom()->GetTid(expPair.second))
                                ->Integral(facePhys);
                    }
                    else if (meshDim == 2 && dim == 1)
                    {
                        Array<OneD, NekDouble> edgePhys;
                        exp->GetTracePhysVals(expPair.second, exp, phys, edgePhys);
                        input =
                            pFields[i]
                                ->GetTrace()
                                ->GetExp(exp->GetGeom()->GetTid(expPair.second))
                                ->Integral(edgePhys);
                    }
                    else
                    {
                        ASSERTL0(false,
                                 "FilterIntegral: Only composite dimensions "
                                 "equal to or one lower than the mesh "
                                 "dimension are supported.")
                    }

                    NekDouble t = sum + input;
                    c += fabs(sum) >= fabs(input) ? (sum - t) + input
                                                  : (input - t) + sum;
                    sum = t;
                }
            }

            // Sum integral values from all ranks
            Array<OneD, NekDouble> sumArray(1, sum + c);
            m_comm->AllReduce(sumArray, LibUtilities::ReduceSum);
            if (m_comm->GetRank() == 0)
            {
                m_outFile << " " << sumArray[0];
            }
        }
    }

    if (m_comm->GetRank() == 0)
    {
        m_outFile << std::endl;
    }
}

void FilterIntegral::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(pFields, time);

    if (m_comm->GetRank() == 0)
    {
        m_outFile.close();
    }
}

bool FilterIntegral::v_IsTimeDependent()
{
    return true;
}
} // namespace SolverUtils
} // namespace Nektar
