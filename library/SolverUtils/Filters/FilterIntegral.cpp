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

#include <boost/core/ignore_unused.hpp>
#include <SolverUtils/Filters/FilterIntegral.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>

#include <boost/algorithm/string.hpp>

namespace Nektar
{
namespace SolverUtils
{
std::string FilterIntegral::className =
    GetFilterFactory().RegisterCreatorFunction("Integral", FilterIntegral::create);

FilterIntegral::FilterIntegral(const LibUtilities::SessionReaderSharedPtr &pSession,
                         const std::weak_ptr<EquationSystem> &pEquation,
                         const ParamMap &pParams)
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
    boost::split(m_splitCompString, it->second,  boost::is_any_of(","));
    for (auto &comp : m_splitCompString)
    {
        size_t first = comp.find_first_of('[') + 1;
        size_t last = comp.find_last_of(']') - 1;
        auto tmpString = comp.substr(first, last - first + 1);

        std::vector<unsigned int> tmpVec;
        bool parseGood = ParseUtils::GenerateSeqVector(tmpString, tmpVec);

        ASSERTL0(parseGood && !tmpVec.empty(),
                 "Unable to read composite regions index range for "
                 "FilterIntegral: " + comp);

        m_compVector.emplace_back(tmpVec);
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

        m_outFile << "Time";

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
    auto expList =  pFields[0]->GetExp();
    for (size_t i = 0; i < expList->size(); ++i)
    {
        auto exp = (*expList)[i];
        m_geomElmtIdToExpId[exp->GetElmtId()] = i;
    }

    std::vector<std::pair<size_t, std::vector<size_t>>>
        pairDimGeomIds(m_compVector.size());


    for (int i = 0; i < m_compVector.size(); ++ i)
    {
        // Get comp list dimension from first composite & element
        auto composites = pFields[0]->GetGraph()->GetComposites();
        size_t dim = composites[m_compVector[0][0]]->m_geomVec[0]->GetShapeDim();


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
            ASSERTL0(dim == compGeom[0]->GetShapeDim(),
                "Differing geometry dimensions specified in FilterIntegral '" +
                    m_splitCompString[i] + "'.");
        }

        if (dim == pFields[0]->GetShapeDimension())
        {
            LocalRegions::ExpansionVector tmpCompExp(compGeomIds.size());

            for (size_t j = 0; j < compGeomIds.size(); ++j)
            {
                tmpCompExp[j] = (*expList)[m_geomElmtIdToExpId[compGeomIds[j]]];
            }

            m_compExpVector.emplace_back(tmpCompExp);
        }
        else if (dim < pFields[0]->GetShapeDimension())
        {
            ASSERTL0(false, "Finding the integral on a composite of dimension"
                            "smaller than the domain is not currently supported.")
        }
    }

    v_Update(pFields, time);
}

void FilterIntegral::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(pFields);

    if (m_index++ % m_outputFrequency > 0)
    {
        return;
    }

    if (m_comm->GetRank() == 0)
    {
        m_outFile << time;
    }

    // Lock equation system pointer
    auto equationSys = m_equ.lock();
    ASSERTL0(equationSys, "Weak pointer expired");

    for (auto &compExp : m_compExpVector)
    {
        for (size_t i = 0; i < m_numVariables; ++i)
        {
            // Evaluate integral using improved Kahan–Babuška summation
            // algorithm to reduce numerical error from adding floating point
            auto phys = pFields[i]->GetPhys();

            NekDouble sum = 0.0;
            NekDouble c   = 0.0;
            for (auto &exp : compExp)
            {
                size_t offset = pFields[0]->GetPhys_Offset(
                    m_geomElmtIdToExpId[exp->GetElmtId()]);

                NekDouble input = exp->Integral(phys + offset);
                NekDouble t = sum + input;
                c += fabs(sum) >= fabs(input) ? (sum - t) + input
                                              : (input - t) + sum;
                sum = t;
            }

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
