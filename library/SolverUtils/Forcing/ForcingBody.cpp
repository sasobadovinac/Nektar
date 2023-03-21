///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingBody.cpp
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
// Description: Body forcing
//
///////////////////////////////////////////////////////////////////////////////
#include <boost/algorithm/string.hpp>

#include <MultiRegions/ExpList.h>
#include <SolverUtils/Forcing/ForcingBody.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

std::string ForcingBody::classNameBody =
    GetForcingFactory().RegisterCreatorFunction("Body", ForcingBody::create,
                                                "Body Forcing");
std::string ForcingBody::classNameField =
    GetForcingFactory().RegisterCreatorFunction("Field", ForcingBody::create,
                                                "Field Forcing");

ForcingBody::ForcingBody(const LibUtilities::SessionReaderSharedPtr &pSession,
                         const std::weak_ptr<EquationSystem> &pEquation)
    : Forcing(pSession, pEquation), m_hasTimeFcnScaling(false)
{
}

void ForcingBody::v_InitObject(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const unsigned int &pNumForcingFields, const TiXmlElement *pForce)
{
    m_NumVariable = pNumForcingFields;

    const TiXmlElement *funcNameElmt = pForce->FirstChildElement("BODYFORCE");
    if (!funcNameElmt)
    {
        funcNameElmt = pForce->FirstChildElement("FIELDFORCE");

        ASSERTL0(funcNameElmt,
                 "Requires BODYFORCE or FIELDFORCE tag "
                 "specifying function name which prescribes body force.");
    }

    m_funcName = funcNameElmt->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName),
             "Function '" + m_funcName + "' not defined.");

    m_homogeneous = pFields[0]->GetExpType() == MultiRegions::e3DH1D ||
                    pFields[0]->GetExpType() == MultiRegions::e3DH2D;

    // Time function is optional
    funcNameElmt = pForce->FirstChildElement("BODYFORCETIMEFCN");
    if (!funcNameElmt)
    {
        funcNameElmt = pForce->FirstChildElement("FIELDFORCETIMEFCN");
    }

    // Load time function if specified
    if (funcNameElmt)
    {
        std::string funcNameTime = funcNameElmt->GetText();

        ASSERTL0(!funcNameTime.empty(),
                 "Expression must be given in BODYFORCETIMEFCN or "
                 "FIELDFORCETIMEFCN.");

        m_session->SubstituteExpressions(funcNameTime);
        m_timeFcnEqn = MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(
            m_session->GetInterpreter(), funcNameTime);

        m_hasTimeFcnScaling = true;
    }

    Array<OneD, Array<OneD, NekDouble>> tmp(m_NumVariable);
    for (int i = 0; i < m_NumVariable; ++i)
    {
        tmp[i] = pFields[i]->GetPhys();
    }

    std::map<std::string, int> varIndex;
    for (int i = 0; i < m_NumVariable; ++i)
    {
        varIndex[m_session->GetVariable(i)] = i;
        if (m_session->DefinesFunction(m_funcName, m_session->GetVariable(i)))
        {
            m_eqnEvars[i] = vector<int>();
        }
    }
    m_hasEvars = false;
    for (auto &it : m_eqnEvars)
    {
        string varStr = m_session->GetVariable(it.first);
        if (LibUtilities::eFunctionTypeExpression ==
            m_session->GetFunctionType(m_funcName, varStr))
        {
            LibUtilities::EquationSharedPtr eqn =
                m_session->GetFunction(m_funcName, varStr);
            string vlist = eqn->GetVlist();
            if (!boost::iequals(vlist, "x y z t"))
            {
                // Coupled forcing
                m_hasEvars = true;
                std::vector<std::string> vars;
                boost::split(vars, vlist, boost::is_any_of(", "));
                for (size_t j = 4; j < vars.size(); ++j)
                {
                    if (vars[j].size() != 0 &&
                        varIndex.find(vars[j]) == varIndex.end())
                    {
                        NEKERROR(ErrorUtil::efatal,
                                 "Variable '" + vars[j] +
                                     "' cannot be used as equation variable in "
                                     "the body force '" +
                                     m_funcName + "'.");
                    }
                    if (vars[j].size() &&
                        varIndex.find(vars[j]) != varIndex.end())
                    {
                        it.second.push_back(varIndex[vars[j]]);
                    }
                }
            }
        }
    }

    m_Forcing = Array<OneD, Array<OneD, NekDouble>>(m_NumVariable);
    for (const auto &it : m_eqnEvars)
    {
        m_Forcing[it.first] =
            Array<OneD, NekDouble>(pFields[0]->GetTotPoints(), 0.0);
    }

    Update(pFields, tmp, 0.0);
}

void ForcingBody::Update(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray, const NekDouble &time)
{
    std::vector<Array<OneD, const NekDouble>> fielddata;
    int nq = pFields[0]->GetNpoints();
    if (m_hasEvars)
    {
        Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq), t(nq, time);
        pFields[0]->GetCoords(xc, yc, zc);
        fielddata.push_back(xc);
        fielddata.push_back(yc);
        fielddata.push_back(zc);
        fielddata.push_back(t);
    }

    for (const auto &it : m_eqnEvars)
    {
        int i = it.first;
        if (it.second.size() > 0)
        {
            // Coupled forcing, reset fielddata for each equation
            fielddata.resize(4);
            for (int j : it.second)
            {
                if (m_homogeneous && pFields[i]->GetWaveSpace())
                {
                    Array<OneD, NekDouble> tmp(nq);
                    pFields[i]->HomogeneousBwdTrans(nq, inarray[j], tmp);
                    fielddata.push_back(tmp);
                }
                else
                {
                    fielddata.push_back(inarray[j]);
                }
            }
            m_session->GetFunction(m_funcName, m_session->GetVariable(i))
                ->Evaluate(fielddata, m_Forcing[i]);
        }
        else
        {
            GetFunction(pFields, m_session, m_funcName, true)
                ->Evaluate(m_session->GetVariable(i), m_Forcing[i], time);
        }

        // If homogeneous expansion is used, transform the forcing term to
        // be in the Fourier space
        if (m_homogeneous)
        {
            pFields[i]->HomogeneousFwdTrans(pFields[i]->GetTotPoints(),
                                            m_Forcing[i], m_Forcing[i]);
        }
    }
}

void ForcingBody::v_Apply(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble &time)
{
    if (m_hasTimeFcnScaling)
    {
        Array<OneD, NekDouble> TimeFcn(1);
        EvaluateTimeFunction(time, m_timeFcnEqn, TimeFcn);

        for (const auto &it : m_eqnEvars)
        {
            int i = it.first;
            Vmath::Svtvp(outarray[i].size(), TimeFcn[0], m_Forcing[i], 1,
                         outarray[i], 1, outarray[i], 1);
        }
    }
    else
    {
        Update(fields, inarray, time);

        for (const auto &it : m_eqnEvars)
        {
            int i = it.first;
            Vmath::Vadd(outarray[i].size(), outarray[i], 1, m_Forcing[i], 1,
                        outarray[i], 1);
        }
    }
}

void ForcingBody::v_ApplyCoeff(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble &time)
{
    int ncoeff = outarray[m_NumVariable - 1].size();
    Array<OneD, NekDouble> tmp(ncoeff, 0.0);

    if (m_hasTimeFcnScaling)
    {
        Array<OneD, NekDouble> TimeFcn(1);
        EvaluateTimeFunction(time, m_timeFcnEqn, TimeFcn);

        for (const auto &it : m_eqnEvars)
        {
            int i = it.first;
            fields[i]->FwdTrans(m_Forcing[i], tmp);
            Vmath::Svtvp(ncoeff, TimeFcn[0], tmp, 1, outarray[i], 1,
                         outarray[i], 1);
        }
    }
    else
    {
        Update(fields, inarray, time);

        for (const auto &it : m_eqnEvars)
        {
            int i = it.first;
            fields[i]->FwdTrans(m_Forcing[i], tmp);
            Vmath::Vadd(ncoeff, outarray[i], 1, tmp, 1, outarray[i], 1);
        }
    }
}

} // namespace SolverUtils
} // namespace Nektar
