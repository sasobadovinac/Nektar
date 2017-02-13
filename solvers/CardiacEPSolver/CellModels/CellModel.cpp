///////////////////////////////////////////////////////////////////////////////
//
// File: CellModel.cpp
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
// Description: Cell model base class.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VmathArray.hpp>

#include <CardiacEPSolver/CellModels/CellModel.h>

#include <StdRegions/StdNodalTriExp.h>
//#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/BasicUtils/Equation.h>

using namespace std;

namespace Nektar
{
CellModelFactory &GetCellModelFactory()
{
    static CellModelFactory instance;
    return instance;
}

/**
 * @class CellModel
 *
 * The CellModel class and derived classes implement a range of cell model
 * ODE systems. A cell model comprises a system of ion concentration
 * variables and zero or more gating variables. Gating variables are
 * time-integrated using the Rush-Larsen method and for each variable y,
 * the corresponding y_inf and tau_y value is computed by Update(). The tau
 * values are stored in separate storage to inarray/outarray, #m_gates_tau.
 */

/**
 * Cell model base class constructor.
 */
CellModel::CellModel(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const MultiRegions::ExpListSharedPtr &pField)
{
    m_session  = pSession;
    m_field    = pField;
    m_lastTime = 0.0;
    m_substeps = pSession->GetParameter("Substeps");
    m_nvar     = 0;
    m_useNodal = false;

    // Number of points in nodal space is the number of coefficients
    // in modified basis
    std::set<enum LibUtilities::ShapeType> s;
    for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
    {
        s.insert(m_field->GetExp(i)->DetShapeType());
    }

    // Use nodal projection if only triangles
    if (s.size() == 1 && (s.count(LibUtilities::eTriangle) == 1 ||
                          s.count(LibUtilities::eTetrahedron) == 1))
    {
        // This is disabled for now as it causes problems at high order.
        // m_useNodal = true;
    }

    // ---------------------------
    // Move to nodal points
    if (m_useNodal)
    {
        m_nq      = pField->GetNcoeffs();
        int order = m_field->GetExp(0)->GetBasis(0)->GetNumModes();

        // Set up a nodal tri
        LibUtilities::BasisKey B0(
            LibUtilities::eModified_A, order,
            LibUtilities::PointsKey(order,
                                    LibUtilities::eGaussLobattoLegendre));
        LibUtilities::BasisKey B1(
            LibUtilities::eModified_B, order,
            LibUtilities::PointsKey(order,
                                    LibUtilities::eGaussRadauMAlpha1Beta0));
        LibUtilities::BasisKey B2(
            LibUtilities::eModified_C, order,
            LibUtilities::PointsKey(order,
                                    LibUtilities::eGaussRadauMAlpha2Beta0));

        m_nodalTri =
            MemoryManager<StdRegions::StdNodalTriExp>::AllocateSharedPtr(
                B0, B1, LibUtilities::eNodalTriEvenlySpaced);
        m_nodalTet =
            MemoryManager<StdRegions::StdNodalTetExp>::AllocateSharedPtr(
                B0, B1, B2, LibUtilities::eNodalTetEvenlySpaced);
    }
    else
    {
        m_nq = pField->GetTotPoints();
    }

    // ----------------------------
    // Load parameter fields from XML/files
    // m_parameters.clear();

    // Scan through CellModel section looking to see if the parameter requested
    // has been given a value.
    TiXmlElement *vParameters =
        m_session->GetElement("Nektar/CellModel/Parameters");
    TiXmlElement *variable = vParameters->FirstChildElement();

    while (variable)
    {
        std::string conditionType = variable->Value();
        std::string variableStr   = variable->Attribute("VAR");

        VariableDefinition varDef;

        boost::to_upper(variableStr);

        WARNINGL0(m_parameters.count(variableStr) == 0,
                  "Parameter " + variableStr +
                      " has multiple entries. Previous entries overwritten.");

        if (conditionType == "E")
        {
            varDef.m_type = eInputTypeExpression;

            // Expression must have a VALUE.
            ASSERTL0(variable->Attribute("VALUE"),
                     "Attribute VALUE expected for variable '" + variableStr +
                         "'.");
            std::string functionStr = variable->Attribute("VALUE");

            ASSERTL0(!functionStr.empty(),
                     (std::string("Expression for var: ") + variableStr +
                      std::string(" must be specified."))
                         .c_str());

            varDef.m_value = functionStr;

            m_parameters[variableStr] = varDef;
        }
        else if (conditionType == "F")
        {
            varDef.m_type = eInputTypeFile;

            // File must have a VALUE.
            ASSERTL0(variable->Attribute("VALUE"),
                     "Attribute VALUE expected for variable '" + variableStr +
                         "'.");
            std::string filenameStr = variable->Attribute("VALUE");

            ASSERTL0(!filenameStr.empty(),
                     "A filename must be specified for the VALUE "
                     "attribute in variable '" +
                         variableStr + "'.");

            varDef.m_value = filenameStr;

            m_parameters[variableStr] = varDef;
        }

        else if (conditionType == "S")
        {

            varDef.m_type = eInputTypeScarMap;

            if (variable->Attribute("VALUE"))
            {
                varDef.m_value = variable->Attribute("VALUE");
            }
            else
            {
                varDef.m_value = "PositiveLinear";
            }

            if (variableStr == "SCARMAP")
            {
                varDef.m_var_params["MIN"] = "0.0";
                varDef.m_var_params["MAX"] = "1.0";
            }

            TiXmlAttribute *pAttrib = variable->FirstAttribute();
            while (pAttrib)
            {
                std::string pAttribName  = pAttrib->Name();
                std::string pAttribValue = pAttrib->Value();

                if (pAttribName != "VAR" && pAttribName != "VALUE")
                {
                    boost::to_upper(pAttribName);
                    varDef.m_var_params[pAttribName] = pAttribValue;
                }

                pAttrib = pAttrib->Next();
            }

            m_parameters[variableStr] = varDef;
        }

        variable = variable->NextSiblingElement();
    }
}

/**
 * Initialise the cell model. Allocate workspace and variable storage.
 */
void CellModel::Initialise()
{
    ASSERTL1(m_nvar > 0, "Cell model must have at least 1 variable.");

    m_cellSol = Array<OneD, Array<OneD, NekDouble>>(m_nvar);
    m_wsp     = Array<OneD, Array<OneD, NekDouble>>(m_nvar);
    for (unsigned int i = 0; i < m_nvar; ++i)
    {
        m_cellSol[i] = Array<OneD, NekDouble>(m_nq);
        m_wsp[i]     = Array<OneD, NekDouble>(m_nq);
    }
    m_gates_tau = Array<OneD, Array<OneD, NekDouble>>(m_gates.size());
    for (unsigned int i = 0; i < m_gates.size(); ++i)
    {
        m_gates_tau[i] = Array<OneD, NekDouble>(m_nq);
    }

    if (m_session->DefinesFunction("CellModelInitialConditions"))
    {
        LoadCellModel();
    }
    else
    {
        v_SetInitialConditions();
    }
}

/**
 * Integrates the cell model for one PDE time-step. Cell model is
 * sub-stepped.
 *
 * Ion concentrations and membrane potential are integrated using forward
 * Euler, while gating variables are integrated using the Rush-Larsen
 * scheme.
 */
void CellModel::TimeIntegrate(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int phys_offset = 0;
    int coef_offset = 0;
    int nvar        = inarray.size();
    Array<OneD, NekDouble> tmp;

    // ---------------------------
    // Check nodal temp array set up
    if (m_useNodal)
    {
        if (!m_nodalTmp.size())
        {
            m_nodalTmp = Array<OneD, Array<OneD, NekDouble>>(nvar);
            for (unsigned int k = 0; k < nvar; ++k)
            {
                m_nodalTmp[k] = Array<OneD, NekDouble>(m_nq);
            }
        }

        // Move to nodal points
        Array<OneD, NekDouble> tmpCoeffs(
            max(m_nodalTri->GetNcoeffs(), m_nodalTet->GetNcoeffs()));

        for (unsigned int k = 0; k < nvar; ++k)
        {
            for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
            {
                phys_offset = m_field->GetPhys_Offset(i);
                coef_offset = m_field->GetCoeff_Offset(i);
                if (m_field->GetExp(0)->DetShapeType() ==
                    LibUtilities::eTriangle)
                {
                    m_field->GetExp(0)->FwdTrans(inarray[k] + phys_offset,
                                                 tmpCoeffs);
                    m_nodalTri->ModalToNodal(tmpCoeffs,
                                             tmp = m_nodalTmp[k] + coef_offset);
                }
                else
                {
                    m_field->GetExp(0)->FwdTrans(inarray[k] + phys_offset,
                                                 tmpCoeffs);
                    m_nodalTet->ModalToNodal(tmpCoeffs,
                                             tmp = m_nodalTmp[k] + coef_offset);
                }
            }
        }
        // Copy new transmembrane potential into cell model
        Vmath::Vcopy(m_nq, m_nodalTmp[0], 1, m_cellSol[0], 1);
    }
    else
    {
        // Copy new transmembrane potential into cell model
        Vmath::Vcopy(m_nq, inarray[0], 1, m_cellSol[0], 1);
    }
    // -------------------------

    NekDouble delta_t = (time - m_lastTime) / m_substeps;

    // Perform substepping
    for (unsigned int i = 0; i < m_substeps - 1; ++i)
    {
        Update(m_cellSol, m_wsp, time);
        // Voltage
        Vmath::Svtvp(m_nq, delta_t, m_wsp[0], 1, m_cellSol[0], 1, m_cellSol[0],
                     1);
        // Ion concentrations
        for (unsigned int j = 0; j < m_concentrations.size(); ++j)
        {
            Vmath::Svtvp(m_nq, delta_t, m_wsp[m_concentrations[j]], 1,
                         m_cellSol[m_concentrations[j]], 1,
                         m_cellSol[m_concentrations[j]], 1);
        }
        // Gating variables: Rush-Larsen scheme
        for (unsigned int j = 0; j < m_gates.size(); ++j)
        {
            Vmath::Sdiv(m_nq, -delta_t, m_gates_tau[j], 1, m_gates_tau[j], 1);
            Vmath::Vexp(m_nq, m_gates_tau[j], 1, m_gates_tau[j], 1);
            Vmath::Vsub(m_nq, m_cellSol[m_gates[j]], 1, m_wsp[m_gates[j]], 1,
                        m_cellSol[m_gates[j]], 1);
            Vmath::Vvtvp(m_nq, m_cellSol[m_gates[j]], 1, m_gates_tau[j], 1,
                         m_wsp[m_gates[j]], 1, m_cellSol[m_gates[j]], 1);
        }
    }

    // Perform final cell model step
    Update(m_cellSol, m_wsp, time);

    // Output dV/dt from last step but integrate remaining cell model vars
    // Transform cell model I_total from nodal to modal space
    if (m_useNodal)
    {
        Array<OneD, NekDouble> tmpCoeffs(
            max(m_nodalTri->GetNcoeffs(), m_nodalTet->GetNcoeffs()));

        for (unsigned int k = 0; k < nvar; ++k)
        {
            for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
            {
                int phys_offset = m_field->GetPhys_Offset(i);
                int coef_offset = m_field->GetCoeff_Offset(i);
                if (m_field->GetExp(0)->DetShapeType() ==
                    LibUtilities::eTriangle)
                {
                    m_nodalTri->NodalToModal(m_wsp[k] + coef_offset, tmpCoeffs);
                    m_field->GetExp(0)->BwdTrans(tmpCoeffs, tmp = outarray[k] +
                                                                  phys_offset);
                }
                else
                {
                    m_nodalTet->NodalToModal(m_wsp[k] + coef_offset, tmpCoeffs);
                    m_field->GetExp(0)->BwdTrans(tmpCoeffs, tmp = outarray[k] +
                                                                  phys_offset);
                }
            }
        }
    }
    else
    {
        Vmath::Vcopy(m_nq, m_wsp[0], 1, outarray[0], 1);
    }

    // Ion concentrations
    for (unsigned int j = 0; j < m_concentrations.size(); ++j)
    {
        Vmath::Svtvp(m_nq, delta_t, m_wsp[m_concentrations[j]], 1,
                     m_cellSol[m_concentrations[j]], 1,
                     m_cellSol[m_concentrations[j]], 1);
    }

    // Gating variables: Rush-Larsen scheme
    for (unsigned int j = 0; j < m_gates.size(); ++j)
    {
        Vmath::Sdiv(m_nq, -delta_t, m_gates_tau[j], 1, m_gates_tau[j], 1);
        Vmath::Vexp(m_nq, m_gates_tau[j], 1, m_gates_tau[j], 1);
        Vmath::Vsub(m_nq, m_cellSol[m_gates[j]], 1, m_wsp[m_gates[j]], 1,
                    m_cellSol[m_gates[j]], 1);
        Vmath::Vvtvp(m_nq, m_cellSol[m_gates[j]], 1, m_gates_tau[j], 1,
                     m_wsp[m_gates[j]], 1, m_cellSol[m_gates[j]], 1);
    }

    m_lastTime = time;
}

Array<OneD, NekDouble> CellModel::GetCellSolutionCoeffs(unsigned int idx)
{
    ASSERTL0(idx < m_nvar, "Index out of range for cell model.");

    Array<OneD, NekDouble> outarray(m_field->GetNcoeffs());
    Array<OneD, NekDouble> tmp;

    if (m_useNodal)
    {
        for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
        {
            int coef_offset = m_field->GetCoeff_Offset(i);
            if (m_field->GetExp(0)->DetShapeType() == LibUtilities::eTriangle)
            {
                m_nodalTri->NodalToModal(m_cellSol[idx] + coef_offset,
                                         tmp = outarray + coef_offset);
            }
            else
            {
                m_nodalTet->NodalToModal(m_cellSol[idx] + coef_offset,
                                         tmp = outarray + coef_offset);
            }
        }
    }
    else
    {
        m_field->FwdTransLocalElmt(m_cellSol[idx], outarray);
    }

    return outarray;
}

Array<OneD, NekDouble> CellModel::GetCellSolution(unsigned int idx)
{
    return m_cellSol[idx];
}

void CellModel::LoadCellModel()
{
    const bool root           = (m_session->GetComm()->GetRank() == 0);
    const std::string fncName = "CellModelInitialConditions";
    const int nvar            = m_cellSol[0].size();
    std::string varName;
    Array<OneD, NekDouble> coeffs(m_field->GetNcoeffs());
    Array<OneD, NekDouble> tmp;
    int j = 0;

    SpatialDomains::MeshGraphSharedPtr vGraph = m_field->GetGraph();

    if (root)
    {
        cout << "Cell model initial conditions: " << endl;
    }

    // First determine all the files we need to load
    std::set<std::string> filelist;
    for (j = 1; j < nvar; ++j)
    {
        // Get the name of the jth variable
        varName = GetCellVarName(j);

        if (m_session->GetFunctionType(fncName, varName) ==
            LibUtilities::eFunctionTypeFile)
        {
            filelist.insert(m_session->GetFunctionFilename(fncName, varName));
        }
    }

    // Read files
    typedef std::vector<LibUtilities::FieldDefinitionsSharedPtr> FDef;
    typedef std::vector<std::vector<NekDouble>> FData;
    std::map<std::string, FDef> FieldDef;
    std::map<std::string, FData> FieldData;
    LibUtilities::FieldMetaDataMap fieldMetaDataMap;

    for (auto &setIt : filelist)
    {
        if (root)
        {
            cout << "  - Reading file: " << setIt << endl;
        }
        FieldDef[setIt]  = FDef(0);
        FieldData[setIt] = FData(0);
        LibUtilities::FieldIOSharedPtr fld =
            LibUtilities::FieldIO::CreateForFile(m_session, setIt);
        fld->Import(setIt, FieldDef[setIt], FieldData[setIt], fieldMetaDataMap);
    }

    // Get time of checkpoint from file if available
    auto iter = fieldMetaDataMap.find("Time");
    if (iter != fieldMetaDataMap.end())
    {
        m_lastTime = boost::lexical_cast<NekDouble>(iter->second);
    }

    // Load each cell model variable
    // j=0 and j=1 are for transmembrane or intra/extra-cellular volt.
    Vmath::Zero(m_nq, m_cellSol[0], 1);
    for (j = 1; j < m_cellSol.size(); ++j)
    {
        // Get the name of the jth variable
        varName = GetCellVarName(j);

        // Check if this variable is defined in a file or analytically
        if (m_session->GetFunctionType(fncName, varName) ==
            LibUtilities::eFunctionTypeFile)
        {
            const std::string file =
                m_session->GetFunctionFilename(fncName, varName);

            if (root)
            {
                cout << "  - Field " << varName << ": from file " << file
                     << endl;
            }

            // Extract the data into the modal coefficients
            for (int i = 0; i < FieldDef[file].size(); ++i)
            {
                m_field->ExtractDataToCoeffs(
                    FieldDef[file][i], FieldData[file][i], varName, coeffs);
            }

            // If using nodal cell model then we do a modal->nodal transform
            // otherwise we do a backward transform onto physical points.
            if (m_useNodal)
            {
                for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
                {
                    int coef_offset = m_field->GetCoeff_Offset(i);
                    if (m_field->GetExp(0)->DetShapeType() ==
                        LibUtilities::eTriangle)
                    {
                        m_nodalTri->ModalToNodal(coeffs + coef_offset,
                                                 tmp = m_cellSol[j] +
                                                       coef_offset);
                    }
                    else
                    {
                        m_nodalTet->ModalToNodal(coeffs + coef_offset,
                                                 tmp = m_cellSol[j] +
                                                       coef_offset);
                    }
                }
            }
            else
            {
                m_field->BwdTrans(coeffs, m_cellSol[j]);
            }
        }
        else if (m_session->GetFunctionType(fncName, varName) ==
                 LibUtilities::eFunctionTypeExpression)
        {
            LibUtilities::EquationSharedPtr equ =
                m_session->GetFunction(fncName, varName);

            if (root)
            {
                cout << "  - Field " << varName << ": " << equ->GetExpression()
                     << endl;
            }

            const unsigned int nphys = m_field->GetNpoints();
            Array<OneD, NekDouble> x0(nphys);
            Array<OneD, NekDouble> x1(nphys);
            Array<OneD, NekDouble> x2(nphys);
            m_field->GetCoords(x0, x1, x2);

            if (m_useNodal)
            {
                Array<OneD, NekDouble> phys(nphys);
                Array<OneD, NekDouble> tmpCoeffs(
                    max(m_nodalTri->GetNcoeffs(), m_nodalTet->GetNcoeffs()));

                equ->Evaluate(x0, x1, x2, phys);
                for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
                {
                    int phys_offset = m_field->GetPhys_Offset(i);
                    int coef_offset = m_field->GetCoeff_Offset(i);
                    if (m_field->GetExp(0)->DetShapeType() ==
                        LibUtilities::eTriangle)
                    {
                        m_field->GetExp(0)->FwdTrans(phys + phys_offset,
                                                     tmpCoeffs);
                        m_nodalTri->ModalToNodal(tmpCoeffs, tmp = m_cellSol[j] +
                                                                  coef_offset);
                    }
                    else
                    {
                        m_field->GetExp(0)->FwdTrans(phys + phys_offset,
                                                     tmpCoeffs);
                        m_nodalTet->ModalToNodal(tmpCoeffs, tmp = m_cellSol[j] +
                                                                  coef_offset);
                    }
                }
            }
            else
            {
                equ->Evaluate(x0, x1, x2, m_cellSol[j]);
            }
        }
    }
}

Array<OneD, NekDouble> CellModel::LoadCellParam(std::string parameter,
                                                NekDouble defaultValue)
{
    const bool root = (m_session->GetComm()->GetRank() == 0);
    boost::to_upper(parameter);

    const unsigned int nphys = m_field->GetNpoints();
    // ASSERTL0(DataOutput.size() == nphys, "Input cell parameter field has
    // incorrect size.")
    Array<OneD, NekDouble> DataOutput(nphys, defaultValue);

    if (m_parameters.count(parameter) > 0)
    {

        VariableDefinition vParameter = m_parameters[parameter];

        if (vParameter.m_type == eInputTypeExpression)
        {
            Array<OneD, NekDouble> x0(nphys);
            Array<OneD, NekDouble> x1(nphys);
            Array<OneD, NekDouble> x2(nphys);
            m_field->GetCoords(x0, x1, x2);

            LibUtilities::Equation paramExpr(m_session->GetInterpreter(),
                                             vParameter.m_value);
            paramExpr.Evaluate(x0, x1, x2, DataOutput);
        }
        else if (vParameter.m_type == eInputTypeFile)
        {

            std::string vValue = vParameter.m_value;

            ReadFromFile(vValue, vValue, DataOutput);
        }
        else if (vParameter.m_type == eInputTypeScarMap)
        {

            if (m_scarmap.size() < nphys)
            {
                std::string vScarMapFilename = m_parameters["SCARMAP"].m_value;
                m_scarmap                    = Array<OneD, NekDouble>(nphys);

                ReadFromFile(vScarMapFilename, "intensity", m_scarmap);
            }

            Array<OneD, NekDouble> vTemp = m_scarmap;
            MappingType vMapType;

            if (vParameter.m_value == "PositiveLinear")
            {
                vMapType = ePosLinThresh;
            }
            else if (vParameter.m_value == "NegativeLinear")
            {
                vMapType = eNegLinThresh;
            }
            else if (vParameter.m_value == "BinaryThreshold")
            {
                vMapType = eBinaryThresh;
            }
            else
            {
                vMapType = eError;
            }

            switch (vMapType)
            {

                case ePosLinThresh:
                {
                    NekDouble f_min, f_max;

                    if (ParamExists(parameter, "CUTLO"))
                    {
                        f_min = GetParamDouble(parameter, "CUTLO");
                    }
                    else
                    {
                        ASSERTL0(m_parameters["SCARMAP"].m_var_params.count(
                                     "CUTLO") > 0,
                                 "Attribute CUTLO required for variable '" +
                                     parameter + "'.");
                        f_min = GetParamDouble("SCARMAP", "CUTLO");
                    }

                    if (ParamExists(parameter, "CUTHI"))
                    {
                        f_max = GetParamDouble(parameter, "CUTHI");
                    }
                    else
                    {
                        ASSERTL0(m_parameters["SCARMAP"].m_var_params.count(
                                     "CUTHI") > 0,
                                 "Attribute CUTHI required for variable '" +
                                     parameter + "'.");
                        f_max = GetParamDouble("SCARMAP", "CUTHI");
                    }
                    ASSERTL0(m_parameters[parameter].m_var_params.count("LO") >
                                 0,
                             "Attribute LO required for variable '" +
                                 parameter + "'.");
                    ASSERTL0(m_parameters[parameter].m_var_params.count("HI") >
                                 0,
                             "Attribute HI required for variable '" +
                                 parameter + "'.");
                    NekDouble scar_min = GetParamDouble(parameter, "LO");
                    NekDouble scar_max = GetParamDouble(parameter, "HI");

                    // Threshold based on d_min, d_max
                    for (int j = 0; j < nphys; ++j)
                    {
                        vTemp[j] = (vTemp[j] < f_min ? f_min : vTemp[j]);
                        vTemp[j] = (vTemp[j] > f_max ? f_max : vTemp[j]);
                    }
                    // Rescale to s \in [0,1] (0 maps to d_min, 1 maps to d_max)
                    Vmath::Sadd(nphys, -f_min, vTemp, 1, vTemp, 1);
                    Vmath::Smul(nphys, 1.0 / (f_max - f_min), vTemp, 1, vTemp,
                                1);
                    Vmath::Smul(nphys, scar_max - scar_min, vTemp, 1, vTemp, 1);
                    Vmath::Sadd(nphys, scar_min, vTemp, 1, vTemp, 1);

                    DataOutput = vTemp;

                    break;
                }

                case eNegLinThresh:
                {
                    NekDouble f_min, f_max;

                    if (ParamExists(parameter, "CUTLO"))
                    {
                        f_min = GetParamDouble(parameter, "CUTLO");
                    }
                    else
                    {
                        ASSERTL0(m_parameters["SCARMAP"].m_var_params.count(
                                     "CUTLO") > 0,
                                 "Attribute CUTLO required for variable '" +
                                     parameter + "'.");
                        f_min = GetParamDouble("SCARMAP", "CUTLO");
                    }

                    if (ParamExists(parameter, "CUTHI"))
                    {
                        f_max = GetParamDouble(parameter, "CUTHI");
                    }
                    else
                    {
                        ASSERTL0(m_parameters["SCARMAP"].m_var_params.count(
                                     "CUTHI") > 0,
                                 "Attribute CUTHI required for variable '" +
                                     parameter + "'.");
                        f_max = GetParamDouble("SCARMAP", "CUTHI");
                    }
                    ASSERTL0(m_parameters[parameter].m_var_params.count("LO") >
                                 0,
                             "Attribute LO required for variable '" +
                                 parameter + "'.");
                    ASSERTL0(m_parameters[parameter].m_var_params.count("HI") >
                                 0,
                             "Attribute HI required for variable '" +
                                 parameter + "'.");
                    NekDouble scar_min = GetParamDouble(parameter, "LO");
                    NekDouble scar_max = GetParamDouble(parameter, "HI");

                    // Threshold based on d_min, d_max
                    for (int j = 0; j < nphys; ++j)
                    {
                        vTemp[j] = (vTemp[j] < f_min ? f_min : vTemp[j]);
                        vTemp[j] = (vTemp[j] > f_max ? f_max : vTemp[j]);
                    }
                    // Rescale to s \in [0,1] (0 maps to d_max, 1 maps to d_min)
                    Vmath::Sadd(nphys, -f_min, vTemp, 1, vTemp, 1);
                    Vmath::Smul(nphys, -1.0 / (f_max - f_min), vTemp, 1, vTemp,
                                1);
                    Vmath::Sadd(nphys, 1.0, vTemp, 1, vTemp, 1);
                    Vmath::Smul(nphys, scar_max - scar_min, vTemp, 1, vTemp, 1);
                    Vmath::Sadd(nphys, scar_min, vTemp, 1, vTemp, 1);

                    DataOutput = vTemp;

                    break;
                }

                case eBinaryThresh:
                {
                    ASSERTL0(vParameter.m_var_params.count("HI") > 0,
                             "Attribute HI required for variable '" +
                                 parameter + "'.");
                    ASSERTL0(vParameter.m_var_params.count("LO") > 0,
                             "Attribute LO required for variable '" +
                                 parameter + "'.");
                    ASSERTL0(vParameter.m_var_params.count("CUTOFF") > 0,
                             "Attribute CUTOFF required for variable '" +
                                 parameter + "'.");

                    NekDouble f_low    = GetParamDouble(parameter, "LO");
                    NekDouble f_high   = GetParamDouble(parameter, "HI");
                    NekDouble f_cutoff = GetParamDouble(parameter, "CUTOFF");

                    for (int j = 0; j < nphys; ++j)
                    {
                        DataOutput[j] = (vTemp[j] < f_cutoff ? f_low : f_high);
                    }

                    break;
                }

                case eError:
                {
                    WARNINGL0(false, "Incorrect VALUE for variable '" +
                                         parameter + "'. Using default value.")
                }
            }
        }
    }
    else
    {
        if (root)
        {
            cout << "Default value used for " << parameter << "." << endl;
        }
    }

    return DataOutput;
}

std::string CellModel::GetParamString(std::string var, std::string param)
{
    ASSERTL0(m_parameters.count(var) > 0 &&
                 m_parameters[var].m_var_params.count(param) > 0,
             "Failed to find variable '" + var + "' for scar parameter: '" +
                 param + "'");

    std::string value = m_parameters[var].m_var_params[param];

    return value;
}

NekDouble CellModel::GetParamDouble(std::string var, std::string param)
{
    std::string strExpression = GetParamString(var, param);

    NekDouble value = 0.0;
    try
    {
        LibUtilities::Equation expession(m_session->GetInterpreter(),
                                         strExpression);
        value = expession.Evaluate();
    }
    catch (const std::runtime_error &)
    {
        ASSERTL0(false, "Error evaluating parameter expression"
                        " '" +
                            strExpression + "' for parameter: '" + param + "'");
    }

    return value;
}

bool CellModel::ParamExists(std::string var, std::string parameter)
{
    ASSERTL0(m_parameters.count(var) > 0,
             "Variable does not exist: '" + var + "'");

    if (m_parameters[var].m_var_params.count(parameter) > 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void CellModel::ReadFromFile(std::string filename, std::string param,
                             Array<OneD, NekDouble> outArray)
{

    LibUtilities::FieldIOSharedPtr fileObj =
        LibUtilities::FieldIO::CreateForFile(m_session, filename);

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FDef(0);
    std::vector<std::vector<NekDouble>> FData(0);
    LibUtilities::FieldMetaDataMap fieldMetaDataMap;

    fileObj->Import(filename, FDef, FData, fieldMetaDataMap);

    const unsigned int ncoeffs = m_field->GetNcoeffs();
    Array<OneD, NekDouble> CoeffOutput(ncoeffs);

    m_field->ExtractDataToCoeffs(FDef[0], FData[0], param, CoeffOutput);

    if (outArray.size() != m_field->GetNpoints())
    {
        outArray = Array<OneD, NekDouble>(m_field->GetNpoints());
    }
    m_field->BwdTrans(CoeffOutput, outArray);
}

} // namespace Nektar
