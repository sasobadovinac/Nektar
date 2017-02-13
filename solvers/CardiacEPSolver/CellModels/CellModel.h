///////////////////////////////////////////////////////////////////////////////
//
// File: CellModel.h
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

#ifndef NEKTAR_SOLVERS_ADRSOLVER_CELLMODELS_CELLMODEL
#define NEKTAR_SOLVERS_ADRSOLVER_CELLMODELS_CELLMODEL

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
//#include <SpatialDomains/SpatialData.h>
#include <LibUtilities/BasicUtils/Equation.h>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/Core/Misc.h>
#include <SolverUtils/UnsteadySystem.h>
#include <StdRegions/StdNodalTetExp.h>
#include <StdRegions/StdNodalTriExp.h>

namespace Nektar
{
// Forward declaration
class CellModel;

typedef std::vector<std::pair<std::string, std::string>> SummaryList;

/// A shared pointer to an EquationSystem object
typedef std::shared_ptr<CellModel> CellModelSharedPtr;
/// Datatype of the NekFactory used to instantiate classes derived from
/// the EquationSystem class.
typedef LibUtilities::NekFactory<std::string, CellModel,
                                 const LibUtilities::SessionReaderSharedPtr &,
                                 const MultiRegions::ExpListSharedPtr &>
    CellModelFactory;
CellModelFactory &GetCellModelFactory();

/// Cell model base class.
class CellModel
{
public:
    CellModel(const LibUtilities::SessionReaderSharedPtr &pSession,
              const MultiRegions::ExpListSharedPtr &pField);

    virtual ~CellModel()
    {
    }

    /// Initialise the cell model storage and set initial conditions
    void Initialise();

    /// Time integrate the cell model by one PDE timestep
    void TimeIntegrate(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                       Array<OneD, Array<OneD, NekDouble>> &outarray,
                       const NekDouble time);

    /// Compute the derivatives of cell model variables
    void Update(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                Array<OneD, Array<OneD, NekDouble>> &outarray,
                const NekDouble time)
    {
        v_Update(inarray, outarray, time);
    }

    /// Print a summary of the cell model
    void GenerateSummary(SummaryList &s)
    {
        v_GenerateSummary(s);
    }

    unsigned int GetNumCellVariables()
    {
        return m_nvar;
    }

    std::string GetCellVarName(unsigned int idx)
    {
        return v_GetCellVarName(idx);
    }

    Array<OneD, NekDouble> GetCellSolutionCoeffs(unsigned int idx);

    Array<OneD, NekDouble> GetCellSolution(unsigned int idx);

    /// Evaluate input expressions
    Array<OneD, NekDouble> LoadCellParam(std::string parameter,
                                         NekDouble defaultValue);

    void ReadFromFile(std::string filename, std::string param,
                      Array<OneD, NekDouble> outArray);

    bool ParamExists(std::string var, std::string parameter);
    std::string GetParamString(std::string strExpression,
                               std::string parameter);
    NekDouble GetParamDouble(std::string strExpression, std::string parameter);

protected:
    /// Session
    LibUtilities::SessionReaderSharedPtr m_session;
    /// Transmembrane potential field from PDE system
    MultiRegions::ExpListSharedPtr m_field;
    /// Number of physical points.
    int m_nq;
    /// Number of variables in cell model (inc. transmembrane voltage)
    int m_nvar;
    /// Timestep for pde model
    NekDouble m_lastTime;
    /// Number of substeps to take
    int m_substeps;

    /// Hold scar map intensity array, if used
    Array<OneD, NekDouble> m_scarmap;

    /// Type of input of given parameter
    enum InputType
    {
        eInputTypeExpression,
        eInputTypeFile,
        eInputTypeScarMap
    };

    /// Type of value-mapping to intensity field (scar map)
    enum MappingType
    {
        ePosLinThresh,
        eNegLinThresh,
        eBinaryThresh,
        eError,
    };

    struct VariableDefinition
    {
        enum InputType m_type;
        std::string m_value;
        /// Hold parameter values for scar map variables
        std::map<std::string, std::string> m_var_params;
    };

    /// Hold parameter functions/filenames
    std::map<std::string, VariableDefinition> m_parameters;

    /// Cell model solution variables
    Array<OneD, Array<OneD, NekDouble>> m_cellSol;
    /// Cell model integration workspace
    Array<OneD, Array<OneD, NekDouble>> m_wsp;

    /// Flag indicating whether nodal projection in use
    bool m_useNodal;
    /// StdNodalTri for cell model calculations
    StdRegions::StdNodalTriExpSharedPtr m_nodalTri;
    StdRegions::StdNodalTetExpSharedPtr m_nodalTet;
    /// Temporary array for nodal projection
    Array<OneD, Array<OneD, NekDouble>> m_nodalTmp;

    /// Indices of cell model variables which are concentrations
    std::vector<int> m_concentrations;
    /// Indices of cell model variables which are gates
    std::vector<int> m_gates;
    /// Storage for gate tau values
    Array<OneD, Array<OneD, NekDouble>> m_gates_tau;

    virtual void v_Update(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const NekDouble time) = 0;

    virtual void v_GenerateSummary(SummaryList &s) = 0;

    virtual std::string v_GetCellVarName(unsigned int idx)
    {
        return "Var" + boost::lexical_cast<std::string>(idx);
    }

    virtual void v_SetInitialConditions() = 0;

    void LoadCellModel();
};

} // namespace Nektar
#endif /* CELLMODEL_H_ */
