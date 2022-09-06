///////////////////////////////////////////////////////////////////////////////
//
// File FilterIntegral.h
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

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERINTEGRAL_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERINTEGRAL_H

#include <SolverUtils/Filters/Filter.h>

namespace Nektar
{
namespace SolverUtils
{
class FilterIntegral : public Filter
{
public:
    friend class MemoryManager<FilterIntegral>;

    /// Creates an instance of this class
    static FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation,
        const std::map<std::string, std::string> &pParams)
    {
        FilterSharedPtr p = MemoryManager<FilterIntegral>::AllocateSharedPtr(
            pSession, pEquation, pParams);
        return p;
    }

    /// Name of the class
    static std::string className;

    /// Constructs the integral filter and parses filter options, opens file
    SOLVER_UTILS_EXPORT FilterIntegral(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation,
        const ParamMap &pParams);

    /// Default destructor
    SOLVER_UTILS_EXPORT virtual ~FilterIntegral() = default;

protected:
    /**
     * Initialises the integral filter and stores the composite expansions in
     * #m_compExpMap for all composites specified in the XML file
     *
     * @param pFields Field data
     * @param time Current time
     */
    virtual void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) final;

    /**
     * Performs the integration on the stored composite expansions and outputs
     * in to the output data file
     *
     * @param pFields Field data
     * @param time Current time
     */
    virtual void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) final;

    /**
     * Closes the output data file
     *
     * @param pFields Field data
     * @param time Current time
     */
    virtual void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) final;

    /// Returns true as filter depends on time
    virtual bool v_IsTimeDependent() final;

private:
    size_t m_index = 0;
    /// Frequency to write to output data file in timesteps
    size_t m_outputFrequency;
    /// Number of fields to perform integral on
    size_t m_numVariables;
    /// Out file
    std::ofstream m_outFile;
    /// Global communicator
    LibUtilities::CommSharedPtr m_comm;
    /// Vector of composite IDs as a single string
    std::vector<std::string> m_splitCompString;
    /// Vector of vector of composites IDs as integers
    std::vector<std::vector<unsigned int>> m_compVector;
    /// Mapping from geometry ID to expansion ID
    std::map<size_t, size_t> m_geomElmtIdToExpId;
    /// Map of composite ID to vector of expansions an face/edge local ID
    std::map<int, std::vector<std::pair<LocalRegions::ExpansionSharedPtr, int>>>
        m_compExpMap;
};
} // namespace SolverUtils
} // namespace Nektar

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERIntegral_H */
