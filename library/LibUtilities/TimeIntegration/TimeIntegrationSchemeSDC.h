///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationSchemeSDC.h
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
// Description: Header file of time integration scheme SDC base class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_SCHEME_SDC
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_SCHEME_SDC

#define LUE LIB_UTILITIES_EXPORT

#include <string>

#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

namespace Nektar
{
namespace LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
/// Class for spectral deferred correction integration.
class TimeIntegrationSchemeSDC : public TimeIntegrationScheme
{
public:
    TimeIntegrationSchemeSDC(std::string variant, unsigned int order,
                             std::vector<NekDouble> freeParams)
        : TimeIntegrationScheme(variant, order, freeParams),
          m_name("SpectralDeferredCorrection")
    {

        ASSERTL0(variant == "Equidistant" ||
                     variant == "GaussLobattoLegendre" ||
                     variant == "GaussRadauLegendre" ||
                     variant == "GaussLobattoChebyshev" ||
                     variant == "GaussRadauChebyshev",
                 "unknow variant (quadrature) type");

        ASSERTL0(1 <= order, "Spectral Deferred Correction Time integration "
                             "scheme bad order numbers (>=1): " +
                                 std::to_string(order));
        ASSERTL0(freeParams.size() == 2,
                 "SDC Time integration scheme invalid number "
                 "of free parameters, expected two "
                 "<theta, number of quadrature>, received " +
                     std::to_string(freeParams.size()));
        ASSERTL0(0.0 <= freeParams[0] && freeParams[0] <= 1.0,
                 "Spectral Deferred Correction Time integration "
                 "scheme bad parameter numbers (0.0 - 1.0): " +
                     std::to_string(freeParams[0]));
        ASSERTL0(1 <= int(freeParams[1]),
                 "Spectral Deferred Correction Time integration "
                 "scheme bad quadrature numbers (>=1): " +
                     std::to_string(int(freeParams[1])));

        m_variant    = variant;
        m_order      = order;
        m_freeParams = freeParams;
        m_theta      = freeParams[0];

        // Compute quadrature points and weights
        if (variant == "Equidistant")
        {
            ASSERTL0(int(m_freeParams[1]) >= order,
                     "Spectral Deferred Correction Time integration "
                     "Maximum order (<= n): " +
                         std::to_string(int(m_freeParams[1])));

            ASSERTL0(2 <= freeParams[1],
                     "Spectral Deferred Correction Time integration "
                     "scheme bad quadrature numbers (>=1): " +
                         std::to_string(freeParams[1]));

            m_nQuadPts = int(m_freeParams[1]);
            m_points   = SingleArray(m_nQuadPts, 0.0);
            m_weights  = SingleArray(m_nQuadPts, 0.0);
            for (int i = 0; i < m_nQuadPts - 1; i++)
            {
                m_points[i + 1] = m_points[i] + 1.0 / (m_nQuadPts - 1);
            }
        }
        else if (variant == "GaussLobattoLegendre")
        {
            ASSERTL0(2 * int(m_freeParams[1]) - 2 >= order,
                     "Spectral Deferred Correction Time integration "
                     "Maximum order (<= 2 * n - 2): " +
                         std::to_string(2 * int(m_freeParams[1]) - 2));

            // m_order    = 2 * int(m_freeParams[1]) - 2;
            m_nQuadPts = int(m_freeParams[1]);
            m_points   = SingleArray(m_nQuadPts, 0.0);
            m_weights  = SingleArray(m_nQuadPts, 0.0);
            Polylib::zwglj(&m_points[0], &m_weights[0], int(m_freeParams[1]),
                           0.0, 0.0);
        }
        else if (variant == "GaussRadauLegendre")
        {
            ASSERTL0(2 * int(m_freeParams[1]) - 1 >= order,
                     "Spectral Deferred Correction Time integration "
                     "Maximum order (<= 2 * n - 1): " +
                         std::to_string(2 * int(m_freeParams[1]) - 1));

            // m_order     = 2 * int(m_freeParams[1]) - 1;
            m_nQuadPts = int(m_freeParams[1]) + 1;
            m_points   = SingleArray(m_nQuadPts, -1.0);
            m_weights  = SingleArray(m_nQuadPts, 0.0);
            if (int(m_freeParams[1]) == 1)
            {
                m_points[1] = 1.0;
            }
            else
            {
                Polylib::zwgrjp(&m_points[1], &m_weights[1],
                                int(m_freeParams[1]), 0.0, 0.0);
            }
        }
        else if (variant == "GaussLobattoChebyshev")
        {
            ASSERTL0(2 * int(m_freeParams[1]) - 2 >= order,
                     "Spectral Deferred Correction Time integration "
                     "Maximum order (<= 2 * n - 2): " +
                         std::to_string(2 * int(m_freeParams[1]) - 2));

            // m_order    = 2 * int(m_freeParams[1]) - 2;
            m_nQuadPts = int(m_freeParams[1]);
            m_points   = SingleArray(m_nQuadPts, 0.0);
            m_weights  = SingleArray(m_nQuadPts, 0.0);
            Polylib::zwglj(&m_points[0], &m_weights[0], int(m_freeParams[1]),
                           -0.5, -0.5);
        }
        else if (variant == "GaussRadauChebyshev")
        {
            ASSERTL0(2 * int(m_freeParams[1]) - 2 >= order,
                     "Spectral Deferred Correction Time integration "
                     "Maximum order (<= 2 * n - 2): " +
                         std::to_string(2 * int(m_freeParams[1]) - 2));

            // m_order     = 2 * int(m_freeParams[1]) - 1;
            m_nQuadPts = int(m_freeParams[1]) + 1;
            m_points   = SingleArray(m_nQuadPts, -1.0);
            m_weights  = SingleArray(m_nQuadPts, 0.0);
            if (int(m_freeParams[1]) == 1)
            {
                m_points[1] = 1.0;
            }
            else
            {
                Polylib::zwgrjp(&m_points[1], &m_weights[1],
                                int(m_freeParams[1]), -0.5, -0.5);
            }
        }
        // Note: m_weights is not used

        // Rescale quadrature points to [0, 1]
        NekDouble denom = m_points[m_nQuadPts - 1] - m_points[0];
        NekDouble x0    = m_points[0];
        for (int i = 0; i < m_nQuadPts; i++)
        {
            m_points[i] = (m_points[i] - x0) / denom;
        }
    }

    /// Destructor
    virtual ~TimeIntegrationSchemeSDC()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(
        std::string variant, unsigned int order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<TimeIntegrationSchemeSDC>::AllocateSharedPtr(
                variant, order, freeParams);

        return p;
    }

    static std::string className;

protected:
    LUE virtual std::string v_GetName() const override;
    LUE virtual std::string v_GetVariant() const override;
    LUE virtual unsigned int v_GetOrder() const override;
    LUE virtual std::vector<NekDouble> v_GetFreeParams() const override;
    LUE virtual TimeIntegrationSchemeType v_GetIntegrationSchemeType()
        const override;
    LUE virtual NekDouble v_GetTimeStability() const override;
    LUE virtual unsigned int v_GetNumIntegrationPhases() const override;

    /**
     * \brief Gets the solution vector of the ODE
     */
    virtual const TripleArray &v_GetSolutionVector() const override;

    /**
     * \brief Sets the solution vector of the ODE
     */
    virtual void v_SetSolutionVector(const int Offset,
                                     const DoubleArray &y) override;

    // The worker methods from the base class that are virtual
    LUE virtual void v_InitializeScheme(
        const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
        const TimeIntegrationSchemeOperators &op) override;

    LUE virtual ConstDoubleArray &v_TimeIntegrate(
        const int timestep, const NekDouble delta_t,
        const TimeIntegrationSchemeOperators &op) override;

    LUE virtual void v_print(std::ostream &os) const override;
    LUE virtual void v_printFull(std::ostream &os) const override;

    // Variables common to all schemes.
    std::string m_name;
    std::string m_variant;
    unsigned int m_order{0}, m_nQuadPts{0};
    bool m_initialized = false;
    std::vector<NekDouble> m_freeParams;
    NekDouble m_time, m_theta = 0.0;

    TimeIntegrationSchemeType m_schemeType{eSpectralDeferredCorrection};

    // Storage of previous states and associated timesteps.
    TripleArray m_Y; /// Array containing the stage values
    TripleArray m_F; /// Array corresponding to the stage Derivatives
    TripleArray m_Fint;
    DoubleArray m_Fn;
    DoubleArray m_tmp;
    DoubleArray m_wMat;
    SingleArray m_coeffs;
    SingleArray m_points;
    SingleArray m_weights;

    int m_nvars{0};   // Number of variables in the integration scheme.
    int m_npoints{0}; // Number of points    in the integration scheme.

    void UpdateIntegratedFlux(const NekDouble &delta_t);
    NekDouble EvaluateInt(NekDouble *coeffs, const int npoints,
                          const NekDouble &t1, const NekDouble &t2) const;
    void EvaluateCoeffs(NekDouble *coeffs, NekDouble *points, const int index,
                        const int npoints);

}; // end class TimeIntegrationSchemeSDC

} // end namespace LibUtilities
} // end namespace Nektar

#endif
