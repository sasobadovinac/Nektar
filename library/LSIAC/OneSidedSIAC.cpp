#include "OneSidedSIAC.h"

namespace Nektar
{
namespace LSIAC
{

OneSidedSIAC::OneSidedSIAC(const int Order,
                           OneSidedSIAC::OneSidedFilterType filter,
                           const int Derivative)
{

    m_filter = filter;
    switch (filter)
    {
        case (OneSidedFilterType::BASIC_SIAC_2kp1):
            m_order  = Order;
            m_nBSpl  = 2 * (m_order - 1) + 1;
            m_coeffs = Array<OneD, NekDouble>((2 * (m_order - 1) + 1), 0.0);
            m_cenBSplinePtr = std::make_shared<CentralBSplines>(Order);
            break;
        case (OneSidedFilterType::VAN_SIAC_4kp1):
            m_order         = Order;
            m_nBSpl         = 4 * (m_order - 1) + 1;
            m_coeffs        = Array<OneD, NekDouble>(m_nBSpl, 0.0);
            m_cenBSplinePtr = std::make_shared<CentralBSplines>(Order);
            break;
        case (OneSidedFilterType::Der_SMOOTH_BASIC_SIAC_2kp1):
            m_order         = Order;
            m_nthDer        = Derivative;
            m_nBSpl         = 2 * (m_order - 1) + 1 + m_nthDer;
            m_coeffs        = Array<OneD, NekDouble>(m_nBSpl, 0.0);
            m_cenBSplinePtr = std::make_shared<CentralBSplines>(Order);
            break;
        case (OneSidedFilterType::Der_BASIC_SIAC_2kp1):
            m_order         = Order;
            m_nthDer        = Derivative;
            m_nBSpl         = 2 * (m_order - 1) + 1;
            m_coeffs        = Array<OneD, NekDouble>(m_nBSpl, 0.0);
            m_cenBSplinePtr = std::make_shared<CentralBSplines>(Order);
            break;
        case (OneSidedFilterType::Der_SMOOTH_BASIC_SIAC_4kp1):
            m_order         = Order;
            m_nthDer        = Derivative;
            m_nBSpl         = 4 * (m_order - 1) + 1 + m_nthDer;
            m_coeffs        = Array<OneD, NekDouble>(m_nBSpl, 0.0);
            m_cenBSplinePtr = std::make_shared<CentralBSplines>(Order);
            break;
        case (OneSidedFilterType::Der_BASIC_SIAC_4kp1):
            m_order         = Order;
            m_nthDer        = Derivative;
            m_nBSpl         = 4 * (m_order - 1) + 1;
            m_coeffs        = Array<OneD, NekDouble>(m_nBSpl, 0.0);
            m_cenBSplinePtr = std::make_shared<CentralBSplines>(Order);
            break;
        case (OneSidedFilterType::XLi_SIAC_2kp2):
            m_order  = Order;
            m_nthDer = 0;
            m_nBSpl  = 2 * (m_order - 1) + 2;
            m_knotMatrix.resize(2 * (m_order - 1) + 2);
            m_coeffs = Array<OneD, NekDouble>(2 * (m_order - 1) + 2, 0.0);
            m_genBSplinePtr = std::make_shared<GeneralBSplines>(Order);
            break;
        case (OneSidedFilterType::Der_XLi_SIAC_2kp2):
            m_order         = Order;
            m_nthDer        = Derivative;
            m_nBSpl         = 2 * (m_order - 1) + 2 + m_nthDer;
            m_coeffs        = Array<OneD, NekDouble>(m_nBSpl, 0.0);
            m_genBSplinePtr = std::make_shared<GeneralBSplines>(Order);
            m_knotMatrix.resize(m_nBSpl); // dont know if this is necessary.
                                          // Might need to adjust accordingly.
            break;
        case (OneSidedFilterType::N_Der_SMOOTH_BASIC_SIAC_2kp1):
            m_order  = Order + Derivative;
            m_nthDer = 0;
            m_nBSpl  = 2 * (Order - 1) + 1;
            m_coeffs = Array<OneD, NekDouble>(m_nBSpl, 0.0);
            m_cenBSplinePtr =
                std::make_shared<CentralBSplines>(Order + Derivative);
        default:
            cout << "OneSidedSIAC::Constructor 1st.";
            break;
    }
    //	EvaluateCoefficients();
}

bool OneSidedSIAC::v_GetBreakPts(const NekDouble scaling,
                                 vector<NekDouble> &valT, const NekDouble shift)
{
    NekDouble tmin, tmax;
    switch (m_filter)
    {
        case (OneSidedFilterType::BASIC_SIAC_2kp1):
        case (OneSidedFilterType::VAN_SIAC_4kp1):
        case (OneSidedFilterType::Der_SMOOTH_BASIC_SIAC_2kp1):
        case (OneSidedFilterType::Der_BASIC_SIAC_2kp1):
        case (OneSidedFilterType::Der_SMOOTH_BASIC_SIAC_4kp1):
        case (OneSidedFilterType::Der_BASIC_SIAC_4kp1):
        case (OneSidedFilterType::N_Der_SMOOTH_BASIC_SIAC_2kp1):
            tmin = -((m_order) / 2.0 + (m_nBSpl - 1.0) / 2.0) * scaling + shift;
            tmax = ((m_order) / 2.0 + (m_nBSpl - 1.0) / 2.0) * scaling + shift;
            valT.clear();
            for (NekDouble t = tmin; t <= tmax; t += scaling)
            {
                valT.push_back(t);
            }
            break;
        case (OneSidedFilterType::XLi_SIAC_2kp2):
            tmin = -((m_order) / 2.0 + (m_order - 1.0)) * scaling + shift;
            tmax = ((m_order) / 2.0 + (m_order - 1.0)) * scaling + shift;
            valT.clear();
            for (NekDouble t = tmin; t <= tmax; t += scaling)
            {
                valT.push_back(t);
            }
            break;
        case (OneSidedFilterType::Der_XLi_SIAC_2kp2):
            tmin = -((m_order) / 2.0 + (m_order - 1.0)) * scaling + shift;
            tmax = ((m_order) / 2.0 + (m_order - 1.0)) * scaling + shift;
            valT.clear();
            if (shift < 0)
            {
                for (NekDouble t = tmin - m_nthDer * scaling; t <= tmax;
                     t += scaling)
                {
                    valT.push_back(t);
                }
            }
            else
            {
                for (NekDouble t = tmin; t <= tmax + m_nthDer * scaling;
                     t += scaling)
                {
                    valT.push_back(t);
                }
            }
            break;
        default:
            assert(false &&
                   "Missed all switch cases. Something is not right..");
    }
    return true;
}

bool OneSidedSIAC::v_GetFilterRange(NekDouble scaling, NekDouble &tmin,
                                    NekDouble &tmax, const NekDouble shift)
{
    switch (m_filter)
    {
        case (OneSidedFilterType::BASIC_SIAC_2kp1):
        case (OneSidedFilterType::VAN_SIAC_4kp1):
        case (OneSidedFilterType::Der_SMOOTH_BASIC_SIAC_2kp1):
        case (OneSidedFilterType::Der_BASIC_SIAC_2kp1):
        case (OneSidedFilterType::Der_SMOOTH_BASIC_SIAC_4kp1):
        case (OneSidedFilterType::Der_BASIC_SIAC_4kp1):
        case (OneSidedFilterType::N_Der_SMOOTH_BASIC_SIAC_2kp1):
            tmin = -((m_order) / 2.0 + (m_nBSpl - 1.0) / 2.0) * scaling + shift;
            tmax = ((m_order) / 2.0 + (m_nBSpl - 1.0) / 2.0) * scaling + shift;
            break;
        case (OneSidedFilterType::XLi_SIAC_2kp2):
            tmin = -1.0 * ((m_order) / 2.0 + (m_order - 1.0)) * scaling + shift;
            tmax = ((m_order) / 2.0 + (m_order - 1.0)) * scaling + shift;
            break;
        case (OneSidedFilterType::Der_XLi_SIAC_2kp2):
            tmin = -1.0 * ((m_order) / 2.0 + (m_order - 1.0)) * scaling + shift;
            tmax = ((m_order) / 2.0 + (m_order - 1.0)) * scaling + shift;
            if (TOLERENCE < std::abs(shift))
            {
                if (shift < 0)
                {
                    tmin = tmin - m_nthDer * scaling;
                }
                else
                {
                    tmax = tmax + m_nthDer * scaling;
                }
            }
            break;
        default:
            assert(false &&
                   "Missed all switch cases. Something is not right..");
    }
    return true;
}

/*
 *
 *
 */
bool OneSidedSIAC::v_EvaluateFilter(const Array<OneD, NekDouble> &x_pos,
                                    Array<OneD, NekDouble> &t_vals,
                                    const NekDouble meshScaling,
                                    const NekDouble meshShift,
                                    const bool evalCoeff)
{
    int nq;
    nq = x_pos.num_elements();
    Vmath::Fill(nq, 0.0, t_vals, 1);
    Array<OneD, NekDouble> t_valTemp(nq);
    if (evalCoeff)
    { // Remeber kernel shift is opposite to mesh shift.
        if (m_filter == OneSidedFilterType::XLi_SIAC_2kp2)
        {
            EvaluateCoefficients(meshShift / meshScaling);
        }
        else
        {
            EvaluateCoefficients(-1.0 * meshShift / meshScaling);
        }
    }
    switch (m_filter)
    {
        case (OneSidedFilterType::BASIC_SIAC_2kp1):
        case (OneSidedFilterType::VAN_SIAC_4kp1):
            for (int i = 0; i < m_nBSpl; i++)
            {
                //  NekDouble it = -m_R*0.5 + i;
                NekDouble it = -(m_nBSpl - 1.0) / 2.0 + i;
                // Some kind of mesh shifting needs to be done while evaluating
                // the Bsplines. Currently this is to be done when necessary.
                // ???
                m_cenBSplinePtr->EvaluateBSplines(
                    x_pos, t_valTemp,
                    (((NekDouble)it) * meshScaling + meshShift), meshScaling);
                Vmath::Smul(nq, m_coeffs[m_nBSpl - 1 - i] / meshScaling,
                            t_valTemp, 1, t_valTemp, 1);
                Vmath::Vadd(nq, t_valTemp, 1, t_vals, 1, t_vals, 1);
            }
            break;
        case (OneSidedFilterType::Der_SMOOTH_BASIC_SIAC_2kp1):
        case (OneSidedFilterType::Der_BASIC_SIAC_2kp1):
        case (OneSidedFilterType::Der_SMOOTH_BASIC_SIAC_4kp1):
        case (OneSidedFilterType::Der_BASIC_SIAC_4kp1):
        case (OneSidedFilterType::N_Der_SMOOTH_BASIC_SIAC_2kp1):
            for (int i = 0; i < m_nBSpl; i++)
            {
                //  NekDouble it = -m_R*0.5 + i;
                NekDouble it = -(m_nBSpl - 1.0) / 2.0 + i;
                // Some kind of mesh shifting needs to be done while evaluating
                // the Bsplines. Currently this is to be done when necessary.
                // ???
                m_cenBSplinePtr->EvaluateBSplines(
                    x_pos, t_valTemp,
                    (((NekDouble)it) * meshScaling + meshShift), meshScaling);
                //		Vmath::Smul(nq,
                // m_coeffs[m_nBSpl-1-i]/meshScaling , t_valTemp,1,t_valTemp,1);
                Vmath::Smul(nq,
                            m_coeffs[m_nBSpl - 1 - i] /
                                std::pow(meshScaling, 1 + m_nthDer),
                            t_valTemp, 1, t_valTemp, 1);
                Vmath::Vadd(nq, t_valTemp, 1, t_vals, 1, t_vals, 1);
            }
            break;
        case (OneSidedFilterType::XLi_SIAC_2kp2):
            m_nBSpl = 2 * (m_order - 1) + 2;
            // loop through all the bsplines.
            for (int i = 0; i < m_nBSpl; i++)
            {
                m_genBSplinePtr->EvaluateBSplines(x_pos, m_knotMatrix[i], 0,
                                                  t_valTemp, 0.0, meshScaling);
                Vmath::Smul(nq, m_coeffs[i] / meshScaling, t_valTemp, 1,
                            t_valTemp, 1);
                Vmath::Vadd(nq, t_valTemp, 1, t_vals, 1, t_vals, 1);
            }
            break;
        case (OneSidedFilterType::Der_XLi_SIAC_2kp2):
            m_nBSpl = m_knotMatrix.size();
            for (int i = 0; i < m_nBSpl; i++)
            {
                m_genBSplinePtr->EvaluateBSplines(x_pos, m_knotMatrix[i], 0,
                                                  t_valTemp, 0.0, meshScaling);
                Vmath::Smul(nq,
                            std::pow(-1.0, m_nthDer) * m_coeffs[i] /
                                std::pow(meshScaling, 1 + m_nthDer),
                            t_valTemp, 1, t_valTemp, 1);
                Vmath::Vadd(nq, t_valTemp, 1, t_vals, 1, t_vals, 1);
            }
            break;
        default:
            assert(false &&
                   "Missed all switch cases. Something is not right..");
            cout << "Assert here, OneSidedFilter::v_EvaluateFilter" << endl;
    }
    return true;
}

/*
 *
 * KernelShift : The kernel shift is supposed to be independent of mesh scaling.
 *
 */
bool OneSidedSIAC::v_EvaluateCoefficients(const NekDouble kernelShift)
{
    bool retValue = true;
    switch (m_filter)
    {
        case (OneSidedFilterType::BASIC_SIAC_2kp1):
        case (OneSidedFilterType::VAN_SIAC_4kp1):
            CalCoeffForWideSymKernel(m_order - 1, m_nBSpl, m_coeffs,
                                     kernelShift);
            // cout << "kernelShift" << kernelShift << endl;
            break;
        case (OneSidedFilterType::Der_SMOOTH_BASIC_SIAC_2kp1):
        case (OneSidedFilterType::Der_BASIC_SIAC_2kp1):
        case (OneSidedFilterType::Der_SMOOTH_BASIC_SIAC_4kp1):
        case (OneSidedFilterType::Der_BASIC_SIAC_4kp1):
            CalCoeffForWideSymKernel(m_order - 1 + m_nthDer, m_nthDer, m_nBSpl,
                                     m_coeffs, kernelShift);
            // printNekArray(m_coeffs,0);
            CalCoeffForCenBSplDerivatives(m_order - 1 + m_nthDer, m_nthDer,
                                          m_nBSpl, m_coeffs);
            // printNekArray(m_coeffs,0);
            // cout << "No of BSplines: " <<m_nBSpl << endl;
            break;
        case (OneSidedFilterType::XLi_SIAC_2kp2):
            //			if (TOLERENCE < std::abs(kernelShifted) )
            {
                int numKnotsPvec = (m_order - 1) + 2;
                int nBSpl        = 2 * (m_order - 1) + 2;
                NekDouble leftknotOfsym =
                    -1.0 * ((m_order) / 2.0 + (m_order - 1.0));
                NekDouble rightknotOfsym = ((m_order) / 2.0 + (m_order - 1.0));
                NekDouble kernelShifted  = kernelShift;
                // create and fill KnotMatrix.
                if (kernelShifted > 0)
                { // Right shift.
                    // The first one is general B-Spline. All the others are
                    // central bsplines.
                    // kernelShifted =-1*leftknotOfsym;
                    m_knotMatrix.clear();
                    for (int tm = 0; tm < nBSpl; tm++)
                    {
                        if (0 == tm)
                        { // for tm=0; i.e. first generalB-Spline.
                            std::vector<NekDouble> genVec(
                                numKnotsPvec, leftknotOfsym + kernelShifted);
                            genVec[numKnotsPvec - 1] =
                                leftknotOfsym + kernelShifted + 1.0;
                            m_knotMatrix.push_back(genVec);
                        }
                        else
                        {
                            std::vector<NekDouble> genVec(numKnotsPvec);
                            for (int tv = 0; tv < numKnotsPvec; tv++)
                            {
                                genVec[tv] = leftknotOfsym + kernelShifted +
                                             tv + (tm - 1);
                            }
                            m_knotMatrix.push_back(genVec);
                        }
                    }
                }
                else // KernelShift < 0
                {    // leftshift
                    m_knotMatrix.clear();
                    for (int tm = 0; tm < nBSpl; tm++)
                    {
                        if ((2 * (m_order - 1) + 1) == tm)
                        { // for tm=0; i.e. first generalB-Spline.
                            std::vector<NekDouble> genVec(
                                numKnotsPvec, rightknotOfsym + kernelShifted);
                            genVec[0] = rightknotOfsym + kernelShifted - 1.0;
                            m_knotMatrix.push_back(genVec);
                        }
                        else
                        {
                            std::vector<NekDouble> genVec(numKnotsPvec);
                            for (int tv = 0; tv < numKnotsPvec; tv++)
                            {
                                genVec[tv] =
                                    leftknotOfsym + kernelShifted + tv + tm;
                            }
                            m_knotMatrix.push_back(genVec);
                        }
                    }
                }
                // Evaluate the coefficents and store in m_coeff array.
                // CalCoeffForKnotMatrixVec(m_order-1, m_knotMatrix, m_coeffs);
                CalCoeffForKnotMatrixVec_Hanieh(m_order - 1, m_knotMatrix,
                                                m_coeffs);
                // CalDerivativesForKnotMatrixVec_Hanieh(m_order-1,m_nthDer,
                // m_knotMatrix, m_coeffs);
            }
            //			else{
            //				assert((std::abs(kernelShifted) >TOLERENCE)
            //&&"you should be using symmetric filter. Somthing is not right.");
            // retValue = false;
            //			}
            break;
        case (OneSidedFilterType::Der_XLi_SIAC_2kp2):
        {
            NekDouble leftknotOfsym =
                -1.0 * ((m_order) / 2.0 + (m_order - 1.0));
            NekDouble rightknotOfsym = ((m_order) / 2.0 + (m_order - 1.0));
            int numKnotsPvec         = (m_order - 1) + 2 + m_nthDer;
            int nBSpl                = 2 * (m_order - 1) + 2;
            // create and fill KnotMatrix.
            if (kernelShift < 0)
            { // Right shift.
                // The first one is general B-Spline. All the others are central
                // bsplines.
                m_knotMatrix.clear();
                for (int tm = 0; tm < nBSpl; tm++)
                {
                    if (0 == tm)
                    { // for tm=0; i.e. first generalB-Spline.
                        std::vector<NekDouble> genVec(numKnotsPvec, -10.0);
                        for (int tv = 0; tv < numKnotsPvec; tv++)
                        {
                            genVec[tv] =
                                (leftknotOfsym - kernelShift) +
                                std::max(0.0, (NekDouble)(tv - (m_order - 1)));
                        }
                        m_knotMatrix.push_back(genVec);
                    }
                    else
                    {
                        std::vector<NekDouble> genVec(numKnotsPvec);
                        for (int tv = 0; tv < numKnotsPvec; tv++)
                        {
                            genVec[tv] =
                                -1.0 + tv + tm + (leftknotOfsym - kernelShift);
                        }
                        m_knotMatrix.push_back(genVec);
                    }
                }
            }
            else // KernelShift > 0
            {    // leftshift
                m_knotMatrix.clear();
                for (int tm = 0; tm < nBSpl; tm++)
                {
                    if ((2 * (m_order - 1) + 1) == tm)
                    { // for tm=0; i.e. first generalB-Spline.
                        std::vector<NekDouble> genVec(numKnotsPvec, -10.0);
                        for (int tv = 0; tv < numKnotsPvec; tv++)
                        {
                            genVec[tv] =
                                (rightknotOfsym - kernelShift) +
                                std::min(0.0, (NekDouble)(tv - m_nthDer - 1));
                        }
                        m_knotMatrix.push_back(genVec);
                    }
                    else
                    {
                        std::vector<NekDouble> genVec(numKnotsPvec);
                        for (int tv = 0; tv < numKnotsPvec; tv++)
                        {
                            genVec[tv] = -3 * (m_order - 1) - 1.0 - m_nthDer +
                                         tv + tm +
                                         (rightknotOfsym - kernelShift);
                        }
                        m_knotMatrix.push_back(genVec);
                    }
                }
            }
            // Evaluate the coefficents and store in m_coeff array.
            // CalCoeffForKnotMatrixVec(m_order-1, m_knotMatrix, m_coeffs);
            CalCoeffForKnotMatrixVec_Hanieh(m_order - 1 + m_nthDer,
                                            m_knotMatrix, m_coeffs);
            // printNekArray(m_coeffs,0);
            CalDerivativesForKnotMatrixVec_Hanieh(m_order - 1, m_nthDer,
                                                  m_knotMatrix, m_coeffs);
            // printNekArray(m_coeffs,0);
        }
        // assert("Need to code this.");
        break;
        case (OneSidedFilterType::N_Der_SMOOTH_BASIC_SIAC_2kp1):
            CalCoeffForWideSymKernel(m_order - 1, m_nBSpl, m_coeffs,
                                     kernelShift);
            break;
            /*
                            case (OneSidedFilterType::Der_XLi_SIAC_2kp2):
                                    {
                                            int numKnotsPvec =
               (m_order-1)+2+m_nthDer; int nBSpl = 2*(m_order-1) +2; NekDouble
               leftknotOfsym = -1.0*( (m_order)/2.0 + (m_order-1.0) +m_nthDer);
                                            NekDouble rightknotOfsym = (
               (m_order)/2.0 + (m_order-1.0) + m_nthDer );
                                            // create and fill KnotMatrix.
                                            if( kernelShift>0 )
                                            { // Right shift.
                                                    // The first one is general
               B-Spline. All the others are central bsplines.
                                                    m_knotMatrix.clear();
                                                    for ( int tm = 0; tm <
               nBSpl;tm++)
                                                    {
                                                            if ( 0 == tm)
                                                            { // for tm=0; i.e.
               first generalB-Spline. std::vector<NekDouble>
               genVec(numKnotsPvec,leftknotOfsym+kernelShift); for (int tv=0; tv
               < numKnotsPvec;tv++)
                                                                    {
                                                                            genVec[tv]
               = kernelShift+ std::max(0.0, leftknotOfsym +kernelShift+tv);
                                                                    }
                                                                    genVec[numKnotsPvec-1]
               = leftknotOfsym+kernelShift +1.0; m_knotMatrix.push_back(genVec);
                                                            }
                                                            else
                                                            {
                                                                    std::vector<NekDouble>
               genVec(numKnotsPvec); for (int tv =0; tv < numKnotsPvec; tv++)
                                                                    {
                                                                            genVec[tv]
               = leftknotOfsym + kernelShift + tv + (tm-1);
                                                                    }
                                                                    m_knotMatrix.push_back(genVec);
                                                            }
                                                    }
                                            }
                                            else // KernelShift < 0
                                            { //leftshift
                                                    m_knotMatrix.clear();
                                                    for ( int tm = 0; tm <
               nBSpl;tm++)
                                                    {
                                                            if (
               (2*(m_order-1)+1 ) == tm ) { // for tm=0; i.e. first
               generalB-Spline. std::vector<NekDouble>
               genVec(numKnotsPvec,rightknotOfsym+kernelShift); genVec[0] =
               rightknotOfsym+kernelShift - 1.0; for (int tv=0; tv <
               numKnotsPvec;tv++)
                                                                    {
                                                                            genVec[tv]
               = std::min(0.0, rightknotOfsym+kernelShift +tv);
                                                                    }
                                                                    m_knotMatrix.push_back(genVec);
                                                            }
                                                            else
                                                            {
                                                                    std::vector<NekDouble>
               genVec(numKnotsPvec); for (int tv =0; tv < numKnotsPvec; tv++)
                                                                    {
                                                                            genVec[tv]
               = leftknotOfsym + kernelShift + tv + tm;
                                                                    }
                                                                    m_knotMatrix.push_back(genVec);
                                                            }
                                                    }
                                            }
                                            // Evaluate the coefficents and
               store in m_coeff array.
                            //CalCoeffForKnotMatrixVec(m_order-1, m_knotMatrix,
               m_coeffs); CalCoeffForKnotMatrixVec_Hanieh(m_order-1+m_nthDer,
               m_knotMatrix, m_coeffs);
                                    }
                                    //assert("Need to code this.");
                                    break;
             */
        default:
            assert("Missed all switch cases. Something is not right..");
            retValue = false;
            break;
    }
    return retValue;
}
} // namespace LSIAC
} // namespace Nektar
