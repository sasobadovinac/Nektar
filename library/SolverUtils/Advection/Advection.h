///////////////////////////////////////////////////////////////////////////////
//
// File: Advection.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Abstract base class for advection.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_ADVECTION
#define NEKTAR_SOLVERUTILS_ADVECTION

#include <string>
#include <functional>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

namespace Nektar
{
namespace SolverUtils
{

/// Defines a callback function which evaluates the flux vector \f$ F(u)
/// \f$ in a conservative advection of the form \f$ \nabla\cdot F(u)
/// \f$.
typedef std::function<void (
    const Array<OneD, Array<OneD, NekDouble> >&,
    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&)>
    AdvectionFluxVecCB;

typedef std::function<void (
    const Array<OneD, Array<OneD, NekDouble> >&,
    const Array<OneD, Array<OneD, NekDouble> >&,
    Array<OneD, Array<OneD, NekDouble> >&)>
    AdvectionFluxVecCBTraceNormal;

/**
 * @brief An abstract base class encapsulating the concept of advection
 * of a vector field.
 *
 * Subclasses override the Advection::v_InitObject function to
 * initialise the object and the Advection::v_Advect function to
 * evaluate the advection of the vector field.
 */
class Advection
{
public:
    /// Interface function to initialise the advection object.
    SOLVER_UTILS_EXPORT void InitObject(
        LibUtilities::SessionReaderSharedPtr               pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);

    /// Interface function to advect the vector field.
    SOLVER_UTILS_EXPORT void Advect(
        const int nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &advVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray,
        const NekDouble                                   &time,
        const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
        const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);
    
    /// Interface function to advect the Volume field.
    SOLVER_UTILS_EXPORT void AdvectVolumeFlux(
        const int                                         nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, Array<OneD, NekDouble>>         &pAdvVel,
        const Array<OneD, Array<OneD, NekDouble>>         &pInarray,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>  &pVolumeFlux,
        const NekDouble                                   &pTime)
    {
        v_AdvectVolumeFlux(nConvectiveFields, pFields, pAdvVel, pInarray, pVolumeFlux,
                   pTime);
    }

    /// Interface function to advect the Trace field.
    SOLVER_UTILS_EXPORT void AdvectTraceFlux(
        const int                                         nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, Array<OneD, NekDouble>>         &pAdvVel,
        const Array<OneD, Array<OneD, NekDouble>>         &pInarray,
        Array<OneD, Array<OneD, NekDouble>>               &pTraceFlux,
        const NekDouble                                   &pTime,
        const Array<OneD, Array<OneD, NekDouble>>         &pFwd = NullNekDoubleArrayofArray,
        const Array<OneD, Array<OneD, NekDouble>>         &pBwd = NullNekDoubleArrayofArray)
    {
        v_AdvectTraceFlux(nConvectiveFields, pFields, pAdvVel, pInarray, 
                  pTraceFlux, pTime, pFwd, pBwd);
    }
    
    SOLVER_UTILS_EXPORT void Advect_coeff(
        const int nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &advVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray,
        const NekDouble                                   &time,
        const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
        const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray,
        const bool                                       &flagFreezeJac=false);


    /**
     * @brief Set the flux vector callback function.
     *
     * This routine is a utility function to avoid the explicit use of
     * std::bind. A function and object can be passed to this function
     * instead.
     */
    template<typename FuncPointerT, typename ObjectPointerT>
    void SetFluxVector(FuncPointerT func, ObjectPointerT obj)
    {
        m_fluxVector = std::bind(
            func, obj, std::placeholders::_1, std::placeholders::_2);
    }

    template<typename FuncPointerT, typename ObjectPointerT>
    void SetFluxVectorMF(FuncPointerT func, ObjectPointerT obj)
    {
        m_fluxVectorMF = std::bind(
            func, obj, std::placeholders::_1, std::placeholders::_2);
    }

    template<typename FuncPointerT, typename ObjectPointerT>
    void SetFluxVectortraceMF(FuncPointerT func, ObjectPointerT obj)
    {
        m_fluxVectortraceMF = std::bind(
            func, obj, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    }

    /**
     * @brief Set a Riemann solver object for this advection object.
     *
     * @param riemann  The RiemannSolver object.
     */
    inline void SetRiemannSolver(RiemannSolverSharedPtr riemann)
    {
        m_riemann = riemann;
    }

    /**
     * @brief Set the flux vector callback function.
     *
     * @param fluxVector  The callback function to override.
     */
    inline void SetFluxVector(AdvectionFluxVecCB fluxVector)
    {
        m_fluxVector = fluxVector;
    }

    /**
     * @brief Set the base flow used for linearised advection objects.
     *
     * @param inarray   Vector to use as baseflow
     */
    inline void SetBaseFlow(
            const Array<OneD, Array<OneD, NekDouble> >& inarray,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
    {
        v_SetBaseFlow(inarray, fields);
    }

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
    SOLVER_UTILS_EXPORT void CalcTraceJac(
        const int                                          nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &advVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >        &pBwd,
        DNekBlkMatSharedPtr &FJac,
        DNekBlkMatSharedPtr &BJac);

    SOLVER_UTILS_EXPORT void AddTraceJacToMat(
            const int                                           nConvectiveFields,
            const int                                           nSpaceDim,
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields,
            const Array<OneD, DNekBlkMatSharedPtr>              &TracePntJacCons,
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >      &gmtxarray,
            const Array<OneD, DNekBlkMatSharedPtr>              &TracePntJacGrad        = NullArrayDNekBlkMatSharedPtr,
            const Array<OneD, Array<OneD, NekDouble> >          &TracePntJacGradSign    = NullNekDoubleArrayofArray)
    {
        v_AddTraceJacToMat(nConvectiveFields,nSpaceDim, pFields, TracePntJacCons, gmtxarray,TracePntJacGrad,TracePntJacGradSign);
    }

    SOLVER_UTILS_EXPORT void AddVolumJacToMat( 
        const Array<OneD, MultiRegions::ExpListSharedPtr>                       &pFields,
        const int                                                               &nConvectiveFields,
        const Array<OneD, const Array<OneD,  Array<OneD, 
            Array<OneD,  Array<OneD,  NekDouble> > > > >                        &ElmtJacArray,
        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >                          &gmtxarray)
    {
        v_AddVolumJacToMat(pFields,nConvectiveFields,ElmtJacArray,gmtxarray);
    }

    SOLVER_UTILS_EXPORT void NumCalRiemFluxJac( 
        const int                                          nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &AdvVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >        &pBwd,
        DNekBlkMatSharedPtr &FJac,
        DNekBlkMatSharedPtr &BJac)
    {
        v_NumCalRiemFluxJac(nConvectiveFields,fields,AdvVel,inarray,pFwd,pBwd,FJac,BJac);
    }

    void CoutBlkMat(
        DNekBlkMatSharedPtr &gmtx, 
        const unsigned int nwidthcolm=12);

    void CoutStandardMat(
        DNekMatSharedPtr &loc_matNvar,
        const unsigned int nwidthcolm=12);

    void Cout1DArrayBlkMat(
        Array<OneD, DNekBlkMatSharedPtr> &gmtxarray,
        const unsigned int nwidthcolm=12);

    void Cout2DArrayBlkMat(
        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,
        const unsigned int nwidthcolm=12);

    void Cout1DArrayStdMat(
        Array<OneD, DNekMatSharedPtr> &gmtxarray,
        const unsigned int nwidthcolm=12);

    void Cout2DArrayStdMat(
        Array<OneD, Array<OneD, DNekMatSharedPtr> > &gmtxarray,
        const unsigned int nwidthcolm=12);

#endif

protected:
    /// Callback function to the flux vector (set when advection is in
    /// conservative form).
    AdvectionFluxVecCB     m_fluxVector;
    AdvectionFluxVecCB     m_fluxVectorMF;
    AdvectionFluxVecCBTraceNormal m_fluxVectortraceMF;
    /// Riemann solver for DG-type schemes.
    RiemannSolverSharedPtr m_riemann;
    /// Storage for space dimension. Used for homogeneous extension.
    int                    m_spaceDim;

    /// Initialises the advection object.
    SOLVER_UTILS_EXPORT virtual void v_InitObject(
        LibUtilities::SessionReaderSharedPtr              pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr>       pFields);

    /// Advects a vector field.
    SOLVER_UTILS_EXPORT virtual void v_Advect(
        const int nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &advVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
              Array<OneD, Array<OneD, NekDouble> >        &outarray,
        const NekDouble                                   &time,
        const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
        const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray)=0;
    
    /// Advects Volume Flux.
    SOLVER_UTILS_EXPORT virtual void v_AdvectVolumeFlux(
        const                                             int nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &advVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>  &pVolumeFlux,
        const NekDouble                                   &time);
    
    /// Advects Trace Flux.
    SOLVER_UTILS_EXPORT virtual void v_AdvectTraceFlux(
        const int nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &advVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble>>               &pTraceFlux,
        const NekDouble                                   &time,
        const Array<OneD, Array<OneD, NekDouble> >        &pFwd = NullNekDoubleArrayofArray,
        const Array<OneD, Array<OneD, NekDouble> >        &pBwd = NullNekDoubleArrayofArray);

    SOLVER_UTILS_EXPORT virtual void v_Advect_coeff(
        const int nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &advVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
              Array<OneD, Array<OneD, NekDouble> >        &outarray,
        const NekDouble                                   &time,
        const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
        const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);

    

    /// Overrides the base flow used during linearised advection
    SOLVER_UTILS_EXPORT virtual void v_SetBaseFlow(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields);
    
#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
    SOLVER_UTILS_EXPORT virtual void v_AddVolumJacToMat( 
        const Array<OneD, MultiRegions::ExpListSharedPtr>                       &pFields,
        const int                                                               &nConvectiveFields,
        const Array<OneD, const Array<OneD,  Array<OneD, 
            Array<OneD,  Array<OneD,  NekDouble> > > > >                        &ElmtJacArray,
        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >                          &gmtxarray);

    SOLVER_UTILS_EXPORT virtual void v_AddTraceJacToMat(
        const int                                           nConvectiveFields,
        const int                                           nSpaceDim,
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields,
        const Array<OneD, DNekBlkMatSharedPtr>              &TracePntJacCons,
        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >      &gmtxarray,
        const Array<OneD, DNekBlkMatSharedPtr>              &TracePntJacGrad        ,
        const Array<OneD, Array<OneD, NekDouble> >          &TracePntJacGradSign    );

    SOLVER_UTILS_EXPORT virtual  void v_NumCalRiemFluxJac( 
        const int                                          nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &AdvVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >        &pBwd,
        DNekBlkMatSharedPtr &FJac,
        DNekBlkMatSharedPtr &BJac);
#endif
};

/// A shared pointer to an Advection object.
typedef std::shared_ptr<Advection> AdvectionSharedPtr;

/// Datatype of the NekFactory used to instantiate classes derived
/// from the Advection class.
typedef LibUtilities::NekFactory<std::string, Advection,
std::string> AdvectionFactory;

/// Gets the factory for initialising advection objects.
SOLVER_UTILS_EXPORT AdvectionFactory& GetAdvectionFactory();

}
}

#endif
