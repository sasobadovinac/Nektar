///////////////////////////////////////////////////////////////////////////////
//
// File: Advection.cpp
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

#include <SolverUtils/Advection/Advection.h>

namespace Nektar
{
namespace SolverUtils
{

/**
 * @returns The advection factory.
 */
AdvectionFactory& GetAdvectionFactory()
{
    static AdvectionFactory instance;
    return instance;
}


/**
 * @param   pSession            Session configuration data.
 * @param   pFields             Array of ExpList objects.
 */
void Advection::InitObject(
    const LibUtilities::SessionReaderSharedPtr        pSession,
    Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
{
    v_InitObject(pSession, pFields);
}


/**
 * @param   nConvectiveFields   Number of velocity components.
 * @param   pFields             Expansion lists for scalar fields.
 * @param   pAdvVel             Advection velocity.
 * @param   pInarray            Scalar data to advect.
 * @param   pOutarray           Advected scalar data.
 * @param   pTime               Simulation time.
 */
void Advection::Advect(
    const int                                          nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble> >        &pAdvVel,
    const Array<OneD, Array<OneD, NekDouble> >        &pInarray,
    Array<OneD, Array<OneD, NekDouble> >              &pOutarray,
    const NekDouble                                   &pTime,
    const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
    const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
{
    v_Advect(nConvectiveFields, pFields, pAdvVel, pInarray,
            pOutarray, pTime, pFwd, pBwd);
}

// Check if the function is supported
// To notice, the const pFwd and pBwd are not initialized, need to be
// initialized in children class
void Advection::v_AdvectVolumeFlux(
    const int nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble>>         &pAdvVel,
    const Array<OneD, Array<OneD, NekDouble>>         &pInarray,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>>  &pVolumeFlux,
    const NekDouble                                   &pTime)
{
    ASSERTL0(false, "Not defined for AdvectVolumeFlux.");
}

// Check if the function is supported
// To notice, the const pFwd and pBwd are not initialized, need to be
// initialized in children class
void Advection::v_AdvectTraceFlux(
    const int nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble>>         &pAdvVel,
    const Array<OneD, Array<OneD, NekDouble>>         &pInarray,
    Array<OneD, Array<OneD, NekDouble>>               &pTraceFlux,
    const NekDouble                                   &pTime,
    const Array<OneD, Array<OneD, NekDouble>>         &pFwd,
    const Array<OneD, Array<OneD, NekDouble>>         &pBwd)
{
    ASSERTL0(false, "Not defined for AdvectTraceFlux.");
}

/**
 * @param   nConvectiveFields   Number of velocity components.
 * @param   pFields             Expansion lists for scalar fields.
 * @param   pAdvVel             Advection velocity.
 * @param   pInarray            Scalar data to advect.
 * @param   pOutarray           Advected scalar data.
 * @param   pTime               Simulation time.
 */
void Advection::Advect_coeff(
    const int                                          nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble> >        &pAdvVel,
    const Array<OneD, Array<OneD, NekDouble> >        &pInarray,
    Array<OneD, Array<OneD, NekDouble> >              &pOutarray,
    const NekDouble                                   &pTime,
    const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
    const Array<OneD, Array<OneD, NekDouble> >        &pBwd,
    const bool                                       flagFreezeJac)
{
    m_flagFreezeJac  = flagFreezeJac;
    v_Advect_coeff(nConvectiveFields, pFields, pAdvVel, pInarray,
            pOutarray, pTime, pFwd, pBwd);
    m_flagFreezeJac  = false;
}

/**
 * This function should be overridden in derived classes to initialise the
 * specific advection data members. However, this base class function should
 * be called as the first statement of the overridden function to ensure the
 * base class is correctly initialised in order.
 *
 * @param   pSession            Session information.
 * @param   pFields             Expansion lists for scalar fields.
 */
void Advection::v_InitObject(
    const LibUtilities::SessionReaderSharedPtr        pSession,
    Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
{
    m_spaceDim = pFields[0]->GetCoordim(0);

    if (pSession->DefinesSolverInfo("HOMOGENEOUS"))
    {
        std::string HomoStr = pSession->GetSolverInfo("HOMOGENEOUS");
        if (HomoStr == "HOMOGENEOUS1D" || HomoStr == "Homogeneous1D" ||
            HomoStr == "1D"            || HomoStr == "Homo1D" ||
            HomoStr == "HOMOGENEOUS2D" || HomoStr == "Homogeneous2D" ||
            HomoStr == "2D"            || HomoStr == "Homo2D")
        {
            m_spaceDim++;
        }
        else
        {
            ASSERTL0(false, "Only 1D homogeneous dimension supported.");
        }
    }
}

/**
 *
 */
void Advection::v_SetBaseFlow(
        const Array<OneD, Array<OneD, NekDouble> >    &inarray,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
{
    ASSERTL0(false,
            "A baseflow is not appropriate for this advection type.");
}

void Advection::v_Advect_coeff(
        const int nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &advVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
              Array<OneD, Array<OneD, NekDouble> >        &outarray,
        const NekDouble                                   &time,
        const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
        {
            ASSERTL0(false, "v_Advect_coeff no defined");
        }

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
    void Advection::v_NumCalRiemFluxJac( 
        const int                                          nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &AdvVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >        &pBwd,
        DNekBlkMatSharedPtr &FJac,
        DNekBlkMatSharedPtr &BJac)
        {
            ASSERTL0(false, "v_NumCalRiemFluxJac no defined");
        }

    void Advection::v_AddTraceJacToMat(
        const int                                           nConvectiveFields,
        const int                                           nSpaceDim,
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields,
        const Array<OneD, DNekBlkMatSharedPtr>              &TracePntJacCons,
        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >      &gmtxarray,
        const Array<OneD, DNekBlkMatSharedPtr>              &TracePntJacGrad        ,
        const Array<OneD, Array<OneD, NekDouble> >          &TracePntJacGradSign    )
    {
        ASSERTL0(false,"v_AddTraceJacToMat NOT SPECIFIED");
        return;
    }
  
    void Advection::v_AddVolumJacToMat( 
        const Array<OneD, MultiRegions::ExpListSharedPtr>                       &pFields,
        const int                                                               &nConvectiveFields,
        const Array<OneD, const Array<OneD,  Array<OneD, 
              Array<OneD,  Array<OneD,  NekDouble> > > > >                      &ElmtJacArray,
        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >                          &gmtxarray)
    {
        ASSERTL0(false,"v_AddVolumJacToMat NOT SPECIFIED");
        return;
    }

    void Advection::CalcTraceJac(
        const int                                          nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &AdvVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >        &pBwd,
        DNekBlkMatSharedPtr &FJac,
        DNekBlkMatSharedPtr &BJac)
    {
        // int nPointsTot      = fields[0]->GetTotPoints();
        // int nCoeffs         = fields[0]->GetNcoeffs();
        int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();
        
        // Store forwards/backwards space along trace space
        Array<OneD, Array<OneD, NekDouble> > Fwd    (nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble> > Bwd    (nConvectiveFields);

        if (pFwd == NullNekDoubleArrayofArray ||
            pBwd == NullNekDoubleArrayofArray)
        {
            for(int i = 0; i < nConvectiveFields; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
        }
        else
        {
            for(int i = 0; i < nConvectiveFields; ++i)
            {
                Fwd[i]     = pFwd[i];
                Bwd[i]     = pBwd[i];
            }
        }

        ASSERTL1(m_riemann,
                    "Riemann solver must be provided for AdvectionWeakDG.");
    
        m_riemann->CalcFluxJacobian(m_spaceDim, Fwd, Bwd, FJac,BJac);
    }


    void Advection::Cout1DArrayBlkMat(Array<OneD, DNekBlkMatSharedPtr> &gmtxarray,const unsigned int nwidthcolm)
    {
        int nvar1 = gmtxarray.num_elements();

        
        for(int i = 0; i < nvar1; i++)
        {
            cout<<endl<<"£$£$£$£$£$£$££$£$£$$$£$££$$£$££$£$$££££$$£$£$£$£$£$£$££$£$$"<<endl<< "Cout2DArrayBlkMat i= "<<i<<endl;
            CoutBlkMat(gmtxarray[i],nwidthcolm);
        }
    }
    
    
    void Advection::Cout2DArrayBlkMat(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,const unsigned int nwidthcolm)
    {
        int nvar1 = gmtxarray.num_elements();
        int nvar2 = gmtxarray[0].num_elements();

        
        for(int i = 0; i < nvar1; i++)
        {
            for(int j = 0; j < nvar2; j++)
            {
                cout<<endl<<"£$£$£$£$£$£$££$£$£$$$£$££$$£$££$£$$££££$$£$£$£$£$£$£$££$£$$"<<endl<< "Cout2DArrayBlkMat i= "<<i<<" j=  "<<j<<endl;
                CoutBlkMat(gmtxarray[i][j],nwidthcolm);
            }
        }
    }
    
    void Advection::CoutBlkMat(DNekBlkMatSharedPtr &gmtx,const unsigned int nwidthcolm)
    {
        DNekMatSharedPtr    loc_matNvar;

        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;
        gmtx->GetBlockSizes(rowSizes,colSizes);

        int nelmts  = rowSizes.num_elements();
        
        // int noffset = 0;
        for(int i = 0; i < nelmts; ++i)
        {
            loc_matNvar =   gmtx->GetBlock(i,i);
            std::cout   <<std::endl<<"*********************************"<<std::endl<<"element :   "<<i<<std::endl;
            CoutStandardMat(loc_matNvar,nwidthcolm);
        }
        return;
    }

    void Advection::Cout1DArrayStdMat(Array<OneD, DNekMatSharedPtr> &gmtxarray,const unsigned int nwidthcolm)
    {
        int nvar1 = gmtxarray.num_elements();

        
        for(int i = 0; i < nvar1; i++)
        {
            cout<<endl<<"£$£$£$£$£$£$££$£$£$$$£$££$$£$££$£$$££££$$£$£$£$£$£$£$££$£$$"<<endl<< "Cout2DArrayBlkMat i= "<<i<<endl;
            CoutStandardMat(gmtxarray[i],nwidthcolm);
        }
    }

    void Advection::Cout2DArrayStdMat(Array<OneD, Array<OneD, DNekMatSharedPtr> > &gmtxarray,const unsigned int nwidthcolm)
    {
        int nvar1 = gmtxarray.num_elements();
        int nvar2 = gmtxarray[0].num_elements();

        
        for(int i = 0; i < nvar1; i++)
        {
            for(int j = 0; j < nvar2; j++)
            {
                cout<<endl<<"£$£$£$£$£$£$££$£$£$$$£$££$$£$££$£$$££££$$£$£$£$£$£$£$££$£$$"<<endl<< "Cout2DArrayBlkMat i= "<<i<<" j=  "<<j<<endl;
                CoutStandardMat(gmtxarray[i][j],nwidthcolm);
            }
        }
    }

    void Advection::CoutStandardMat(DNekMatSharedPtr &loc_matNvar,const unsigned int nwidthcolm)
    {
        int nrows = loc_matNvar->GetRows();
        int ncols = loc_matNvar->GetColumns();
        NekDouble tmp=0.0;
        std::cout   <<"ROW="<<std::setw(3)<<-1<<" ";
        for(int k = 0; k < ncols; k++)
        {
            std::cout   <<"   COL="<<std::setw(nwidthcolm-7)<<k;
        }
        std::cout   << endl;

        for(int j = 0; j < nrows; j++)
        {
            std::cout   <<"ROW="<<std::setw(3)<<j<<" ";
            for(int k = 0; k < ncols; k++)
            {
                tmp =   (*loc_matNvar)(j,k);
                std::cout   <<std::scientific<<std::setw(nwidthcolm)<<std::setprecision(nwidthcolm-8)<<tmp;
            }
            std::cout   << endl;
        }
    }

#endif

}
}
