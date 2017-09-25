///////////////////////////////////////////////////////////////////////////////
//
// File ExecuteSharpy.h
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
// Description: *****.
//
///////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_EXECUTESHARPY_H
#define NEKTAR_EXECUTESHARPY_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <LibUtilities/SHARPy/SHARPy.hpp>


namespace Nektar
{
namespace Sharpy
{
    class Beam;


    class Beam
    {
        public:

            Beam(){};

            ~Beam(){};

            void InitialiseStatic(const LibUtilities::SessionReaderSharedPtr& pSession);

            void InitialiseDynamic(const LibUtilities::SessionReaderSharedPtr& pSession,
                    const Array<OneD, Array<OneD, NekDouble> >& HydCoeffs);

            void SolvenlnStatic(const LibUtilities::SessionReaderSharedPtr& pSession,
                    const Array<OneD, Array<OneD, NekDouble> > &HydroForces);

            void SolvenlnDynamic(const LibUtilities::SessionReaderSharedPtr& pSession,
                    const Array<OneD, Array<OneD, NekDouble> > &HydroForces,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &motions);


        private:

            void Input_elem(const LibUtilities::SessionReaderSharedPtr& pSession);

            void Input_node(const LibUtilities::SessionReaderSharedPtr& pSession);

            void AddGravityLoad(const LibUtilities::SessionReaderSharedPtr& pSession);

            void Skew(const Array<OneD, NekDouble> &Vector,
                    Array<OneD, Array<OneD, NekDouble> > &SkewMat);

            void RotCRV(const Array<OneD, NekDouble> &Psi,
                    Array<OneD, Array<OneD, NekDouble> > &C);

            ///SHARPy Variables
            bool m_FollowerForce;        // Forces follow local deflections
            bool m_FollowerForceRig;     // Forces follow the body-fixed frame
            bool m_PrintInfo;            // ???
            bool m_OutInBframe;          // Print velocities in local-frame (if not, use body-attached frame)
            bool m_OutInaframe;          // Print velocities in body-fixed-frame (if not, use inertial frame)

            int m_MaxElNod;              // Max number of nodes per element.
            int m_NumElems;              // Elements discreising the beam
            int m_ElemProj;              // Element info computed in (1) global frame (2) fixed element frame (3) moving element frame.
            int m_MaxIterations;         // Newton-Raphson iterations for nonlinear solution
            int m_NumLoadSteps;          // Number of load increments
            int m_NumGauss;              // Gauss points in element integration
            int m_Solution;              // 902/912: cbeam3 linear/nonlinear flexible-body dynamic
            int m_NumNodesElem;
            int m_NumNodes_tot;
            int m_NumDof;
            int m_DimMat;
            int m_ms;
            int m_cs;
            int m_ks;
            int m_fs;
            int m_as;
            int m_NumDof2;
            int m_NumSteps;
            int m_iStep;
            int m_NForcedDisp;

            NekDouble m_t0;
            NekDouble m_tf;
            NekDouble m_dt;
            NekDouble m_gamma;
            NekDouble m_beta;
            NekDouble m_DeltaCurved;
            NekDouble m_MinDelta;
            NekDouble m_NewmarkDamp;
            NekDouble m_BeamLength;
            NekDouble m_D;

            NekDouble m_U_inf;
            NekDouble m_rho_f;

            Array<OneD, int> m_NumNodes;
            Array<OneD, int> m_MemNo;
            Array<OneD, int> m_Conn;
            Array<OneD, int> m_Master_Array;
            Array<OneD, int> m_BoundConds;
            Array<OneD, int> m_Nod_Master;
            Array<OneD, int> m_Nod_Vdof;
            Array<OneD, int> m_Nod_Fdof;
            Array<OneD, int> m_Nod_Sflag;
            Array<OneD, int> m_Master;
            Array<OneD, int> m_Vdof;
            Array<OneD, int> m_Fdof;
            Array<OneD, int> m_ListIN;
            Array<OneD, int> m_NodeForcedDisp_Array;

            Array<OneD, char>m_OutFile;


            Array<OneD, NekDouble> m_Length;
            Array<OneD, NekDouble> m_PreCurv;
            Array<OneD, NekDouble> m_Psi;
            Array<OneD, NekDouble> m_Vector;
            Array<OneD, NekDouble> m_Mass_Array;
            Array<OneD, NekDouble> m_Stiff_Array;
            Array<OneD, NekDouble> m_InvStiff_Array;
            Array<OneD, NekDouble> m_RBMass_Array;
            Array<OneD, NekDouble> m_PosIni_Array;
            Array<OneD, NekDouble> m_ForceStatic_Array;
            Array<OneD, NekDouble> m_PhiNodes;
            Array<OneD, NekDouble> m_PsiIni_Array;
            Array<OneD, NekDouble> m_AppForces_Array;
            Array<OneD, NekDouble> m_PosDefor_Array;
            Array<OneD, NekDouble> m_PsiDefor_Array;
            Array<OneD, NekDouble> m_PosDeforDot_Array;
            Array<OneD, NekDouble> m_PsiDeforDot_Array;
            Array<OneD, NekDouble> m_X;
            Array<OneD, NekDouble> m_DX;
            Array<OneD, NekDouble> m_DXdt;
            Array<OneD, NekDouble> m_DXddt;
            Array<OneD, NekDouble> m_Force_Array;
            Array<OneD, NekDouble> m_Vrel;
            Array<OneD, NekDouble> m_VrelDot;
            Array<OneD, NekDouble> m_Mglobal_Array;
            Array<OneD, NekDouble> m_Mvel_Array;
            Array<OneD, NekDouble> m_Cglobal_Array;
            Array<OneD, NekDouble> m_Cvel_Array;
            Array<OneD, NekDouble> m_Kglobal_Array;
            Array<OneD, NekDouble> m_Fglobal_Array;
            Array<OneD, NekDouble> m_Qglobal;
            Array<OneD, NekDouble> m_PosDeforDDot_Array;
            Array<OneD, NekDouble> m_PsiDeforDDot_Array;
            Array<OneD, NekDouble> m_Quat_Array;
            Array<OneD, NekDouble> m_Cao_Array;
            Array<OneD, NekDouble> m_Asys_Array;
            Array<OneD, NekDouble> m_F0;
            Array<OneD, NekDouble> m_Fa;
            Array<OneD, NekDouble> m_FTime;
            Array<OneD, NekDouble> m_Time;
            Array<OneD, NekDouble> m_ForcedVel;
            Array<OneD, NekDouble> m_ForcedVelDot;
            Array<OneD, NekDouble> m_PosForcedDisp_Array;
            Array<OneD, NekDouble> m_PsiA_G_Array;
    };
}
}

#endif

