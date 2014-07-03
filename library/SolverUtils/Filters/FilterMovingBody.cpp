///////////////////////////////////////////////////////////////////////////////
//
// File FilterAeroForces.cpp
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
// Description: Output values of aerodynamic forces during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <iomanip>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <SolverUtils/Filters/FilterMovingBody.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string FilterMovingBody::className = GetFilterFactory().RegisterCreatorFunction("MovingBody", FilterMovingBody::create);

        /**
         *
         */
        FilterMovingBody::FilterMovingBody(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::map<std::string, std::string> &pParams) :
            Filter(pSession)
        {
            if (pParams.find("OutputFile") == pParams.end())
            {
                m_outputFile_fce = m_session->GetSessionName();
                m_outputFile_mot = m_session->GetSessionName();
            }
            else
            {
                ASSERTL0(!(pParams.find("OutputFile")->second.empty()),"Missing parameter 'OutputFile'.");
                
                m_outputFile_fce = pParams.find("OutputFile")->second;
                m_outputFile_mot = pParams.find("OutputFile")->second;
            }
            if (!(m_outputFile_fce.length() >= 4 && m_outputFile_fce.substr(m_outputFile_fce.length() - 4) == ".fce"))
            {
                m_outputFile_fce += ".fce";
            }

            if (!(m_outputFile_mot.length() >= 4 && m_outputFile_mot.substr(m_outputFile_mot.length() - 4) == ".mot"))
            {
                m_outputFile_mot += ".mot";
            }
            if (pParams.find("OutputFrequency") == pParams.end())
            {
                m_outputFrequency = 1;
            }
            else
            {
                m_outputFrequency = atoi(pParams.find("OutputFrequency")->second.c_str());
            }


            m_session->MatchSolverInfo("Homogeneous", "1D", m_isHomogeneous1D, false);
            
            ASSERTL0(m_isHomogeneous1D,"Moving Body implemented just for 3D Homogeneous 1D discetisations.");

            //specify the boundary to calculate the forces
            if (pParams.find("Boundary") == pParams.end())
            {
                ASSERTL0(false, "Missing parameter 'Boundary'.");
            }
            else
            {
                ASSERTL0(!(pParams.find("Boundary")->second.empty()),
                         "Missing parameter 'Boundary'.");
                m_BoundaryString = pParams.find("Boundary")->second;
            }
        }

        /**
         *
         */
        FilterMovingBody::~FilterMovingBody()
        {

        }


        /**
         *
         */
        void FilterMovingBody::v_Initialise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time)
        {
			m_index = 0;
			m_outputStream =  Array<OneD, std::ofstream>(2);
            // Parse the boundary regions into a list.
            std::string::size_type FirstInd = m_BoundaryString.find_first_of('[') + 1;
            std::string::size_type LastInd  = m_BoundaryString.find_last_of(']') - 1;

            ASSERTL0(FirstInd <= LastInd,
                    (std::string("Error reading boundary region definition:") +
                     m_BoundaryString).c_str());

            std::string IndString = m_BoundaryString.substr(FirstInd, LastInd - FirstInd + 1);
            
            bool parseGood = ParseUtils::GenerateSeqVector(IndString.c_str(),m_boundaryRegionsIdList);
            
            ASSERTL0(parseGood && !m_boundaryRegionsIdList.empty(),
                     (std::string("Unable to read boundary regions index "
                      "range for FilterAeroForces: ") + IndString).c_str());

            // determine what boundary regions need to be considered
            int cnt;
            
            unsigned int numBoundaryRegions = pFields[0]->GetBndConditions().num_elements();
            
            m_boundaryRegionIsInList.insert(m_boundaryRegionIsInList.end(),numBoundaryRegions, 0);

            SpatialDomains::BoundaryConditions bcs(m_session,pFields[0]->GetGraph());
            
            const SpatialDomains::BoundaryRegionCollection &bregions = bcs.GetBoundaryRegions();
            
            SpatialDomains::BoundaryRegionCollection::const_iterator it;

            for (cnt = 0, it = bregions.begin(); it != bregions.end();
                    ++it, cnt++)
            {
                if ( std::find(m_boundaryRegionsIdList.begin(),
                               m_boundaryRegionsIdList.end(), it->first) !=
                        m_boundaryRegionsIdList.end() )
                {
                    m_boundaryRegionIsInList[cnt] = 1;
                }
            }

            LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

            if (vComm->GetRank() == 0)
            {
                // Open output stream for cable forces
                m_outputStream[0].open(m_outputFile_fce.c_str());
                m_outputStream[0] << "#";
                m_outputStream[0].width(7);
                m_outputStream[0] << "Time";
                m_outputStream[0].width(25);
                m_outputStream[0] << "z";
                m_outputStream[0].width(25);
                m_outputStream[0] << "Fx (press)";
                m_outputStream[0].width(25);
                m_outputStream[0] << "Fx (visc)";
                m_outputStream[0].width(25);
                m_outputStream[0] << "Fx (tot)";
                m_outputStream[0].width(25);
                m_outputStream[0] << "Fy (press)";
                m_outputStream[0].width(25);
                m_outputStream[0] << "Fy (visc)";
                m_outputStream[0].width(25);
                m_outputStream[0] << "Fy (tot)";
                m_outputStream[0] << endl;

                // Open output stream for cable motions
                m_outputStream[1].open(m_outputFile_mot.c_str());
                m_outputStream[1] << "#";
                m_outputStream[1].width(7);
                m_outputStream[1] << "Time";
                m_outputStream[1].width(25);
                m_outputStream[1] << "z";
                m_outputStream[1].width(25);
                m_outputStream[1] << "Disp_x";
                m_outputStream[1].width(25);
                m_outputStream[1] << "Vel_x";
                m_outputStream[1].width(25);
                m_outputStream[1] << "Acel_x";
                m_outputStream[1].width(25);
                m_outputStream[1] << "Disp_y";
                m_outputStream[1].width(25);
                m_outputStream[1] << "Vel_y";
                m_outputStream[1].width(25);
                m_outputStream[1] << "Acel_y";
                m_outputStream[1] << endl;
            }

            v_Update(pFields, time);
        }


        /**
         *
         */
        void FilterMovingBody::v_Update(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time)
        {
            int n, cnt, elmtid, nq, offset, nt, boundary;
            nt = pFields[0]->GetNpoints();
            int dim = pFields.num_elements()-1;

            StdRegions::StdExpansionSharedPtr elmt;
            Array<OneD, int> BoundarytoElmtID;
            Array<OneD, int> BoundarytoTraceID;
            Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;

            Array<OneD, const NekDouble> P(nt);
            Array<OneD, const NekDouble> U(nt);
            Array<OneD, const NekDouble> V(nt);
            Array<OneD, const NekDouble> W(nt);

            Array<OneD, Array<OneD, NekDouble> > gradU(dim);
            Array<OneD, Array<OneD, NekDouble> > gradV(dim);
            Array<OneD, Array<OneD, NekDouble> > gradW(dim);

            Array<OneD, Array<OneD, NekDouble> > fgradU(dim);
            Array<OneD, Array<OneD, NekDouble> > fgradV(dim);
            Array<OneD, Array<OneD, NekDouble> > fgradW(dim);
            
            LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
            
            // set up storage space for forces on all the planes (z-positions)
            // on each processors some of them will remain empty as we may
            // have just few planes per processor
            int Num_z_pos = pFields[0]->GetHomogeneousBasis()->GetNumModes();
            Array<OneD, NekDouble> Fx(Num_z_pos,0.0);
            Array<OneD, NekDouble> Fxp(Num_z_pos,0.0);
            Array<OneD, NekDouble> Fxv(Num_z_pos,0.0);
            Array<OneD, NekDouble> Fy(Num_z_pos,0.0);
            Array<OneD, NekDouble> Fyp(Num_z_pos,0.0);
            Array<OneD, NekDouble> Fyv(Num_z_pos,0.0);
            

            NekDouble rho = (m_session->DefinesParameter("rho"))
                    ? (m_session->GetParameter("rho"))
                    : 1;
            NekDouble mu = rho*m_session->GetParameter("Kinvis");
            
            for(int i = 0; i < pFields.num_elements(); ++i)
            {
                pFields[i]->SetWaveSpace(false);
                pFields[i]->BwdTrans(pFields[i]->GetCoeffs(),
                                     pFields[i]->UpdatePhys());
                pFields[i]->SetPhysState(true);
            }

            // Get the number of local planes on the process and their IDs 
            // to properly locate the forces in the Fx, Fy etc. vectors. 
            Array<OneD, unsigned int> ZIDs;
            ZIDs = pFields[0]->GetZIDs();
            int local_planes = ZIDs.num_elements();
            
            // Homogeneous 1D case  Compute forces on all WALL boundaries
            // This only has to be done on the zero (mean) Fourier mode.
            for(int plane = 0 ; plane < local_planes; plane++)
            {
                pFields[0]->GetPlane(plane)->GetBoundaryToElmtMap(BoundarytoElmtID,BoundarytoTraceID);
                BndExp = pFields[0]->GetPlane(plane)->GetBndCondExpansions();
                StdRegions::StdExpansion1DSharedPtr bc;
                
                // loop over the types of boundary conditions
                for(cnt = n = 0; n < BndExp.num_elements(); ++n)
                {
                    if(m_boundaryRegionIsInList[n] == 1)
                    {
                        for(int i = 0; i <  BndExp[n]->GetExpSize(); ++i, cnt++)
                        {
                            // find element of this expansion.
                            elmtid = BoundarytoElmtID[cnt];
                            elmt   = pFields[0]->GetPlane(plane)->GetExp(elmtid);
                            nq     = elmt->GetTotPoints();
                            offset = pFields[0]->GetPlane(plane)->GetPhys_Offset(elmtid);
                            
                            // Initialise local arrays for the velocity
                            // gradients size of total number of quadrature
                            // points for each element (hence local).
                            for(int j = 0; j < dim; ++j)
                            {
                                gradU[j] = Array<OneD, NekDouble>(nq,0.0);
                                gradV[j] = Array<OneD, NekDouble>(nq,0.0);
                                gradW[j] = Array<OneD, NekDouble>(nq,0.0);
                            }

                            // identify boundary of element
                            boundary = BoundarytoTraceID[cnt];

                            // Extract  fields
                            U = pFields[0]->GetPlane(plane)->GetPhys() + offset;
                            V = pFields[1]->GetPlane(plane)->GetPhys() + offset;
                            P = pFields[3]->GetPlane(plane)->GetPhys() + offset;

                            // compute the gradients
                            elmt->PhysDeriv(U,gradU[0],gradU[1]);
                            elmt->PhysDeriv(V,gradV[0],gradV[1]);

                            // Get face 1D expansion from element expansion
                            bc =  boost::dynamic_pointer_cast<LocalRegions::Expansion1D> (BndExp[n]->GetExp(i));

                            // number of points on the boundary
                            int nbc = bc->GetTotPoints();

                            // several vectors for computing the forces
                            Array<OneD, NekDouble> Pb(nbc,0.0);

                            for(int j = 0; j < dim; ++j)
                            {
                                fgradU[j] = Array<OneD, NekDouble>(nbc,0.0);
                                fgradV[j] = Array<OneD, NekDouble>(nbc,0.0);
                            }

                            Array<OneD, NekDouble>  drag_t(nbc,0.0);
                            Array<OneD, NekDouble>  lift_t(nbc,0.0);
                            Array<OneD, NekDouble>  drag_p(nbc,0.0);
                            Array<OneD, NekDouble>  lift_p(nbc,0.0);
                            Array<OneD, NekDouble>  temp(nbc,0.0);
                            Array<OneD, NekDouble>  temp2(nbc,0.0);

                            // identify boundary of element .
                            boundary = BoundarytoTraceID[cnt];

                            // extraction of the pressure and wss on the
                            // boundary of the element
                            elmt->GetEdgePhysVals(boundary,bc,P,Pb);

                            for(int j = 0; j < dim; ++j)
                            {
                                elmt->GetEdgePhysVals(boundary,bc,gradU[j],fgradU[j]);
                                elmt->GetEdgePhysVals(boundary,bc,gradV[j],fgradV[j]);
                            }

                            //normals of the element
                            const Array<OneD, Array<OneD, NekDouble> > &normals= elmt->GetEdgeNormal(boundary);

                            //
                            // Compute viscous tractive forces on wall from
                            //
                            //  t_i  = - T_ij * n_j  (minus sign for force 
                            //                        exerted BY fluid ON wall),
                            //
                            // where
                            //
                            //  T_ij = viscous stress tensor (here in Cartesian
                            //         coords)
                            //                          dU_i    dU_j
                            //       = RHO * KINVIS * ( ----  + ---- ) .
                            //                          dx_j    dx_i

                            //a) DRAG TERMS
                            //-rho*kinvis*(2*du/dx*nx+(du/dy+dv/dx)*ny)

                            Vmath::Vadd(nbc,fgradU[1],1,fgradV[0],1,drag_t,1);
                            Vmath::Vmul(nbc,drag_t,1,normals[1],1,drag_t,1);

                            Vmath::Smul(nbc,2.0,fgradU[0],1,fgradU[0],1);
                            Vmath::Vmul(nbc,fgradU[0],1,normals[0],1,temp2,1);
                            Vmath::Smul(nbc,0.5,fgradU[0],1,fgradU[0],1);

                            Vmath::Vadd(nbc,temp2,1,drag_t,1,drag_t,1);
                            Vmath::Smul(nbc,-mu,drag_t,1,drag_t,1);

                            //zero temporary storage vector
                            Vmath::Zero(nbc,temp,0);
                            Vmath::Zero(nbc,temp2,0);


                            //b) LIFT TERMS
                            //-rho*kinvis*(2*dv/dy*ny+(du/dy+dv/dx)*nx)

                            Vmath::Vadd(nbc,fgradU[1],1,fgradV[0],1,lift_t,1);
                            Vmath::Vmul(nbc,lift_t,1,normals[0],1,lift_t,1);

                            Vmath::Smul(nbc,2.0,fgradV[1],1,fgradV[1],1);
                            Vmath::Vmul(nbc,fgradV[1],1,normals[1],1,temp2,1);
                            Vmath::Smul(nbc,-0.5,fgradV[1],1,fgradV[1],1);


                            Vmath::Vadd(nbc,temp2,1,lift_t,1,lift_t,1);
                            Vmath::Smul(nbc,-mu,lift_t,1,lift_t,1);

                            // Compute normal tractive forces on all WALL
                            // boundaries
  
                            Vmath::Vvtvp(nbc,Pb,1,normals[0],1,drag_p,1,drag_p, 1);
                            Vmath::Vvtvp(nbc,Pb,1,normals[1],1,lift_p,1,lift_p,1);

                            //integration over the boundary
                            Fxv[ZIDs[plane]] += bc->Integral(drag_t);
                            Fyv[ZIDs[plane]] += bc->Integral(lift_t);

                            Fxp[ZIDs[plane]] += bc->Integral(drag_p);
                            Fyp[ZIDs[plane]] += bc->Integral(lift_p);
                        }
                    }
                    else
                    {
                        cnt += BndExp[n]->GetExpSize();
                    }
                }
            }
            
            for(int i = 0; i < pFields.num_elements(); ++i)
            {
                pFields[i]->SetWaveSpace(true);
                pFields[i]->BwdTrans(pFields[i]->GetCoeffs(),
                                     pFields[i]->UpdatePhys());
                pFields[i]->SetPhysState(false);
            }
            
            // In case we are using an hybrid parallelisation we need
            // to reduce the forces on the same plane coming from
            // different mesh partitions.
            // It is quite an expensive communication, therefore
            // we check to make sure it is actually required.
            if(vComm->GetRowComm()->GetSize() > 0)
            {
                // NOTE 1: We can eventually sum the viscous and pressure
                // component before doing the communication, thus
                // reducing by a factor 2 the communication.
                // NOTE 2: We may want to set up in the Comm class an AllReduce
                // routine wich can handle arrays more efficiently
                for(int plane = 0 ; plane < local_planes; plane++)
                {
                    vComm->GetRowComm()->AllReduce(Fxp[ZIDs[plane]], LibUtilities::ReduceSum);
                    vComm->GetRowComm()->AllReduce(Fxv[ZIDs[plane]], LibUtilities::ReduceSum);
                    vComm->GetRowComm()->AllReduce(Fyp[ZIDs[plane]], LibUtilities::ReduceSum);
                    vComm->GetRowComm()->AllReduce(Fyv[ZIDs[plane]], LibUtilities::ReduceSum);
                }
            }
            
            // At this point on rank (0,n) of the Mesh partion communicator we have
            // the total areo forces on the planes which are on the same column communicator.
            // Since the planes are scattered on different processors some of the entries
            // of the vector Fxp, Fxp etc. are still zero.
            // Now we need to reduce the values on a single vector on rank (0,0) of the
            // global communicator.
            if(vComm->GetRowComm()->GetRank() == 0)
            {
                for(int z = 0 ; z < Num_z_pos; z++)
                {
                    vComm->GetColumnComm()->AllReduce(Fxp[z], LibUtilities::ReduceSum);
                    vComm->GetColumnComm()->AllReduce(Fxv[z], LibUtilities::ReduceSum);
                    vComm->GetColumnComm()->AllReduce(Fyp[z], LibUtilities::ReduceSum);
                    vComm->GetColumnComm()->AllReduce(Fyv[z], LibUtilities::ReduceSum);
                }
            }

            //set the forces imparted on the cable's wall
            Array<OneD, NekDouble> fces(2);
            for(int plane = 0 ; plane < local_planes; plane++)
            {
                fces[0] = Fxp[ZIDs[plane]] + Fxv[ZIDs[plane]];
                fces[1] = Fyp[ZIDs[plane]] + Fyv[ZIDs[plane]];
                pFields[0]->GetPlane(plane)->SetMovBodyForces(fces[0],fces[1]);
            }
           
			//get and output moving body variables
			Array<OneD, Array<OneD, NekDouble> > Motion_x;
			Array<OneD, Array<OneD, NekDouble> > Motion_y;
			int nStrVars = 3;
        	Motion_x = Array<OneD, Array<OneD, NekDouble> > (nStrVars);
        	Motion_y = Array<OneD, Array<OneD, NekDouble> > (nStrVars);
        	for(int i = 0; i < Motion_x.num_elements(); i++)
        	{
            	Motion_x[i] = Array<OneD, NekDouble>(local_planes,0.0);
            	Motion_y[i] = Array<OneD, NekDouble>(local_planes,0.0);
        	}
        	for(int plane = 0; plane < local_planes; plane++)
        	{
            	Array<OneD, NekDouble> tmp0(nStrVars);
            	Array<OneD, NekDouble> tmp1(nStrVars);

            	pFields[0]->GetPlane(plane)->GetMovBodyMotionVars(tmp0);
            	pFields[1]->GetPlane(plane)->GetMovBodyMotionVars(tmp1);
                for (int var = 0; var < nStrVars; var++)
                {
                    Motion_x[var][plane] = tmp0[var];
                    Motion_y[var][plane] = tmp1[var];
                }
        	}
            // Only output every m_outputFrequency.
            if ((m_index++) % m_outputFrequency)
            {
                return;
            }
 
            // At thi point in rank (0,0) we have the full vectors
            // containing Fxp,Fxv,Fyp and Fyv where different positions
            // in the vectors correspond to different planes.
            // Here we write it to file. We do it just on one porcess

            Array<OneD, NekDouble> z_coords(Num_z_pos,0.0);
            Array<OneD, const NekDouble> pts = pFields[0]->GetHomogeneousBasis()->GetZ();

            NekDouble LZ;
            m_session->LoadParameter("LZ", LZ);
            Vmath::Smul(Num_z_pos,LZ/2.0,pts,1,z_coords,1);
            Vmath::Sadd(Num_z_pos,LZ/2.0,z_coords,1,z_coords,1);
            if (vComm->GetRank() == 0)
            {
                
                Vmath::Vadd(Num_z_pos,Fxp,1,Fxv,1,Fx,1);
                Vmath::Vadd(Num_z_pos,Fyp,1,Fyv,1,Fy,1);
                
                for(int i = 0 ; i < Num_z_pos; i++)
                {
                    m_outputStream[0].width(8);
                    m_outputStream[0] << setprecision(6) << time;
                    
                    m_outputStream[0].width(25);
                    m_outputStream[0] << setprecision(6) << z_coords[i];
                    
                    m_outputStream[0].width(25);
                    m_outputStream[0] << setprecision(8) << Fxp[i];
                    m_outputStream[0].width(25);
                    m_outputStream[0] << setprecision(8) << Fxv[i];
                    m_outputStream[0].width(25);
                    m_outputStream[0] << setprecision(16) << Fx[i];
                    
                    m_outputStream[0].width(25);
                    m_outputStream[0] << setprecision(8) << Fyp[i];
                    m_outputStream[0].width(25);
                    m_outputStream[0] << setprecision(8) << Fyv[i];
                    m_outputStream[0].width(25);
                    m_outputStream[0] << setprecision(8) << Fy[i];
                    m_outputStream[0] << endl;
                }
            }

			//
        	Array <OneD, NekDouble> CableAccelX;
        	Array <OneD, NekDouble> CableVelocX;
        	Array <OneD, NekDouble> CableDisplX;
        	Array <OneD, NekDouble> CableAccelY;
        	Array <OneD, NekDouble> CableVelocY;
        	Array <OneD, NekDouble> CableDisplY;

        	int npoints = Motion_x[0].num_elements();
        	CableAccelX = Array <OneD, NekDouble>(npoints);
        	CableVelocX = Array <OneD, NekDouble>(npoints);
        	CableDisplX = Array <OneD, NekDouble>(npoints);
        	CableAccelY = Array <OneD, NekDouble>(npoints);
        	CableVelocY = Array <OneD, NekDouble>(npoints);
        	CableDisplY = Array <OneD, NekDouble>(npoints);

        	Vmath::Vcopy(npoints,Motion_x[0],1,CableDisplX,1);
        	Vmath::Vcopy(npoints,Motion_x[1],1,CableVelocX,1);
        	Vmath::Vcopy(npoints,Motion_x[2],1,CableAccelX,1);
        	Vmath::Vcopy(npoints,Motion_y[0],1,CableDisplY,1);
        	Vmath::Vcopy(npoints,Motion_y[1],1,CableVelocY,1);
        	Vmath::Vcopy(npoints,Motion_y[2],1,CableAccelY,1);

        	int colrank = vComm->GetColumnComm()->GetRank();
        	int nproc   = vComm->GetColumnComm()->GetSize();
        	// Send to root process.
        	if (colrank == 0)
        	{
            	for (int j = 0; j <Motion_x[0].num_elements(); j++)
            	{
                	m_outputStream[1].width(8);
                	m_outputStream[1] << setprecision(6) << time;
                	m_outputStream[1].width(25);
                	m_outputStream[1] << setprecision(6) << z_coords[j];
                	m_outputStream[1].width(25);
                	m_outputStream[1] << setprecision(8) << CableDisplX[j];
                	m_outputStream[1].width(25);
                	m_outputStream[1] << setprecision(8) << CableVelocX[j];
                	m_outputStream[1].width(25);
                	m_outputStream[1] << setprecision(8) << CableAccelX[j];
                	m_outputStream[1].width(25);
                	m_outputStream[1] << setprecision(8) << CableDisplY[j];
                	m_outputStream[1].width(25);
                	m_outputStream[1] << setprecision(8) << CableVelocY[j];
                	m_outputStream[1].width(25);
                	m_outputStream[1] << setprecision(8) << CableAccelY[j];
                	m_outputStream[1] << endl;
            	}

            	for (int i = 1; i < nproc; ++i)
            	{
                 	vComm->GetColumnComm()->Recv(i, CableAccelX);
                 	vComm->GetColumnComm()->Recv(i, CableVelocX);
                 	vComm->GetColumnComm()->Recv(i, CableDisplX);
                 	vComm->GetColumnComm()->Recv(i, CableAccelY);
                 	vComm->GetColumnComm()->Recv(i, CableVelocY);
                 	vComm->GetColumnComm()->Recv(i, CableDisplY);

                	for (int j = 0; j < Motion_x[0].num_elements(); ++j)
                	{
                    	m_outputStream[1].width(8);
                    	m_outputStream[1] << setprecision(6) << time;
                    	m_outputStream[1].width(25);
                    	m_outputStream[1] << setprecision(6) << z_coords[2 * i + j];
                    	m_outputStream[1].width(25);
                    	m_outputStream[1] << setprecision(8) << CableDisplX[j];
                    	m_outputStream[1].width(25);
                    	m_outputStream[1] << setprecision(8) << CableVelocX[j];
                    	m_outputStream[1].width(25);
                    	m_outputStream[1] << setprecision(8) << CableAccelX[j];
                    	m_outputStream[1].width(25);
                    	m_outputStream[1] << setprecision(8) << CableDisplY[j];
                    	m_outputStream[1].width(25);
                    	m_outputStream[1] << setprecision(8) << CableVelocY[j];
                    	m_outputStream[1].width(25);
                    	m_outputStream[1] << setprecision(8) << CableAccelY[j];
                    	m_outputStream[1] << endl;
                	}
            	}
        	}
        	else
        	{
            	vComm->GetColumnComm()->Send(0, CableAccelX);
            	vComm->GetColumnComm()->Send(0, CableVelocX);
            	vComm->GetColumnComm()->Send(0, CableDisplX);
            	vComm->GetColumnComm()->Send(0, CableAccelY);
            	vComm->GetColumnComm()->Send(0, CableVelocY);
            	vComm->GetColumnComm()->Send(0, CableDisplY);
        	}
        }


        /**
         *
         */
        void FilterMovingBody::v_Finalise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, 
            const NekDouble &time)
        {
            if (pFields[0]->GetComm()->GetRank() == 0)
            {
                m_outputStream[0].close();
                m_outputStream[1].close();
            }
        }


        /**
         *
         */
        bool FilterMovingBody::v_IsTimeDependent()
        {
            return true;
        }
    }
}
