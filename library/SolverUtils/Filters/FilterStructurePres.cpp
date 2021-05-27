///////////////////////////////////////////////////////////////////////////////
//
// File FilterStructurePres.cpp
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
// Description: Output Structure preserved coeffs.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>
//#include <LibUtilities/Polylib/Polylib.h>
#include <SolverUtils/Filters/FilterStructurePres.h>

//#include <SolverUtils/Filters/FilterInterfaces.hpp>
//#include <LibUtilities/Memory/NekMemoryManager.hpp>


namespace Nektar
{
  namespace SolverUtils
  {
    
    std::string FilterStructurePres::className = GetFilterFactory().
      RegisterCreatorFunction("StructurePres", FilterStructurePres::create);

    FilterStructurePres::FilterStructurePres(
					     const LibUtilities::SessionReaderSharedPtr &pSession,
					     const std::weak_ptr<EquationSystem>        &pEquation,
					     const ParamMap &pParams)
      : Filter        (pSession, pEquation)
    {
      //setup filter
      auto it = pParams.find("c_gd");
      if (it != pParams.end())
	{
	  demo.Setchold(std::stod( it->second));
	}

      it = pParams.find("gam_gd");
      if (it != pParams.end())
        {
	  demo.Setgamhold(std::stod(it->second));
	}
      
      // OutputFrequency
      it = pParams.find("OutputFrequency");
      if (it == pParams.end())
	{
	  m_outputFrequency = 1;
	}
      
    }

    FilterStructurePres::~FilterStructurePres()
    {
      
    }
    
    void FilterStructurePres::v_Initialise(
					   const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
					   const NekDouble &time)
    {
      // ASSUMPTION: All different kinds of elements (2D + 3D) supported
      // iff the pointstype, basistype, order, npts across all 2D elements remain the same
      // and the same for all 3D elements
      // across all the elements, including the 2D sides of 3D elements

      nfields = pFields.size(); // assumption: only 1 variable

      //loop this for every variable and initialize only the ptrs that are still null
      nelmt = pFields[0]->GetExpSize();

      MultiRegions::ExpansionType exptype =  pFields[0]->GetExpType();
      // for(int k = 0; k< nelmt; k++)
      // 	{
      // 	  LocalRegions::ExpansionSharedPtr exp1 = pFields[0]->GetExp(k);
      // 	  cout<<"\n "<<exp1->DetShapeType()<<" "<< exp1->GetBasis(0)->GetNumModes();
      // 	}

      if (E3seg == nullptr)
	{
	  Array<OneD, int> nmodesperel = pFields[0]->EvalBasisNumModesMaxPerExp();

	  int orderseg3 = 3*Vmath::Vmax( nmodesperel.size(), nmodesperel, 1);
	  call_setup_seg(pFields[0]->GetExp(0), orderseg3);

	}

      C = demo.formCompanion(E3seg->GetBasis(0)->GetNumModes());
      if(exptype > 1)
	{
	  for(int k = 0; k < nfields ; k++)
	    {
	      // find out what is the max num of quad points in x-dir
	      //	      int segP = 0;
	      for(int i = 0; i < nelmt; i++)
		{
		  
		  LocalRegions::ExpansionSharedPtr exp = pFields[k]->GetExp(i);
		  switch(exp->DetShapeType())
		    {
		    case LibUtilities::eHexahedron:
		      if(Equad == nullptr)
			{
			  call_setup_quad(exp);
			}
		      
		      if(Ehex == nullptr)
		        {
			  call_setup_hex(exp);
			}
	
		      break;
		    case LibUtilities::eTetrahedron:

		      if(Etri == nullptr)
			{
			  call_setup_tri(exp);
			}

		      if(Etet == nullptr)
			{
			  call_setup_tet(exp);

			}
		      break;
		    case LibUtilities::ePyramid:
		      if(Etri == nullptr)
			{
			  call_setup_tri(exp);
			}
		      if(Equad == nullptr)
			{
			  call_setup_quad(exp);
			}
		      if(Epyr == nullptr)
			{
			  call_setup_pyr(exp);
			}
		      break;
		    case LibUtilities::eQuadrilateral:
		      if(Equad == nullptr)
			{
			  call_setup_quad(exp);
			}
	              break;
		      
		    case  LibUtilities::eTriangle:
		      
		      if(Etri == nullptr)
			{
			  call_setup_tri(exp);
			}
		      
		      break;
		      
		    default:
		      cout<<"\n unsupported element for this filter!\n";
		      exit(0);
		      break;
		    }
		}
	      // LibUtilities::PointsType pointsTypeCheb = LibUtilities::eGaussGaussChebyshev;
	      // LibUtilities::PointsKey pkeycheb(segP, pointsTypeCheb);
	      // LibUtilities::BasisKey bkeycheb(LibUtilities::eChebyshev, 3*(segN),  pkeycheb);
	      // E3seg = new StdSegExp(bkeycheb);
	    }
	}
      v_Update(pFields, time);

    }

    void FilterStructurePres::call_setup_tri(LocalRegions::ExpansionSharedPtr exp)
    {
      int nmodes0 = exp->GetBasis(0)->GetNumModes();
      
      int nmodes1 = exp->GetBasis(1)->GetNumModes();
      int npts0 = exp->GetBasis(0)->GetNumPoints();
      int npts1 = exp->GetBasis(1)->GetNumPoints();
      // LibUtilities::PointsType ptypeA =  LibUtilities::eGaussLobattoLegendre;
      // LibUtilities::PointsType ptypeB =  LibUtilities::eGaussGaussLegendre;
      PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
      PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
      BasisKey  b0(LibUtilities::eOrtho_A,  nmodes0, p0);
      BasisKey  b1(LibUtilities::eOrtho_B,  nmodes1, p1);
      Etri = new StdTriExp(b0, b1);

      demo.storage2dt =  Etri->GetPhysEvaluateStorage(); 
      demo.coordtri = demo.GetCoords(Etri);
      demo.coordmidtri = demo.GetQuadratureMidCoords(demo.coordtri);
      demo.coordlatticetri = demo.GetLatticeCoords(demo.coordtri, demo.coordmidtri);

      demo.midptevaltri = Etri->PhysEvaluateBasis(demo.coordmidtri, demo.storage2dt, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

      Array<OneD, Array<OneD, NekDouble> > edgeptsin(2), edgexy(2);
      for(int k = 0; k <2; k++)
	{
	  edgexy[k] = Array<OneD, NekDouble>(E3seg->GetBasis(0)->GetNumPoints());
	  edgeptsin[k] =   Array<OneD, NekDouble>(edgexy[k]);
	}
      Array<OneD, NekDouble> edgexytemp = E3seg->GetBasis(0)->GetZ();
      int totszedges = edgexytemp.size()*(Etri->GetNcoeffs());
      
      // left (x = -1)
      Vxm1t = Array<OneD, NekDouble>(totszedges);
      Vdyxm1t = Array<OneD, NekDouble>(totszedges);
      Vdxxm1t = Array<OneD, NekDouble>(totszedges);
      
      // left x = -1
      edgexy[1] = edgexytemp;
      edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
      Vxm1t = Etri->PhysEvaluateBasis(edgexy, demo.storage2dt, Vdxxm1t, Vdyxm1t, NullNekDouble1DArray);
      
      // bot (y = -1)
      Vym1t  = Array<OneD, NekDouble>(totszedges);
      Vdxym1t  = Array<OneD, NekDouble>(totszedges);
      Vdyym1t  = Array<OneD, NekDouble>(totszedges);
      
      // bot y = -1
      edgexy[0] = edgexytemp;             
      edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
      
      Vym1t = Etri->PhysEvaluateBasis(edgexy,  demo.storage2dt, Vdxym1t, Vdyym1t, NullNekDouble1DArray);
      
      // hypt tri (y = -x)
      Vxyhypt = Array<OneD, NekDouble>(totszedges);
      Vdxxyhypt = Array<OneD, NekDouble>(totszedges);
      Vdyxyhypt = Array<OneD, NekDouble>(totszedges);
      
      edgexy[0] = edgexytemp;
      Vmath::Smul(edgexy[0].size(), -1.0, edgexy[0], 1, edgexy[1] , 1);
      Vxyhypt = Etri->PhysEvaluateBasis(edgexy, demo.storage2dt,  Vdxxyhypt, Vdyxyhypt, NullNekDouble1DArray);
      
    }

    void FilterStructurePres::call_setup_seg(LocalRegions::ExpansionSharedPtr exp, int orderseg3)
    {

      LibUtilities::PointsKey pkeycheb(exp->GetBasis(0)->GetNumPoints()  , LibUtilities::eGaussLobattoChebyshev);
      LibUtilities::BasisKey bkeyorth(LibUtilities::eOrtho_A, orderseg3,  pkeycheb);
      E3seg = new StdSegExp(bkeyorth);
      
    }
    
    void FilterStructurePres::call_setup_quad(LocalRegions::ExpansionSharedPtr exp)
    {
      // get # of quad pts in curr exp
      // get order of current exp
      // get ptypes and btypes
      int nmodes0 = exp->GetBasis(0)->GetNumModes();
      int nmodes1 = exp->GetBasis(1)->GetNumModes();
      int npts0 = exp->GetBasis(0)->GetNumPoints();
      int npts1 = exp->GetBasis(1)->GetNumPoints();
      
      PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
      PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
      BasisKey  b0(LibUtilities::eOrtho_A,  nmodes0, p0);
      BasisKey  b1(LibUtilities::eOrtho_A,  nmodes1, p1);
      Equad = new StdQuadExp(b0, b1);
      demo.storage2dq =  Equad->GetPhysEvaluateStorage();
      demo.coordquad = demo.GetCoords(Equad);
      demo.coordmidquad = demo.GetQuadratureMidCoords(demo.coordquad);
      demo.coordlatticequad = demo.GetLatticeCoords(demo.coordquad, demo.coordmidquad);

      demo.midptevalquad = Equad->PhysEvaluateBasis(demo.coordmidquad, demo.storage2dq, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
      

      Array<OneD, Array<OneD, NekDouble> > edgeptsin(2), edgexy(2);
      for(int k = 0; k <2; k++)
	{
	  edgexy[k] = Array<OneD, NekDouble>(E3seg->GetBasis(0)->GetNumPoints());
	  edgeptsin[k] =   Array<OneD, NekDouble>(edgexy[k]);
	}
      Array<OneD, NekDouble> edgexytemp =   E3seg->GetBasis(0)->GetZ();
      int totszedges = edgexytemp.size()*(Equad->GetNcoeffs());
      
      // left (x = -1)
      Vxm1q = Array<OneD, NekDouble>(totszedges);
      Vdyxm1q = Array<OneD, NekDouble>(totszedges);
      Vdxxm1q = Array<OneD, NekDouble>(totszedges);
      
      // left x = -1
      edgexy[1] = edgexytemp;
      edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
      Vxm1q = Equad->PhysEvaluateBasis(edgexy, demo.storage2dq, Vdxxm1q, Vdyxm1q, NullNekDouble1DArray);
      
      // bot (y = -1)
      Vym1q  = Array<OneD, NekDouble>(totszedges);
      Vdxym1q  = Array<OneD, NekDouble>(totszedges);
      Vdyym1q  = Array<OneD, NekDouble>(totszedges);
      
      // bot y = -1
      edgexy[0] = edgexytemp;             
      edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
      
      Vym1q = Equad->PhysEvaluateBasis(edgexy,  demo.storage2dq, Vdxym1q, Vdyym1q, NullNekDouble1DArray);

      	  
      // right quad (x = 1)
      Vx1q = Array<OneD, NekDouble>(totszedges);
      Vdxx1q = Array<OneD, NekDouble>(totszedges);
      Vdyx1q = Array<OneD, NekDouble>(totszedges);
      
      // top quad (y = 1)
      Vdxy1q = Array<OneD, NekDouble>(totszedges);
      Vdyy1q = Array<OneD, NekDouble>(totszedges);
      Vy1q = Array<OneD, NekDouble>(totszedges);
      
      // right x = 1
      edgexy[1] = edgexytemp; 
      edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), 1.0); 
      Vx1q = Equad->PhysEvaluateBasis(edgexy, demo.storage2dq, Vdxx1q, Vdyx1q, NullNekDouble1DArray);
      
      //top y = 1
      edgexy[0] = edgexytemp; 
      edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), 1.0); 
      Vy1q = Equad->PhysEvaluateBasis(edgexy, demo.storage2dq, Vdxy1q, Vdyy1q, NullNekDouble1DArray);
      
    }
    
    void FilterStructurePres::call_setup_tet(LocalRegions::ExpansionSharedPtr exp)
    {
      int nmodes0 = exp->GetBasis(0)->GetNumModes();
      int nmodes1 = exp->GetBasis(1)->GetNumModes();
      int nmodes2 = exp->GetBasis(2)->GetNumModes();
      int npts0 = exp->GetBasis(0)->GetNumPoints();	
      int npts1 = exp->GetBasis(1)->GetNumPoints();
      int npts2 = exp->GetBasis(2)->GetNumPoints();
      PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
      PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
      PointsKey p2(npts2, exp->GetBasis(2)->GetPointsType());
      BasisKey  b0(LibUtilities::eOrtho_A,  nmodes0, p0);
      BasisKey  b1(LibUtilities::eOrtho_B,  nmodes1, p1);
      BasisKey  b2(LibUtilities::eOrtho_C,  nmodes2, p2);
      Etet = new StdTetExp(b0, b1, b2);
      demo.storage3dtet =  Etet->GetPhysEvaluateStorage();
      demo.coordtet = demo.GetCoords(Etet);
      demo.coordmidtet = demo.GetQuadratureMidCoords(demo.coordtet);
      demo.coordlatticetri = demo.GetLatticeCoords(demo.coordtet, demo.coordmidtet);

      demo.midptevaltet = Etet->PhysEvaluateBasis(demo.coordmidtet, demo.storage3dtet, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

      Array<OneD, NekDouble> edgexyztemp =   E3seg->GetBasis(0)->GetZ();
      int totszedges1d =  edgexyztemp.size()*(exp->GetNcoeffs());
      Array<OneD, Array<OneD, NekDouble> > edgeptsin(dimension);
      for(int p = 0; p < dimension; p++)
	{
	  edgeptsin[p] =   Array<OneD, NekDouble>(edgexyztemp);
	}
      
      // edge front left (AD) (x = -1) (y = -1)                                                
      Vxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
      Vdxxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
      Vdyxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
      Vdzxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
      edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
      edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

      Vxm1ym1ztet = Etet->PhysEvaluateBasis(edgeptsin, demo.storage3dtet, Vdxxm1ym1ztet, Vdyxm1ym1ztet, Vdzxm1ym1ztet);

      //edge front hypt (DB) (y = -1) (z = -x)                                                 
      Vym1xmztet = Array<OneD, NekDouble>(totszedges1d);
      Vdxym1xmztet = Array<OneD, NekDouble>(totszedges1d);
      Vdyym1xmztet = Array<OneD, NekDouble>(totszedges1d);
      Vdzym1xmztet = Array<OneD, NekDouble>(totszedges1d);

      Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[2][0], 1, &edgeptsin[0][0], 1);

      Vym1xmztet = Etet->PhysEvaluateBasis(edgeptsin,demo.storage3dtet, Vdxym1xmztet, Vdyym1xmztet, Vdzym1xmztet);

      //edge front bot (AB) (y = -1) (z = -1)                                                
      Vym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
      Vdxym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
      Vdyym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
      Vdzym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
      edgeptsin[0] = edgexyztemp;
      edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

      Vym1xzm1tet = Etet->PhysEvaluateBasis(edgeptsin, demo.storage3dtet, Vdxym1xzm1tet, Vdyym1xzm1tet, Vdzym1xzm1tet);

      //edge left hypt (DC) ( x = -1) (z = -y)                                               
      Vxm1ymztet = Array<OneD, NekDouble>(totszedges1d);
      Vdxxm1ymztet= Array<OneD, NekDouble>(totszedges1d);
      Vdyxm1ymztet = Array<OneD, NekDouble>(totszedges1d);
      Vdzxm1ymztet= Array<OneD, NekDouble>(totszedges1d);
      edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
      edgeptsin[1] = edgexyztemp;
      Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[1][0], 1, &edgeptsin[2][0], 1);
  
      Vxm1ymztet = Etet->PhysEvaluateBasis(edgeptsin,demo.storage3dtet, Vdxxm1ymztet, Vdyxm1ymztet, Vdzxm1ymztet);

      // edge bot diag (BC) (z = -1) (y = -x)                                                
      Vxmyzm1tet = Array<OneD, NekDouble>(totszedges1d);
      Vdxxmyzm1tet = Array<OneD, NekDouble>(totszedges1d)   ;
      Vdyxmyzm1tet = Array<OneD, NekDouble>(totszedges1d)    ;
      Vdzxmyzm1tet = Array<OneD, NekDouble>(totszedges1d)    ;
      edgeptsin[0] = edgexyztemp;
      Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[0][0], 1, &edgeptsin[1][0], 1);
      edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[2].size(), -1.0);
      Vxmyzm1tet = Etet->PhysEvaluateBasis(edgeptsin, demo.storage3dtet, Vdxxmyzm1tet, Vdyxmyzm1tet, Vdzxmyzm1tet);
  
      //edge CA bot left (x = -1) (z = -1)                                                   
      Vxm1yzm1tet   = Array<OneD, NekDouble>(totszedges1d);
      Vdxxm1yzm1tet = Array<OneD, NekDouble>(totszedges1d)   ;
      Vdyxm1yzm1tet = Array<OneD, NekDouble>(totszedges1d)    ;
      Vdzxm1yzm1tet = Array<OneD, NekDouble>(totszedges1d)    ;
      edgeptsin[1] = edgexyztemp;
      edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[2].size(), -1.0);
  
      Vxm1yzm1tet = Etet->PhysEvaluateBasis(edgeptsin, demo.storage3dtet, Vdxxm1yzm1tet, Vdyxm1yzm1tet, Vdzxm1yzm1tet);
      int totsurf2d =  (demo.coordtri[0].size())*Etet->GetNcoeffs();
      Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);
  
      for(int k = 0; k < dimension-1; k++)
	{
	  surfptsin[k] = Array<OneD, NekDouble>(demo.coordtri[k]);
	  surfptsintemp[k] = Array<OneD, NekDouble>(demo.coordtri[k]);
	}
      surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.coordtri[0].size());
      surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.coordtri[0].size());

      int totpt = surfptsin[0].size();

      //surface bot z = -1, (ABC) 
      Vxyzm1tet = Array<OneD, NekDouble>(totsurf2d);
      surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
      surfptsintemp[0] = surfptsin[0];
      surfptsintemp[1] = surfptsin[1];
      Vxyzm1tet = Etet->PhysEvaluateBasis(surfptsintemp, demo.storage3dtet, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

      //surface left x = -1  (DAC)  
      Vxm1ypz0tet =  Array<OneD, NekDouble>(totsurf2d);
      surfptsintemp[0] = Array<OneD, NekDouble> (totpt, -1.0);
      surfptsintemp[1] = surfptsin[0];
      surfptsintemp[2] = surfptsin[1];
  
      Vxm1ypz0tet = Etet->PhysEvaluateBasis(surfptsintemp, demo.storage3dtet, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

      //surf front y = -1,  (DAB) 
      Vxpz0ym1tet =  Array<OneD, NekDouble>(totsurf2d);
      surfptsintemp[1] = Array<OneD, NekDouble> (totpt, -1.0);
      surfptsintemp[0] = surfptsin[0];
      surfptsintemp[2] = surfptsin[1];
      Vxpz0ym1tet = Etet->PhysEvaluateBasis(surfptsintemp, demo.storage3dtet, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

      //surf DCB (x + y + z = -1), 
      Vxpypzm1tet =  Array<OneD, NekDouble>(totsurf2d);
      surfptsintemp[1] = surfptsin[1];;

      Vmath::Vadd(totpt, &surfptsin[0][0], 1, &surfptsin[1][0], 1, &surfptsintemp[2][0], 1);
      Vmath::Smul(totpt, -1.0, &surfptsintemp[2][0], 1, &surfptsintemp[2][0], 1);
      Vmath::Sadd(totpt, -1.0, &surfptsintemp[2][0], 1, &surfptsintemp[2][0], 1);
      Vxpypzm1tet = Etet->PhysEvaluateBasis(surfptsintemp,demo.storage3dtet , NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

      
    }

    void FilterStructurePres::call_setup_hex(LocalRegions::ExpansionSharedPtr exp)
    {
      int nmodes0 = exp->GetBasis(0)->GetNumModes();
      int nmodes1 = exp->GetBasis(1)->GetNumModes();
      int nmodes2 = exp->GetBasis(2)->GetNumModes();
      int npts0 = exp->GetBasis(0)->GetNumPoints();
      int npts1 = exp->GetBasis(1)->GetNumPoints();
      int npts2 = exp->GetBasis(2)->GetNumPoints(); 
      PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
      PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
      PointsKey p2(npts2, exp->GetBasis(2)->GetPointsType());
      BasisKey  b0(LibUtilities::eOrtho_A,  nmodes0, p0);                
      BasisKey  b1(LibUtilities::eOrtho_A,  nmodes1, p1);
      BasisKey  b2(LibUtilities::eOrtho_A,  nmodes2, p2);                        
      Ehex = new StdHexExp(b0, b1, b2);
      demo.storage3dhex =  Ehex->GetPhysEvaluateStorage();
      demo.coordhex = demo.GetCoords(Ehex);
      demo.coordmidhex = demo.GetQuadratureMidCoords(demo.coordhex);
      demo.coordlatticehex = demo.GetLatticeCoords(demo.coordhex, demo.coordmidhex);

      demo.midptevalhex = Ehex->PhysEvaluateBasis(demo.coordmidhex, demo.storage3dhex, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
      
      Array<OneD, NekDouble> edgexyztemp =   E3seg->GetBasis(0)->GetZ();
      int totszedges1d = edgexyztemp.size()*(Ehex->GetNcoeffs());
      Array<OneD, Array<OneD, NekDouble> > edgeptsin (dimension);    
      for(int p = 0; p < dimension; p++)
	{
	  edgeptsin[p] = Array<OneD, NekDouble>(edgexyztemp);
	} 

      // edge front left (x = -1) (y = -1)
      Vxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
      Vdxxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
      Vdyxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
      Vdzxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
      edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
      edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

            
      Vxm1ym1z = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);  
                
      //edge front right (x = 1) (y = -1)
      Vx1ym1z = Array<OneD, NekDouble>(totszedges1d);
      Vdxx1ym1z = Array<OneD, NekDouble>(totszedges1d);
      Vdyx1ym1z = Array<OneD, NekDouble>(totszedges1d);
      Vdzx1ym1z = Array<OneD, NekDouble>(totszedges1d);
      edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);       
      Vx1ym1z = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxx1ym1z, Vdyx1ym1z, Vdzx1ym1z);  


      //edge front top (y = -1) (z = 1)
      Vym1xz1 = Array<OneD, NekDouble>(totszedges1d);
      Vdxym1xz1 = Array<OneD, NekDouble>(totszedges1d);
      Vdyym1xz1 = Array<OneD, NekDouble>(totszedges1d);
      Vdzym1xz1 = Array<OneD, NekDouble>(totszedges1d);
      edgeptsin[0] = edgexyztemp; 
      edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
      Vym1xz1 = Ehex->PhysEvaluateBasis(edgeptsin,demo.storage3dhex, Vdxym1xz1, Vdyym1xz1, Vdzym1xz1);  

                
      //edge front bot (y = -1) (z = -1)
      Vym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
      Vdxym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
      Vdyym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
      Vdzym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
      edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
      Vym1xzm1 =             Ehex->PhysEvaluateBasis( edgeptsin, demo.storage3dhex, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);  
	    
	    
      // edge back left (y = 1), (x = -1)
      Vxm1y1z = Array<OneD, NekDouble>(totszedges1d);
      Vdxxm1y1z = Array<OneD, NekDouble>(totszedges1d);
      Vdyxm1y1z = Array<OneD, NekDouble>(totszedges1d);
      Vdzxm1y1z = Array<OneD, NekDouble>(totszedges1d);
      edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
      edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);      
      edgeptsin[2] =    edgexyztemp;


      Vxm1y1z = Ehex->PhysEvaluateBasis(edgeptsin,demo.storage3dhex, Vdxxm1y1z, Vdyxm1y1z, Vdzxm1y1z);  

      //edge back right (x = 1), (y = 1))
      Vx1y1z = Array<OneD, NekDouble>(totszedges1d);
      Vdxx1y1z = Array<OneD, NekDouble>(totszedges1d);
      Vdyx1y1z = Array<OneD, NekDouble>(totszedges1d);
      Vdzx1y1z = Array<OneD, NekDouble>(totszedges1d);
      edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);                
         
      Vx1y1z  = Ehex->PhysEvaluateBasis(edgeptsin,demo.storage3dhex, Vdxx1y1z, Vdyx1y1z, Vdzx1y1z);  

                
      //edge back top ( y = 1) (z = 1)
      Vy1xz1 = Array<OneD, NekDouble>(totszedges1d);
      Vdxy1xz1 = Array<OneD, NekDouble>(totszedges1d);
      Vdyy1xz1 = Array<OneD, NekDouble>(totszedges1d);
      Vdzy1xz1 = Array<OneD, NekDouble>(totszedges1d);
      edgeptsin[0] = edgexyztemp;
      edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
      Vy1xz1 = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxy1xz1, Vdyy1xz1, Vdzy1xz1 );  

      
      //edge back bot (y = 1) (z = -1))
      Vy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
      Vdxy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
      Vdyy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
      Vdzy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
      edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
            
      Vy1xzm1 = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);  

      // edge left bot (z = -1), (x = -1)
      Vxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
      Vdxxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
      Vdyxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
      Vdzxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
      edgeptsin[1] = edgexyztemp;
      edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

      Vxm1yzm1 = Ehex->PhysEvaluateBasis(edgeptsin,demo.storage3dhex, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);  
                
      //edge left top (x = -1), (z = 1))
      Vxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
      Vdxxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
      Vdyxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
      Vdzxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
                
      edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);           
      Vxm1yz1 = Ehex->PhysEvaluateBasis(edgeptsin,demo.storage3dhex, Vdxxm1yz1, Vdyxm1yz1, Vdzxm1yz1);  

      //edge right bot ( z = -1) (x = 1)
      Vx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
      Vdxx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
      Vdyx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
      Vdzx1yzm1 = Array<OneD, NekDouble>(totszedges1d);

      edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
      edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);           
      Vx1yzm1 = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);   

      //edge right top (z  1) (x  1))
      Vx1yz1 = Array<OneD, NekDouble>(totszedges1d);
      Vdxx1yz1 = Array<OneD, NekDouble>(totszedges1d);
      Vdyx1yz1 = Array<OneD, NekDouble>(totszedges1d);
      Vdzx1yz1 = Array<OneD, NekDouble>(totszedges1d);
                
      edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);

      Vx1yz1 = Ehex->PhysEvaluateBasis(edgeptsin,demo.storage3dhex, Vdxx1yz1, Vdyx1yz1, Vdzx1yz1);   
      int totszsurf2d = (demo.coordquad[0].size())*Ehex->GetNcoeffs();
            
      Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);
      for(int k = 0; k < dimension-1; k++)
	{
	  surfptsin[k] = Array<OneD, NekDouble>(demo.coordquad[k]);
	  surfptsintemp[k] = Array<OneD, NekDouble>(demo.coordquad[k]);
	}
            
      surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.coordquad[0].size()); 
      surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.coordquad[0].size()); 
            
      //surface bot z = -1
      Vxyzm1 =  Array<OneD, NekDouble>(totszsurf2d);
      surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
      surfptsintemp[0] = surfptsin[0];
      surfptsintemp[1] = surfptsin[1];
      Vxyzm1 =             Ehex->PhysEvaluateBasis(surfptsintemp,demo.storage3dhex, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

      //surface right x = 1
      Vx1yz =  Array<OneD, NekDouble>(totszsurf2d);
      surfptsintemp[0] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
      surfptsintemp[1] =  surfptsin[0];
      surfptsintemp[2] =  surfptsin[1];
      Vx1yz = Ehex->PhysEvaluateBasis(surfptsintemp,demo.storage3dhex,NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   


      //surface top z = 1
      Vxyz1 =  Array<OneD, NekDouble>(totszsurf2d);
      surfptsintemp[0] = surfptsin[0];
      surfptsintemp[1] = surfptsin[1];
      surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
      Vxyz1 = Ehex->PhysEvaluateBasis(surfptsintemp, demo.storage3dhex, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

      //surface left x = -1
      Vxm1yz =  Array<OneD, NekDouble>(totszsurf2d);
      surfptsintemp[0] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
      surfptsintemp[1] = surfptsin[0];
      surfptsintemp[2] = surfptsin[1];
      Vxm1yz = Ehex->PhysEvaluateBasis(surfptsintemp,demo.storage3dhex, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

      //surface front y = -1
      Vxym1z =  Array<OneD, NekDouble>(totszsurf2d);
      surfptsintemp[1] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
      surfptsintemp[0] = surfptsin[0];
      surfptsintemp[2] = surfptsin[1];
      Vxym1z = Ehex->PhysEvaluateBasis(surfptsintemp,demo.storage3dhex, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   
      
      //surface back y = 1
      Vxy1z =  Array<OneD, NekDouble>(totszsurf2d);
      surfptsintemp[1] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
      surfptsintemp[0] = surfptsin[0];
      surfptsintemp[2] = surfptsin[1];
      Vxy1z = Ehex->PhysEvaluateBasis(surfptsintemp,demo.storage3dhex, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
    }

    void FilterStructurePres::call_setup_pyr(LocalRegions::ExpansionSharedPtr exp)
    {
      int nmodes0 = exp->GetBasis(0)->GetNumModes();
      int nmodes1 = exp->GetBasis(1)->GetNumModes();
      int nmodes2 = exp->GetBasis(2)->GetNumModes();
      int npts0 = exp->GetBasis(0)->GetNumPoints();
      int npts1 = exp->GetBasis(1)->GetNumPoints();
      int npts2 = exp->GetBasis(2)->GetNumPoints();
      PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
      PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
      PointsKey p2(npts2, exp->GetBasis(2)->GetPointsType());
      BasisKey  b0(LibUtilities::eOrtho_A,  nmodes0, p0);
      BasisKey  b1(LibUtilities::eOrtho_A,  nmodes1, p1);
      BasisKey  b2(LibUtilities::eOrthoPyr_C,  nmodes2, p2);
      Epyr = new StdPyrExp(b0, b1, b2);
      demo.storage3dpyr =  Epyr->GetPhysEvaluateStorage();
      demo.coordpyr = demo.GetCoords(Epyr);
      demo.coordmidpyr = demo.GetQuadratureMidCoords(demo.coordpyr);
      demo.coordlatticepyr = demo.GetLatticeCoords(demo.coordpyr, demo.coordmidpyr);

      demo.midptevalpyr = Ehex->PhysEvaluateBasis(demo.coordmidpyr, demo.storage3dpyr, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

      //we need qzin + qzinmid as input points                                                         
      Array<OneD, NekDouble> edgexyztemp =   E3seg->GetBasis(0)->GetZ();
      int totszedges1d =  edgexyztemp.size()*(Epyr->GetNcoeffs());      
      Array<OneD, Array<OneD, NekDouble> >edgeptsin(dimension);
      for(int p = 0; p < dimension; p++)
	{
	  edgeptsin[p] = Array<OneD, NekDouble>(edgexyztemp);
	}
		
      // edge front left EA (x = -1) (y = -1)                 
      Vxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);
      Vdxxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);
      Vdyxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);
      Vdzxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);

      edgeptsin[0]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
      edgeptsin[1]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
      Vxm1ym1zpyr   = Epyr->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxxm1ym1zpyr, Vdyxm1ym1zpyr, Vdzxm1ym1zpyr);

      //edge front hypt EB (y = -1) (z + x = 0)
      Vym1xmzpyr    = Array<OneD, NekDouble>(totszedges1d);
      Vdxym1xmzpyr  = Array<OneD, NekDouble>(totszedges1d);
      Vdyym1xmzpyr  = Array<OneD, NekDouble>(totszedges1d);
      Vdzym1xmzpyr  = Array<OneD, NekDouble>(totszedges1d);
      Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[2][0], 1, &edgeptsin[0][0], 1);

      Vxm1ym1zpyr =  Epyr->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxym1xmzpyr, Vdxym1xmzpyr, Vdzym1xmzpyr);
      //edge back hypt (z = -x  or x = -z) (x = y) EC
      Vxeyxmzpyr    = Array<OneD, NekDouble>(totszedges1d);
      Vdxxeyxmzpyr  = Array<OneD, NekDouble>(totszedges1d);
      Vdyxeyxmzpyr  = Array<OneD, NekDouble>(totszedges1d);
      Vdzxeyxmzpyr  = Array<OneD, NekDouble>(totszedges1d);

      edgeptsin[1]  =  Array<OneD, NekDouble>(edgexyztemp);
      Vxeyxmzpyr    =  Epyr->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxxeyxmzpyr, Vdyxeyxmzpyr, Vdzxeyxmzpyr);

      //edge back left (y = -z) (x = -1) ED
      Vx1ymzpyr     = Array<OneD, NekDouble>(totszedges1d);
      Vdxx1ymzpyr   = Array<OneD, NekDouble>(totszedges1d);
      Vdyx1ymzpyr   = Array<OneD, NekDouble>(totszedges1d);
      Vdzx1ymzpyr   = Array<OneD, NekDouble>(totszedges1d); ;
      Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[2][0], 1, &edgeptsin[1][0], 1);

      edgeptsin[0]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
      Vx1ymzpyr     = Epyr->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxx1ymzpyr, Vdyx1ymzpyr, Vdzx1ymzpyr);
  
      //edge front bot (y = -1) (z = -1) AB
      Vym1xzm1pyr   = Array<OneD, NekDouble>(totszedges1d);
      Vdxym1xzm1pyr = Array<OneD, NekDouble>(totszedges1d)   ;
      Vdyym1xzm1pyr = Array<OneD, NekDouble>(totszedges1d)    ;
      Vdzym1xzm1pyr = Array<OneD, NekDouble>(totszedges1d)    ;

      edgeptsin[1]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
      edgeptsin[2]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
      edgeptsin[0]  = Array<OneD, NekDouble>(edgexyztemp);
      Vym1xzm1pyr   = Epyr->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxym1xzm1pyr, Vdyym1xzm1pyr,Vdzym1xzm1pyr);

  
      //edge back bot (y = 1) (z = -1)) DC
      Vy1xzm1pyr    = Array<OneD, NekDouble>(totszedges1d);
      Vdyy1xzm1pyr  = Array<OneD, NekDouble>(totszedges1d)   ;
      Vdxy1xzm1pyr  = Array<OneD, NekDouble>(totszedges1d);
      Vdzy1xzm1pyr  = Array<OneD, NekDouble>  (totszedges1d) ;

      edgeptsin[0]  = edgexyztemp;
      edgeptsin[1]  = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
      Vy1xzm1pyr   = Epyr->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxy1xzm1pyr, Vdyy1xzm1pyr, Vdzy1xzm1pyr);

  
      // edge left bot (z = -1), (x = -1) AD
      Vxm1yzm1pyr   = Array<OneD, NekDouble>(totszedges1d);
      Vdyxm1yzm1pyr = Array<OneD, NekDouble>(totszedges1d)   ;
      Vdxxm1yzm1pyr = Array<OneD, NekDouble>(totszedges1d)   ;
      Vdzxm1yzm1pyr = Array<OneD, NekDouble>(totszedges1d)   ;

      edgeptsin[1]  = edgexyztemp;
      edgeptsin[0]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
      edgeptsin[2]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
      Vxm1yzm1pyr   =  Epyr->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxxm1yzm1pyr, Vdyxm1yzm1pyr, Vdzxm1yzm1pyr);

  
      //edge right bot ( z = -1) (x = 1) BC
      Vx1yzm1pyr    = Array<OneD, NekDouble>(totszedges1d);
      Vdyx1yzm1pyr  = Array<OneD, NekDouble>(totszedges1d);
      Vdxx1yzm1pyr  = Array<OneD, NekDouble>(totszedges1d);
      Vdzx1yzm1pyr  = Array<OneD, NekDouble>(totszedges1d);

      edgeptsin[1]  = edgexyztemp;
      edgeptsin[0]  = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
      edgeptsin[2]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
      Vx1yzm1pyr    =  Epyr->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxx1yzm1pyr, Vdyx1yzm1pyr, Vdzx1yzm1pyr);

      int totsurf2d = (demo.coordquad[0].size())*Epyr->GetNcoeffs();
            
      Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);
      for(int k = 0; k < dimension-1; k++)
	{
	  surfptsin[k] = Array<OneD, NekDouble>(demo.coordquad[k]);
	  surfptsintemp[k] = Array<OneD, NekDouble>(demo.coordquad[k]);
	}
            
      surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.coordquad[0].size()); 
      surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.coordquad[0].size()); 
            
      //surface bot z = -1, (ABC)  
      Vxyzm1pyr = Array<OneD, NekDouble>(totsurf2d);
      surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
      surfptsintemp[0] = surfptsin[0];
      surfptsintemp[1] = surfptsin[1];
      Vxyzm1pyr = Epyr->PhysEvaluateBasis(surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  
      totsurf2d =  (demo.coordtri[0].size())*Epyr->GetNcoeffs();
      //Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);

      for(int k = 0; k < dimension-1; k++)
	{
	  surfptsin[k] = Array<OneD, NekDouble>(demo.coordtri[k]);
	  surfptsintemp[k] = Array<OneD, NekDouble>(demo.coordtri[k]);
	}
      surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.coordtri[0].size());
      surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.coordtri[0].size());

      int totpt = surfptsin[0].size();
  
      //surface hypt (tri) x+z = 0 && y + z = 0   
      Vxmzypyr = Array<OneD, NekDouble>(totsurf2d);
      Vmath::Smul(totpt, -1.0, &surfptsintemp[0][0], 1, &surfptsintemp[2][0], 1);
      Vmath::Vcopy(totpt,  &surfptsintemp[0][0], 1,  &surfptsintemp[1][0], 1);
      Vxmzypyr = Epyr->PhysEvaluateBasis(surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

      //surface left (tri) x = -1   
      Vxm1yzpyr = Array<OneD, NekDouble>(totsurf2d);
      surfptsintemp[1] = surfptsin[1];
      surfptsintemp[2] = surfptsin[2];
      surfptsintemp[0] = Array<OneD, NekDouble>(totpt, -1.0);
      Vxm1yzpyr = Epyr->PhysEvaluateBasis(surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
  
      //surface front (tri) y = -1     
      Vxym1zpyr = Array<OneD, NekDouble>(totsurf2d);
      surfptsintemp[1] = Array<OneD, NekDouble>(totpt, -1.0);
      surfptsintemp[0] = surfptsin[0];
      Vxym1zpyr = Epyr->PhysEvaluateBasis(surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

      //surface back (tri) y = -z
      Vxymzpyr = Array<OneD, NekDouble>(totsurf2d);
      Vmath::Smul(totpt, -1.0, &surfptsintemp[2][0], 1, &surfptsintemp[1][0], 1);
      Vxymzpyr =  Epyr->PhysEvaluateBasis(surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

    }

    // strategy 1 (current impl:): Orthonormalize all elements and then  flag them using stdelements
    // predef using ortho basis
    
    // strategy 2 or,create stdelements with appropriate basis "non-ortho" and evaluate phys on lattice
    // then orthonormalize only flagged elements
    
    Array<OneD, Array<OneD, NekDouble> > FilterStructurePres::FindIds(
								      const MultiRegions::ExpListSharedPtr &pFields)
    {

      
      vector<Array<OneD, NekDouble>> minvalandpt; 
      Array<OneD, Array<OneD, NekDouble> >retarr;
      Array<OneD, NekDouble> physallel, coeffsallel;

      coeffsallel = pFields->GetCoeffs();
      int nelmt = pFields->GetExpSize();
      int ctp = 0, coeffs_offset = 0, phys_offset = 0;
      int dim = pFields->GetShapeDimension();
      Array<OneD, Array<OneD, NekDouble> > coords(dim);

      for(int p = 0; p < dim; p++)
	{
	  coords[p] = Array<OneD, NekDouble>( pFields->GetPhys().size());
	  
	}
      if(dim >2)
	pFields->GetCoords(coords[0], coords[1], coords[2]);
      else if(dim >1)
	pFields->GetCoords(coords[0], coords[1], NullNekDouble1DArray);
      //to-do:handle seg case here
      for(int k = 0; k < nelmt; k++)
      {
	coeffs_offset = pFields->GetCoeff_Offset(k);
	phys_offset = pFields->GetPhys_Offset(k);

	int numlocalphys = pFields->GetTotPoints(k);

	Array<OneD, NekDouble> uhatslocal(pFields->GetNcoeffs(k));
	Vmath::Vcopy(uhatslocal.size(), &coeffsallel[0]+coeffs_offset, 1, &uhatslocal[0], 1);
	Array<OneD, NekDouble> physlocal(numlocalphys);

	LocalRegions::ExpansionSharedPtr ptr = pFields->GetExp(k); 
	ptr->BwdTrans(uhatslocal, physlocal);
	//Vmath::Vcopy(numlocalphys, &physallel[0]+phys_offset, 1, &physlocal[0], 1);
	NekDouble localmin = Vmath::Vmin(numlocalphys, &physlocal[0], 1);
	if(localmin < 0 && abs(localmin) > 1e-8)
	  {
	    int idxmin = Vmath::Imin(numlocalphys, &physlocal[0], 1); 

	    // val, coords, eleid
	    Array<OneD, NekDouble> temprow(2+dim);
	    temprow[0] = localmin;
	    for(int i = 0; i < dim; i++)
	      {
		temprow[i+1] = coords[i][phys_offset+idxmin];
	      }
	    temprow[1+dim] = k;
	    minvalandpt.push_back(temprow);
	  }
	ctp += numlocalphys;
      }
      
      // int ctrc = 0;
      // vector<Array<OneD, NekDouble>> minvalandpt;	
      // for(int k = 0; k < nelmt; k++)
      // 	{
      // 	  LocalRegions::ExpansionSharedPtr ptr = pFields->GetExp(k);
      // 	  int nc = pFields->GetNcoeffs(k);
      // 	  Array<OneD, NekDouble> uhatslocal(nc);
      // 	  Vmath::Vcopy(nc, &uhatsallel[ctrc], 1, &uhatslocal[0], 1);
      // 	  ctrc = ctrc + nc;
	
      // 	  Array<OneD, NekDouble> minvalandpthold;

	  // // check on mid-quadrature pts of staggered lattice
	  // switch(ptr->DetShapeType())
	  //   {
	  //   case LibUtilities::eQuadrilateral:
	  //     minvalandpthold = demo.FindLatticeEval(uhatslocal, demo.storage2dq, demo.coordquad, demo.midptevalquad, demo.coordmidquad);
	  //     break;
	  //   case LibUtilities::eTriangle:
	  //     minvalandpthold = demo.FindLatticeEval(uhatslocal, demo.storage2dt, demo.coordtri, demo.midptevaltri, demo.coordmidtri);
	  //     break;
	  //   case LibUtilities::eHexahedron:
	  //     minvalandpthold = demo.FindLatticeEval(uhatslocal, demo.storage3dhex, demo.coordhex, demo.midptevalhex, demo.coordmidhex);
	  //     break;
	  //   case LibUtilities::eTetrahedron:
	  //     minvalandpthold = demo.FindLatticeEval(uhatslocal, demo.storage3dtet, demo.coordtet, demo.midptevaltet, demo.coordmidtet);
	  //     break;
	  //   case LibUtilities::ePyramid:
	  //     minvalandpthold = demo.FindLatticeEval(uhatslocal, demo.storage3dpyr, demo.coordpyr, demo.midptevalpyr, demo.coordmidpyr);
	  //     break;
	  //   default:
	  //     cout<<"\n not implemented for this shape type yet!";
	  //     exit(0);
	    
	  //   }
      // 	  int sz = minvalandpthold.size();
      // 	  if(minvalandpthold[sz-1] < 0 && abs(minvalandpthold[sz-1])>1e-8)
      // 	    {
      // 	      Array<OneD, NekDouble> tmp(sz+1);
      // 	      for(int p = 0; p < sz; p++)
      // 		{
      // 		  tmp[p] = minvalandpthold[p];
      // 		}
      // 	      tmp[sz] = k; // add ele id at the end
      // 	      minvalandpt.push_back(tmp);
      // 	    }	    
	
      // 	}
       retarr = Array<OneD, Array<OneD, NekDouble> > (minvalandpt.size());

       for(int i = 0; i < retarr.size(); i++)
       	{
       	  retarr[i] = Array<OneD, NekDouble>(minvalandpt[i].size(), minvalandpt[i].data());
       	}
      return retarr;
    }
  
    
    void FilterStructurePres::v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
				       const NekDouble &time)
    {
      boost::ignore_unused(time);

      for(int i = 0; i < nfields; i++)
	{
	  //	  bool pstate = pFields[i]->GetPhysState();
	  Array<OneD, Array<OneD, NekDouble> >flaggedele = FindIds(pFields[i]);
	  
	  if(flaggedele.size() > 0)
	    {

	      // if (m_session->GetComm()->GetRank() == 0 &&
	      // 	  m_session->DefinesCmdLineArgument("verbose"))
	      // 	{
		  cout<<"\nfound neg at: \n val coord_localxyz \t id \n ";

		  for(int k = 0; k < flaggedele.size(); k++)
		    {
		      for(int p = 0; p < flaggedele[k].size(); p++)
			{
			  cout<< flaggedele[k][p]<<"\t";
			}
		      cout<<"\n";
		      //		    }
		}
		  
	      Array<OneD, NekDouble> uhatsallelhold, uhatsallel = pFields[i]->GetCoeffs();
	      Array<OneD, NekDouble>physallel (pFields[i]->GetPhys().size());
	      uhatsallelhold = Array<OneD, NekDouble>(uhatsallel.size());
	     pFields[i]->BwdTrans(uhatsallel,physallel);    
	      Vmath::Vcopy(uhatsallel.size(), uhatsallel, 1, uhatsallelhold, 1);
	      int dim = flaggedele[0].size() - 2;
	      int coeffs_offset = 0, phys_offset = 0, ct = 0, ctc = 0, ctp = 0;
	      int nelmt = pFields[0]->GetExpSize();  
	      Array<OneD, Array<OneD, NekDouble> > coords(dim);
	      for(int p = 0; p < dim; p++)
		{
		  coords[p] = Array<OneD, NekDouble>( physallel.size());
		  
		}
	      if(dim >2)
		pFields[i]->GetCoords(coords[0], coords[1], coords[2]);
	      else if(dim >1)
		pFields[i]->GetCoords(coords[0], coords[1], NullNekDouble1DArray);
	      //to-do:handle seg case here

	      //&& ct < flaggedele.size()
	      for(int k = 0;  k < nelmt; k++)
		{
		  phys_offset = pFields[i]->GetPhys_Offset(k); 
		  coeffs_offset = pFields[i]->GetCoeff_Offset(k); 
		  LocalRegions::ExpansionSharedPtr ptr = pFields[i]->GetExp(k); 
		  
		  Array<OneD, NekDouble> uhatslocal(pFields[i]->GetNcoeffs(k));
		  int numlocalphys = pFields[i]->GetTotPoints(k);
		  Array<OneD, NekDouble> physlocal(numlocalphys);
		  Vmath::Vcopy(numlocalphys, &physallel[0]+phys_offset, 1, &physlocal[0], 1);
		  Vmath::Vcopy(uhatslocal.size(), &uhatsallel[0]+coeffs_offset, 1, &uhatslocal[0], 1);
		  
		  if(ct < flaggedele.size() && flaggedele[ct][dim+1] == k)
		    {
		      cout<<"\n el = "<<k<<" ";
		      Optimize(ptr, uhatslocal, physlocal, flaggedele[ct]);

		      // cout<<"\n AFTER:for el = "<<ct<<" \ncoord local and  phys local=\n";
		      // for(int m = 0; m < physlocal.size(); m++)
		      // 	cout<<" \n"<<" "<<coords[0][m]<<" "<<coords[1][m]<<" ="<<physlocal[m];
		      ptr->FwdTrans(physlocal, uhatslocal);
				  
		      ct++;		     
		      // if(ct == flaggedele.size())
		      // 	{
		      // 	  break;
		      // 	}
		    }
		  
		  Vmath::Vcopy(uhatslocal.size(), &uhatslocal[0], 1, &uhatsallelhold[ctc], 1);

		  ctc += pFields[i]->GetNcoeffs(k);
		  ctp += pFields[i]->GetTotPoints(k);

		}
	      
	      Vmath::Vcopy(uhatsallel.size(), uhatsallelhold, 1,  pFields[i]->UpdateCoeffs(), 1);
	      //	      pFields[i]->SetCoeffsArray(uhatsallelhold);
	      pFields[i]->BwdTrans(pFields[i]->GetCoeffs(), pFields[i]->UpdatePhys());
	      pFields[i]->SetPhysState(false);
	      // cout<<"\n diff b/w coeffs:\n";
	      // for(int k = 0; k <uhatsallel.size(); k++)
	      // 	{
	      // 	  cout<<"\n "<<uhatsallelhold[k]<<" "<<uhatsallel[k]<<" "<<uhatsallelhold[k] - uhatsallel[k]<<" ";
	      // 	}
	      Array<OneD, Array<OneD, NekDouble> >flaggedeleverify = FindIds(pFields[i]);
	      // if (m_session->GetComm()->GetRank() == 0 &&
	      //           m_session->DefinesCmdLineArgument("verbose"))
	      //         {
	      cout<<"\nfound neg after  at: \ncoord_localxyz \t val\t id \n ";
	      if(flaggedeleverify.size() > 0)
		{
		  for(int k = 0; k < flaggedeleverify.size(); k++)
		    {
		      for(int p = 0; p < flaggedeleverify[k].size(); p++)
			{
			  cout<< flaggedeleverify[k][p]<<"\t";
			}
		      cout<<"\n";
		    }
		  exit(0);
		}
	      
	     
	    } //end if negative val found
	  else
	    {
	      pFields[i]->BwdTrans(pFields[i]->GetCoeffs(), pFields[i]->UpdatePhys());
	      //pFields[i]->SetPhysState(true);
	    }
	}// fields loop
    }

    void FilterStructurePres::DoOrthonormalize(const MultiRegions::ExpListSharedPtr &pFields,
						 int flag)
    {
	
      Array<OneD, NekDouble> inarray = pFields->GetCoeffs();
      Array<OneD, NekDouble> cfs(inarray.size()),cfs2(inarray.size());
      
	Vmath::Vcopy(cfs.size(), inarray, 1, cfs, 1);

	int nelmt = pFields->GetExpSize();
        LibUtilities::BasisType btype0, btype1;
	LocalRegions::ExpansionSharedPtr exp;

	int ctc = 0;

	for (int i = 0; i < nelmt; i++)
          {
	    exp = pFields->GetExp(i);
            int nmodes0  = exp->GetBasis(0)->GetNumModes();
            int nmodes1  = exp->GetBasis(1)->GetNumModes();

	    int n_coeffs = pFields->GetNcoeffs(i);
            Array<OneD, NekDouble> coeffPerEl(n_coeffs,0.0);

	    if(exp->DetShapeType() == LibUtilities::eQuadrilateral)
	      {
		btype0 = LibUtilities::eOrtho_A;
		btype1 = LibUtilities::eOrtho_A;
	      }
	    else if(exp->DetShapeType() == LibUtilities::eTriangle)
	      {
		btype0 = LibUtilities::eOrtho_A;
                btype1 = LibUtilities::eOrtho_B;
	      }
	    
	    LibUtilities::BasisKey bkey0 (
					  exp->GetBasis(0)->GetBasisType(),
					  nmodes0,
					  exp->GetBasis(0)->GetPointsKey());
	    LibUtilities::BasisKey bkey1 (
					  exp->GetBasis(1)->GetBasisType(),
					  nmodes1,
					  exp->GetBasis(1)->GetPointsKey());
	    
            // Storage for orthogonal coefficients                                                 
            Array<OneD, NekDouble> coeffsOrth(n_coeffs,0.0);

            Vmath::Vcopy(n_coeffs, &cfs[ctc], 1, &coeffPerEl[0],1);

            LibUtilities::BasisKey bkeyOrth0(
                                             btype0, nmodes0, exp->GetBasis(0)->GetPointsKey());

            LibUtilities::BasisKey bkeyOrth1(
                                             btype1, nmodes1,   exp->GetBasis(1)->GetPointsKey());


            // Set uo basis key for orthogonal basis
	    if(flag == 1)
              {
		
                // Project from coeffs -> orthogonal coeffs                                        
                LibUtilities::InterpCoeff2D(bkey0, bkey1, coeffPerEl, bkeyOrth0, bkeyOrth1, coeffsOrth);
              }
            else //todo: asset flag is either 1 or 2 always                                        
              {
                LibUtilities::InterpCoeff2D( bkeyOrth0, bkeyOrth1, coeffPerEl, bkey0, bkey1, coeffsOrth);
        
              }

            Vmath::Vcopy(n_coeffs, &coeffsOrth[0], 1, &cfs2[ctc], 1);

	    ctc += n_coeffs;
          }
        //cout<<"\n cfs2:";                                                                        
	// for(int uu = 0; uu<cfs2.size(); uu++)                                                   
        //      cout<<cfs2[uu]<<" ";                                                               

	//        Vmath::Vcopy(inarray.size(), cfs2, 1, inarray, 1);
	Vmath::Vcopy(cfs2.size(), cfs2, 1, pFields->UpdateCoeffs(), 1);
	// for(int k = 0; k < inarray.size(); k++)
	//   {
	    
	//   }
      }


    
    void FilterStructurePres::Optimize(LocalRegions::ExpansionSharedPtr exp, Array<OneD, NekDouble> &coeffs, Array<OneD, NekDouble> &phys, Array<OneD, NekDouble> flaggedelecoord)
    {
      boost::ignore_unused(phys);
      StdExpansion *E;
      int ne, ns, dim;
      int stype = exp->DetShapeType();
      Array<OneD, Array<OneD, NekDouble> > surfaceuhats, storage, Pf, tmpcoord;
      int N1 = coeffs.size(),  Nseg = E3seg->GetNcoeffs();
      
      switch(stype)
	{
	case LibUtilities::eQuadrilateral:
	  E =  Equad;
	  ne = 4;
	  ns = 0;
	  dim = 2;
	  storage = demo.storage2dq;
	  Pf = Array<OneD, Array<OneD, NekDouble> > (ne);
	  for(int k= 0; k < ne; k++)
	    {
	      Pf[k] = Array<OneD, NekDouble>(Nseg);
	    }
	  
	  break;
	case LibUtilities::eTriangle:
	  E = Etri;
	  ne = 3;
	  ns = 0;
	  dim = 2;
	  storage = demo.storage2dt;
	  Pf = Array<OneD, Array<OneD, NekDouble> > (ne);
	  for(int k= 0; k < ne; k++)
	    {
	      Pf[k] = Array<OneD, NekDouble>(Nseg);
	    }
	  break;
	case LibUtilities::eTetrahedron:
	  E = Etet;
	  ne = 6;
	  ns = 4;
	  dim = 3;
	  storage = demo.storage3dtet;
	  surfaceuhats =  Array<OneD, Array<OneD, NekDouble> >(ns);
	  for(int k = 0; k < ns; k++)
	    {
	      surfaceuhats[k] = Array<OneD, NekDouble>(Etri->GetNcoeffs()); 
	    }
	  Pf = Array<OneD, Array<OneD, NekDouble> > (ne);
          for(int k= 0; k < ne; k++)
            {
              Pf[k] = Array<OneD, NekDouble>(Nseg);
            }

	  break;
	case LibUtilities::eHexahedron:
	  E = Ehex;
	  ne = 12;
	  ns = 6;
	  dim = 3;
	  storage = demo.storage3dhex;
	  surfaceuhats =  Array<OneD, Array<OneD, NekDouble> >(ns);
	  for(int k = 0; k < ns; k++)
	    {
	      surfaceuhats[k] =	Array<OneD, NekDouble>(Equad->GetNcoeffs());
	    }
	  Pf = Array<OneD, Array<OneD, NekDouble> > (ne);
          for(int k= 0; k < ne; k++)
            {
              Pf[k] = Array<OneD, NekDouble>(Nseg);
            }
	  break;
	case LibUtilities::ePyramid:
	  E = Epyr;
	  ne = 8;
	  ns = 5;
	  dim = 3;
	  storage = demo.storage3dpyr;
	  surfaceuhats =  Array<OneD, Array<OneD, NekDouble> >(ns);
	  surfaceuhats[0] = Array<OneD, NekDouble>(Equad->GetNcoeffs());
	  for(int k = 1; k < ns; k++)
	    {
	      surfaceuhats[k] =	Array<OneD, NekDouble>(Etri->GetNcoeffs());
	    }
	  Pf = Array<OneD, Array<OneD, NekDouble> > (ne);
          for(int k= 0; k < ne; k++)
            {
              Pf[k] = Array<OneD, NekDouble>(Nseg);
            }
	  break;
	default:
	  cout<<"\n element type not supported yet";
	  exit(0);
	  
	}
      Array<OneD, BasisType> btorth(dim);
    
      for(int k = 0; k < dim; k++)
	{
	  btorth[k] = E->GetBasis(k)->GetBasisType(); 
	}	  
      E->FwdTrans(phys , coeffs);  
      // ALWAYS transform to orthonormal:
      //      demo.OrthoNormalize(exp, coeffs, btorth, 0);
      //E->FwdTrans(phys , coeffs);
      double inf = numeric_limits<double>::infinity();

      int niter = 1e3, counter = 0;
      NekDouble tol = 1e-11, minv = inf;         // constraint specific tolerances
      
      Array<OneD, NekDouble> optima(dim), pqvalcoords(dim+1), xastarr(dim), utemp(N1), wsp1(1), Vtmp(N1);

      tmpcoord = Array<OneD, Array<OneD, NekDouble> > (dim);

      NekDouble pqval, timeprojectsurf = 0.0, timeprojectedges = 0.0;

      NekDouble avgiterGD = 0.0, avgiterhold = 0.0, roots1dtimehold = 0.0, roots2dtimehold = 0.0, roots3dtimehold = 0.0;

      NekDouble startcoordx = flaggedelecoord[0], startcoordy, startcoordz;;
      if(dim > 1)
	{
	  startcoordy = flaggedelecoord[1];
	}
      if(dim > 2)
	{
	  startcoordz = flaggedelecoord[2];
	}
      Timer t;
      while (counter <= niter)
	{
	  NekDouble roots1dtime = 0.0, roots2dtime = 0.0, roots3dtime = 0.0 ;
	  pqval = inf;
	  utemp = coeffs;
	  if (counter > -1)
	    {
	      t.Start();
	      //pq = @(xx) cfun(xx) .* (constraints(1,jj).*Vc(xx)*uhats);
	      project_edges(utemp, Pf, E);
 	      t.Stop();
	      timeprojectedges += t.TimePerTest(1);
	      if(dim > 2)
		{
		  t.Start();
		  project_surfaces(utemp, surfaceuhats, E);
		  t.Stop();
		  timeprojectsurf += t.TimePerTest(1);
		}

	      optima = call_find_roots(utemp, avgiterhold, Pf, surfaceuhats,  minv, E, roots1dtime, roots2dtime, roots3dtime);

	      //   cout<<"\n optima = "<<optima[0]<<" ,"<<optima[1]<<" minv = "<<minv<<"\n";
	      roots1dtimehold += roots1dtime;
	      roots2dtimehold += roots2dtime;
	      roots3dtimehold += roots3dtime;
	      avgiterGD += avgiterhold;
	      
	    }// end if(counter > 0)                                                                   
	  else
	    {
	      t.Start();
	      //pq = @(xx) cfun(xx) .* (constraints(1,jj).*Vc(xx)*uhats);                              

	      // project_edges(utemp, Pf, E);

	      // t.Stop();
	      // timeprojectedges += t.TimePerTest(1);
	      // if(dim > 2)
	      // 	{
	      // 	  t.Start();
	      // 	  project_surfaces( utemp, surfaceuhats, E);
	      // 	  t.Stop();
	      // 	  timeprojectsurf += t.TimePerTest(1);
	      // 	}
	      // optima = call_find_roots(utemp, avgiterhold, Pf, surfaceuhats,  minv, E, roots1dtime, roots2dtime, roots3dtime);

	      // roots1dtimehold += roots1dtime;
	      // roots2dtimehold += roots2dtime;
	      // roots3dtimehold += roots3dtime;
	      // avgiterGD += avgiterhold;

	      // Vtmp is evaluation of basis fn at optima  
	      for(int k = 0; k < dim; k++)
		{
		  tmpcoord[k] = Array<OneD, NekDouble>(1);
		}

	      tmpcoord[0][0] = startcoordx;
	      tmpcoord[1][0] = startcoordy;
	      if(dim > 2)
		{
		  tmpcoord[2][0] = startcoordz;
		}
	      optima[0] = startcoordx;
	      optima[1] = startcoordy;
	      if(dim > 2)
		{
		  optima[2] = startcoordz;
		}
	      // if(dim > 2)
	      // 	{
	      Vtmp = E->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
	      // 	}
	      // else if(dim > 1)
	      // 	{
	      // 	  Vtmp = E->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
	      // 	}
	      demo.pq(utemp, tmpcoord, Vtmp, pqvalcoords, wsp1);
	      minv = pqvalcoords[0];
	    }
	  for(int k = 0; k < dim; k++)
	    {
	      tmpcoord[k] = Array<OneD, NekDouble>(1, optima[k]);
	    }
	  Vtmp = E->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
	  demo.pq(utemp, tmpcoord, Vtmp, pqvalcoords, wsp1);
	  //	  cout<<"\n optima = "<<tmpcoord[0][0]<<" ,"<<tmpcoord[1][0]<<" minv = "<<minv<<"pqvalcoords[0 = "<<pqvalcoords[0]<<"\n";

	  if (pqvalcoords[0] < pqval)
	    {
	      xastarr = optima;
	      pqval = minv;//pqvalcoords[0];
	    }
	  //cout<<"\n pqval = "<<pqval<<" at "<< optima[0]<<" "<<optima[1];
	  // If minimum is non-negative, we're done                                                    
	  if (pqval >= -tol)
	    {
	      break;
	    }

	  Array<OneD, NekDouble>  Vastsq (N1);
	  Array<OneD, NekDouble>  Vast (N1);
	  Array<OneD, Array<OneD, NekDouble> > xastarrofarr(dim);
	  for(int k = 0; k < dim; k++)
	    {
	      xastarrofarr[k] = Array<OneD, NekDouble>(1, xastarr[k]);  
	    }
	  Array<OneD, NekDouble> tmp;
	  NekDouble vastsqsum;

	  for( int i = 0; i < N1; i++)
	    {
	      // if(dim > 2)
	      // 	{
	      tmp = E->PhysEvaluateBasis(xastarrofarr, storage, i);
	      //		}
	      // else if(dim > 1)
	      // 	{
	      // 	  tmp = E->PhysEvaluateBasis(xastarrofarr, demo.storage2d, i);
	      // 	}
	      Vast[i] = tmp[0];
	      Vastsq[i] = (Vast[i]*Vast[i]);

	    }

	  vastsqsum = Vmath::Vsum(N1, &Vastsq[0], 1);

	  Array<OneD, NekDouble>  qast(N1);

	  for(int i = 0; i<N1; i++)
	    {
	      qast[i] = ((1/sqrt(vastsqsum))*(Vast[i]));
	    }
	  Vmath::Smul(N1, pqval, &qast[0], 1, &qast[0], 1);

	  Vmath::Vsub(utemp.size(), &utemp[0], 1, &qast[0], 1, &qast[0], 1);
	  coeffs = qast;

	  counter = counter + 1;


	}
      //doorthonormalize back to original basis
      //demo.OrthoNormalize(exp, coeffs, btorth, 1);
      E->BwdTrans(coeffs, phys);
      //      E->Bwd
      roots1dtimehold = roots1dtimehold/(counter-1);
      roots2dtimehold = roots2dtimehold/(counter-1);
      roots3dtimehold = roots3dtimehold/(counter-1);
      timeprojectedges = timeprojectedges/(counter);
      timeprojectsurf = timeprojectsurf/(counter);
      //      NekDouble itersGD2 = (avgiterGD)/(counter);
      int iterstaken = counter;
cout<<"\nsphere_rotation took "<<iterstaken<<"iterations";
	  
      if (m_session->GetComm()->GetRank() == 0 &&
          m_session->DefinesCmdLineArgument("verbose"))
        {
	  cout<<"\nAvg times per iter (1drootfinder,2drootfinder,3drootfinder) = : "  << roots1dtimehold<<     ", "<<roots2dtimehold<<", "<<roots3dtimehold ;//<< timeprojectedges<<", "<<timeprojectsurf<<",     return coeffprev;
	}
    }
    
    
    void FilterStructurePres::project_edges( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret, StdExpansion *E)
    {
      if(E->DetShapeType()  == LibUtilities::eQuadrilateral) // quad
	{
	  // bot edge
	  edgederpquhats(uhats, ret[0], E->GetNcoeffs(),  Vym1q, Vdxym1q, Vdyym1q, NullNekDouble1DArray);
	  // right edge
	  edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vx1q, Vdxx1q, Vdyx1q, NullNekDouble1DArray);
	  // top edge
	  edgederpquhats(uhats, ret[2], E->GetNcoeffs(), Vy1q, Vdxy1q, Vdyy1q, NullNekDouble1DArray);
	  // left edge
	  edgederpquhats(uhats, ret[3], E->GetNcoeffs(),  Vxm1q, Vdxxm1q, Vdyxm1q, NullNekDouble1DArray);
	  
	}
      if(E->DetShapeType() == LibUtilities::eTriangle) // tri
	{	     
	  // bot edge  
	  edgederpquhats(uhats, ret[0], E->GetNcoeffs(),  Vym1t, Vdxym1t, Vdyym1t,NullNekDouble1DArray);
	  
	  // left edge
	  edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vxm1t, Vdxxm1t, Vdyxm1t, NullNekDouble1DArray);
	  
	  // hypto
	  edgederpquhats(uhats, ret[2], E->GetNcoeffs(),  Vxyhypt, Vdxxyhypt, Vdyxyhypt, NullNekDouble1DArray);
	        
	  return;
	}
      
      if(E->DetShapeType() == LibUtilities::eHexahedron) // hex
	{

	  // edge front left (x = -1) (y = -1)                                                         
	  edgederpquhats(uhats, ret[0], E->GetNcoeffs(), Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);
	  //edge front right (x = 1) (y = -1)                                                          
	  edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vx1ym1z, Vdxx1ym1z, Vdyx1ym1z, Vdzx1ym1z);

	  //edge front top (y = -1) (z = 1)                                                            
	  edgederpquhats(uhats, ret[2], E->GetNcoeffs(),  Vym1xz1, Vdxym1xz1, Vdyym1xz1, Vdzym1xz1);

	  //edge front bot (y = -1) (z = -1)                                                           
	  edgederpquhats(uhats, ret[3], E->GetNcoeffs(), Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);

	  // edge back left (y = 1), (x = -1)                                                          
	  edgederpquhats(uhats, ret[4], E->GetNcoeffs(), Vxm1y1z, Vdxxm1y1z, Vdyxm1y1z, Vdzxm1y1z);

	  //edge back right (x = 1), (y = 1))                                                          
	  edgederpquhats(uhats, ret[5], E->GetNcoeffs(),  Vx1y1z, Vdxx1y1z, Vdyx1y1z, Vdzx1y1z);

	  //edge back top ( y = 1) (z = 1)                                                             
	  edgederpquhats(uhats, ret[6], E->GetNcoeffs(),  Vy1xz1, Vdxy1xz1, Vdyy1xz1, Vdzy1xz1);

	  //edge back bot (y = 1) (z = -1))                                                            
	  edgederpquhats(uhats, ret[7], E->GetNcoeffs(),  Vy1xzm1, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);

	  // edge left bot (z = -1), (x = -1)                                                          
	  edgederpquhats(uhats, ret[8], E->GetNcoeffs(), Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);
  
	  //edge left top (x = -1), (z = 1))                                                           
	  edgederpquhats(uhats, ret[9], E->GetNcoeffs(),   Vxm1yz1,Vdxxm1yz1, Vdyxm1yz1, Vdzxm1yz1);
  
	  //edge right bot ( z = -1) (x = 1)                                                           
	  edgederpquhats(uhats, ret[10], E->GetNcoeffs(),   Vx1yzm1, Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);


	  //edge right top (z  1) (x  1))                                                              
	  edgederpquhats(uhats, ret[11], E->GetNcoeffs(),Vx1yz1, Vdxx1yz1, Vdyx1yz1, Vdzx1yz1);

  
	}
      else if(E->DetShapeType() == LibUtilities::eTetrahedron) //tet                                                                     
	{
	  // edge front left (AD) (x = -1) (y = -1)                                                    
	  edgederpquhats(uhats, ret[0], E->GetNcoeffs(),  Vxm1ym1ztet, Vdxxm1ym1ztet, Vdyxm1ym1ztet, Vdzxm1ym1ztet);

	  //edge front hypt (DB) (y = -1) (z = -x)                                                     
	  edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vym1xmztet, Vdxym1xmztet, Vdyym1xmztet, Vdzym1xmztet);

	  //edge front bot (AB) (y = -1) (z = -1)                                                      
	  edgederpquhats(uhats, ret[2], E->GetNcoeffs(), Vym1xzm1tet, Vdxym1xzm1tet, Vdyym1xzm1tet, Vdzym1xzm1tet );

	  //edge left hypt (DC) ( x = -1) (z = -y)                                                     
	  edgederpquhats(uhats, ret[3], E->GetNcoeffs(), Vxm1ymztet, Vdxxm1ymztet, Vdyxm1ymztet, Vdzxm1ymztet);

	  // edge bot diag (BC) (z = -1) (y = -x)                                                      
	  edgederpquhats(uhats, ret[4], E->GetNcoeffs(),  Vxmyzm1tet, Vdxxmyzm1tet, Vdyxmyzm1tet, Vdzxmyzm1tet);

	  //edge CA bot left (x = -1) (z = -1)                                                         
	  edgederpquhats(uhats, ret[5], E->GetNcoeffs(),Vxm1yzm1tet, Vdxxm1yzm1tet,  Vdyxm1yzm1tet,  Vdzxm1yzm1tet);


	}
      else if(E->DetShapeType() == LibUtilities::ePyramid) //pyr                                                                     
	{

	  // edge front left (x = -1) (y = -1) EA                                                      
	  edgederpquhats(uhats, ret[0], E->GetNcoeffs(), Vxm1ym1zpyr, Vdxxm1ym1zpyr, Vdyxm1ym1zpyr, Vdzxm1ym1zpyr);

	  //edge front hypt (y = -1) (z = -x) EB                                                       
	  edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vym1xmzpyr, Vdxym1xmzpyr, Vdyym1xmzpyr, Vdzym1xmzpyr);

	  //edge back hypt (z = -x  or x = -z) (x = y) EC                                              
	  edgederpquhats(uhats, ret[2], E->GetNcoeffs(), Vxeyxmzpyr, Vdxxeyxmzpyr, Vdyxeyxmzpyr, Vdzxeyxmzpyr);

	  //edge back left (y = -z) (x = -1) ED                                                        
	  edgederpquhats(uhats, ret[3], E->GetNcoeffs(), Vxeyxmzpyr, Vdxxeyxmzpyr, Vdyxeyxmzpyr, Vdzxeyxmzpyr);

	  //edge front bot (y = -1) (z = -1) AB                                                        
	  edgederpquhats(uhats, ret[4], E->GetNcoeffs(), Vym1xzm1pyr, Vdxym1xzm1pyr, Vdyym1xzm1pyr, Vdzym1xzm1pyr);

	  //edge back bot (y = 1) (z = -1)) DC                                                         
	  edgederpquhats(uhats, ret[5], E->GetNcoeffs(), Vy1xzm1pyr, Vdxy1xzm1pyr, Vdyy1xzm1pyr, Vdzy1xzm1pyr);

	  // edge left bot (z = -1), (x = -1) AD                                                       
	  edgederpquhats(uhats, ret[6], E->GetNcoeffs(), Vxm1yzm1pyr, Vdxxm1yzm1pyr, Vdyxm1yzm1pyr, Vdzxm1yzm1pyr);

	  //edge right bot ( z = -1) (x = 1) BC                                                        
	  edgederpquhats(uhats, ret[7], E->GetNcoeffs(), Vx1yzm1pyr, Vdxx1yzm1pyr, Vdyx1yzm1pyr, Vdzx1yzm1pyr);
	}
    }

    void FilterStructurePres::edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  int modes,  Array <OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2)
    {
      int uhatstot = uhats.size();
      int totpts =  E3seg->GetBasis(0)->GetZ().size();

      Array<OneD, NekDouble> temp(totpts), temp2(modes), temp3(uhatstot);
      Array<OneD, NekDouble> pqeval(totpts);
      NekDouble v1, v2;

      for(int i = 0; i<totpts; i++)
	{
	  Vmath::Vmul(modes, &Vxy[0]+i, totpts, &Vxyd0[0]+i, totpts, &temp2[0], 1);
	  v1  = Vmath::Vsum(modes, temp2, 1);
	  Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);
	  v2  = Vmath::Vsum(uhatstot, temp3, 1);

	  v1 = v2*v1;

	  // At this point,                                                                            
	  // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point                                          

	  // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)                                                    
	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &temp2[0], 1);

	  v2  = Vmath::Vsum(uhatstot, temp2, 1);

	  Vmath::Vmul(uhats.size(), &Vxyd0[i], totpts, &uhats[0], 1, &temp2[0], 1);

	  v1= v2*Vmath::Vsum(uhats.size(), temp2, 1)- v1;

	  pqeval[i] = v1;

	  if(dimension > 1)
	    {
	      Vmath::Vmul(modes, &Vxy[0]+i, totpts, &Vxyd1[0]+i, totpts, &temp2[0], 1);
	      v1  = Vmath::Vsum(modes, temp2, 1);
	      Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);

	      v2  = Vmath::Vsum(uhatstot, temp3, 1);
	      v1 = v2*v1;

	      // At this point,                                                                        
	      // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point                                      

	      // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)                                                
	      Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &temp2[0], 1);



	      v2  = Vmath::Vsum(uhatstot, temp2, 1);
	      Vmath::Vmul(uhats.size(), &Vxyd1[i], totpts, &uhats[0], 1, &temp2[0], 1);

	      v1= v2*Vmath::Vsum(uhats.size(), temp2, 1) - v1;
	      pqeval[i] += v1;
  
	    }
	  if(dimension == 3)
	    {
	      Vmath::Vmul(modes, &Vxy[0]+i, totpts, &Vxyd2[0]+i, totpts, &temp2[0], 1);
	      v1  = Vmath::Vsum(modes, temp2, 1);
	      Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);

	      v2  = Vmath::Vsum(uhatstot, temp3, 1);
	      v1 = v2*v1;

	      // At this point,                                                                        
	      // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point                                      

	      // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)                                                
	      Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &temp2[0], 1);

	      v2  = Vmath::Vsum(uhatstot, temp2, 1);
	      Vmath::Vmul(uhats.size(), &Vxyd2[i], totpts, &uhats[0], 1, &temp2[0], 1);

	      v1= v2*Vmath::Vsum(uhats.size(), temp2, 1) - v1;
	      pqeval[i] += v1;

	    }
	}

      E3seg->FwdTrans(pqeval, ret);
   
    }

    void FilterStructurePres::project_surfaces(Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret,  StdExpansion *E)
    {

      if(E->DetShapeType() == LibUtilities::eHexahedron) //hex                                                                       
	{
	  //bot surface  z = -1                                                                        
	  surfuhats(uhats, ret[0], Vxyzm1, Equad);

	  //right surface x = 1                                                                        
	  surfuhats(uhats, ret[1], Vx1yz, Equad);

	  //surface top z = 1                                                                          
	  surfuhats(uhats, ret[2], Vxyz1, Equad);

	  //left surface x = -1                                                                        
	  surfuhats(uhats, ret[3], Vxm1yz, Equad);

	  //surface front y = -1                                                                       
	  surfuhats(uhats, ret[4], Vxym1z, Equad);

	  //surface back y = 1                                                                         
	  surfuhats(uhats, ret[5], Vxy1z, Equad);
	}
      else if(E->DetShapeType() == LibUtilities::eTetrahedron) //tets                                                                 
	{
	  //surface bot z = -1 restriction on GD: x+y = 0 (ABC)                                        
	  surfuhats(uhats, ret[0],  Vxyzm1tet, Etri);

	  //surface left x = -1 restriction on GD: y + z = 0 (DAC)                                     
	  surfuhats(uhats, ret[1], Vxm1ypz0tet, Etri);

	  //surf front y = -1   restriction on GD: x + z = 0 (DAB)                                     
	  surfuhats(uhats, ret[2], Vxpz0ym1tet, Etri);

	  //surf DCB restriction  on GD: (x + y + z = -1)                                              
	  surfuhats(uhats, ret[3], Vxpypzm1tet, Etri);


	}
      else if(E->DetShapeType() == LibUtilities::ePyramid)
	{
	  //surface bot z = -1                                                                         
	  surfuhats(uhats, ret[0], Vxyzm1pyr, Equad);

	  //surface hypt x+z = 0                                                                       
	  surfuhats(uhats, ret[1], Vxmzypyr, Etri);

	  //surface left x = -1                                                                      
	  surfuhats(uhats, ret[2], Vxm1yzpyr, Etri);

	  //surface front y = -1                                                                     
	  surfuhats(uhats, ret[3], Vxym1zpyr, Etri);

	  //surface back y +z =  0                                                                   
	  surfuhats(uhats, ret[4], Vxymzpyr, Etri);
	}
    }
    
    void FilterStructurePres::surfuhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble>  Vxyz, StdExpansion *Etemp)
    {
      int modes = uhats.size();

      Array<OneD, NekDouble> temp(uhats.size());
      int totpts = Etemp->GetTotPoints();
      Array<OneD, NekDouble> vals(totpts);
      //Vxyz*uhats -> project to -> E3tri or Equad                                                     
      for(int k = 0; k < totpts; k++)
	{

	  Vmath::Vmul(uhats.size(), &Vxyz[k], totpts, &uhats[0], 1, &temp[0], 1);

	  vals[k]  = Vmath::Vsum(modes, temp, 1);
  
	}
      //    vals = demo.blasmatvec(Vxyz, uhats, totpts, modes);                                        
  
      Etemp->FwdTrans(vals, ret);
    }

    Array<OneD, NekDouble>  FilterStructurePres::call_find_roots(Array<OneD,  NekDouble> &uhats , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats, NekDouble &minv, StdExpansion *E, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime)
    {

      //      NekDouble dum;
      int dimension = 0;
      Array<OneD, NekDouble> temp(uhats.size());
      roots1dtime = 0.0;
      roots2dtime = 0.0;
      roots3dtime = 0.0;
      NekDouble  inf = numeric_limits<double>::infinity();
      minv = inf;

      NekDouble avgiterGDhold = 0;
      Timer t;
      Array<OneD, Array<OneD, NekDouble > >  rethold, storage;
      vector<vector<NekDouble> > retall;

      // EDGES
      if(E->DetShapeType() == LibUtilities::eQuadrilateral)
	{
	  dimension = 2;
	  retall = vector<vector<NekDouble> >(dimension); 
	  int numedges = 4;
	  //	  NekDouble dum;
	  for(int i = 0; i < numedges; i++)
	    {
	      storage = demo.storage2dq;
	      t.Start();
	      rethold  = demo.call_companion_rf(uhatsedges[i], C);//(demo.find_roots(uhatsedges[i], E, demo.storage2dq, dum,  0, 0, 0, 0)) ;
	      t.Stop();
	      roots1dtime +=  t.TimePerTest(1);
	      for(int p = 0; p < rethold[0].size(); p++)
		{

		  switch(i)
		    {
		    case 0:
		      //if(i == 0) // edge bot (y = -1)                                               
		      retall[1].push_back(   (-1.0));
		      retall[0].push_back( (rethold[0][p]));
		      break;
		    case 1:
		      //if(i == 1) // edge right (x = 1)                                              
		      retall[0].push_back(   (1.0));
		      retall[1].push_back( (rethold[0][p]));
		      break;
		    case 2:
		      //if(i == 2) // edge top (y = 1)                                                
		      retall[1].push_back(   (1.0));
		      retall[0].push_back( (rethold[0][p]));
		      break;
		    case 3:
		      //if(i == 3) // edge left (x = -1)                                              
		      retall[0].push_back(   (-1.0));
		      retall[1].push_back( (rethold[0][p]));
		      break;
		    default:
		      break;
		    }
		}
	    }
	  // add 4 corners                                                                             
	  retall[0].push_back(-1);
	  retall[1].push_back(-1);
	  
	  retall[0].push_back(1);
	  retall[1].push_back(1);
	  
	  retall[0].push_back(-1);
	  retall[1].push_back(1);
	  
	  retall[0].push_back(1);
	  retall[1].push_back(-1);

	  rethold = NullNekDoubleArrayOfArray;

	  t.Start();
	  demo.steepestgradient_descent2Dquad(uhats, Equad, rethold, avgiterGDhold); 
	  roots2dtime +=  t.TimePerTest(1);
	  avgiterGD += avgiterGDhold;
	  if(rethold != NullNekDoubleArrayOfArray)
	    {
	      retall[0].push_back(rethold[0][0]);
	      retall[1].push_back(rethold[1][0]);
	      
	    }
		  
	}
      else if(E->DetShapeType() == LibUtilities::eTriangle)
	{
	  dimension = 2;

	  retall = vector<vector<NekDouble> >(dimension);
	  storage = demo.storage2dt;
	    
	  int numedges = 3;
	  //	  NekDouble dum = 0;

	  for(int i = 0; i < numedges; i++)
	    {
	      // cout<<"\ni = "<<i<<"\n\n";      
	      // for(int k = 0; k < uhatsedges[i].size(); k++)
	      // 	cout<<" "<<uhatsedges[i][k];
	      // cout<<"\n";
	      t.Start();
	      rethold  = demo.call_companion_rf(uhatsedges[i], C);
	      t.Stop();
	      roots1dtime +=  t.TimePerTest(1);
	      //	      cout<<"\n done roots edge:"<<i<<" ";
	      if(rethold != NullNekDoubleArrayOfArray)
		{
		  for(int p = 0;  p < rethold[0].size(); p++)
		    {
		      switch(i)
			{
			case 0:
			  //if(i == 0) // edge bot (y = -1)
			  //cout<<"\n in case 0\n";
			  retall[1].push_back((-1.0));
			  retall[0].push_back( (rethold[0][p]));
			  break;
			case 1:
			  //if(i == 1) // edge left (x = -1)                                              
			  retall[0].push_back((-1.0));
			  retall[1].push_back( (rethold[0][p]));

			  break;
			case 2:
			  //if(i == 2) // edge hypt (y = -x)                                              
			  retall[1].push_back((-rethold[0][p]));
			  retall[0].push_back((rethold[0][p]));

			  break;
			default:
			  break;
			}
		    }

		}
	    }
          // add 3 corners                                                                             
	  retall[0].push_back(-1);
	  retall[1].push_back(-1);

	  retall[0].push_back(-1);
	  retall[1].push_back(1);

	  retall[0].push_back(1);
	  retall[1].push_back(-1);

	  rethold = NullNekDoubleArrayOfArray;
	  t.Start();

	  demo.steepestgradient_descent2Dtri(uhats, Etri, rethold, avgiterGDhold); 
	  roots2dtime +=  t.TimePerTest(1);
	  avgiterGD += avgiterGDhold;
	  if(rethold != NullNekDoubleArrayOfArray)
	    {
	      retall[0].push_back(rethold[0][0]);
	      retall[1].push_back(rethold[1][0]);
	      
	    }
	}
      else if(E->DetShapeType() == LibUtilities::eTetrahedron) // tet
	{
	  dimension = 3;

	  retall = vector<vector<NekDouble> >(dimension);
	  int numedges = 6;
	  storage = demo.storage3dtet;
	  for(int i = 0; i < numedges; i++)
	    {
	      t.Start();
	      rethold  = demo.call_companion_rf(uhatsedges[i], C);//(demo.find_roots(uhatsedges[i], E, NullNekDoubleArrayOfArray, dum,  0, 0, 0, 0));
	      t.Stop();

	      roots1dtime +=  t.TimePerTest(1);
	      for(int p = 0; p < rethold[0].size(); p++)
		{
		  switch(i)
		    {
		    case 0:

		      retall[0].push_back(-1);
		      retall[1].push_back(-1);
		      retall[2].push_back(rethold[0][p]);
		      break;
		    case 1:  //edge front hypt (DB) (y = -1) (z = -x)                                  
		      retall[0].push_back(rethold[0][p]);
		      retall[2].push_back(-rethold[0][p]);
		      retall[1].push_back(-1);

		      retall[0].push_back(-rethold[0][p]);
		      retall[2].push_back(rethold[0][p]);
		      retall[1].push_back(-1);
		      break;
		    case 2: //edge front bot (AB) (y = -1) (z = -1)                                    


		      retall[0].push_back(rethold[0][p]);
		      retall[1].push_back(-1);
		      retall[2].push_back(-1);

  
		      break;
		    case 3:    //edge left hypt (DC) ( x = -1) (z = -y)                                
  
		      retall[0].push_back(-1);
		      retall[1].push_back(rethold[0][p]);
		      retall[2].push_back(-rethold[0][p]);
  
		      retall[0].push_back(-1);
		      retall[1].push_back(-rethold[0][p]);
		      retall[2].push_back(rethold[0][p]);
		      break;
		    case 4:   // edge bot diag (BC) (z = -1) (y = -x)                                  

		      retall[2].push_back(-1);
		      retall[1].push_back(-rethold[0][p]);
		      retall[0].push_back(rethold[0][p]);

		      retall[2].push_back(-1);
		      retall[1].push_back(rethold[0][p]);
		      retall[0].push_back(-rethold[0][p]);

		      break;
		    case 5:  //edge CA bot left (x = -1) (z = -1)                                      

		      retall[1].push_back(rethold[0][p]);
		      retall[0].push_back(-1);
		      retall[2].push_back(-1);
		      break;
		    default:
		      cout<<"\n wrong edge id\n";
		      exit(0);

		    }

		}


	    }
	}
      else if(E->DetShapeType() == LibUtilities::eHexahedron) //hex
	{
	  dimension = 3;

	  retall = vector<vector<NekDouble> >(dimension);
	  storage = demo.storage3dhex;
	  int numedges = 12;
	  for(int i = 0; i < numedges; i++)
	    {

	      t.Start();
	      rethold  = demo.call_companion_rf(uhatsedges[i], C);//(demo.find_roots(uhatsedges[i], E, NullNekDoubleArrayOfArray, dum,  0, 0, 0, 0));
	      t.Stop();

	      roots1dtime +=  t.TimePerTest(1);

	      for(int p = 0; p < rethold[0].size(); p++)
		{
		  switch(i)
		    {
		    case 0: // edge front left (x = -1) (y = -1)                                       
		      {
			retall[0].push_back(   (-1.0));
			retall[1].push_back( (-1.0));
			retall[2].push_back( (rethold[0][p]));
		      }
		      break;
		    case 1: //edge front right (x = 1) (y = -1)                                        
		      {
			retall[0].push_back(  (1.0));
			retall[1].push_back( (-1.0));
			retall[2].push_back( (rethold[0][p]));

		      }
		      break;
		    case 2: //edge front top (y = -1) (z = 1)                                          
		      {
			retall[0].push_back(  (rethold[0][p]));
			retall[1].push_back( (-1.0));
			retall[2].push_back( (1.0));
		      }
		      break;
		    case 3: //edge front bot (y = -1) (z = -1)                                         
		      {
			retall[0].push_back( (rethold[0][p]));
			retall[1].push_back( (-1.0));
			retall[2].push_back( (-1.0));
		      }
		      break;
		    case 4: //edge back left (y = 1), (x = -1)                                         
		      {
			retall[0].push_back(  (-1.0));
			retall[1].push_back( (1.0));
			retall[2].push_back( (rethold[0][p]));
		      }
		      break;
		    case 5: //edge back right (x = 1), (y = 1)                                         
		      {
			retall[0].push_back(  (1.0));
			retall[1].push_back( (1.0));
			retall[2].push_back( (rethold[0][p]));
		      }
		      break;
		    case 6: //edge back top ( y = 1) (z = 1)                                           
		      {
			retall[0].push_back(  (rethold[0][p]));
			retall[1].push_back( (1.0));
			retall[2].push_back( (1.0));

		      }
		      break;
		    case 7: //edge back bot (y = 1) (z = -1)                                           
		      {
			retall[0].push_back( (rethold[0][p]));
			retall[1].push_back( (1.0));
			retall[2].push_back( (-1.0));
		      }
		      break;
		    case 8: //edge left bot (z = -1), (x = -1)                                         
		      {
			retall[0].push_back(  (-1.0));
			retall[1].push_back( (rethold[0][p]));
			retall[2].push_back( (-1.0));
		      }
		      break;
		    case 9: //edge left top (x = -1), (z = 1)                                          
		      {
			retall[0].push_back(  (-1.0));
			retall[1].push_back( (rethold[0][p]));
			retall[2].push_back( (1.0));

		      }
		      break;
		    case 10: //edge right bot (z = -1) (x = 1)                                         
		      {
			retall[0].push_back(  (1.0));
			retall[1].push_back( (rethold[0][p]));
			retall[2].push_back( (-1.0));
		      }
		      break;
		    case 11: //edge right top (z  1) (x  1)
		      
		      {
			retall[0].push_back(  (1.0));
			retall[1].push_back( (rethold[0][p]));
			retall[2].push_back( (1.0));;
		      }
		      break;
		    default:
		      {
			cout<<"\n edge number not valid!\n\n";
			exit(0);
		      }
		    }
		}

	    }
	}
      else if(E->DetShapeType() == LibUtilities::ePyramid) //pyr                                                                     
	{
	  dimension = 3;

	  retall = vector<vector<NekDouble> >(dimension);
	  storage = demo.storage3dpyr;
	  int numedges = 8;
	  for(int i = 0; i < numedges; i++)
	    {
	      t.Start();
	      rethold  = demo.call_companion_rf(uhatsedges[i], C);//(demo.find_roots(uhatsedges[i], E, NullNekDoubleArrayOfArray, dum,  0, 0, 0,  0));
	      t.Stop();

	      roots1dtime +=  t.TimePerTest(1);
	      for(int p = 0; p < rethold[0].size(); p++)
		{
		  if(i == 0)
		    {
		      retall[0].push_back(-1);
		      retall[1].push_back(-1);
		      retall[2].push_back(rethold[0][p]);
		    }
		  else if(i == 1) //edge front hypt (y = -1) (z = -x)                   
		    {
		      retall[0].push_back(rethold[0][p]);
		      retall[2].push_back(-rethold[0][p]);
		      retall[1].push_back(-1);

		      retall[0].push_back(-rethold[0][p]);
		      retall[2].push_back(rethold[0][p]);
		      retall[1].push_back(-1);
		    }
                  else if(i == 2) //edge back hypt ( y = -z) (z +x = 0)                           
		    {
		      retall[0].push_back(rethold[0][p]);
		      retall[1].push_back(rethold[0][p]);
		      retall[2].push_back(-rethold[0][p]);

		      retall[0].push_back(-rethold[0][p]);
		      retall[1].push_back(-rethold[0][p]);
		      retall[2].push_back(rethold[0][p]);

		    }
		  else if(i == 3)  //edge back left (y = -z) (x = -1) ED                         
		    {
		      retall[0].push_back(-1);
		      retall[1].push_back(-rethold[0][p]);
		      retall[2].push_back(rethold[0][p]);

		      retall[0].push_back(-1);
		      retall[1].push_back(rethold[0][p]);
		      retall[2].push_back(-rethold[0][p]);
		    }
		  else if(i == 4)  //edge front bot (y = -1) (z = -1)                               
		    {
		      retall[1].push_back(-1);
		      retall[0].push_back(rethold[0][p]);
		      retall[2].push_back(-1);
		    }
		  else if(i == 5) //edge back bot (y = 1) (z = -1)                                  
		    {
		      retall[1].push_back(1);
		      retall[0].push_back(rethold[0][p]);
		      retall[2].push_back(-1);
		    }
		  else if(i == 6)//edge left bot (z = -1), (x = -1)                                 
		    {
		      retall[0].push_back(-1);
		      retall[1].push_back(rethold[0][p]);
		      retall[2].push_back(-1);
		    }
		  else if(i == 7) //edge right bot ( z = -1) (x = 1)                                
		    {
		      retall[2].push_back(-1);
		      retall[1].push_back(rethold[0][p]);
		      retall[0].push_back(1);
		    }

		}
	    }
}	

      //SURFACES:                                                                                      
      if(E->DetShapeType()  == LibUtilities::eTetrahedron) //tet
	{
	  int numsurfaces = 4;
	  NekDouble avgiterGDhold;

	  // call 2D rootfinder on each surface:                                                       
	  for(int i = 0; i < numsurfaces; i++)
	    {
	      t.Start();
	      demo.steepestgradient_descent2Dtri(surfaceuhats[i], Etri,  rethold, avgiterGDhold);//demo.find_roots(surfaceuhats[i], Etri, demo.demo.storage2dt,  avgiterGDhold, 0, 1, 0, 1);
	      // last arg = 1 coz all 2d ele are triangles                                             

	      t.Stop();
	      roots2dtime +=  t.TimePerTest(1);
	      avgiterGD += avgiterGDhold;

	      switch(i)
		{
		case 0:
		  //surface bot z = -1, x + y = 0 (ABC)                                                
		  retall[0].push_back( (rethold[0][0]));
		  retall[1].push_back( (rethold[1][0]));
		  retall[2].push_back( (-1.0));
		  break;
		case 1:
		  //surface left x = -1 y + z = 0 (DAC)                                                
		  retall[0].push_back(-1);
		  retall[2].push_back(rethold[1][0]);
		  retall[1].push_back(rethold[0][0]);
		  break;
		case 2:
		  //surf front y = -1, x + z = 0 (DAB)                                                 
		  retall[2].push_back(rethold[1][0]);
		  retall[1].push_back(-1);
		  retall[0].push_back(rethold[0][0]);
		  break;
		case 3:
		  //surf DCB (x + y + z = -1)                                                          
		  retall[2].push_back(-1.0 - rethold[1][0] - rethold[0][0]);
		  retall[1].push_back(rethold[1][0]);
		  retall[0].push_back(rethold[0][0]);
		  break;
		default: cout<<"\n invalid surface id\n";
		  exit(0);
		}
	    }
	  //3d gd:
	  t.Start();
	  demo.steepestgradientdescent3D(uhats, Etet, demo.storage3dtet, demo.coordtet, demo.coordmidtet, demo.midptevaltet, rethold, avgiterGDhold);//(demo.find_roots(uhats, E, demo.storage3d,  avgiterGDhold, 0, 0, 1, 0)) ;
	  t.Stop();
	  roots3dtime +=  t.TimePerTest(1);

	  avgiterGD += avgiterGDhold;
  
	  for(int k = 0; k < dimension; k++)
	    {
	      if(rethold[k][0] <inf)
		{
		  retall[k].push_back(rethold[k][0]);
		}
	    }

	}
      else if(E->DetShapeType()  == LibUtilities::eHexahedron) //hex 
	{
	  int numsurfaces = 6;
	  NekDouble avgiterGDhold;
	  for(int i = 0; i < numsurfaces; i++)
	    {
	      t.Start();

	      demo.steepestgradient_descent2Dquad(surfaceuhats[i], Equad, rethold, avgiterGDhold);
	      // rethold = demo.find_roots(surfaceuhats[i], Equad, demo.storage2dq,  avgiterGDhold, 0, 1, 0 , 0);
	      t.Stop();

	      roots2dtime +=  t.TimePerTest(1);
	      avgiterGD += avgiterGDhold;

	      if(i == 0) // surf bot (z = -1)                                                         
		{
		  retall[0].push_back( (rethold[0][0]));
		  retall[1].push_back( (rethold[1][0]));
		  retall[2].push_back( (-1.0));
		}
	      else if(i == 1) //surf right (x = 1)                                                    
		{
		  retall[0].push_back  (1.0);
		  retall[1].push_back (rethold[0][0]);
		  retall[2].push_back (rethold[1][0]);

		}
	      else if(i == 2) //surf top (z = 1)                                                      
		{
		  retall[0].push_back  (rethold[0][0]);
		  retall[1].push_back (rethold[1][0]);
		  retall[2].push_back (1.0);
		}
	      else if(i == 3) //surf left (x = -1)                                                    
		{
		  retall[0].push_back  (-1.0);
		  retall[1].push_back (rethold[0][0]);
		  retall[2].push_back (rethold[1][0]);
		}
	      else if(i == 4) //surf front (y = -1)                                                   
		{
		  retall[0].push_back  (rethold[0][0]);
		  retall[1].push_back (-1.0);
		  retall[2].push_back (rethold[1][0]);
		}
	      else if(i == 5) //surf back (y = 1)                                                     
		{
		  retall[0].push_back  (rethold[0][0]);
		  retall[1].push_back (1.0);
		  retall[2].push_back (rethold[1][0]);
		}

	    }

	  //3d gd:
	  t.Start();
	  demo.steepestgradientdescent3D(uhats, Ehex, demo.storage3dhex, demo.coordhex, demo.coordmidhex, demo.midptevalhex, rethold, avgiterGDhold);//(demo.find_roots(uhats, E, storage3d,  avgiterGDhold, 0, 0, 1, 0)) ;
	  t.Stop();
	  roots3dtime +=  t.TimePerTest(1);

	  avgiterGD += avgiterGDhold;
	  
	  for(int k = 0; k < dimension; k++)
	    {
	      if(rethold[k][0] <inf)
		{
		  retall[k].push_back(rethold[k][0]);
		}
	    }

	}
      else if(E->DetShapeType()  == LibUtilities::ePyramid) //pyr
	{
	  int numsurfaces = 5;
	  NekDouble avgiterGDhold;
	  for(int i = 0; i < numsurfaces; i++)
	    {
	      if( i == 0) // bot surface: quad, z = -1                                                
		{
		  t.Start();
		  demo.steepestgradient_descent2Dquad(surfaceuhats[i], Equad, rethold, avgiterGDhold);

		    //demo.find_roots(surfaceuhats[i], Equad, demo.storage2dq,  avgiterGDhold, 0, 1, 0 , 0);
		  t.Stop();

		  roots2dtime +=  t.TimePerTest(1);
		  avgiterGD += avgiterGDhold;
 
		  retall[2].push_back(-1);
		  retall[1].push_back(rethold[1][0]);
		  retall[0].push_back(rethold[0][0]);
		  
		}
	      else  // rest surfaces: tri                                                              
		{
		  t.Start();
		  demo.steepestgradient_descent2Dtri(surfaceuhats[i], Etri, rethold, avgiterGDhold);

		  // rethold = demo.find_roots(surfaceuhats[i], Etri, demo.storage2dt,  avgiterGDhold, 0, 1, 0 , 0);
		  t.Stop();
		  roots2dtime +=  t.TimePerTest(1);
		  avgiterGD += avgiterGDhold;
		  
		  switch(i)
		    {
		    case 1: // x+z = 0                                                                 
		      retall[0].push_back(rethold[0][0]);
		      retall[1].push_back(rethold[1][0]);
		      retall[2].push_back(-rethold[0][0]);
		      break;
		    case 2: // x = -1                                                                  
		      retall[2].push_back(rethold[1][0]);
		      retall[1].push_back(rethold[0][0]);
		      retall[0].push_back(-1.0);
		      break;
		    case 3: //y = -1                                                                   
		      retall[2].push_back(rethold[1][0]);
		      retall[1].push_back(-1.0);
		      retall[0].push_back(rethold[0][0]);
		      break;
		    case 4: // y+z = 0                                                                 
		      retall[2].push_back(-rethold[1][0]);
		      retall[1].push_back(rethold[1][0]);
		      retall[0].push_back(rethold[0][0]);

		    }
		}
	    }

	  //3d gd:
	  t.Start();
	  demo.steepestgradientdescent3D(uhats, Epyr, demo.storage3dpyr, demo.coordpyr, demo.coordmidpyr, demo.midptevalpyr, rethold, avgiterGDhold);//(demo.find_roots(uhats, E, demo.storage3d,  avgiterGDhold, 0, 0, 1, 0)) ;
	  t.Stop();
	  roots3dtime +=  t.TimePerTest(1);


	  avgiterGD += avgiterGDhold;
	  
	  for(int k = 0; k < dimension; k++)
	    {
	      if(rethold[k][0] <inf)
		{
		  retall[k].push_back(rethold[k][0]);
		}
	    }
	}
      
      // cout<<"\n 3d rootfinder ret: "<<rethold[0][0]<<" "<<rethold[1][0]<<" "<<rethold[2][0];
      Array<OneD, Array<OneD, NekDouble> > tmpcoord(dimension);
      for(int p = 0; p < dimension; p++)
	{
	  tmpcoord[p] = Array<OneD, NekDouble>(retall[0].size());
	}
      
      for(int p = 0; p < dimension; p++)
	{
	  for(int q = 0; q < retall[0].size(); q++)
	    {
	      tmpcoord[p][q] = retall[p][q];


	    }
	}
    NekDouble tempmin = inf;
    Array<OneD, NekDouble> evalroots, ret(dimension);
    
    // if(dimension > 2)
    //   {    cout<<"\n before eval\n";
    evalroots = E->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
	 //      }
    // else if(dimension > 1)
    //   {
    // 	if(E->DetShapeType() == LibUtilities::eQuadrilateral)
    // 	  {
    // 	    evalroots = E->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
    // 	  }
    // 	else //tri
    // 	  {
    // 	     evalroots = E->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
    // 	  }
    //   }
    Array<OneD, NekDouble> minvandx(dimension+1);
    Array<OneD, NekDouble> tmparr(retall[0].size());
    demo.pq(uhats, tmpcoord, evalroots, minvandx, tmparr);
    // cout<<"\n vals at roots:\n";                                                                  
    // for(int k = 0; k < tmparr.size(); k++)                                                        
    //   cout<<" "<<tmparr[k]<<" ";                                                                  
    // cout<<"\n";                                                                                   
    
    tempmin = minvandx[0];
    //    cout<<"\n out of pq minvandxsz = "<<minvandx.size()<<":\n";
    for(int k = 0; k < dimension; k++)
      {
	ret[k] = minvandx[k+1];
      } 
    minv = tempmin;
    return ret;
    
  } 

  void FilterStructurePres::v_Finalise(
					 const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
					 const NekDouble &time)
    {
      boost::ignore_unused(pFields, time);
    
    }
  
    bool FilterStructurePres::v_IsTimeDependent()
    {
      return false;
    }
  }
}
//     Array<OneD, Array<OneD, NekDouble> > FilterStructurePres::formConf(NekDouble N)
//     {
//       int evalPts = N+1; //+1 for luck
//       int i, j;
//       Nektar::Array<Nektar::OneD, Nektar::NekDouble> ab1(evalPts);
//       Nektar::Array<Nektar::OneD, Nektar::NekDouble> ab2(evalPts);
//       Polylib::RecCoeff(evalPts, &ab1[0], &ab2[0], -0.5, -0.5);
//       // Form confederate matrix
//       // a = 2*a
//       // b = 2*b
//       Vmath::Smul(ab1.size(), 2.0, ab1, 1, ab1, 1);
//       Vmath::Smul(ab2.size(), 2.0, ab2, 1, ab2, 1);

//       //  J = full(spdiags([[b(3:n);0.5;0] a(1:n) b(1:n)], -1:1, n, n));
//       vector<vector<NekDouble> > J;//, Array<OneD, NekDouble>(N,0.0));

//       vector<NekDouble> row(N);
//       NekDouble tt = ab1[0];
//       row[0] = tt;
//       tt = ab2[1];
//       row[1] = tt;
//       J.push_back(row);
//       for(i = 1; i < N-1; ++i)
// 	{
// 	  vector<NekDouble> row(N,0.0);
// 	  for(j = 1; j < N-1; ++j)
// 	    {
// 	      if( i == j)
// 		{
// 		  NekDouble t1 = ab1[i];
// 		  NekDouble t2 = ab2[i+1];
// 		  NekDouble t3 = ab2[i+1];
// 		  row[j] = t1;
// 		  row[j+1] = t2;
// 		  row[j-1] = t3;
// 		}
// 	    }
// 	  J.push_back(row);
// 	}
//       vector<NekDouble> rowN(N);
//       tt = ab1[N-1];
//       rowN[N-1] = tt;
//       tt = ab2[N-1];
//       rowN[N-2] = tt;
//       J.push_back(rowN);

//       Array<OneD, Array<OneD, NekDouble> > C(J[0].size());
//       for(i = 0; i < J[0].size(); i++)
// 	{
// 	  C[i] = Array<OneD, NekDouble>(J.size());
// 	}
//       for(int k = 0; k < C.size(); k++)
// 	{
// 	  for(int p = 0; p < C[0].size(); p++)
// 	    {
// 	      C[k][p] = J[p][k];
// 	    }
// 	}
//       return C;
//     }
    
//   } // namespace SolverUtils
// } // namespace Nektar



    
