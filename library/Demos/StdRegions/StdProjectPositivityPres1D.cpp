///////////////////////////////////////////////////////////////////////////////
//
// File: StdProjectPositivityPres1D.cpp
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
// Description: Demo for 1D Positivity pres filter
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "StdDemoSupport.hpp"
#include <LibUtilities/BasicUtils/Timer.h>

namespace po = boost::program_options;

//declare filter call
void Do_optimize(StdExpansion *E, Array<OneD, Array<OneD, NekDouble> > &storage, Array<OneD, NekDouble> &uhats);

Array<OneD,NekDouble>  call_find_roots(Array<OneD,  NekDouble> &uhats1 , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats ,Array<OneD, Array<OneD, NekDouble> >&storage,  StdExpansion *E, NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime);

// uhatpqd stuff: 
void edgederpquhats(Array<OneD, NekDouble> &uhats1, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> &wsp,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2);


string funcdef;

int iterstaken;
Array<OneD, Array<OneD, NekDouble> > edgexy(2);
Array<OneD, Array<OneD, NekDouble> > storage, storageorth;
// Array<OneD, Array<OneD, NekDouble> > storage3;

int numedges, numsurfaces;
StdExpansion *E, *Eorth;

Array<OneD, NekDouble> Vx, Vdx, startcoord;
int dimension ;
            
DemoSupport demo;
LibUtilities::ShapeType stype;
Array<OneD, Array<OneD, NekDouble> > edgeptsin;

int main(int argc, char *argv[])
{
  demo.GetOptions().add_options()("optm,z", "positivity preserving optimizer");
  demo.ParseArguments(argc, argv);

  po::variables_map vm = demo.GetVariableMap();

  E = demo.CreateStdExpansion();
  storage = E->GetPhysEvaluateStorage(); 
  
  if (E == nullptr)
    {
      return 1;
    }
  dimension = E->GetShapeDimension();

  if(dimension != 1 )
    {
      cout<<"\n dimension should be 1, try using StdProjectPositivityPres1D or StdProjectPositivityPres3DD for other dimensions\n\n";
      exit(0);
    }
  startcoord = Array<OneD, NekDouble>(dimension);
  std::vector<int> order;
  std::vector<BasisType> btype(3, eNoBasisType);
  stype = E->DetShapeType();
  for (int i = 0; i < dimension; ++i)
    {
      btype[i] = E->GetBasisType(i);
      order.push_back(E->GetBasisNumModes(i));
    }

    
  switch(E->DetShapeType())
    {
    case LibUtilities::eSegment:
       demo.testcoord1dpts = demo.GetCoords(E);
       demo.testcoord1dmidptspts = demo.GetQuadratureMidCoords(E, demo.testcoord1dpts);
       demo.testcoord1dlattice = demo.GetLatticeCoords(demo.testcoord1dpts, demo.testcoord1dmidptspts);
       demo.interioreval1dmidpts = E->PhysEvaluateBasis(demo.testcoord1dmidptspts, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
       numedges = 1;
       numsurfaces = 0;
       break;
	
    
    default:
      cout<<"\n This demo works only for seg,  try using StdProjectPositivityPres2D or StdProjectPositivityPres3D for other shapes\n\n";
      exit(0);

    }
  const auto totPoints = (unsigned) E->GetTotPoints();

  Array<OneD, NekDouble> x = Array<OneD, NekDouble>(totPoints);
  Array<OneD, NekDouble> y = Array<OneD, NekDouble>(totPoints);
  Array<OneD, NekDouble> z = Array<OneD, NekDouble>(totPoints);
  Array<OneD, NekDouble> dx = Array<OneD, NekDouble>(totPoints);
  Array<OneD, NekDouble> dy = Array<OneD, NekDouble>(totPoints);
  Array<OneD, NekDouble> dz = Array<OneD, NekDouble>(totPoints);
  Array<OneD, NekDouble> sol = Array<OneD, NekDouble>(totPoints);


  switch (dimension)
    {
    case 1:
      {
	E->GetCoords(x);
	break;
      }

    case 2:
      {
	E->GetCoords(x, y);
	break;
      }

    case 3:
      {
	E->GetCoords(x, y, z);
	break;
      }
    default:
      break;
    }

  //get solution array
  for (int i = 0; i < totPoints; ++i)
    {

      if(dimension ==1)
        {

	  sol[i] =2*pow(x[i],3) -0.5*pow(x[i]-0.1, 2)-x[i];
	  //pow(x[i]-0.1,2) - 0.09;
	  //(floor(0.5*(x[i]<=-0.4)+0.5*(x[i]>=-0.8)));
	   
	  //pow(x[i]-0.1,2) -0.1;
	  funcdef = "pow(x[i]-0.1,2) - 0.09;";        
	}
      else
	{
	  cout<<"\n This demo only works for Segments\n";
	  exit(0);
	}
    }

  Array<OneD, NekDouble> phys(totPoints);
  Array<OneD, NekDouble> coeffs((unsigned) E->GetNcoeffs());
 
  //Project onto expansion
  E->FwdTrans(sol, coeffs);

  //Backward transform solution to get projected values
  E->BwdTrans(coeffs, phys);
  LibUtilities::Timer     timer;
  NekDouble elapsed       = 0.0;


  // check for -ve values and apply opt if necessary
  if (vm.count("optm"))
    {
      timer.Start();

      Array<OneD, BasisType> btypearr(dimension);
      btypearr[0] = LibUtilities::eOrtho_A;
      // storage3 = E3seg->GetPhysEvaluateStorage(); 
      edgeptsin = Array<OneD, Array<OneD, NekDouble> >(dimension);
      Array<OneD, NekDouble> edgexytemp =   demo.E3seg->GetBasis(0)->GetZ();
	  
      for(int k = 0; k < dimension; k++)
	{
	  edgexy[k] = edgexytemp;
	  edgeptsin[k] =   Array<OneD, NekDouble>(edgexy[k]);

	}
      edgexytemp =   demo.E3seg->GetBasis(0)->GetZ();
      int totszedges = edgexytemp.size()*(E->GetNcoeffs());

      Vx = Array<OneD, NekDouble>(totszedges);
      Vdx = Array<OneD, NekDouble>(totszedges);
      if(demo.Eorthseg == nullptr)
	{
	  Eorth = E;
	}
      else
	{
	  Eorth = demo.Eorthseg;
	}
      storageorth = Eorth->GetPhysEvaluateStorage();
      Vx = Eorth->PhysEvaluateBasis(edgexy, storageorth, Vdx,  NullNekDouble1DArray, NullNekDouble1DArray);

      timer.Stop();
      elapsed  = timer.TimePerTest(1);
      cout<<"\n setup phase took "<<elapsed<<"s\n\n";

      timer.Start();

      int found_negb4 = 0;
      NekDouble  inf = numeric_limits<double>::infinity();
      NekDouble minall = inf;
      // evaluate fn on quad pts
      Array<OneD, NekDouble> tmp1(coeffs.size());
      NekDouble tempv;
      for(int k = 0; k < phys.size(); k++)
	{
	  if(phys[k] < 0 && abs(phys[k]) >1e-9)
	    {
	      if(minall > phys[k])
		{
		  minall = phys[k];
		  startcoord[0] = demo.testcoord1dpts[0][k];
		}
	      cout<<"\n phys neg = "<<phys[k]<<" at "<<demo.testcoord1dpts[0][k];
	      found_negb4 = 1;
	    }
	}
      
      //cout<<"\n vals at other pts of lattice:\n";
      for(int k = 0; k < demo.testcoord1dmidptspts[0].size(); ++k)
	{
	  
	  Vmath::Vmul(coeffs.size(), &demo.interioreval1dmidpts[k], demo.testcoord1dmidptspts[0].size(), &coeffs[0], 1, &tmp1[0], 1);
	  tempv = Vmath::Vsum(tmp1.size(), tmp1, 1);
	  if(tempv < 0 &&abs(tempv) > 1e-9)
	    {
		  if(minall >  tempv)
		    {
		      cout<<"\n val = "<<tempv<<" at "<<demo.testcoord1dmidptspts[0][k];
		      minall = tempv;
		      startcoord[0] = demo.testcoord1dmidptspts[0][k];
		      found_negb4 = 1;
		    }
		}
	}
      
      timer.Stop();
      elapsed += timer.TimePerTest(1);
      cout<<"\n checking for -vity took "<<timer.TimePerTest(1)<<"s\n\n";
      
      if(found_negb4)
	{
	  cout<<"\n need optimization\n\n";
	  timer.Start();
	  int orthoflag = (E->GetBasisType(0) != LibUtilities::eOrtho_A); 
	  if(orthoflag)//demo.Eorthseg != nullptr)
	    {
	      demo.OrthoNormalize(E, coeffs, btypearr, 0);
	    }
	  Do_optimize(Eorth, storageorth,  coeffs);
	  
	  if(orthoflag)
	    {
	      demo.OrthoNormalize(E, coeffs, btypearr, 1);

	    }
	  timer.Stop();
	  elapsed += timer.TimePerTest(1);
	  cout<<"\n optimizertook "<<timer.TimePerTest(1)<<"s\n\n";

	}
      else{
	elapsed += timer.TimePerTest(1);
	  
	cout<<"\n optimizer no need\n\n";
      }
      //Backward transform solution to get projected values
      E->BwdTrans(coeffs, phys);

      int found_neg = 0;
    
      // evaluate fn on quad pts
      for(int k = 0; k < phys.size(); k++)
	{
	  if(phys[k] < 0 && abs(phys[k]) >1e-9)
	    {
	      cout<<"\n phys neg = "<<phys[k]<<" at "<<demo.testcoord1dpts[0][k]<<".";
	      found_neg = 1;
	    }
	  
	  for(int k = 0; k < demo.testcoord1dmidptspts[0].size(); ++k)
	    {
	      Vmath::Vmul(coeffs.size(), &demo.interioreval1dmidpts[k], demo.testcoord1dmidptspts[0].size(), &coeffs[0], 1, &tmp1[0], 1);
	      tempv = Vmath::Vsum(tmp1.size(), tmp1, 1);
	      if(tempv < 0 &&abs(tempv) > 1e-9)
		{
		  cout<<"\n tempv = "<<tempv<< " at "<<demo.testcoord1dmidptspts[0][k]<<"\n";;
		  found_neg = 1;
		}
	    }
	}
      
      if(found_neg)
	{
	  timer.Stop();
	  elapsed += timer.TimePerTest(1);
	  cout<<"\n Verification took "<<timer.TimePerTest(1)<<"s\n\n";
	  
	  cout<<" fail\n";
	  cout<<"\n func="<<funcdef<<"\n";
	  if(numedges == 3)
	    cout<<"\ntri\n";
	  else
	    cout<<"\n quad\n";
	}
      else
	{
	  timer.Stop();
	  elapsed += timer.TimePerTest(1);
	  cout<<"\n Verification took "<<timer.TimePerTest(1)<<"s\n\n";
	      
	  cout<<"pass\n";
	      
	  cout<<"\n func="<<funcdef<<"\n";;
	  if(numedges == 3)
	    cout<<"\ntri\n";
	  else
	    cout<<"\n quad\n";
	      
	  phys = sol;
	}
	  
    }


  //Calculate L_inf & L_2 error
  cout << "L infinity error: \t" << E->Linf(phys, sol) << endl;
  if (stype != ePoint)
    {
      cout << "L 2 error: \t \t \t" << E->L2(phys, sol) << endl;
    }

  return 0;
    
}


void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> &wsp,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2)
{
  int uhatstot = uhats.size();
  int totpts = edgeptsin[0].size();

  Array<OneD, NekDouble> pqeval(totpts);
  NekDouble v1, v2;
 
  for(int i = 0; i<totpts; i++)
    {

      Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd0[0]+i, totpts, &wsp[0], 1);
      v1  = Vmath::Vsum(uhatstot, wsp, 1); 

      Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &wsp[0], 1);
      v2  = Vmath::Vsum(uhatstot, wsp, 1);  
      v1 = v2*v1;
      
      // At this point,                 
      // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point
      // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
      Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &wsp[0], 1);
      v2  = Vmath::Vsum(uhatstot, wsp, 1);  
 
      Vmath::Vmul(uhats.size(), &Vxyd0[i], totpts, &uhats[0], 1, &wsp[0], 1);
      v1= v1-v2*Vmath::Vsum(uhats.size(), wsp, 1);
  
      pqeval[i] = v1;
 
      if(dimension > 1)
	{
	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd1[0]+i, totpts, &wsp[0], 1);
	  v1  = Vmath::Vsum(uhatstot, wsp, 1); 

	  Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &wsp[0], 1);
	  v2  = Vmath::Vsum(uhatstot, wsp, 1);  
	  v1 = v2*v1;

	  // At this point,                 
	  // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point
	  // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &wsp[0], 1);
	  v2  = Vmath::Vsum(uhatstot, wsp, 1);  

 	  Vmath::Vmul(uhats.size(), &Vxyd1[i], totpts, &uhats[0], 1, &wsp[0], 1);
	  v1=  v2*Vmath::Vsum(uhats.size(), wsp, 1) - v1;
	  pqeval[i] += v1;
 
	}
      if(dimension == 3)
	{
	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd2[0]+i, totpts, &wsp[0], 1);
	  v1  = Vmath::Vsum(uhatstot, wsp, 1); 

	  Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &wsp[0], 1);
	  v2  = Vmath::Vsum(uhatstot, wsp, 1);  
	  v1 = v2*v1;

	  // At this point,                 
	  // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point
	  // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &wsp[0], 1);
	  v2  = Vmath::Vsum(uhatstot, wsp, 1);  
 
	  Vmath::Vmul(uhats.size(), &Vxyd2[i], totpts, &uhats[0], 1, &wsp[0], 1);
	  v1= v2*Vmath::Vsum(uhats.size(), wsp, 1) - v1;
 
	  pqeval[i] += v1;
	}
    }
  demo.E3seg->FwdTrans(pqeval, ret);
}

Array<OneD,  NekDouble> call_find_roots(Array<OneD,  NekDouble> &uhats , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats , Array<OneD, Array<OneD, NekDouble> >&storage, StdExpansion *Etmp,  NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime)
{
  
  boost::ignore_unused(surfaceuhats, roots2dtime, roots3dtime, roots1dtime, avgiterGD);
  int dimension = E->GetShapeDimension(); 
  NekDouble dummy;
  Array<OneD, NekDouble> temp(uhats.size());

  NekDouble  inf = numeric_limits<double>::infinity();
  //    NekDouble avgiterGDhold;
  Timer t;
  Array<OneD, Array<OneD, NekDouble > >  rethold(dimension);
  vector<vector<NekDouble> > retall(dimension);
  Array<OneD, Array<OneD, NekDouble> > tmpcoord(dimension);

  Array<OneD, NekDouble> ret(dimension);
  if(numedges == 1) // segment
    {
      t.Start();
      rethold  = (demo.find_roots(uhatsedges[0], nullptr, storage, dummy,  0, 0, 0, 0)) ;
      t.Stop();
      cout<<"\n roots returned:\n";
      for(int p = 0; p < rethold[0].size(); p++)    
	{
	  cout<<" "<<rethold[0][p];
	  retall[0].push_back(rethold[0][p]);
	}
    }
  for(int p = 0; p < dimension; p++)
    {
      tmpcoord[p] = Array<OneD, NekDouble>(retall[0].size());
    }
  for(int p = 0; p < retall[0].size(); p++)
    {
      for(int j = 0; j < dimension; j++)
	{
	  tmpcoord[j][p] = retall[j][p];
	}
      // 	  tmpcoord[1][p] = retall[1][p];
    }

	
  NekDouble tempmin = inf;
  Array<OneD, NekDouble> evalroots;
  evalroots = Etmp->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
	    

  Array<OneD, NekDouble> minvandx(dimension+1);
  Array<OneD, NekDouble> tmparr(retall[0].size());
  demo.pq(uhats, tmpcoord, evalroots, minvandx, tmparr);
  tempmin = minvandx[0];
  for(int k = 0; k < dimension; k++)
    ret[k] = minvandx[k+1];
  
  minv = tempmin;
  return ret;

}


void Do_optimize(StdExpansion *E, Array<OneD, Array<OneD, NekDouble> > &storage, Array<OneD, NekDouble> &uhats)
{

  int dim = E->GetShapeDimension();
  double inf = numeric_limits<double>::infinity();
    
  //assert(size(constraints, 1) == N+2);
    
  int N1 = uhats.size();
  vector<Array<OneD,NekDouble> > d;
 
  d.push_back(uhats);
  int counter = 0;

  // constraint specific tolerances
  vector<double> tols;
  tols.push_back(1e-12);
    
  int niter = 1e3;
  NekDouble minv;

  //number of constraints, only positivity for now
  //    int NC = 1; 

  Array<OneD, NekDouble>  pqvalcoords(dim+1), xastarr(dim), utemp(N1), wsp1(1),
    optima(dim), wsp(N1), Vtmp(E->GetNcoeffs()), tmp,  Vastsq (N1), Vast (N1);;

  Array<OneD, Array<OneD, NekDouble> >  tmpcoord(dim), Pf(numedges), xastarrofarr(dim);
  for(int k= 0; k < numedges; k++)
    {
      Pf[k] = Array<OneD, NekDouble>(demo.E3seg->GetNcoeffs()); 
    }
  cout<<"\n uhats sz = "<<uhats.size()<<" Pf[0].sz = "<<Pf[0].size();
  NekDouble pqval, roots1dtimehold = 0.0, timeprojectedges = 0.0, dum, vastsqsum;
  Timer t;
  while (counter <= niter)
    {
      NekDouble roots1dtime = 0.0 ;
      pqval = inf;
      utemp = d.back();
      if (counter > 0)
	{
	  t.Start();
	  
	  //pq = @(xx) cfun(xx) .* (constraints(1,jj).*Vc(xx)*uhats);  
	  edgederpquhats(utemp, Pf[0], wsp, Vx, Vdx, NullNekDouble1DArray,
			 NullNekDouble1DArray);
	  //Vmath::Vcopy(uhats.size(), uhats, 1, Pf[0], 1);
	  t.Stop();
	  timeprojectedges += t.TimePerTest(1);
	  optima = call_find_roots(utemp, dum, Pf, NullNekDoubleArrayOfArray, storage, E, minv, roots1dtime, dum, dum);
	  roots1dtimehold += roots1dtime;
	  
	  // Vtmp is evaluation of basis fn at optima
	  for(int k = 0; k < dim; k++)
	    {
	      tmpcoord[k] = Array<OneD, NekDouble>(1, optima[k]);
	    }
	  
	  Vtmp = E->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
	}
      else // end of counter > 0
	{
	  for(int k = 0; k < dim; k++)
	    {
	      tmpcoord[k] = Array<OneD, NekDouble>(1, startcoord[k]);
	    }
	  Vtmp = E->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
	}

      demo.pq(utemp, tmpcoord, Vtmp, pqvalcoords, wsp1);
      if (pqvalcoords[0] < pqval)
	{
	  xastarr = optima;
	  pqval = pqvalcoords[0];
	}
      
      cout<<"\n pqval = "<<pqval<<" at "<<xastarr[0]<<"\n";        
      // If minimum is non-negative, we're done
      if (pqval >= -tols.at(0))
	{
	  break;
	}
      
      Array<OneD, NekDouble>  qast(N1);
      for(int j = 0; j < dim; j++)
	{
	  xastarrofarr[j] = Array<OneD, NekDouble>(1, xastarr[j]);
	}
      
      Vast = E->PhysEvaluateBasis(xastarrofarr, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
      
      for( int i = 0; i < N1; i++)
	{
	  Vastsq[i] = (Vast[i]*Vast[i]);
	}

      vastsqsum = Vmath::Vsum(N1, &Vastsq[0], 1);

      
      for(int i = 0; i<N1; i++)
	{
	  qast[i] = ((1/sqrt(vastsqsum))*(Vast[i]));
	}

      Vmath::Smul(N1, pqval, &qast[0], 1, &qast[0], 1);
      Vmath::Vsub(utemp.size(), &utemp[0], 1, &qast[0], 1, &qast[0], 1);
      d.push_back(qast);

      counter = counter + 1;
    
    
    }
  roots1dtimehold = roots1dtimehold/(counter-1);
  timeprojectedges = timeprojectedges/(counter);
  iterstaken = counter;
  cout<<"sphere_rotation took "<<counter<<"iterations.";
  uhats = d.back();
}

