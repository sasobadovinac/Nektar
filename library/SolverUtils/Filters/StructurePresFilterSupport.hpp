//////////////////////////////////////////////////////////////////////////////
//
// File: StructurePresFilterSupport.cpp
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
// Description: Demo for testing functionality of StdProject
//
///////////////////////////////////////////////////////////////////////////////

#ifndef DEMOS_STDREGIONS_STRUCTUREPRESFILTERSUPPORT_HPP
#define DEMOS_STDREGIONS_STRUCTUREPRESFILTERSUPPORT_HPP

#include <string>
#include <vector>

#include <StdRegions/StdPointExp.h>
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdHexExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdNodalPrismExp.h>
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdNodalTetExp.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/Foundations/InterpCoeff.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <LibUtilities/BasicUtils/Timer.h>        
using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;

class StructurePresFilterSupport
{
public:

  
  Array<OneD, Array<OneD, NekDouble > > storage3dhex, storage3dtet, storage3dpyr, storage3dpri, storage2dt, storage2dq, coordmidhex, coordmidtet, coordmidquad,  coordmidtri, coordmidpri, coordmidpyr, coordhex, coordtet, coordpyr, coordpri, coordtri, coordquad, coordlatticequad, coordlatticetri, coordlatticehex, coordlatticetet, coordlatticepyr, coordlatticepri;
  
  Array<OneD, NekDouble> midptevalhex, midptevaltet, midptevalpyr, midptevalquad, midptevaltri, midptevalpri;
  NekDouble ub = numeric_limits<double>::infinity();
  int details = 0;
  Array<OneD, Array<OneD, NekDouble> > formCompanion(NekDouble N)
  {
    // returns monomial-connection matrix as a base
    // over which we can construct companion matrix
    
    // recurrence coeff:
    Nektar::Array<Nektar::OneD, Nektar::NekDouble> ab1;
    Nektar::Array<Nektar::OneD, Nektar::NekDouble> ab2;

    ab1 = Array<OneD, NekDouble>(N);
    ab2 = Array<OneD, NekDouble>(N);
    Polylib::RecCoeff(N, &ab1[0], &ab2[0], 0, 0);   

    Array<OneD, Array<OneD, NekDouble> > R(N);
    //b = sqrt(b)
    transform(ab2.begin(), ab2.end(), ab2.begin(), (double(*)(double)) sqrt);

    //b = 1./b
    Array<OneD, NekDouble> r(ab2.size(), 1.0);

    Vmath::Vdiv(r.size(), r, 1, ab2, 1, r, 1);

    //r = cumprod(1./b)
    partial_sum (&r[0], &r[0]+r.size(), &r[0], multiplies<double>());


    //Build the connection matrix for this dimension and then distribute it accordingly

    for(int i = 0; i<N; i++)
      {
	R[i] = Array<OneD, NekDouble>(N);
	for(int j = 0; j<N; j++)
	  {
	    if( i == j)
	      R[i][j] = r[i];
	    else
	      R[i][j] = 0;
	  }
      }
    for(int i = 1; i<N; i++)
      {
	//"side" conditions (i.e. those w/o left/right boundary points)
	if(i<2)
	  {
	    R[i][0] = (1/ab2[i])*(-ab1[i-1]*R[i-1][0]);
	  }
	else
	  {
	    R[i][0] = (1/ab2[i])*(-ab1[i-1]*R[i-1][0] - ab2[i-1]*R[i-2][0]);

	    for(int j = 1; j<i-1; j++)
	      {
		R[i][j] = R[i-1][j-1] - ab1[i-1]*R[i-1][j] - ab2[i-1]*R[i-2][j];
		R[i][j] = R[i][j]/ab2[i];
	      }

	  }
      }
    return R;
  }
  
  Array<OneD, Array<OneD, NekDouble> > formConf(NekDouble N)
  {
    
    // recurrence coeff:
    Nektar::Array<Nektar::OneD, Nektar::NekDouble> ab1;
    Nektar::Array<Nektar::OneD, Nektar::NekDouble> ab2;
    
    ab1 = Array<OneD, NekDouble>(N);
    ab2 = Array<OneD, NekDouble>(N);

    Polylib::RecCoeff(N, &ab1[0], &ab2[0], -0.5, -0.5);   
    int i, j;
    // Form confederate matrix
    // a = 2*a
    // b = 2*b
    Vmath::Smul(ab1.size(), 2.0, ab1, 1, ab1, 1);
    Vmath::Smul(ab2.size(), 2.0, ab2, 1, ab2, 1);
    
    //  J = full(spdiags([[b(3:n);0.5;0] a(1:n) b(1:n)], -1:1, n, n));
    vector<vector<NekDouble> > J;
    vector<NekDouble> row(N);
    NekDouble tt = ab1[0];
    row[0] = tt;
    tt = ab2[1];
    row[1] = tt;
    J.push_back(row);

    for(i = 1; i < N-1; ++i)
      {
	vector<NekDouble> row(N,0.0);
	for(j = 1; j < N-1; ++j)
	  {
	    if( i == j)
	      {
		NekDouble t1 = ab1[i];
		NekDouble t2 = ab2[i+1];
		NekDouble t3 = ab2[i+1];
		row[j] = t1;
		row[j+1] = t2;
		row[j-1] = t3;
	      }
	  }
	J.push_back(row);
      }
    vector<NekDouble> rowN(N);
    tt = ab1[N-1];
    rowN[N-1] = tt;
    tt = ab2[N-1];
    rowN[N-2] = tt;
    J.push_back(rowN);

    Array<OneD, Array<OneD, NekDouble> > C(J[0].size());
    for(i = 0; i < J[0].size(); i++)
      {
	C[i] = Array<OneD, NekDouble>(J.size());
      }
    for(int k = 0; k < C.size(); k++)
      {
	for(int p = 0; p < C[0].size(); p++)
	  {
	    C[k][p] = J[p][k];
	  }
      }
    return C;
  }
  
  Array<OneD, Array<OneD, NekDouble>  >FindLatticeEval(Array<OneD, NekDouble> coeffs,  Array<OneD, Array< OneD, NekDouble> >storage,   Array<OneD, Array<OneD, NekDouble> > quadcoord , Array<OneD, NekDouble> evalbasisatmid, Array<OneD, Array<OneD, NekDouble>> quadmidcoord )
  {
    boost::ignore_unused(quadmidcoord, evalbasisatmid);
    int uhatstot = coeffs.size();
    Array<OneD, NekDouble> temp(uhatstot); 
    double inf = numeric_limits<double>::infinity();
    NekDouble minv = inf;
    NekDouble maxv = -1;

    int dimension = quadcoord.size();
    Array<OneD, Array<OneD, NekDouble > >ret(2);
    ret[0] = Array<OneD, NekDouble>(dimension+1);
    ret[1] = Array<OneD, NekDouble>(dimension+1);
    
    int totpts = quadcoord[0].size();
    
    for(int i = 0; i<totpts; i++)
      {
        Vmath::Vmul(uhatstot, &storage[0][0]+i, totpts, &coeffs[0], 1, &temp[0], 1 );
        NekDouble holdv = Vmath::Vsum(uhatstot, temp, 1);
		
	if(holdv < minv)
          {
            minv = holdv;
	  
            for(int k = 0; k < dimension; k++)
              {
                ret[0][k] = quadcoord[k][i];
	      }
            ret[0][dimension] = minv;
          }
	if(holdv > maxv)
	  {
	    maxv = holdv;
	    for(int k = 0; k < dimension; k++)
	      {
		ret[1][k] = quadcoord[k][i];
	      }
	    ret[1][dimension] = maxv;
	  }
	
      }
    
    
    return ret;
  }
  

  
  //if flag = 0 -> Non-ortho to Ortho
  //if flag = 1 -> Ortho to Non-ortho
  void OrthoNormalize(LocalRegions::ExpansionSharedPtr E, Array<OneD, NekDouble>  &coeffs,  StdExpansion *Eorth, int flag)
  {
    vector<LibUtilities::BasisKey> bkeyarr, bkeyortharr;
    for(int k = 0; k < E->GetShapeDimension(); k++)
      {
	int nmodes0 = E->GetBasis(k)->GetNumModes();
	LibUtilities::BasisKey orthoBasisKey0 (
					       Eorth->GetBasis(k)->GetBasisType(),
					       nmodes0,
					       E->GetBasis(k)->GetPointsKey());
	LibUtilities::BasisKey nonOrthoBasisKey0 (
						  E->GetBasis(k)->GetBasisType(),
						  nmodes0,
						  E->GetBasis(k)->GetPointsKey());
	bkeyarr.push_back(nonOrthoBasisKey0);
	bkeyortharr.push_back(orthoBasisKey0);
      }
    Array<OneD,NekDouble> temp1(coeffs.size());
    
    if(bkeyarr.size() == 1) //dimension = 1
      {
	
	if(flag == 0) //Non-ortho to Ortho 
	  { 
	    LibUtilities::InterpCoeff1D(bkeyarr[0],  coeffs, bkeyortharr[0], temp1);
	  }
	else  if(flag == 1) // Ortho to Non-ortho
	  {
	    LibUtilities::InterpCoeff1D(bkeyortharr[0],  coeffs,  bkeyarr[0], temp1);
	  }
      }
    else if(bkeyarr.size() == 2) //dimension = 2
      {
	if(flag == 0) //Non-ortho to Ortho 
	  {
	    LibUtilities::InterpCoeff2D(bkeyarr[0],  bkeyarr[1],  coeffs, bkeyortharr[0], bkeyortharr[1], temp1);
	  }
	else if(flag == 1) // Ortho to Non-ortho
	  {

	    LibUtilities::InterpCoeff2D(bkeyortharr[0], bkeyortharr[1],  coeffs,  bkeyarr[0], bkeyarr[1], temp1);
	  }

      }
    Vmath::Vcopy(coeffs.size(), temp1, 1, coeffs, 1);
	
  }

  

  void FindEigenval(Array<OneD, Array<OneD, NekDouble> > CM,
  		    Array<OneD,NekDouble> &EIG_R,
  		    Array<OneD,NekDouble> &EIG_I)
  {

    boost::ignore_unused(EIG_I);
    const unsigned int n = CM.size();
    const int sz = (n)*(n);
    Array<OneD, NekDouble> A(sz,0.0);

    for(int i = 0; i<n; i++)
      {
	Vmath::Vcopy(n, &CM[i][0], 1, &A[0]+i*n, 1);
      }
    Array<OneD, NekDouble>  wr(n), wi(n);
    Nektar::FullMatrixFuncs::EigenSolve(n, A, wr, wi, NullNekDouble1DArray);
    vector<NekDouble>EIG_Rvec;
    for(int k = 0; k < wr.size(); k++)
      {
	if(abs(wi[k])<tol_gd && abs(wr[k])<=1.0)
  	  {
  	    EIG_Rvec.push_back(wr[k]);
  	  }
      }

    EIG_R = Array<OneD, NekDouble>(EIG_Rvec.size());
    for(int k = 0; k < EIG_Rvec.size(); k++)
      {
  	EIG_R[k] = EIG_Rvec[k];
      }
    
  }

  // nu = 1 -> lower bound, else -1 for upper bound
  void pq(
	  Array<OneD,NekDouble> uhats,
	  Array<OneD, Array<OneD,  NekDouble> > roots,
	  Array<OneD, NekDouble> V1,
	  Array<OneD,NekDouble> &pqevalxast,
	  Array<OneD,NekDouble> &fvals,
	  NekDouble nu)
  {
    int N = uhats.size();
    Array<OneD,NekDouble> w1(N);
    Array<OneD, NekDouble> Vsumsq(roots[0].size());
    fvals = blasmatvec(V1, uhats, roots[0].size(), uhats.size());
    Vmath::Smul(fvals.size(), nu, fvals, 1, fvals, 1);
    for( int i = 0; i < roots[0].size(); i++)
      {
	Vmath::Vmul(N, &V1[i], roots[0].size(), &V1[i], roots[0].size(), &w1[0], 1);

	Vsumsq[i] = (pow(Vmath::Vsum(N, &w1[0], 1),-0.5));
      }
    Vmath::Vmul(roots[0].size(), &Vsumsq[0], 1, &fvals[0], 1, &fvals[0], 1);
    if(pqevalxast.size() > 0)
      {
	int minidx ;
	minidx = Vmath::Imin(fvals.size(), &fvals[0], 1);
	
	pqevalxast[0] = fvals[minidx];
	for(int k = 0; k < pqevalxast.size()-1; k++)
	  {
	    pqevalxast[k+1] = roots[k][minidx];
	  }
      }

  }

    
  void derpq(Array<OneD, NekDouble> &uhats, NekDouble &ret,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd)
  {
    int N = uhats.size();
    Array<OneD, NekDouble> temp(N);

    //\sum(phi^2)
    Vmath::Vmul(N, Vxy, 1, Vxy, 1, temp, 1);
    NekDouble tmp = Vmath::Vsum(N, temp, 1);

    //\sum(u*phi')
    Vmath::Vmul(N, uhats, 1,Vxyd, 1, temp, 1);
    NekDouble v1 = Vmath::Vsum(N, temp, 1);

    v1 = v1*tmp;

    //\sum(u*phi)
    Vmath::Vmul(N, uhats, 1,Vxy, 1, temp, 1);
    NekDouble v2 = Vmath::Vsum(N, temp, 1);

    // //\sum(phi*\phi')
    Vmath::Vmul(N, Vxyd, 1, Vxy, 1, temp, 1);
    tmp = Vmath::Vsum(N, temp, 1);

    v2 = v2*tmp;

    ret = v1-v2;

  }

  //nq1 = no. of rows of mat
  //nq2 = no. of cols (size of vec which is = cols of mat)
  Array<OneD, NekDouble> blasmatvec(Array<OneD, NekDouble> M, Array<OneD, NekDouble> vec, int \
				    nq1, int nq2 )
  {

    Array<OneD, NekDouble> hold(nq1);
    // we want fn eval on lattice
    //    for(int k = 0; k < GetAvgNum(); k++)
    //{
	double alpha= 1.0, beta= 0.0;
	char no= 'N';
	int m = nq1, n = nq2, lda = m, incx = 1, incy = 1;
	double *A = &M[0];
	//t1.Start();
	Blas::dgemv_(no,m,n,alpha,A,lda,&vec[0],incx,beta,&hold[0],incy);
	//t1.Stop();
	//}
    return hold;
  }

  void steepestgradient_descent2D(Array<OneD, NekDouble> &uhats, StdExpansion *E, Array<OneD, Array<OneD, NekDouble> >&retarr, NekDouble &avgiterGD, NekDouble &totalcallstobasis, NekDouble tol, NekDouble nu)
  {
    boost::ignore_unused(tol);
    Timer tmr;
    tmr.Start();
    NekDouble totalbasiscalls = 0;  
    int dim = 2;
    Array<OneD, NekDouble> g(dim);
    double inf = numeric_limits<double>::infinity();
    retarr = Array<OneD, Array<OneD, NekDouble> >(dim);
    for(int p = 0; p < dim; p++)
      {
	retarr[p] = Array<OneD,NekDouble>(1);
      }
    // flag that determines if GD is called from GDBackTrackingTest.cpp
    int stype = E->DetShapeType();
    Array<OneD, Array<OneD, NekDouble> > storage, coords, coordsmid;
    Array<OneD, NekDouble> midpteval;
    if(stype == 3)
      {
	storage = storage2dt;
	coords = coordtri;	
	midpteval = midptevaltri;
	coordsmid = coordmidtri;
      }
    else
      {
	storage = storage2dq;
	coords = coordquad;
	midpteval = midptevalquad;
	coordsmid = coordmidquad;
      }
    Array<OneD, NekDouble> interioreval = midpteval;

    Array<OneD,NekDouble> xnew(dim) ,x(dim);
    NekDouble gprev = inf;
    Array<OneD, NekDouble> nullarr(0), temp2(coords[0].size()), temp3, alltemp;
    NekDouble gnew = inf;

    int idxgprev;
        
    int sz = coords[0].size();
    Array<OneD, Array<OneD, NekDouble > > tempeval(4);
    Array<OneD, NekDouble> xstart(dim);
 
    for(int i = 0; i < 4; i++)
      {
        tempeval[i] = Array<OneD, NekDouble>(E->GetNcoeffs());
      }

    pq(uhats, coords , storage[0], nullarr, temp2, nu);
    
    gprev = Vmath::Vmin(sz, temp2, 1);
    idxgprev = Vmath::Imin(temp2.size(), temp2, 1);
    for(int p = 0; p < dim; p++)
      {
	xnew[p] = coords[p][idxgprev];
	xstart[p] = xnew[p];  
      }
    
    pq(uhats, coordsmid, midpteval, nullarr, temp3, nu );

    NekDouble gprevmid =  Vmath::Vmin(temp3.size(), temp3, 1);

    if(gprev > gprevmid)
      {
        gprev = gprevmid;
        idxgprev = Vmath::Imin(temp3.size(), temp3, 1);
	for(int p = 0; p < dim; p++)
	  {
	    xnew[p] = coordsmid[p][idxgprev];  
	    xstart[p] = xnew[p]; 
	  }
      }
    if(startarr.size() > 1)
      {
	cout<<"\n expected root roughly = "<< startarr[0]<<" , "<<startarr[1]<<"\n";
	
      }
    
    Array<OneD, NekDouble> dereval(dim);
    if(gprev < 0 && abs(gprev)>tol_gd)
      {
	
	Array<OneD, Array<OneD, NekDouble > > xastaa(2);
	xastaa[0] = Array<OneD, NekDouble> (1, xnew[0]);
	xastaa[1] = Array<OneD, NekDouble> (1, xnew[1]);
	
	NekDouble c = chold;
	tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2], tempeval[3]);
	//tot_interp_points += xastaa[0].size();

	totalbasiscalls++;
	for(int p = 0; p < dim; p++)
	  {
	    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
	    g[p] = dereval[p];
	  }
	gnew = pow((pow(g[0],2) + pow(g[1],2)) ,0.5);
	if(abs(g[0]) > 1 || abs(g[1])>1)
	  {
	    g[0] = g[0]/gnew;
	    g[1] = g[1]/gnew;
	  }

	Array<OneD, NekDouble> t1(1);
	pq(uhats, xastaa, tempeval[0], NullNekDouble1DArray, t1, nu);
	Array<OneD, NekDouble > holdxandval(dim+1), saveholdxandval(dim+1), savesavehold(dim+1);
	holdxandval[0] = xastaa[0][0];
	holdxandval[1] = xastaa[1][0];
	holdxandval[2] = t1[0];
	savesavehold[2] = inf;
	int truth_val = abs(g[0]+g[1]) > tol;
	int ctr = 0, ct  = 0;
	NekDouble savevalprev0 = 0, savevalprev1 = 0, saveholdval = 0, fac = 1,  gnew0 = gnew, gnew1 = gnew, gnew2 = gnew;
    		    	
	if(truth_val)
	  {
	    NekDouble iter = secarg;
	    int counter_min_change = 0, max_backtracksteps = 5;	

	    while(ctr < iter && ((abs(g[0]) > tol_gd) || abs(g[1]) > tol_gd) && fac > 1e-6  && counter_min_change < max_backtracksteps )
	      {
		if(ctr > 2 && abs(gnew0 - gnew2) < tol_gd)
		  {
		    break;
		  }
		if(ctr > 1)
		  {
		    fac = gamhold;
		  }
		saveholdxandval[0] = xastaa[0][0];
		saveholdxandval[1] = xastaa[1][0];
		saveholdxandval[2] = t1[0];
	    
		xastaa[0][0] = saveholdxandval[0]  - fac*g[0];
		xastaa[1][0] = saveholdxandval[1]  - fac*g[1];
	    
		bool chkval = abs(xastaa[0][0] )>1 ||abs( xastaa[1][0]) > 1;
		if(stype == 3)
		  {
		    chkval = chkval ||  abs(xastaa[0][0] + xastaa[1][0]) > 1;
		  }
		while(chkval)
		  {
		    xastaa[0][0] = holdxandval[0]  - fac*g[0];
		    xastaa[1][0] = holdxandval[1]  - fac*g[1];
		    fac = fac*gamhold;
		    if(fac < 1e-6)
		       break;
		    chkval = abs(xastaa[0][0] )>1 ||abs( xastaa[1][0]) > 1;
		    if(stype == 3)
		      {
			chkval = chkval ||  abs(xastaa[0][0] + xastaa[1][0]) > 1;
		      }
		  }
		
		ctr++;
		//tmr.Start();
		tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2], tempeval[3]);
		//tot_interp_points += xastaa[0].size();

		totalbasiscalls++;
		//tmr.Stop();
		//cout<<"\n 2D eval basis 2 call takes = "<<tmr.TimePerTest(1);
		pq(uhats, xastaa, tempeval[0], NullNekDouble1DArray, t1, nu);
		for(int p = 0; p < dim; p++)
		  {
		    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
		    g[p] = dereval[p];
		  }
	    
		if(savesavehold[dim] > t1[0])
		  {
		    for(int j = 0; j < dim; j++)
		      {
			savesavehold[j] = xastaa[j][0];
		      }
		    savesavehold[dim] = t1[0];
		  }
	    
		gnew2 = gnew1;
		gnew1 = gnew0;
		gnew = pow((pow(g[0],2) + pow(g[1],2)) ,0.5);
		gnew0 = gnew;
		if(abs(g[0]) > 1 || abs(g[1])>1)
		  {
		    g[0] = g[0]/gnew;
		    g[1] = g[1]/gnew;
		  }
		
		savevalprev0 = savevalprev1;
		savevalprev1 = saveholdval;
		saveholdval = t1[0];
		NekDouble holdval = t1[0];
		Array<OneD, NekDouble> gsave(2);
		gsave[0] = g[0];
		gsave[1] = g[1];
		ct = 0;
		while(holdval > saveholdxandval[2] - c*fac*(gsave[0] + gsave[1]) && ct < iter &&
		      fac > 1e-6)
		  {
		    ct++;
		    holdxandval[0] = xastaa[0][0];
		    holdxandval[1] = xastaa[1][0];
		    holdxandval[2] = holdval;
		
		    fac = fac*gamhold;
		    xastaa[0][0] = holdxandval[0]  - fac*g[0];
		    xastaa[1][0] = holdxandval[1]  - fac*g[1];
		    chkval = abs(xastaa[0][0] )>1 ||abs( xastaa[1][0]) > 1;
		    if(stype == 3)
		      {
			chkval = chkval ||  abs(xastaa[0][0] + xastaa[1][0]) > 1;
		      }
		    while(chkval)
		      {
			xastaa[0][0] = holdxandval[0]  - fac*g[0];
			xastaa[1][0] = holdxandval[1]  - fac*g[1];
			fac = fac*gamhold;
			if(fac < 1e-6)
			  break;
			chkval = abs(xastaa[0][0] )>1 ||abs( xastaa[1][0]) > 1;
			if(stype == 3)
			  {
			    chkval = chkval ||  abs(xastaa[0][0] + xastaa[1][0]) > 1;
			  }
		      }
		    // tmr.Start();
		    tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2]
						       , tempeval[3]);
		    //tot_interp_points += xastaa[0].size();
		    totalbasiscalls++;
		    //tmr.Stop();
		    //cout<<"\n 2D eval basis 3 call takes: "<<tmr.TimePerTest(1);
		    for(int p = 0; p < dim; p++)
		      {
			derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
			g[p] = nu*dereval[p];
		      }
		    gnew2 = gnew1;
		    gnew1 = gnew0;
		    gnew = pow((pow(g[0],2) + pow(g[1],2)) ,0.5);
		    gnew0 = gnew;
		
		    if(abs(g[0]) > 1 || abs(g[1])>1)
		      {
			g[0] = g[0]/gnew;
			g[1] = g[1]/gnew;
		      }
		    pq(uhats, xastaa, tempeval[0], NullNekDouble1DArray, t1, nu);
		    if(savesavehold[dim] > t1[0])
		      {
			for(int j = 0; j < dim; j++)
			  savesavehold[j] = xastaa[j][0];
			savesavehold[dim] = t1[0];
		      }
		
		    if(abs(holdval -t1[0]) > 1e-5)
		      holdval = t1[0];
		    else
		      break;
		  }
		avgiterGD = ctr;

		if(saveholdval < holdval)
		  {
		    
		    if((saveholdxandval[2]<savesavehold[2]) && abs(saveholdxandval[2] - savesavehold[2]) > tol_gd)
		      {
			savesavehold[0] = saveholdxandval[0];
			savesavehold[1] = saveholdxandval[1];
			savesavehold[2] = saveholdxandval[2];
			counter_min_change = 0;
		      }
		    else
		      {
			counter_min_change++;
		      }
		  
		    xastaa[0][0] = saveholdxandval[0];
		    xastaa[1][0] = saveholdxandval[1];
		
		    t1[0] = saveholdval;
		  }
		else
		  {
		   
		    if((saveholdxandval[2]<savesavehold[2]) && abs(saveholdxandval[2] - savesavehold[2]) > tol_gd)
		      {
			savesavehold[0] = saveholdxandval[0];
			savesavehold[1] = saveholdxandval[1];
			savesavehold[2] = saveholdxandval[2];
			counter_min_change = 0;
		      }
		    else
		      {
		    	counter_min_change++;
		      }
		
		    t1[0] = holdval;
		    if(abs(xastaa[0][0])>1 || abs(xastaa[1][0]) > 1)
		      {
			xastaa[0][0] = savesavehold[0];
			xastaa[1][0] = savesavehold[1];
		      }
		    fac = 1;
		    
		    if(ctr > 2 && abs(savevalprev0 - saveholdval) < 1e-5)
		      {
		    	break;
		      }
		  }
	      }
	  }
	    
	else
	  {
	    savesavehold[0] = xastaa[0][0];
	    savesavehold[1] = xastaa[1][0];
	  }

	retarr = Array<OneD, Array<OneD, NekDouble> >(2);
	retarr[0] = Array<OneD, NekDouble>(1, savesavehold[0]);
	retarr[1] = Array<OneD, NekDouble>(1, savesavehold[1]);
	avgiterGD = ctr;
	tmr.Stop();
	//	totalbasiscalls = tmr.TimePerTest(1);
	//	cout<<"\n iters taken by 2D GD during 1 call = "<<ctr;
	totalcallstobasis = totalbasiscalls; 
	return;	
      }
    else
      {
	//	cout<<"\n returning NULL after ctr = "<<avgiterGD;
	tmr.Stop();
	//totalbasiscalls = tmr.TimePerTest(1);
	//	cout<<"\n iters taken by 2D GD during 1 call = "<<ctr;
	totalcallstobasis = totalbasiscalls; 
	
	retarr = NullNekDoubleArrayOfArray;
	return;
      }
    
  }

  

  //GD with backtracking
  // void steepestgradientdescent3D(Array<OneD, NekDouble> &uhats, StdExpansion *E, Array<OneD, Array<OneD, NekDouble> > &storage, Array<OneD, Array<OneD, NekDouble> > &coordpts, Array<OneD, Array<OneD, NekDouble> > &coordmidpts, Array<OneD, NekDouble> &interioreval, Array<OneD, Array<OneD, NekDouble> >&retarr, NekDouble &avgiterGD, NekDouble &totalcallstobasis, NekDouble tol, NekDouble nu, int   &tot_interp_points)
  // {
  //   Timer tbeg;
  //   tbeg.Start();
  //   NekDouble totalbasiscalls = 0;
  //   int dim = 3;
  //   Array<OneD, NekDouble> g(dim);
  //   double inf = numeric_limits<double>::infinity();
  //   retarr = Array<OneD, Array<OneD, NekDouble> >(dim);

  //   for(int p = 0; p < dim; p++)
  //     {
  // 	retarr[p] = Array<OneD,NekDouble>(1);
  //     }

    
  //   Array<OneD,NekDouble> xnew(dim) ,x(dim);
  //   NekDouble gprev = inf;
  //   Array<OneD, NekDouble>  temp3, temp2(coordpts[0].size()), xstart(dim);
  //   NekDouble gnew = inf;

  //   int idxgprev;
  //   int sz = coordpts[0].size();
  //   Array<OneD, Array<OneD, NekDouble > > tempeval(4);
  //   for(int k = 0; k < tempeval.size(); k++)
  //     {
  // 	tempeval[k] = Array<OneD, NekDouble>(uhats.size());
  //     }
  //   //NekDouble tmpq = 0;
  //   //tbeg.Start();
  //   pq(uhats, coordpts, storage[0], NullNekDouble1DArray, temp2, nu );
  //   //tbeg.Stop();
  //   //tmpq += tbeg.TimePerTest(1);
    
  //   gprev = Vmath::Vmin(sz, temp2, 1);
  //   idxgprev = Vmath::Imin(temp2.size(), temp2, 1);
    
  //   for(int p = 0; p < dim; p++)
  //     {
  // 	xnew[p] = coordpts[p][idxgprev];
  // 	xstart[p] = xnew[p];   	    
  //     }
    
  //   //tbeg.Start();     
  //   pq(uhats, coordmidpts, interioreval, NullNekDouble1DArray, temp3, nu );
  //   //tbeg.Stop();
  //   // tmpq += tbeg.TimePerTest(1);   
  //   NekDouble gprevmid =  Vmath::Vmin(temp3.size(), temp3, 1);
  //   if(gprev > gprevmid)
  //     {
  //       gprev = gprevmid;
  //       idxgprev = Vmath::Imin(temp3.size(), temp3, 1);
  // 	for(int p = 0; p < dim; p++)
  // 	  {
  // 	    xnew[p] = coordmidpts[p][idxgprev];
  // 	    xstart[p] = xnew[p];   	    
    
  // 	  }

  //     }

  //   // if(startarr.size() > 0)
  //     //{
  // 	//	    	cout<<"\nexpected root roughly = "<< startarr[0]<<" , "<<startarr[1]<<", "<<startarr[2]<<" val = "<<startarr[3]<< " GD starts at: "<<xnew[0]<<","<<xnew[1]<<","<<xnew[2];

  // 	//	flagTest = 0;
  // 	//}
  //   // tbeg.Stop();
  //   // cout<<"\n (re) eval pq on lattice takes "<<tbeg.TimePerTest(1);
  //   Array<OneD, NekDouble> dereval(dim);
  //   //    cout<<"\n gprev = "<<gprev<<" tol_gd = "<<tol_gd<<"\n";
  //   if(gprev < 0 && abs(gprev)>tol_gd) 
  //     {
  //       Array<OneD, Array<OneD, NekDouble > > xastaa(dim);
  // 	for(int p = 0; p < dim; p++)
  // 	  {
  // 	    xastaa[p] = Array<OneD, NekDouble> (1, xnew[p]);
  // 	  }
	
  // 	NekDouble c = chold;
  // 	tbeg.Start();
  // 	tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2], tempeval[3]);
  // 	tot_interp_points += xastaa[0].size();

  // 	//	totalbasiscalls++;
  // 	tbeg.Stop();
  // 	totalbasiscalls += tbeg.TimePerTest(1);

  // 	//tbeg.Stop();
  // 	//cout<<"\n 3D eval basis 1 call takes "<<tbeg.TimePerTest(1)<<"\n\n";
  // 	//exit(0);
  // 	for(int p = 0; p < dim; p++)
  // 	  {
  // 	    //tbeg.Start();
  // 	    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
  // 	    //tbeg.Stop();
  // 	    //tmpq += tbeg.TimePerTest(1);
  // 	    g[p] = nu*dereval[p];
  // 	  }
  // 	gnew = 0;
  // 	for(int p = 0; p < dim; p++)
  // 	  {
  // 	    gnew += pow(g[p],2);
  // 	  }
	    
  // 	gnew = pow( gnew, 0.5);
  // 	if(abs(g[0]) > 1 || abs(g[1]) > 1 || abs(g[2]) > 1)
  //         {
  //           g[0] = g[0]/gnew;
  //           g[1] = g[1]/gnew;
  // 	    g[2] = g[2]/gnew;
  // 	  }
  // 	Array<OneD, NekDouble> t1(1);
  // 	Timer tmr;
  // 	pq(uhats, xastaa, tempeval[0], NullNekDouble1DArray, t1, nu);
  // 	Array<OneD, NekDouble > holdxandval(dim+1), saveholdxandval(dim+1), savesavehold(dim+1), minall(dim+1);

  // 	for(int p = 0; p < dim; p++)
  // 	  {
  // 	    holdxandval[p] = xastaa[p][0];
  // 	  }
  // 	holdxandval[dim] = t1[0];
  // 	savesavehold[dim] = inf;//t1[0];
  // 	NekDouble allg;
  // 	for(int p = 0; p < dim; p++)
  // 	  {
  // 	    allg = allg + g[p];
  // 	  }
  // 	int truth_val = abs(allg) > tol;
  // 	int ctr = 0, ct  = 0, max_backtracksteps = 1;
  // 	NekDouble savevalprev0 = 0, savevalprev1 = 0, saveholdval = 0, fac = 1,  gnew0 = gnew, gnew1 = gnew, gnew2 = gnew ;

  // 	if(truth_val)
  // 	  {
  // 	    NekDouble iter = secarg;
  // 	    int counter_min_change = 0; 
  // 	    while(ctr < iter  && fac > 1e-3 && ((abs(g[0]) > tol_gd) || abs(g[1]) > tol_gd ||  abs(g[2]) > tol_gd) && counter_min_change <max_backtracksteps )
  // 	      {
		
  // 		if(ctr > 2 && abs(gnew0 - gnew2) < tol_gd)// && abs(gnew0 - gnew1) < tol_gd && abs(gnew1- gnew2) < tol_gd)
  // 		  {
  // 		    break;
  // 		  }
		
  // 		if(ctr > 1)
  // 		  {
  // 		    fac = gamhold;
  // 		  }
  // 		for(int p = 0; p < dim; p++)
  // 		  {
  // 		    saveholdxandval[p] = xastaa[p][0];
  // 		  }
  // 		saveholdxandval[dim] = t1[0];
  // 		for(int p = 0; p < dim; p++)
  // 		  {
  // 		    xastaa[p][0] = saveholdxandval[p]  - fac*g[p];
  // 		  }
	        
  // 		// switch(E->DetShapeType())
  // 		//   {
  // 		//   case LibUtilities::eHexahedron:
  // 		    while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1)
  // 		      {
  // 			for(int p = 0; p < dim; p++)
  // 			  {
  // 			    xastaa[p][0] = holdxandval[p]  - fac*g[p];
  // 			  }
  // 			fac = fac*gamhold;
  // 			if(fac < 1e-4)
  // 			  {
  // 			    break;
  // 			  }
  // 		      }
  // 		  //   break;
  // 		  // case LibUtilities::eTetrahedron:
  // 		  //   // while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[2][0] + xastaa[1][0] > 1e-8 || xastaa[0][0] + xastaa[2][0] > 1e-8 || xastaa[1][0] + xastaa[2][0] > 1e-8)
  // 		  //   while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[2][0] + xastaa[1][0] > 0 || xastaa[0][0] + xastaa[2][0] > 0 || xastaa[1][0] + xastaa[2][0] > 0)
  // 		  //     {
  // 		  // 	for(int p = 0; p < dim; p++)
  // 		  // 	  {
  // 		  // 	    xastaa[p][0] = holdxandval[p]  - fac*g[p];
  // 		  // 	  }
  // 		  // 	fac = fac*gamhold;
  // 		  // 	if(fac <1e-7)
  // 		  // 	  break;
  // 		  //     }
  // 		  //   break;
  // 		  // case LibUtilities::ePyramid:
  // 		  //   while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[0][0] + xastaa[1][0] > 0 ||  xastaa[0][0] + xastaa[2][0] > 0 ) // x+y+z > 0??
  // 		  //     {
  //                 //       for(int p = 0; p < dim; p++)
  //                 //         {
  //                 //           xastaa[p][0] = holdxandval[p]  - fac*g[p];
  //                 //         }
  //                 //       fac = fac*gamhold;
  // 		  // 	if(fac < 1e-7)
  // 		  // 	break;
  // 		  //     }
  // 		  //   break;
  // 		  // case LibUtilities::ePrism:
  // 		  //   while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[0][0] + xastaa[2][0] > 1e-8 )
  // 		  //     {
  //                 //       for(int p = 0; p < dim; p++)
  //                 //         {
  //                 //           xastaa[p][0] = holdxandval[p]  - fac*g[p];
  //                 //         }
  //                 //       fac = fac*gamhold;
  // 		  // 	if(fac < 1e-7)
  // 		  // 	break;
  // 		  //     }

  // 		  //   break;
  // 		  // default: cout<<"\n invalid element, not yet supported!";
  // 		  //   exit(0);
  // 		  // }
  // 		ctr++;
  // 		tmr.Start();
  // 		tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2], tempeval[3]);
  // 		tot_interp_points += xastaa[0].size();

  // 		//totalbasiscalls++;
  // 		// tmr.Stop();
  // 		    tmr.Stop();
  // 		    totalbasiscalls += tmr.TimePerTest(1);

  // 		// cout<<"\n 3D eval basis 2 call takes: "<<tmr.TimePerTest(1);
  // 		//tbeg.Start();
  // 		pq(uhats, xastaa, tempeval[0], NullNekDouble1DArray, t1, nu);
  // 		//	tbeg.Stop();
  // 		//tmpq += tbeg.TimePerTest(1);
  // 		if(savesavehold[dim] > t1[0])
  // 		  {
  // 		    minall[dim] = savesavehold[dim];
  // 		    for(int j = 0; j < dim; j++)
  // 		      {
  // 			minall[j] = savesavehold[j];
  // 			savesavehold[j] = xastaa[j][0];
  // 		      }
  // 		    savesavehold[dim] = t1[0];
  // 		  }
  // 		for(int p = 0; p < dim; p++)
  // 		  {
  // 		    //tbeg.Start();
  // 		    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
  // 		    //tbeg.Stop();
  // 		    // tmpq += tbeg.TimePerTest(1);
  // 		    g[p] = dereval[p];
  // 		  }

  // 		gnew2 = gnew1;
  // 		gnew1 = gnew0;
  // 		gnew = 0;
  // 		for(int p = 0; p < dim; p++)
  // 		  {
  // 		    gnew += pow(g[p],2);
  // 		  }

  // 		gnew = pow( gnew, 0.5);
  // 		gnew0 = gnew;
	        
  // 		if(abs(g[0]) > 1 || abs(g[1])>1 || abs(g[2])>1)
  // 		  {
  // 		    for(int p = 0; p < dim; p++)
  // 		      {
  // 			g[p] = g[p]/gnew;
  // 		      }
  // 		  }
		
  // 		savevalprev0 = savevalprev1;
  // 		savevalprev1 = saveholdval;
  // 		saveholdval = t1[0];
  // 		NekDouble holdval = t1[0];
  // 	        Array<OneD, NekDouble> gsave(dim);

  // 		for(int p = 0; p < dim; p++)
  // 		  {
  // 		    gsave[p] = g[p];
  // 		  }
  // 		ct = 0;
  // 		while(holdval > saveholdxandval[dim] - c*fac*(gsave[0] + gsave[1] + gsave[2]) && ct < iter && fac > 1e-3)
  // 		  {
  // 		    ct++;
  // 		    for(int p = 0; p < dim; p++)
  // 		      {
  // 			holdxandval[p] = xastaa[p][0];
  // 		      }
  // 		    holdxandval[dim] = holdval;

  // 		    fac = fac*gamhold;
  // 		    for(int p = 0; p < dim; p++)
  // 		      {
  // 			xastaa[p][0] = holdxandval[p]  - fac*g[p];
  // 		      }
		    
  // 		    // switch(E->DetShapeType())
  // 		    //   {
  // 		    //		      case LibUtilities::eHexahedron:
  // 			while((abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1) )
  // 			  {
  // 			    for(int p = 0; p < dim; p++)
  // 			      {
  // 				xastaa[p][0] = holdxandval[p]  - fac*g[p];
				
  // 			      }
  // 			    fac = fac*gamhold;
  // 			    if(fac < 1e-4)
  // 			      {
  // 				break;
  // 			      }
  // 			  }
  // 		      // 	break;
  // 		      // case LibUtilities::eTetrahedron:
		        
  // 		      // 	while((abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[0][0] + xastaa[1][0] > 0 || xastaa[0][0] + xastaa[2][0] > 0 || xastaa[1][0] + xastaa[2][0] > 0)) 
  // 		      // 	  {
  // 		      // 	    for(int p = 0; p < dim; p++)
  // 		      // 	      {
  // 		      // 		xastaa[p][0] = holdxandval[p]  - fac*g[p];
  // 		      // 	      }
			
  // 		      // 	    fac = fac*gamhold;
  // 		      // 	if(fac < 1e-7)
  // 		      // 	  break;
  // 		      // 	  }
  // 		      // 	break;
  // 		      // 	case LibUtilities::ePyramid:
  // 		      // 	  while((abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1  ||  xastaa[0][0] + xastaa[2][0] > 0 || xastaa[0][0] + xastaa[1][0] > 0) )
  // 		      // 	    {
  // 		      // 	      for(int p = 0; p < dim; p++)
  // 		      // 		{
  // 		      // 		  xastaa[p][0] = holdxandval[p]  - fac*g[p];
				  
  // 		      // 		}
  // 		      // 	      fac = fac*gamhold;
  // 		      // 	      if(fac < 1e-7)
  // 		      // 		break;
				  
  // 		      // 	    }
			    
  // 		      // 	  break;
  // 		      // case LibUtilities::ePrism:

  // 		      // 	while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[0][0] + xastaa[2][0] > 0 )
  // 		      // 	  {
  // 		      // 	    for(int p = 0; p < dim; p++)
  // 		      // 	      {
  // 		      // 		xastaa[p][0] = holdxandval[p]  - fac*g[p];
  // 		      // 	      }
  // 		      // 	    fac = fac*gamhold;
  // 		      // 	    if(fac < 1e-7)
  // 		      // 	      break;
  // 		      // 	  }
  // 		      // 	break;
  // 		      // default:
  // 		      // 	cout<<"\n invalid element!\n";
  // 		      // 	exit(0);
  // 		      // }
  // 		    //    if(fac < 1e-6)
  // 		    //  break;
		    
  // 			tmr.Start();
  // 		    tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2]
  // 						       , tempeval[3]);
  // 		    tot_interp_points += xastaa[0].size();

  // 		    //totalbasiscalls++;
  // 		    tmr.Stop();
  // 		    totalbasiscalls += tmr.TimePerTest(1);
		    
  // 		    // cout<<"\n 3D eval basis 3 call takes: "<<tmr.TimePerTest(1); 

		    
  // 		    for(int p = 0; p < dim; p++)
  // 		      {
  // 			//	tbeg.Start();
  // 			derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
  // 			//tbeg.Stop();
  // 			//tmpq += tbeg.TimePerTest(1);
  // 			g[p] = nu*dereval[p];
  // 		      }
  // 		    gnew2 = gnew1;
  // 		    gnew1 = gnew0;

  // 		    gnew = 0;
  // 		    for(int p = 0; p < dim; p++)
  // 		      {
  // 			gnew += pow(g[p],2);
  // 		      }
		    
  // 		    gnew = pow(gnew ,0.5);
  // 		    gnew0 = gnew;

		    
  // 		    if(abs(g[2])>1 || abs(g[0]) > 1 || abs(g[1])>1)
  // 		      {
  // 			for(int p = 0; p < dim; p++)
  // 			  {
  // 			    g[p] = g[p]/gnew;
  // 			  }
  // 		      }

  // 		    //tbeg.Start(); 
  // 		    pq(uhats, xastaa, tempeval[0], NullNekDouble1DArray, t1, nu);
  // 		    //tbeg.Stop();
  // 		    // tmpq += tbeg.TimePerTest(1);  
  // 		    if(savesavehold[dim] > t1[0])
  // 		      {
  // 			minall[dim] = savesavehold[dim];
  // 			for(int j = 0; j < dim; j++)
  // 			  {
  // 			    minall[j] = savesavehold[j];
  // 			    savesavehold[j] = xastaa[j][0];
  // 			  }
  // 			savesavehold[dim] = t1[0];
  // 		      }
  // 	            if(abs(holdval -t1[0]) > 1e-5)
  // 		      {
  // 			holdval = t1[0];
  // 		      }
  // 		    else
  // 		      break;
  // 		  }
		
  // 		avgiterGD = ctr;
  // 		// if(saveholdval < holdval)
  // 		//   {
  // 		//     if((saveholdxandval[dim]<savesavehold[dim]) && abs(saveholdxandval[dim] - savesavehold[dim] >tol_gd))//ctr > 1 && minall[dim] <  savesavehold[dim] && abs(minall[dim] -  savesavehold[dim]) >tol_gd)
  // 		//       {
  // 		// 	for(int p = 0; p < dim; p++)
  // 		// 	  { 
  // 		// 	    savesavehold[p] = saveholdxandval[p];
  // 		// 	  }
  // 		// 	//savesavehold[dim] = minall[dim];
  // 		// 	counter_min_change = 0;
  // 		//       }
  // 		//     else
  // 		//       {
  // 		// 	counter_min_change++;
  // 		//       }
  // 		//     for(int p = 0; p < dim; p++)
  // 		//       {
  // 		// 	xastaa[p][0] = saveholdxandval[p];
  // 		//       }
  // 		//     // t1[0] = holdval;
  // 		//     t1[0] = saveholdval;
  // 		//   }
  // 		// else
  // 		//   {
  // 		//     if((saveholdxandval[dim]<savesavehold[dim]) && abs(saveholdxandval[dim] - savesavehold[dim])>tol_gd)
  // 		//       {
  // 		// 	for(int p = 0; p < dim+1; p++)
  // 		// 	  {
  // 		// 	    savesavehold[p] = saveholdxandval[p]; 
  // 		// 	  }
  // 		// 	counter_min_change = 0;      
  // 		//       }
  // 		//     else
  // 		//       {
  // 		// 	counter_min_change++; 
  // 		//       }
  // 		//     t1[0] = holdval;
		    
		    
  // 		//     fac = 1;
		    
  // 		//     if(ctr > 2 && abs(savevalprev0 - saveholdval) < 1e-5)
  // 		//       {
  // 		// 	break;
  // 		//       }
		
  // 		//   }

  // 	      }
  // 	    		    int opt = 0;
  // 			opt = abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0])>1;
  // 			if(opt)
  // 		      {
			
  // 			cout<<"\n opt = "<<opt<<"\n";
  // 		cout<<"\n retarr =  ("<<retarr[0][0]<<","<<retarr[1][0]<<","<<retarr[2][0]<<")";
  // 			cout<<"\n savesavehold =  ("<<savesavehold[0]<<","<<savesavehold[1]<<","<<savesavehold[2][0]<<")";
	
  // 			for(int p = 0; p < dim; p++)
  // 			  {
  // 			    xastaa[p][0] = savesavehold[p];
  // 			  }
  // 		      }


  // 	  }
  // 	else
  // 	  {
  // 	    cout<<"\n in 3d GD else\n";
  // 		    int opt = 0;
  // 		    // switch(E->DetShapeType())
  // 		    //   {
  // 		    //  case LibUtilities::eHexahedron:
  // 			opt = abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0])>1;
  // 			//break;
			
  // 		      // case LibUtilities::eTetrahedron:
  // 		      // 	opt = abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1|| xastaa[0][0] + xastaa[1][0] > 0 || xastaa[0][0] + xastaa[2][0] > 0 || xastaa[1][0] + xastaa[2][0] > 0;
  // 		      // 	break;
			
  // 		      // case LibUtilities::ePyramid:
  // 		      // 	opt = abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1|| xastaa[0][0] + xastaa[2][0] > 0 || xastaa[2][0] + xastaa[1][0] > 0;  
  // 		      // 	break;
  // 		      // case LibUtilities::ePrism:
  // 		      // 	opt = abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1|| xastaa[0][0] + xastaa[2][0] > 0;
  // 		      // 	break;
  // 		      // default: cout<<"\n invalid element type!\n";
  // 		      // 	exit(0);
  // 		      // }
		    
  // 		    if(opt)
  // 		      {
  // 			cout<<"\n opt = "<<opt<<"\n";
  // 			cout<<"\n retarr =  ("<<retarr[0][0]<<","<<retarr[1][0]<<","<<retarr[2][0]<<")";
  // 			cout<<"\n savesavehold =  ("<<savesavehold[0][0]<<","<<savesavehold[1][0]<<","<<savesavehold[2][0]<<")";
	
  // 			for(int p = 0; p < dim; p++)
  // 			  {
  // 			    xastaa[p][0] = savesavehold[p];
  // 			  }
  // 		      }

  // 	  //   for(int p = 0; p < dim; p++)
  // 	  //     savesavehold[p] = xastaa[p][0];
  // 	   }
  // 	//	cout<<"\n at ctr = "<<ctr<<"  returning: "; 
  // 	for(int p = 0; p < dim; p++)
  // 	  {
  // 	    retarr[p] = Array<OneD, NekDouble>(1, savesavehold[p]);
  // 	    //	    cout<<" "<<savesavehold[p]<<" ";
  // 	  }
  // 	//	cout<<"\n iters taken by GD during 1 call = "<<ctr <<" calls to basis eval = "<<totalbasiscalls<<" returning ("<<retarr[0][0]<<","<<retarr[1][0]<<","<<retarr[2][0]<<")";
  // 	//	totalcallstobasis = tbeg.TimePerTest(1);//totalbasiscalls;
  // 	avgiterGD = ctr;
  // 	//cout<<"\n time spent inside pq and derpq="<<tmpq<<"\n";
  // 	//exit(0);
  // 	tbeg.Stop();
  // 	totalcallstobasis = totalbasiscalls;
	
  // 		cout<<"\n number of iteration taken by 3D GD ="<<ctr;
  // 	return;
  //     }
  //   else
  //     {
  // 	retarr = NullNekDoubleArrayOfArray;
  // 		cout<<"\n returning NULL after ctr = "<<avgiterGD;
  // 	//cout<<"\n time spent inside pq and derpq="<<tmpq<<"\n";
  // 	tbeg.Stop();
  // 	totalcallstobasis = totalbasiscalls;
  // 		cout<<"\n number of iteration taken by 3D GD =0";
  // 	//exit(0);
  // 	return;
  //     }
    
  // }
    // void steepestgradientdescent3D(Array<OneD, NekDouble> &uhats, StdExpansion *E, Array<OneD, Array<OneD, NekDouble> > &storage, Array<OneD, Array<OneD, NekDouble> > &coordpts, Array<OneD, Array<OneD, NekDouble> > &coordmidpts, Array<OneD, NekDouble> &interioreval, Array<OneD, Array<OneD, NekDouble> >&retarr, NekDouble &avgiterGD, NekDouble &totalcallstobasis, NekDouble tol, NekDouble nu, )
  //GD with backtracking
  void steepestgradientdescent3D(Array<OneD, NekDouble> &uhats, StdExpansion *E, Array<OneD, Array<OneD, NekDouble> > &storage, Array<OneD, Array<OneD, NekDouble> > &coordpts, Array<OneD, Array<OneD, NekDouble> > &coordmidpts, Array<OneD, NekDouble> &interioreval, Array<OneD, Array<OneD, NekDouble> >&retarr, NekDouble &avgiterGD, NekDouble &totalcallstobasis, NekDouble tol, NekDouble nu, int   &tot_interp_points)
  {
    boost::ignore_unused(tot_interp_points);
    Timer tbeg;
    // tbeg.Start();
    int totalbasiscalls = 0;
    int dim = 3;
    Array<OneD, NekDouble> g(dim);
    double inf = numeric_limits<double>::infinity();
    retarr = Array<OneD, Array<OneD, NekDouble> >(dim);

    for(int p = 0; p < dim; p++)
      {
	retarr[p] = Array<OneD,NekDouble>(1);
      }


    Array<OneD,NekDouble> xnew(dim) ,x(dim);
    NekDouble gprev = inf;
    Array<OneD, NekDouble>  temp3, temp2(coordpts[0].size()), xstart(dim);
    NekDouble gnew = inf;

    int idxgprev;
    int sz = coordpts[0].size();
    Array<OneD, Array<OneD, NekDouble > > tempeval(4);
    for(int k = 0; k < tempeval.size(); k++)
      {
	tempeval[k] = Array<OneD, NekDouble>(uhats.size());
      }

    pq(uhats, coordpts, storage[0], NullNekDouble1DArray, temp2, nu );

    gprev = Vmath::Vmin(sz, temp2, 1);
    idxgprev = Vmath::Imin(temp2.size(), temp2, 1);

    for(int p = 0; p < dim; p++)
      {
	xnew[p] = coordpts[p][idxgprev];
	xstart[p] = xnew[p];
      }


    pq(uhats, coordmidpts, interioreval, NullNekDouble1DArray, temp3, nu );

    NekDouble gprevmid =  Vmath::Vmin(temp3.size(), temp3, 1);
    if(gprev > gprevmid)
      {
	gprev = gprevmid;
	idxgprev = Vmath::Imin(temp3.size(), temp3, 1);
	for(int p = 0; p < dim; p++)
	  {
	    xnew[p] = coordmidpts[p][idxgprev];
	    xstart[p] = xnew[p];

	  }

      }

    if(startarr.size() > 0)
      {
	//    cout<<"\nexpected root roughly = "<< startarr[0]<<" , "<<startarr[1]<<", "<<startarr[2]<<" val = "<<startarr[3]<< " GD starts at: "<<xnew[0]<<","<<xnew[1]<<","<<xnew[2];

	//flagTest = 0;
      }
    // tbeg.Stop();
    // cout<<"\n (re) eval pq on lattice takes "<<tbeg.TimePerTest(1);
    Array<OneD, NekDouble> dereval(dim);
    //    cout<<"\n gprev = "<<gprev<<" tol_gd = "<<tol_gd<<"\n";
    if(gprev < 0 && abs(gprev)>tol_gd)
      {
	Array<OneD, Array<OneD, NekDouble > > xastaa(dim);
	for(int p = 0; p < dim; p++)
	  {
	    xastaa[p] = Array<OneD, NekDouble> (1, xnew[p]);
	  }

	NekDouble c = chold;
	//tbeg.Start();
	tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2], tempeval[3]);
	totalbasiscalls++;
	//tbeg.Stop();
	//cout<<"\n 3D eval basis 1 call takes "<<tbeg.TimePerTest(1)<<"\n\n";
	//exit(0);
	for(int p = 0; p < dim; p++)
	  {
	    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
	    g[p] = nu*dereval[p];
	  }
	gnew = 0;
	for(int p = 0; p < dim; p++)
	  {
	    gnew += pow(g[p],2);
	  }

	gnew = pow( gnew, 0.5);
	if(abs(g[0]) > 1 || abs(g[1]) > 1 || abs(g[2]) > 1)
	  {
	    g[0] = g[0]/gnew;
	    g[1] = g[1]/gnew;
	    g[2] = g[2]/gnew;
	  }
	Array<OneD, NekDouble> t1(1);
	Timer tmr;
	pq(uhats, xastaa, tempeval[0], NullNekDouble1DArray, t1, nu);
	Array<OneD, NekDouble > holdxandval(dim+1), saveholdxandval(dim+1), savesavehold(dim+1), minall(dim+1);

	for(int p = 0; p < dim; p++)
	  {
	    holdxandval[p] = xastaa[p][0];
	  }
	holdxandval[dim] = t1[0];
	savesavehold[dim] = inf;//t1[0];
	NekDouble allg = 0;
	for(int p = 0; p < dim; p++)
	  {
	    allg = allg + g[p];
	  }
	int truth_val = abs(allg) > tol;
	int ctr = 0, ct  = 0, max_backtracksteps = 1;
	NekDouble savevalprev0 = 0, savevalprev1 = 0, saveholdval = 0, fac = 1,  gnew0 = gnew, gnew1 = gnew, gnew2 = gnew ;
	if(truth_val)
	  {
	    NekDouble iter = secarg;
	    int counter_min_change = 0;

	    while(ctr < iter  && fac > 1e-7 && ((abs(g[0]) > tol_gd) || abs(g[1]) > tol_gd ||  abs(g[2]) > tol_gd) && counter_min_change <max_backtracksteps )
	      {

		if(ctr > 2 && abs(gnew0 - gnew2) < tol_gd)// && abs(gnew0 - gnew1) < tol_gd && abs(gnew1- gnew2) < tol_gd)
		  {
		    break;
		  }

		if(ctr > 1)
		  {
		    fac = gamhold;
		  }
		for(int p = 0; p < dim; p++)
		  {
		    saveholdxandval[p] = xastaa[p][0];
		  }
		saveholdxandval[dim] = t1[0];
		for(int p = 0; p < dim; p++)
		  {
		    xastaa[p][0] = saveholdxandval[p]  - fac*g[p];
		  }

		switch(E->DetShapeType())
		  {
		  case LibUtilities::eHexahedron:
		    while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1)
		      {
			for(int p = 0; p < dim; p++)
			  {
			    xastaa[p][0] = holdxandval[p]  - fac*g[p];
			  }
			fac = fac*gamhold;
			if(fac < 1e-7)
			  {
			    break;
			  }
		      }
		    break;
		  case LibUtilities::eTetrahedron:
		    // while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[2][0] + xastaa[1][0] > 1e-8 || xastaa[0][0] + xastaa[2][0] > 1e-8 || xastaa[1][0] + xastaa[2][0] > 1e-8)
		    while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[2][0] + xastaa[1][0] > 0 || xastaa[0][0] + xastaa[2][0] > 0 || xastaa[1][0] + xastaa[2][0] > 0)
		      {
			for(int p = 0; p < dim; p++)
			  {
			    xastaa[p][0] = holdxandval[p]  - fac*g[p];
			  }
			fac = fac*gamhold;
			if(fac <1e-7)
			  break;
		      }
		    break;
		  case LibUtilities::ePyramid:
		    while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[0][0] + xastaa[1][0] > 0 ||  xastaa[0][0] + xastaa[2][0] > 0 ) // x+y+z > 0??
		      {
			for(int p = 0; p < dim; p++)
			  {
			    xastaa[p][0] = holdxandval[p]  - fac*g[p];
			  }
			fac = fac*gamhold;
			if(fac < 1e-7)
			  break;
		      }
		    break;
		  case LibUtilities::ePrism:
		    while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[0][0] + xastaa[2][0] > 1e-8 )
		      {
			for(int p = 0; p < dim; p++)
			  {
			    xastaa[p][0] = holdxandval[p]  - fac*g[p];
			  }
			fac = fac*gamhold;
			if(fac < 1e-7)
			  break;
		      }

		    break;
		  default: cout<<"\n invalid element, not yet supported!";
		    exit(0);
		  }
		ctr++;
		//tmr.Start();
		tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2], tempeval[3]);
		totalbasiscalls++;
		// tmr.Stop();
		// cout<<"\n 3D eval basis 2 call takes: "<<tmr.TimePerTest(1);
		pq(uhats, xastaa, tempeval[0], NullNekDouble1DArray, t1, nu);

		if(savesavehold[dim] > t1[0])
		  {
		    minall[dim] = savesavehold[dim];
		    for(int j = 0; j < dim; j++)
		      {
			minall[j] = savesavehold[j];
			savesavehold[j] = xastaa[j][0];
		      }
		    savesavehold[dim] = t1[0];
		  }
		for(int p = 0; p < dim; p++)
		  {
		    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
		    g[p] = dereval[p];
		  }

		gnew2 = gnew1;
		gnew1 = gnew0;
		gnew = 0;
		for(int p = 0; p < dim; p++)
		  {
		    gnew += pow(g[p],2);
		  }

		gnew = pow( gnew, 0.5);
		gnew0 = gnew;

		if(abs(g[0]) > 1 || abs(g[1])>1 || abs(g[2])>1)
		  {
		    for(int p = 0; p < dim; p++)
		      {
			g[p] = g[p]/gnew;
		      }
		  }

		savevalprev0 = savevalprev1;
		savevalprev1 = saveholdval;
		saveholdval = t1[0];
		NekDouble holdval = t1[0];
		Array<OneD, NekDouble> gsave(dim);

		for(int p = 0; p < dim; p++)
		  {
		    gsave[p] = g[p];
		  }
		ct = 0;
		while(holdval > saveholdxandval[dim] - c*fac*(gsave[0] + gsave[1] + gsave[2]) && ct < iter && fac > 1e-7)
		  {
		    ct++;
		    for(int p = 0; p < dim; p++)
		      {
			holdxandval[p] = xastaa[p][0];
		      }
		    holdxandval[dim] = holdval;

		    fac = fac*gamhold;
		    for(int p = 0; p < dim; p++)
		      {
			xastaa[p][0] = holdxandval[p]  - fac*g[p];
		      }

		    switch(E->DetShapeType())
		      {
		      case LibUtilities::eHexahedron:
			while((abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1) )
			  {
			    for(int p = 0; p < dim; p++)
			      {
				xastaa[p][0] = holdxandval[p]  - fac*g[p];

			      }
			    fac = fac*gamhold;
			    if(fac < 1e-7)
			      {
				break;
			      }
			  }
			break;
		      case LibUtilities::eTetrahedron:

			while((abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[0][0] + xastaa[1][0] > 0 || xastaa[0][0] + xastaa[2][0] > 0 || xastaa[1][0] + xastaa[2][0] > 0))
			  {
			    for(int p = 0; p < dim; p++)
			      {
				xastaa[p][0] = holdxandval[p]  - fac*g[p];
			      }

			    fac = fac*gamhold;
			    if(fac < 1e-7)
			      break;
			  }
			break;
		      case LibUtilities::ePyramid:
			while((abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1  ||  xastaa[0][0] + xastaa[2][0] > 0 || xastaa[0][0] + xastaa[1][0] > 0) )
			  {
			    for(int p = 0; p < dim; p++)
			      {
				xastaa[p][0] = holdxandval[p]  - fac*g[p];

			      }
			    fac = fac*gamhold;
			    if(fac < 1e-7)
			      break;

			  }

			break;
		      case LibUtilities::ePrism:

			while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[0][0] + xastaa[2][0] > 0 )
			  {
			    for(int p = 0; p < dim; p++)
			      {
				xastaa[p][0] = holdxandval[p]  - fac*g[p];
			      }
			    fac = fac*gamhold;
			    if(fac < 1e-7)
			      break;
			  }
			break;
		      default:
			cout<<"\n invalid element!\n";
			exit(0);
		      }
		    //    if(fac < 1e-6)
		    //  break;

		    //    tmr.Start();
		    tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2]
						       , tempeval[3]);
		    totalbasiscalls++;
		    // tmr.Stop();
		    // cout<<"\n 3D eval basis 3 call takes: "<<tmr.TimePerTest(1);


		    for(int p = 0; p < dim; p++)
		      {
			derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
			g[p] = nu*dereval[p];
		      }
		    gnew2 = gnew1;
		    gnew1 = gnew0;

		    gnew = 0;
		    for(int p = 0; p < dim; p++)
		      {
			gnew += pow(g[p],2);
		      }

		    gnew = pow(gnew ,0.5);
		    gnew0 = gnew;


		    if(abs(g[2])>1 || abs(g[0]) > 1 || abs(g[1])>1)
		      {
			for(int p = 0; p < dim; p++)
			  {
			    g[p] = g[p]/gnew;
			  }
		      }


		    pq(uhats, xastaa, tempeval[0], NullNekDouble1DArray, t1, nu);

		    if(savesavehold[dim] > t1[0])
		      {
			minall[dim] = savesavehold[dim];
			for(int j = 0; j < dim; j++)
			  {
			    minall[j] = savesavehold[j];
			    savesavehold[j] = xastaa[j][0];
			  }
			savesavehold[dim] = t1[0];
		      }
		    if(abs(holdval -t1[0]) > 1e-5)
		      {
			holdval = t1[0];
		      }
		    else
		      break;
		  }
		avgiterGD = ctr;
		if(saveholdval < holdval)
		  {
		    if((saveholdxandval[dim]<savesavehold[dim]) && abs(saveholdxandval[dim] - savesavehold[dim] >tol_gd))//ctr > 1 && minall[dim] <  savesavehold[dim] && abs(minall[dim] -  savesavehold[dim]) >tol_gd)
		      {
			for(int p = 0; p < dim; p++)
			  {
			    savesavehold[p] = saveholdxandval[p];
			  }
			//savesavehold[dim] = minall[dim];
			counter_min_change = 0;
		      }
		    else
		      {
			counter_min_change++;
		      }
		    for(int p = 0; p < dim; p++)
		      {
			xastaa[p][0] = saveholdxandval[p];
		      }
		    // t1[0] = holdval;
		    t1[0] = saveholdval;
		  }
		else
		  {
		    if((saveholdxandval[dim]<savesavehold[dim]) && abs(saveholdxandval[dim] - savesavehold[dim])>tol_gd)
		      {
			for(int p = 0; p < dim+1; p++)
			  {
			    savesavehold[p] = saveholdxandval[p];
			  }
			counter_min_change = 0;
		      }
		    else
		      {
			counter_min_change++;
		      }
		    t1[0] = holdval;

		    int opt = 0;
		    switch(E->DetShapeType())
		      {
		      case LibUtilities::eHexahedron:
			opt = abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0])>1;
			break;

		      case LibUtilities::eTetrahedron:
			opt = abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1|| xastaa[0][0] + xastaa[1][0] > 0 || xastaa[0][0] + xastaa[2][0] > 0 || xastaa[1][0] + xastaa[2][0] > 0;
			break;

		      case LibUtilities::ePyramid:
			opt = abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1|| xastaa[0][0] + xastaa[2][0] > 0 || xastaa[2][0] + xastaa[1][0] > 0;
			break;
		      case LibUtilities::ePrism:
			opt = abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1|| xastaa[0][0] + xastaa[2][0] > 0;
			break;
		      default: cout<<"\n invalid element type!\n";
			exit(0);
		      }

		    if(opt)
		      {
			for(int p = 0; p < dim; p++)
			  {
			    xastaa[p][0] = savesavehold[p];
			  }
		      }

		    fac = 1;

		    if(ctr > 2 && abs(savevalprev0 - saveholdval) < 1e-5)
		      {
			break;
		      }

		  }
	      }

	  }
	else
	  {
	    for(int p = 0; p < dim; p++)
	      savesavehold[p] = xastaa[p][0];
	  }
	//cout<<"\n at ctr = "<<ctr<<"  returning: ";
	for(int p = 0; p < dim; p++)
	  {
	    retarr[p] = Array<OneD, NekDouble>(1, savesavehold[p]);
	    //    cout<<" "<<savesavehold[p]<<" ";
	  }
	//cout<<"\n iters taken by GD during 1 call = "<<ctr <<" calls to basis eval = "<<totalbasiscalls<<" returning ("<<retarr[0][0]<<","<<retarr[1][0]<<","<<retarr[2][0]<<")";
	totalcallstobasis = totalbasiscalls;
	avgiterGD = ctr;
	return;
      }
    else
      {
	retarr = NullNekDoubleArrayOfArray;
	//cout<<"\n returning NULL after ctr = "<<avgiterGD;
	return;
      }
  }
  
  Array<OneD, Array<OneD, NekDouble> > call_companion_rf(Array<OneD, NekDouble> uhats, Array<OneD, Array<OneD, NekDouble> > &C)
  {
    Array<OneD, NekDouble> uhatshold;
    Array<OneD, Array<OneD, NekDouble> > Cmat;
    Array<OneD, Array<OneD, NekDouble> > retarr(1);
    int N = uhats.size();
    vector<NekDouble>ret;
    while(abs(uhats[N-1])<1e-8)
      {
	N = N-1;
	if(N == 0 || N == 1)
	  {
	    ret.push_back(-1.0);
	    ret.push_back(1.0);
	    retarr[0] = Array<OneD, NekDouble>(ret.size(), ret.data());
	    return retarr;
	  }
      }
    uhatshold = Array<OneD, NekDouble>(N);;
    while(abs(uhatshold[N-1])<1e-3)
      {

	for(int k = 0; k < N; k++)
	  {
	    
	    for(int i = 0; i < N; i ++)
	      {
		// C[i][k] since transpose
		uhatshold[k] += C[i][k]*uhats[i];
	      }
	  }
	N = N-1;
	if(N == 1)
	  {
	    ret.push_back(-1.0);
	    ret.push_back(1.0);
	    
	    retarr[0] = Array<OneD, NekDouble> (ret.size(), ret.data());
	    return retarr;
	  }
      }
    Vmath::Smul(N, 1/uhatshold[N-1], uhatshold, 1, uhatshold, 1); //last ele = 1
    Cmat =  Array<OneD, Array<OneD, NekDouble> >(N-1);
    for(int i = 0; i<N-1; i++)
      {

	Array<OneD, NekDouble> row(N-1, 0.0);
		
	if(i > 0)
	  {
	    row[i-1] = 1;
	  } 
	row[N-2] = -uhatshold[i];
		
	Cmat[i] = row;
		
      }
    Array<OneD, NekDouble> EIG_R , EIG_I;
    FindEigenval(Cmat, EIG_R, EIG_I);
    for(int kk = 0; kk <EIG_R.size(); kk++)
      {
	ret.push_back( EIG_R[kk] );
      }
    ret.push_back(-1.0);
    ret.push_back(1.0);
    retarr[0] = Array<OneD, NekDouble> (ret.size(), ret.data());
    return retarr;
    
  }
  
  Array<OneD, Array<OneD, NekDouble>> GetQuadratureMidCoords(Array<OneD, Array<OneD, NekDouble>> &coords)
  {

    int dimension = coords.size();
    int totPoints = coords[0].size()-1;

    // formlattice by adding midpoints between each pair of quad points
    Array<OneD, Array<OneD, NekDouble> > quadraturemidcoords(dimension);
    for(int k = 0 ; k < dimension; k++)
      {
	quadraturemidcoords[k] = Array<OneD, NekDouble> (totPoints);

	for(int p = 0; p < quadraturemidcoords[k].size(); p++)
	  {
	    quadraturemidcoords[k][p] = (coords[k][p]+coords[k][p+1])/2;
	  }
      }

    return quadraturemidcoords;
  }

  
  Array<OneD, Array<OneD, NekDouble>> GetLatticeCoords( Array<OneD, Array<OneD, NekDouble>> &coords, Array<OneD, Array<OneD, NekDouble>> &midcoords)

  {
    Array<OneD, Array<OneD, NekDouble>> latticecoords(coords.size());

    for(int k = 0; k < coords.size(); k++)
      {
	latticecoords[k] = Array<OneD, NekDouble>(coords[k].size() + midcoords[k].size());
	Vmath::Vcopy(coords[k].size(), &coords[k][0], 1, &latticecoords[k][0], 1);
	Vmath::Vcopy(midcoords[k].size(), &midcoords[k][0], 1, &latticecoords[k][coords[k].size()], 1);
      }
  
    return latticecoords;
  }

  
  
  Array<OneD, Array<OneD, NekDouble>> GetCoords(StdExpansion *E)
  {
    int dimension = E->GetShapeDimension();
    const auto totPoints = (unsigned) E->GetTotPoints();

    Array<OneD, NekDouble> x(totPoints), y(totPoints), z(totPoints);
    Array<OneD, Array<OneD, NekDouble>> coords(dimension);

    switch (dimension)
      {
      case 1:
	{
	  E->GetCoords(x);
	  coords[0] = x;
	  break;
	}

      case 2:
	{
	  E->GetCoords(x, y);
	  coords[0] = x;
	  coords[1] = y;
	  break;
	}

      case 3:
	{
	  E->GetCoords(x, y, z);
	  coords[0] = x;
	  coords[1] = y;
	  coords[2] = z;
	  break;
	}
      default:
	break;
      }

    return coords;
  }
    
  string funcdef;

  void SetStartArr(Array<OneD, NekDouble> &st)
  {
    this->startarr = st;
  }
  void Setgamhold(NekDouble g)
  {
    gamhold = g;
  }
  
  void Setchold(NekDouble c)
  {
    chold = c;
  }

  void SettolGD(NekDouble c)
  {
    tol_gd = c;
  }
  NekDouble GettolGD()
  {
    return    tol_gd; 
  }

protected:

  NekDouble iterGD = 1e-2;
  NekDouble tol_gd = 1e-6;
  NekDouble secarg = 1e4;
  NekDouble eps = -0.1;

  int avgnum = 1;
  //defaults
  NekDouble chold = 0.2;
  NekDouble gamhold = 0.2;
  Array<OneD, NekDouble> startarr;  
  Array<OneD, NekDouble> bcheb;  
  //  int verbose = 0;
};

#endif
