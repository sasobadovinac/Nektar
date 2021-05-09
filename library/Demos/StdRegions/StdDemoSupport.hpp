//////////////////////////////////////////////////////////////////////////////
//
// File: StdDemoSupport.cpp
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

#ifndef DEMOS_STDREGIONS_STDDEMOSUPPORT_HPP
#define DEMOS_STDREGIONS_STDDEMOSUPPORT_HPP

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
namespace po = boost::program_options;

class DemoSupport
{
public:
  Array<OneD, Array<OneD, NekDouble> > testcoord3dqpts;
  Array<OneD, Array<OneD, NekDouble> > testcoord3dqmidpts;
  Array<OneD, Array<OneD, NekDouble> > testcoord2dqqpts;
  Array<OneD, Array<OneD, NekDouble> > testcoord2dtqpts;
  Array<OneD, Array<OneD, NekDouble> > testcoord2dqqmidpts;
  Array<OneD, Array<OneD, NekDouble> > testcoord2dtqmidpts;

  Array<OneD, Array<OneD, NekDouble> > testcoord3dlattice;
  Array<OneD, Array<OneD, NekDouble> > testcoord2dqlattice;
  Array<OneD, Array<OneD, NekDouble> > testcoord2dtlattice;

  Array<OneD, Array<OneD, NekDouble> > testcoord1dpts;
  Array<OneD, Array<OneD, NekDouble> > testcoord1dmidptspts;
  Array<OneD, Array<OneD, NekDouble> > testcoord1dlattice;
  
  Array<OneD, NekDouble > interioreval1dmidpts;
  Array<OneD, NekDouble > interioreval3dqmidpts;
  Array<OneD, NekDouble > interioreval2dqqmidpts;
  Array<OneD, NekDouble > interioreval2dtqmidpts;
  
  DemoSupport() : m_desc("Available options")
  {
    m_desc.add_options()
      ("help,h",
       "Produce this help message and list basis and shape types.")
      ("nodal,n",
       po::value<string>(&m_ntype),
       "Optional nodal type, autofills shape and basis choices.")
      ("shape,s",
       po::value<string>(&m_shape),
       "Region shape to project function on.")
      ("basis,b",
       po::value<vector<string>>(&m_basis)->multitoken(),
       "Basis type, separate by spaces for higher dimensions.")
      ("order,o",
       po::value<vector<int>>(&m_order)->multitoken()->required(),
       "Order of basis sets, separate by spaces for higher dimensions.")
      ("points,p",
       po::value<vector<int>>(&m_points)->multitoken()->required(),
       "Number of quadrature points, separate by spaces for "
       "higher dimensions.")
      ("iterGD,t",
       po::value<NekDouble>(&iterGD)->multitoken(),
       "Tol of gradient descent.")
      ("eps,e",
       po::value<NekDouble>(&eps)->multitoken(),
       "Negativity factor.")
      ("avgnum,a",
       po::value<int>(&avgnum)->multitoken(),
       "Average number.")
      ("secarg,i",
       po::value<NekDouble>(&secarg)->multitoken(),
       "Max iters of gradient descent.")
      ("gamhold,g",
       po::value<NekDouble>(&gamhold)->multitoken(),
       "gamma value for steepest descent with line-search.")
      ("chold,c",
       po::value<NekDouble>(&chold)->multitoken(),
       "c value for steepest descent with line-search.")

      ("pointstype,P",
       po::value<vector<string>>(&m_pointstype)->multitoken(),
       "Optional points type, separate by spaces for higher dimensions.");                      
  }

  void ParseArguments(int argc, char *argv[])
  {
    try
      {
	po::store(po::parse_command_line(argc, argv, m_desc), m_vm);
	if (m_vm.count("help"))
	  {
	    cout << m_desc;
	    cout << endl << "All nodal types, -n [ --nodal ], are:" << endl;
	    for (int i = 22; i < SIZE_PointsType; ++i)
	      {
		cout << kPointsTypeStr[i] << endl;
	      };
	    cout << endl << "All shape types, -s [ --shape ], are:" << endl;
	    for (int i = 1; i < SIZE_ShapeType; ++i)
	      {
		cout << ShapeTypeMap[i] << endl;
	      };
	    cout << endl << "All basis types, -b [ --basis ], are:" << endl;
	    for (int i = 1; i < SIZE_BasisType; ++i)
	      {
		cout << BasisTypeMap[i] << endl;
	      };
	    cout << endl << "All points types, -P [ --pointstype ], are:"
		 << endl;
	    for (int i = 1; i < SIZE_PointsType; ++i)
	      {
		cout << kPointsTypeStr[i] << endl;
	      };
	    exit(0);
	  }
	po::notify(m_vm);
      }
    catch (const exception &e)
      {
	cerr << "Error: " << e.what() << endl << m_desc;
	exit(1);
      }
  }

  Array<OneD, Array<OneD, NekDouble> > formConf(NekDouble N)
  {
    int evalPts = N+1; //+1 for luck
    int i, j;
    Nektar::Array<Nektar::OneD, Nektar::NekDouble> ab1(evalPts);
    Nektar::Array<Nektar::OneD, Nektar::NekDouble> ab2(evalPts);
    Polylib::RecCoeff(evalPts, &ab1[0], &ab2[0], -0.5, -0.5);
    // Form confederate matrix
    // a = 2*a
    // b = 2*b
    Vmath::Smul(ab1.size(), 2.0, ab1, 1, ab1, 1);
    Vmath::Smul(ab2.size(), 2.0, ab2, 1, ab2, 1);
    
    //  J = full(spdiags([[b(3:n);0.5;0] a(1:n) b(1:n)], -1:1, n, n));
    vector<vector<NekDouble> > J;//, Array<OneD, NekDouble>(N,0.0));

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
  
  StdExpansion *CreateStdExpansion()
  {
    vector<PointsType> ptype;
    if (m_vm.count("pointstype"))
      {
	for (auto &p : m_pointstype)
	  {
	    PointsType tmp = eNoPointsType;

	    // starts at nodal points
	    for (int i = 1; i < SIZE_PointsType; ++i)
	      {
		if (boost::iequals(kPointsTypeStr[i], p))
		  {
		    tmp = static_cast<PointsType>(i);
		    break;
		  }
		ASSERTL0(i != SIZE_PointsType - 1,
			 "The points type '" + p + "' does not exist");
	      }

	    ptype.push_back(tmp);
	  }
      }

    // Convert string input argument to nodal type
    PointsType nodaltype = eNoPointsType;
    ShapeType stype = eNoShapeType;
    vector<BasisType> btype(3, eNoBasisType);
    if (m_vm.count("nodal"))
      {
	for (int i = 22; i < SIZE_PointsType; ++i) // starts at nodal points
	  {
	    if (boost::iequals(kPointsTypeStr[i], m_ntype))
	      {
		nodaltype = static_cast<PointsType>(i);
		break;
	      }
	    ASSERTL0(i != SIZE_PointsType - 1,
		     "The nodal type '" + m_ntype + "' does not exist");
	  }
	switch (nodaltype)
	  {
	  case eNodalTriElec:
	  case eNodalTriFekete:
	  case eNodalTriSPI:
	  case eNodalTriEvenlySpaced:
	    btype[0] = eOrtho_A;
	    btype[1] = eOrtho_B;
	    stype = eTriangle;
	    break;
	  case eNodalQuadElec:
	    btype[0] = eOrtho_A;
	    btype[1] = eOrtho_B;
	    stype = eQuadrilateral;
	    break;
	  case eNodalTetElec:
	  case eNodalTetSPI:
	  case eNodalTetEvenlySpaced:
	    btype[0] = eOrtho_A;
	    btype[1] = eOrtho_B;
	    btype[2] = eOrtho_C;
	    stype = eTetrahedron;
	    break;
	  case eNodalPrismElec:
	  case eNodalPrismSPI:
	  case eNodalPrismEvenlySpaced:
	    btype[0] = eOrtho_A;
	    btype[1] = eOrtho_A;
	    btype[2] = eOrtho_B;
	    stype = ePrism;
	    break;
	  case eNodalHexElec:
	    btype[0] = eOrtho_A;
	    btype[1] = eOrtho_A;
	    btype[2] = eOrtho_A;
	    stype = eHexahedron;
	    break;
	  default:
	    ASSERTL0(!nodaltype, ("The nodal type '" + m_ntype +
				  "' is invalid for StdProject."));
	    break;
	  }
      }
    //Convert string input argument to shape type
    if (stype == eNoShapeType)
      {
	for (int i = 1; i < SIZE_ShapeType; ++i)
	  {
	    if (boost::iequals(ShapeTypeMap[i], m_shape))
	      {
		stype = static_cast<ShapeType>(i);
		break;
	      }
	    ASSERTL0(i != SIZE_ShapeType - 1,
		     "The shape type '" + m_shape + "' does not exist");
	  }
      }

    // Check arguments supplied equals dimension
    const int dimension = (stype == ePoint) ? 1 : ShapeTypeDimMap[stype];
    ASSERTL0(m_order.size() == dimension,
	     "Number of orders supplied should match shape dimension");
    ASSERTL0(m_points.size() == dimension,
	     "Number of points supplied should match shape dimension");
    ASSERTL0(ptype.size() == dimension || ptype.size() == 0,
	     "Number of points types should match shape dimension if "
	     "supplied.");
    if (!m_vm.count("nodal"))
      {
	ASSERTL0(m_basis.size() == dimension,
		 "Number of bases supplied should match shape dimension");

	// Convert string input argument to basis types
	for (int i = 0; i < dimension; ++i)
	  {
	    for (int j = 1; j < SIZE_BasisType; ++j)
	      {
		if (boost::iequals(BasisTypeMap[j], m_basis[i]))
		  {
		    btype[i] = static_cast<BasisType>(j);
		    break;
		  }
		ASSERTL0(j != SIZE_BasisType - 1,
			 "Basis type '" + m_basis[i] + "' does not exist");
	      }
	  }
      }

    //check basis selection is permitted for chosen shape
    map<ShapeType, vector<vector<BasisType>>> allowableBasis;
    allowableBasis[ePoint] = {
      {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
       eLegendre, eChebyshev, eMonomial, eFourierSingleMode,
       eFourierHalfModeRe, eFourierHalfModeIm}
    };
    allowableBasis[eSegment] = { allowableBasis[ePoint][0] };
    allowableBasis[eTriangle] = {
      {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
      {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange}
    };
    allowableBasis[eQuadrilateral] = {
      allowableBasis[eSegment][0], allowableBasis[eSegment][0]
    };
    allowableBasis[eTetrahedron] = {
      {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
      {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange},
      {eOrtho_C, eModified_C, eGLL_Lagrange, eGauss_Lagrange}
    };
    allowableBasis[ePyramid] = {
      {eOrtho_A,    eModified_A,    eGLL_Lagrange, eGauss_Lagrange},
      {eOrtho_A,    eModified_A,    eGLL_Lagrange, eGauss_Lagrange},
      {eOrthoPyr_C, eModifiedPyr_C, eGLL_Lagrange, eGauss_Lagrange}
    };
    allowableBasis[ePrism] = {
      {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
      {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
      {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange}
    };
    allowableBasis[eHexahedron] = {
      {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
       eLegendre, eChebyshev, eMonomial},
      {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
       eLegendre, eChebyshev, eMonomial},
      {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
       eLegendre, eChebyshev, eMonomial}
    };

    for (int i = 0; i < dimension; ++i)
      {
	const unsigned int basisListLength = allowableBasis[stype][i].size();
	for (int j = 0; j < basisListLength; ++j)
	  {
	    if (allowableBasis[stype][i][j] == btype[i])
	      {
		break;
	      }
	    ASSERTL0(j != basisListLength - 1,
		     ("The basis type '" +
		      static_cast<string>(BasisTypeMap[btype[i]]) +
		      "' is invalid for basis argument " + to_string(i + 1) +
		      " for shape '" + ShapeTypeMap[stype] + "'."))
	      }
      }

    //Declaration of other variables needed
    StdExpansion *E = nullptr;

    // Assign points type according to basis type selection, if not already
    // assigned.
    if (ptype.size() == 0)
      {
	ptype.resize(dimension);
	for (int i = 0; i < dimension; ++i)
	  {
	    if (btype[i] == eFourier)
	      {
		ptype[i] = eFourierEvenlySpaced;
	      }
	    else if (btype[i] == eFourierSingleMode ||
		     btype[i] == eFourierHalfModeRe ||
		     btype[i] == eFourierHalfModeIm)
	      {
		ptype[i] = eFourierSingleModeSpaced;
	      }
	    else
	      {
		if (i == 1 && (stype == eTriangle || stype == eTetrahedron))
		  {
		    ptype[i] = eGaussRadauMAlpha1Beta0;
		  }
		else if (i == 2 && (stype == eTetrahedron || stype == ePyramid))
		  {
		    ptype[i] = eGaussRadauMAlpha2Beta0;
		  }
		else if (i == 2 && stype == ePrism)
		  {
		    ptype[i] = eGaussRadauMAlpha1Beta0;
		  }
		else
		  {
		    ptype[i] = eGaussLobattoLegendre;
		  }
	      }
	  }
      }

    vector<PointsKey> pkey;
    vector<BasisKey> bkey;
    for (int i = 0; i < dimension; ++i)
      {
	pkey.emplace_back(PointsKey(m_points[i], ptype[i]));
	bkey.emplace_back(BasisKey(btype[i], m_order[i], pkey[i]));
      }

    switch (stype)
      {
      case ePoint:
	{
	  E = new StdPointExp(bkey[0]);
	  break;
	}
      case eSegment:
	{
	  E = new StdSegExp(bkey[0]);
	  break;
	}
      case eTriangle:
	{
	  E = nodaltype != eNoPointsType ? new StdNodalTriExp(bkey[0],
							      bkey[1],
							      nodaltype)
	    : new StdTriExp(bkey[0], bkey[1]);
	  break;
	}
      case eQuadrilateral:
	{
	  E = new StdQuadExp(bkey[0], bkey[1]);
	  break;
	}
      case eTetrahedron:
	{
	  E = nodaltype != eNoPointsType ? new StdNodalTetExp(bkey[0],
							      bkey[1],
							      bkey[2],
							      nodaltype)
	    : new StdTetExp(bkey[0], bkey[1],
			    bkey[2]);
	  break;
	}
      case ePyramid:
	{
	  E = new StdPyrExp(bkey[0], bkey[1], bkey[2]);
	  break;
	}
      case ePrism:
	{
	  E = nodaltype != eNoPointsType ? new StdNodalPrismExp(bkey[0],
								bkey[1],
								bkey[2],
								nodaltype)
	    : new StdPrismExp(bkey[0], bkey[1],
			      bkey[2]);
	  break;
	}
      case eHexahedron:
	{
	  E = new StdHexExp(bkey[0], bkey[1], bkey[2]);
	  break;
	}
      default:
	break;
      }

    m_bkey = bkey;
    m_pkey = pkey;
    C = formConf((dimension)*pow(3*m_order[0]+1,2));
    stypeglo = stype;
    return E;
  }

  void OrthoNormalize(Array<OneD, NekDouble>  &coeffs,  LibUtilities::BasisKey orthoBasisKey0, LibUtilities::BasisKey orthoBasisKey1, int flag)
  {
    Array<OneD,NekDouble> temp1(coeffs.size());

    if(flag == 0)
      {
	LibUtilities::InterpCoeff2D( (m_bkey[0]), (m_bkey[1]), coeffs,  orthoBasisKey0, orthoBasisKey1, temp1);
	Vmath::Vcopy(coeffs.size(), temp1, 1, coeffs, 1);
      }
    if(flag == 1)
      {
	LibUtilities::InterpCoeff2D( orthoBasisKey0, orthoBasisKey1, coeffs,  (m_bkey[0]), (m_bkey[1]), temp1);
	Vmath::Vcopy(coeffs.size(), temp1, 1, coeffs, 1);

      }

  }


  void FindEigenval(Array<OneD, Array<OneD, NekDouble> > CM,
		    Array<OneD,NekDouble> &EIG_R,
		    Array<OneD,NekDouble> &EIG_I)
  {
    const int sz = CM.size();
    Array<OneD, NekDouble> Mat(sz*sz,0.0);
    int info = 0;

    const int  lwork = 3*(sz);
    Array<OneD, NekDouble> work1(lwork);
    for(int i = 0; i<sz; i++)

      {
	Vmath::Vcopy(sz, &CM[i][0], 1, &Mat[0]+i, sz);
      }

    char jobvl = 'N';
    char jobvr = 'N';

    NekDouble dum1, dum2;

    Lapack::Dgeev(jobvl,jobvr,sz,Mat.get(),sz,EIG_R.get(),EIG_I.get(),&dum1,1,&dum2,1,&work1[0],lwork,info);
    
    vector<NekDouble> tmp;
    for(int k = 0; k < EIG_R.size(); k++)
      {
	//	cout<<"\n i = "<<EIG_I[k]<<" r = "<<EIG_R[k];
	if(abs(EIG_I[k]) < 1e-5 && abs(EIG_R[k]) < 1)
	  tmp.push_back(EIG_R[k]);
      }
    EIG_R = Array<OneD, NekDouble>(tmp.size());
    for(int k = 0; k < tmp.size(); k++)
      {
	EIG_R[k] = tmp[k];
      }

  }
  
  void pq(
	  Array<OneD,NekDouble> uhats,
	  Array<OneD, Array<OneD,  NekDouble> > roots,
	  Array<OneD, NekDouble> V1,
	  Array<OneD,NekDouble> &pqevalxast,
	  Array<OneD,NekDouble> &fvals
	  )
  {
    int N = uhats.size();
    Array<OneD,NekDouble> w1(N);
    Array<OneD, NekDouble> Vsumsq(roots[0].size());
    fvals = blasmatvec(V1, uhats, roots[0].size(), uhats.size());

    for( int i = 0; i < roots[0].size(); i++)
      {
	Vmath::Vmul(N, &V1[i], roots[0].size(), &V1[i], roots[0].size(), &w1[0], 1);

	Vsumsq[i] = (pow(Vmath::Vsum(N, &w1[0], 1),-0.5));
      }
    Vmath::Vmul(roots[0].size(), &Vsumsq[0], 1, &fvals[0], 1, &fvals[0], 1);

    if(pqevalxast.size() > 0)
      {
	int minidx = Vmath::Imin(fvals.size(), &fvals[0], 1);
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
    for(int k = 0; k < GetAvgNum(); k++)
      {
	double alpha= 1.0, beta= 0.0;
	char no= 'N';
	int m = nq1, n = nq2, lda = m, incx = 1, incy = 1;
	double *A = &M[0];
	//t1.Start();
	Blas::dgemv_(no,m,n,alpha,A,lda,&vec[0],incx,beta,&hold[0],incy);
	//t1.Stop();
      }
    return hold;
  }
  
  
  // With backtracking new
  void steepestgradient_descent2Dquad(Array<OneD, NekDouble> &uhats, StdExpansion *E, Array<OneD, Array<OneD, NekDouble> > &storage, Array<OneD, Array<OneD, NekDouble> >&retarr, NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> > &testcoord2d, Array<OneD, Array<OneD, NekDouble> > &testcoord2dqmid, Array<OneD, NekDouble> &interioreval2dqmid, Array<OneD, Array<OneD, NekDouble> > &testcoord2dlattice)
  {
    int dim = 2;
    Array<OneD, NekDouble> g(dim);
    double inf = numeric_limits<double>::infinity();
    retarr = Array<OneD, Array<OneD, NekDouble> >(2);
    retarr[0] = Array<OneD,NekDouble>(1);
    retarr[1] = Array<OneD,NekDouble>(1);

    // flag that determines if GD is called from GDBackTrackingTest.cpp
    int flagTest = 1;
    
    Array<OneD,NekDouble> xnew(dim) ,x(dim);
    NekDouble gprev = inf;
    Array<OneD, NekDouble> nullarr(0), temp2(testcoord2d[0].size()), temp3, alltemp;
    NekDouble gnew = inf;

    int idxgprev;
        
    int sz = testcoord2d[0].size();
    Array<OneD, Array<OneD, NekDouble > > tempeval(4);

    for(int i = 0; i < 4; i++)
      {
        tempeval[i] = Array<OneD, NekDouble>(E->GetNcoeffs());
      }

    pq(uhats, testcoord2d, storage[0], nullarr, temp2 );

    gprev = Vmath::Vmin(sz, temp2, 1);
    idxgprev = Vmath::Imin(temp2.size(), temp2, 1);
    pq(uhats, testcoord2dqmid, interioreval2dqmid, nullarr, temp3 );

    NekDouble gprevmid =  Vmath::Vmin(temp3.size(), temp3, 1);

    if(gprev > gprevmid)
      {
        gprev = gprevmid;
        idxgprev = sz+Vmath::Imin(temp3.size(), temp3, 1);

      }
    if(startarr.size() > 1)
      {
	cout<<"\n expected root roughly = "<< startarr[0]<<" , "<<startarr[1]<<"\n";
	flagTest = 0;
      }
    xnew[0] = testcoord2dlattice[0][idxgprev];
    xnew[1] = testcoord2dlattice[1][idxgprev];

    Array<OneD, NekDouble> xstart(2);
    xstart[0] = xnew[0];
    xstart[1] = xnew[1];
    Array<OneD, NekDouble> dereval(dim);
    //    cout<<"    starting pt = "<<xnew[0]<<","<<xnew[1]<<" N "<< uhats.size();
    // fstream fio;
    // fio.open("dumquad.txt", ios::app | ios::out | ios::in);
    // fio<<"\n func = "<<funcdef<<"; starting pt = "<<xnew[0]<<","<<xnew[1]<<" N = 5";

    int call_GD;
    if(flagTest == 0)
      call_GD = 1;
    else
      call_GD = gprev < 0 && abs(gprev)>1e-13; 
    if(call_GD)
      {

	Array<OneD, Array<OneD, NekDouble > > tempeval(4);
	for(int k = 0; k < tempeval.size(); k++)
	  {
	    tempeval[k] = Array<OneD, NekDouble>(uhats.size());
	  }
	Array<OneD, Array<OneD, NekDouble > > xastaa(2);
	xastaa[0] = Array<OneD, NekDouble> (1, xnew[0]);
	xastaa[1] = Array<OneD, NekDouble> (1, xnew[1]);

	NekDouble c = chold;
	tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2], tempeval[3]);
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
	pq(uhats, xastaa, tempeval[0], nullarr, t1);
	Array<OneD, NekDouble > holdxandval(3), saveholdxandval(3), savesavehold(3) ;//dim+1
	holdxandval[0] = xastaa[0][0];
	holdxandval[1] = xastaa[1][0];
	holdxandval[2] = t1[0];
	savesavehold[2] = inf;
	int truth_val = abs(g[0]+g[1]) > 1e-9;
	int ctr = 0, ct  = 0;
	NekDouble fac = 1,  gnew0 = gnew, gnew1 = gnew, gnew2 = gnew;

	if(truth_val)
	  {
	    NekDouble iter = secarg;
	    int counter_min_change = 0;

	    while(ctr < iter && ((abs(g[0]) > 1e-9) || abs(g[1]) > 1e-9) && fac > 1e-6 && counter_min_change < 5 )
	      {
		if(ctr > 2 && abs(gnew0 - gnew2) < 1e-8)
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


		while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1)
		  {
		    xastaa[0][0] = holdxandval[0]  - fac*g[0];
		    xastaa[1][0] = holdxandval[1]  - fac*g[1];
		    fac = fac*gamhold;
		  }

		ctr++;
		tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2], tempeval[3]);
		pq(uhats, xastaa, tempeval[0], nullarr, t1);
		for(int p = 0; p < dim; p++)
		  {
		    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
		    g[p] = dereval[p];
		  }
		
		if(savesavehold[dim] > t1[0])
		  {
		    for(int j = 0; j < dim; j++)
		      savesavehold[j] = xastaa[j][0];
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

		NekDouble holdval = t1[0], saveholdval = t1[0];
		Array<OneD, NekDouble> gsave(2);
		gsave[0] = g[0];
		gsave[1] = g[1];
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

		    while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1)
		      {
			xastaa[0][0] = holdxandval[0]  - fac*g[0];
			xastaa[1][0] = holdxandval[1]  - fac*g[1];
			fac = fac*gamhold;
		      }
		    tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2]
						       , tempeval[3]);
		    for(int p = 0; p < dim; p++)
		      {
			derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
			g[p] = dereval[p];
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
		    pq(uhats, xastaa, tempeval[0], nullarr, t1);
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
		    if((saveholdxandval[2]<savesavehold[2]) && abs(saveholdxandval[2] - savesavehold[2]) > 1e-6)
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
		    //fio<<" "<<savesavehold[2];

		    xastaa[0][0] = saveholdxandval[0];
		    xastaa[1][0] = saveholdxandval[1];

		    t1[0] = saveholdval;

		  }
		else
		  {
 		    
		    
		    if((saveholdxandval[2]<savesavehold[2]) && abs(saveholdxandval[2] - savesavehold[2]) > 1e-6)
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

		  }
	      }
	  }
	retarr = Array<OneD, Array<OneD, NekDouble> >(2);
	retarr[0] = Array<OneD, NekDouble>(1, savesavehold[0]);
	retarr[1] = Array<OneD, NekDouble>(1, savesavehold[1]);
	avgiterGD = ctr;

	return;
      }

  }

  void steepestgradient_descent2Dtri(Array<OneD, NekDouble> &uhats, StdExpansion *E, Array<OneD,Array<OneD, NekDouble> > &storage, Array<OneD, Array<OneD, NekDouble> >&retarr, NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> > &testcoord2d, Array<OneD, Array<OneD, NekDouble> > &testcoord2dqmid, Array<OneD, NekDouble> &interioreval2dqmid, Array<OneD, Array<OneD, NekDouble> > &testcoord2dlattice)
  {
    int dim = 2;
    Array<OneD, NekDouble> g(dim);
    double inf = numeric_limits<double>::infinity();
    retarr = Array<OneD, Array<OneD, NekDouble> >(2);
    retarr[0] = Array<OneD,NekDouble>(1);
    retarr[1] = Array<OneD,NekDouble>(1);

    // flag that determines if GD is called from GDBackTrackingTest.cpp
    int flagTest = 1;
    Array<OneD,NekDouble> xnew(dim) ,x(dim);
    NekDouble gprev = inf;
    Array<OneD, NekDouble> nullarr(0), temp2(testcoord2d[0].size()), temp3, alltemp;

    NekDouble gnew = inf;
    int idxgprev;

    int sz = testcoord2d[0].size();
    Array<OneD, Array<OneD, NekDouble > > tempeval(4);

    for(int i = 0; i < 4; i++)
      {
	tempeval[i] = Array<OneD, NekDouble>(E->GetNcoeffs());
      }

    pq(uhats, testcoord2d, storage[0], nullarr, temp2 );

    gprev = Vmath::Vmin(sz, temp2, 1);
    idxgprev = Vmath::Imin(temp2.size(), temp2, 1);
    pq(uhats, testcoord2dqmid, interioreval2dqmid, nullarr, temp3 );

    NekDouble gprevmid =  Vmath::Vmin(temp3.size(), temp3, 1);

    if(gprev > gprevmid)
      {
	gprev = gprevmid;
	idxgprev = sz+Vmath::Imin(temp3.size(), temp3, 1);

      }
    xnew[0] = testcoord2dlattice[0][idxgprev];
    xnew[1] = testcoord2dlattice[1][idxgprev];

    //expected min (rougly, using dense grid)  
    if(startarr.size() > 1)
      {
	cout<<"\n expected root roughly = "<< startarr[0]<<" , "<<startarr[1]<<"\n";
	flagTest = 0;
      }
    Array<OneD, NekDouble> xstart(2);
    xstart[0] = xnew[0];
    xstart[1] = xnew[1];
    Array<OneD, NekDouble> dereval(dim);

    int call_GD;
    if(flagTest == 0)
      call_GD = 1;
    else
      call_GD = gprev < 0 && abs(gprev)>1e-13;
	
    if(call_GD)
      {

	Array<OneD, Array<OneD, NekDouble > > tempeval(4);
	for(int k = 0; k < tempeval.size(); k++)
	  {
	    tempeval[k] = Array<OneD, NekDouble>(uhats.size());
	  }
	Array<OneD, Array<OneD, NekDouble > > xastaa(2);
	xastaa[0] = Array<OneD, NekDouble> (1, xnew[0]);
	xastaa[1] = Array<OneD, NekDouble> (1, xnew[1]);


	NekDouble c = chold;
	tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2], tempeval[3]);
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
	pq(uhats, xastaa, tempeval[0], nullarr, t1);
	Array<OneD, NekDouble > holdxandval(3), saveholdxandval(3), savesavehold(3) ;

	holdxandval[0] = xastaa[0][0];
	holdxandval[1] = xastaa[1][0];
	holdxandval[2] = t1[0];
	savesavehold[2] = inf;
	int truth_val = abs(g[0]+g[1]) > 1e-9;
	int ctr = 0, ct  = 0;
	NekDouble fac = 1, gnew0 = gnew, gnew1 = gnew, gnew2 = gnew;
	if(truth_val)
	  {
	    NekDouble iter = secarg;
	    int counter_min_change = 0;
	    while(ctr < iter && ((abs(g[0]) > 1e-9) || abs(g[1]) > 1e-9) && fac > 1e-6 && counter_min_change <5)
	      {
		if(ctr > 2 && abs(gnew0 - gnew2) < 1e-8)
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

		while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || (xastaa[0][0] + xastaa[1][0])>1)
		  {
		    xastaa[0][0] = holdxandval[0]  - fac*g[0];
		    xastaa[1][0] = holdxandval[1]  - fac*g[1];
		    fac = fac*gamhold;
		  }

		ctr++;
		tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2], tempeval[3]);
		pq(uhats, xastaa, tempeval[0], nullarr, t1);
		if(savesavehold[dim] > t1[0])
		  {
		    for(int j = 0; j < dim; j++)
		      savesavehold[j] = xastaa[j][0];
		    savesavehold[dim] = t1[0];
		  }
		
		for(int p = 0; p < dim; p++)
		  {
		    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
		    g[p] = dereval[p];
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

		NekDouble holdval = t1[0], saveholdval = t1[0];
		Array<OneD, NekDouble> gsave(2);
		gsave[0] = g[0];
		gsave[1] = g[1];
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

		    while((abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || (xastaa[0][0] + xastaa[1][0])>1 ) && fac > 1e-7)
		      {
			xastaa[0][0] = holdxandval[0]  - fac*g[0];
			xastaa[1][0] = holdxandval[1]  - fac*g[1];
			fac = fac*gamhold;
		      }
		    tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2]
						       , tempeval[3]);

		    for(int p = 0; p < dim; p++)
		      {
			derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
			g[p] = dereval[p];
		      }
		    gnew2 = gnew1;
		    gnew1 = gnew0;
		    gnew = pow((pow(g[0],2) + pow(g[1],2)) ,0.5);
		    if(abs(g[0]) > 1 || abs(g[1])>1)
		      {
			g[0] = g[0]/gnew;
			g[1] = g[1]/gnew;
		      }
		    gnew0 = gnew;
		    pq(uhats, xastaa, tempeval[0], nullarr, t1);
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
		//fio<<"\n ctr = "<< ctr<<" "<<xastaa[0][0]<<" "<<xastaa[1][0]<<" "<<t1[0]<<" ct = "<<ct;
		avgiterGD = ctr;
		if(saveholdval < holdval )
		  {
		    if((saveholdxandval[2]<savesavehold[2]) && abs(saveholdxandval[2] - savesavehold[2]) > 1e-6)
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
		    //fio<<" "<<savesavehold[2];
		    xastaa[0][0] = saveholdxandval[0];
		    xastaa[1][0] = saveholdxandval[1];
		    
		    t1[0] = saveholdval;
		    
		  }
		else
		  {

		    
		    if((saveholdxandval[2]<savesavehold[2]) && abs(saveholdxandval[2] - savesavehold[2]) > 1e-6)
		      {
			savesavehold[0] = saveholdxandval[0];
			savesavehold[1] = saveholdxandval[1];
			savesavehold[2] = saveholdxandval[2];
			//reset counter_min_change
			counter_min_change = 0;
		      }
		    else
		      {
			counter_min_change++;
		      }
		    //fio<<" "<<savesavehold[2];
		    t1[0] = holdval;
		    if(abs(xastaa[0][0])>1 || abs(xastaa[1][0]) > 1 || (xastaa[0][0] + xastaa[1][0])>1)
		      {
			xastaa[0][0] = savesavehold[0];
			xastaa[1][0] = savesavehold[1];
		      }

		    fac = 1;
		  }

	      }
	  }
	retarr = Array<OneD, Array<OneD, NekDouble> >(2);
	retarr[0] = Array<OneD, NekDouble>(1, savesavehold[0]);
	retarr[1] = Array<OneD, NekDouble>(1, savesavehold[1]);
	avgiterGD = ctr;

	return;
      }

  }

  //GD with backtracking
  void steepestgradientdescent3D(Array<OneD, NekDouble> &uhats, StdExpansion *E, Array<OneD, Array<OneD, NekDouble> > &storage, Array<OneD, Array<OneD, NekDouble> >&retarr, NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> > &testcoord3d, Array<OneD, Array<OneD, NekDouble> > &testcoord3dqmid, Array<OneD, NekDouble> &interioreval3dqmid, Array<OneD, Array<OneD, NekDouble> > &testcoord3dlattice)
  {
    int dim = 3;
    Array<OneD, NekDouble> g(dim);
    double inf = numeric_limits<double>::infinity();
    retarr = Array<OneD, Array<OneD, NekDouble> >(dim);
    for(int p = 0; p < dim; p++)
      {
	retarr[p] = Array<OneD,NekDouble>(1);
      }

    // flag that determines if GD is called from GDBackTrackingTest.cpp
    int flagTest = 1;
    
    Array<OneD,NekDouble> xnew(dim) ,x(dim);
    NekDouble gprev = inf;
    Array<OneD, NekDouble> nullarr(0), temp2(testcoord3d[0].size()), temp3, alltemp;
    NekDouble gnew = inf;

    int idxgprev;
    
    int sz = testcoord3d[0].size();

    pq(uhats, testcoord3d, storage[0], nullarr, temp2 );
    gprev = Vmath::Vmin(sz, temp2, 1);
    
    idxgprev = Vmath::Imin(temp2.size(), temp2, 1);
    pq(uhats, testcoord3dqmid, interioreval3dqmid, nullarr, temp3 );

    NekDouble gprevmid =  Vmath::Vmin(temp3.size(), temp3, 1);

    if(gprev > gprevmid)
      {
        gprev = gprevmid;
        idxgprev = sz+Vmath::Imin(temp3.size(), temp3, 1);

      }
    for(int p = 0; p < dim; p++)
      {
	xnew[p] = testcoord3dlattice[p][idxgprev];
      }

    if(startarr.size() > 0)
      {
	cout<<"\nexpected root roughly = "<< startarr[0]<<" , "<<startarr[1]<<", "<<startarr[2]<<" val = "<<startarr[3];
	flagTest = 0;
      }
    Array<OneD, NekDouble> xstart(dim);
    for(int p = 0; p < dim; p++)
      {
	xstart[p] = xnew[p];
      }
    Array<OneD, NekDouble> dereval(dim);

    // fstream fio;
    // fio.open("dumquad.txt", ios::app | ios::out | ios::in);

    int call_GD;
    if(flagTest == 0)
      call_GD = 1;
    else
      call_GD = gprev < 0 && abs(gprev)>1e-13;

    if(call_GD)
      {
        Array<OneD, Array<OneD, NekDouble > > tempeval(4);
	for(int k = 0; k < tempeval.size(); k++)
	  {
	    tempeval[k] = Array<OneD, NekDouble>(uhats.size());
	  }
	Array<OneD, Array<OneD, NekDouble > > xastaa(dim);
	for(int p = 0; p < dim; p++)
	  {
	    xastaa[p] = Array<OneD, NekDouble> (1, xnew[p]);
	  }

	NekDouble c = chold;
	tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2], tempeval[3]);
	for(int p = 0; p < dim; p++)
	  {
	    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
	    g[p] = dereval[p];
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
	pq(uhats, xastaa, tempeval[0], nullarr, t1);
	Array<OneD, NekDouble > holdxandval(dim+1), saveholdxandval(dim+1), savesavehold(dim+1);

	for(int p = 0; p < dim; p++)
	  {
	    holdxandval[p] = xastaa[p][0];
	  }
	holdxandval[dim] = t1[0];
	savesavehold[dim] = inf;
	NekDouble allg;
	for(int p = 0; p < dim; p++)
	  {
	    allg = allg + g[p];
	  }
	int truth_val = abs(allg) > 1e-8;
	int ctr = 0, ct  = 0;
	NekDouble fac = 1,  gnew0 = gnew, gnew1 = gnew, gnew2 = gnew;

	if(truth_val)
	  {
	    NekDouble iter = secarg;
	    int counter_min_change = 0;
	    while(ctr < iter  && fac > 1e-6 && counter_min_change < 5 && ((abs(g[0]) > 1e-9) || abs(g[1]) > 1e-9 ||  abs(g[2]) > 1e-9) )
	      {
		if(ctr > 2 && abs(gnew0 - gnew2) < 1e-6)
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
		      }
		    break;
		  case LibUtilities::eTetrahedron:
		    while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[2][0] + xastaa[1][0] > 1e-9 || xastaa[0][0] + xastaa[2][0] > 1e-9 || xastaa[1][0] + xastaa[2][0] > 1e-9)
		      {
			for(int p = 0; p < dim; p++)
			  {
			    xastaa[p][0] = holdxandval[p]  - fac*g[p];
			  }
			fac = fac*gamhold;
		      }
		    break;
		  case LibUtilities::ePyramid:
		    while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[0][0] + xastaa[1][0] > 1e-9 ||  xastaa[0][0] + xastaa[2][0] > 1e-9 )
		      {
                        for(int p = 0; p < dim; p++)
                          {
                            xastaa[p][0] = holdxandval[p]  - fac*g[p];
                          }
                        fac = fac*gamhold;
                      }
		    break;
		  default: cout<<"\n invalid element, not yet supported!";
		    exit(0);
		  }

		ctr++;
		tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2], tempeval[3]);
		pq(uhats, xastaa, tempeval[0], nullarr, t1);
		if(savesavehold[dim] > t1[0])
		  {
		    for(int j = 0; j < dim; j++)
		      savesavehold[j] = xastaa[j][0];
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
		NekDouble holdval = t1[0], saveholdval = t1[0];
		Array<OneD, NekDouble> gsave(dim);

		for(int p = 0; p < dim; p++)
		  {
		    gsave[p] = g[p];
		  }

		ct = 0;
		while(holdval > saveholdxandval[dim] - c*fac*(gsave[0] + gsave[1] + gsave[2]) && ct < iter && fac > 1e-6)
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
			while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1)
			  {
			    for(int p = 0; p < dim; p++)
			      {
				xastaa[p][0] = holdxandval[p]  - fac*g[p];
				
			      }
			    fac = fac*gamhold;
			  }
			break;
		      case LibUtilities::eTetrahedron:
			while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[0][0] + xastaa[1][0] > 1e-9 || xastaa[0][0] + xastaa[2][0] > 1e-9 || xastaa[1][0] + xastaa[2][0] > 1e-9) 
			  {
			    for(int p = 0; p < dim; p++)
			      {
				xastaa[p][0] = holdxandval[p]  - fac*g[p];

			      }
			    fac = fac*gamhold;
			  }
			break;
			case LibUtilities::ePyramid:
			  while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[0][0] + xastaa[1][0] > 1e-9 ||  xastaa[0][0] + xastaa[2][0] > 1e-9 )
			    {
			      for(int p = 0; p < dim; p++)
				{
				  xastaa[p][0] = holdxandval[p]  - fac*g[p];
				  
				}
			      fac = fac*gamhold;
			    }
			  break;
		      default:
			cout<<"\n invalid element!\n";
			exit(0);
		      }
		
		    
		    tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2]
						       , tempeval[3]);

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
		    
		    gnew = pow(gnew ,0.5);
		    gnew0 = gnew;

		    
		    if(abs(g[2])>1 || abs(g[0]) > 1 || abs(g[1])>1)
		      {
			for(int p = 0; p < dim; p++)
			  {
			    g[p] = g[p]/gnew;
			  }
		      }

		        
		    pq(uhats, xastaa, tempeval[0], nullarr, t1);
		    
		    if(savesavehold[dim] > t1[0])
		      {
			for(int j = 0; j < dim; j++)
			  savesavehold[j] = xastaa[j][0];
			savesavehold[dim] = t1[0];
		      }
		    
		    
	            if(abs(holdval -t1[0]) > 1e-6)
		      {
			holdval = t1[0];
		      }
		    else
		      break;
		  }
		//cout<<"\n ctr = "<< ctr<<" ct = "<<ct;
		avgiterGD = ctr;

		if(saveholdval < holdval)
		  {
		    if((saveholdxandval[dim]<savesavehold[dim]) && abs(saveholdxandval[dim] - savesavehold[dim]) > 1e-5)
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
		    for(int p = 0; p < dim; p++)
		      { 
			xastaa[p][0] = saveholdxandval[p];
		      }
		    
		    t1[0] = saveholdval;
		    
		  }
		else
		  {
		    
		    
		    if((saveholdxandval[dim]<savesavehold[dim]) && abs(saveholdxandval[dim] - savesavehold[dim]) > 1e-5)
		      {
			for(int p = 0; p < dim; p++)
			  { 
			    savesavehold[p] = saveholdxandval[p];
			  }
			savesavehold[dim] = saveholdxandval[dim];
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
			opt = abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1|| xastaa[0][0] + xastaa[1][0] > 1e-9 || xastaa[0][0] + xastaa[2][0] > 1e-9 || xastaa[1][0] + xastaa[2][0] > 1e-9;
			break;

		      case LibUtilities::ePyramid:
			opt = abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1|| xastaa[2][0] + xastaa[1][0] > 1e-9 || xastaa[0][0] + xastaa[2][0] > 1e-9;
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
		    
		  }
	      }
	    
	  }
	for(int p = 0; p < dim; p++)
	  {
	    retarr[p] = Array<OneD, NekDouble>(1, savesavehold[p]);
	  }
	avgiterGD = ctr;
	
	return;
      }
    
  }
  
  
  // upon return, coords will have the only point with min value out of all vals
  Array<OneD, Array<OneD, NekDouble> >find_roots( Array<OneD, NekDouble> uhats, StdExpansion *E , Array<OneD, Array<OneD, NekDouble> > &storage, NekDouble &avgiterGD, int d , int surfflag, int volflag, int surfid)
  {
    boost::ignore_unused(d);
    int dimension;
    if(surfflag == 0 && volflag == 0)
      dimension = 1;
    else if(volflag == 0)
      dimension = 2;
    else
      dimension = 3;

    Array<OneD, Array<OneD, NekDouble> > coords(dimension);
    boost::ignore_unused(surfid);

    //Confederate matrix approach
    if(surfflag == 0 && volflag == 0)
      {
	Array<OneD, Array<OneD, NekDouble> > retarr(1);

	vector<NekDouble>ret;
        while(true)
	  {
	    int N = uhats.size();
	    //	    cout<<"\n uhats:\n";
	    // for(int k = 0; k < uhats.size(); k++)
	    //   cout<<" "<<uhats[k];
	    
	    while(abs(uhats[N-1])<1e-10 && N > 1)
	      {
		N = N-1;

		if(N == 1)
		  {

		    ret.push_back(-1.0);
		    ret.push_back(1.0);
		    break;
		  }

	      }
	    if(N == 0)
	      break;
	    Array<OneD, Array<OneD,NekDouble> > Cnew(N-1);
	    for(int ii = 0; ii<N-1; ii++)
	      {
		Cnew[ii] = Array<OneD, NekDouble>(N-1);
		if(ii < N-1)
		  Vmath::Vcopy(N-1, &C[ii][0], 1, &Cnew[ii][0],1);
		  
	      }
	      
	    vector< NekDouble> temp2;
	    //cout<<"\n uhats =\n";
	    for(int k = 0; k < N-1; k++)
	      {
		uhats[k] = 0.5*(uhats[k])/uhats[N-1];
		//cout<<" "<<uhats[k];
	      }
	    for(int k = 0; k < N-1; k++)
	      {
		for(int jj= 0 ; jj < N-1; jj++)
		  {
		    if(k == N-2)
		      {
			Cnew[jj][k] = Cnew[jj][k] - uhats[jj];
		      }
		    
		      
		    temp2.push_back(  Cnew[jj][k]);
		  }
	      }
	    if( N>2)
	      {
		Cnew[N-3][N-2] =  - uhats[N-3];
	      }
	    else
	      {
		return retarr;	
	      }
	    Array<OneD, NekDouble>	      EIG_R   (N);
	    Array<OneD, NekDouble>	      EIG_I   (N);
	    FindEigenval(Cnew, EIG_R, EIG_I);

	    for(int kk = 0; kk <EIG_R.size(); kk++)
	      {
		ret.push_back( EIG_R[kk] );
	      }
	    ret.push_back(-1.0);
	    ret.push_back(1.0);
	    break;
	      
	  }
	retarr[0] = Array<OneD, NekDouble> (ret.size(), ret.data());

	return retarr;
      }

    else if(dimension > 1 && surfflag == 1 && volflag == 0 )
      {
	Array<OneD, Array<OneD, NekDouble> > coords(dimension);
	for(int k = 0; k < dimension; k++)
	  {
	    coords[k] = Array<OneD, NekDouble> (1);
	  }
	if(surfid == 0)
	  {
	    steepestgradient_descent2Dquad(uhats, E, storage, coords,  avgiterGD, testcoord2dqqpts, testcoord2dqqmidpts, interioreval2dqqmidpts, testcoord2dqlattice);

	  }
	else //tri
	  {
	    steepestgradient_descent2Dtri(uhats, E, storage, coords,  avgiterGD, testcoord2dtqpts, testcoord2dtqmidpts, interioreval2dtqmidpts, testcoord2dtlattice);
	
	  }
	return coords;

      }
    else if(volflag == 1)
      {
	for(int k = 0; k < dimension; k++)
	  {
	    coords[k] = Array<OneD, NekDouble> (1);
	  }
	
	steepestgradientdescent3D(uhats, E, storage, coords,  avgiterGD, testcoord3dqpts, testcoord3dqmidpts, interioreval3dqmidpts, testcoord3dlattice);
      }

    return coords;
  
  }
  
  
  po::options_description &GetOptions()
  {
    return m_desc;
  }

  po::variables_map &GetVariableMap()
  {
    return m_vm;
  }

  std::vector<string> &GetPointsType()
  {
    return m_pointstype;
  }
  vector<int> &GetPoints()
  {
    return m_points;
  }
  NekDouble &GetTol()
  {
    return iterGD;
  }
  NekDouble &GetEps()
  {
    return eps;
  }
  int &GetAvgNum()
  {
    return avgnum;
  }
  NekDouble &GetMaxIter()
  {
    return secarg;
  }
  NekDouble &GetGamHold()
  {
    return gamhold;
  }
  void SetC(Array<OneD, Array<OneD, NekDouble> >&newC)
  {
    this->C = newC;
  }
  NekDouble &GetCHold()
  {
    return chold;
  }

  vector<PointsKey> &GetPointsKey()
  {
    return m_pkey;
  }

  vector<BasisKey>  &GetBasisKey()
  {
    return m_bkey;
  }
  LibUtilities::ShapeType  &GetStype()
  {
    return stypeglo;
  }
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
  
  Array<OneD, Array<OneD, NekDouble>> GetQuadratureMidCoords(StdExpansion *E, Array<OneD, Array<OneD, NekDouble>> &coords)

  {

    int dimension = E->GetShapeDimension();
    //   const auto totPoints = (unsigned) E->GetTotPoints();
    int totPoints = coords[0].size()-1;//pow(E->GetBasis(0)->GetZ().size()-1, dimension);

    // form lattice by adding midpoints between each pair of quad points
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
    
protected:
  po::options_description m_desc;
  po::variables_map m_vm;
  vector<PointsKey> m_pkey;
  vector<BasisKey> m_bkey;
  NekDouble iterGD = 1e-2;
  NekDouble secarg = 1e3;
  NekDouble eps = -0.1;
  int avgnum = 1;
  NekDouble chold = 0.4;
  NekDouble gamhold = 0.7;
  std::string    m_shape;
  std::string    m_ntype;
  vector<string> m_basis{3, "NoBasisType"};
  vector<string> m_pointstype{3, "NoPointsType"};
  vector<int>    m_order;
  vector<int>    m_points;
  Array<OneD, Array<OneD, NekDouble> > C;
  LibUtilities::ShapeType stypeglo;
  Array<OneD, NekDouble> startarr;  
  Array<OneD, NekDouble> bcheb;  

};

#endif
