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
  Array<OneD, Array<OneD, NekDouble> > storage2dquad;
  Array<OneD, Array<OneD, NekDouble> > storage2dtri;

  Array<OneD, Array<OneD, NekDouble> > testcoord3dqpts;
  Array<OneD, Array<OneD, NekDouble> > testcoord3dqmidpts;
  Array<OneD, Array<OneD, NekDouble> > testcoord2dqqpts;
  Array<OneD, Array<OneD, NekDouble> > testcoord2dtqpts;
  Array<OneD, Array<OneD, NekDouble> > testcoord2dqqmidpts;
  Array<OneD, Array<OneD, NekDouble> > testcoord2dtqmidpts;

  Array<OneD, Array<OneD, NekDouble> > testcoord3dlattice;
  Array<OneD, Array<OneD, NekDouble> > testcoord2dqlattice;
  Array<OneD, Array<OneD, NekDouble> > testcoord2dtlattice;

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
		NekDouble t3 = ab2[i];
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

    //print J:

    for(int yy = 0; yy<J.size(); yy++)
      {
	for(int uu  = 0 ; uu<J[0].size(); uu++)
	  {
	    J[yy][uu] = J[yy][uu]/2;
	  }
      }

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


  vector< NekDouble> FindEigenval( vector<NekDouble> &uhatsdiff,
				   int N)
  {
    const unsigned int n = N-1;
    const int sz = (n)*(n);
    Array<OneD, NekDouble> CMdiff(sz,0.0);

    int sizemat = uhatsdiff.size();
    for(int k = 0; k < sizemat-2; k++)
      {
	CMdiff[(sizemat-1)*k+k+1] = 1.0;
      }
    
    Vmath::Vcopy(uhatsdiff.size()-1, &uhatsdiff[0], 1, &CMdiff[sz-(sizemat-1)], 1);
    Vmath::Smul(uhatsdiff.size()-1, -1.0,  &CMdiff[sz-(sizemat-1)], 1, &CMdiff[sz-(sizemat-1)], 1);

    const Array<OneD, const double>& A = CMdiff;

    Array<OneD, NekDouble>  wr(n), wi(n);
    Nektar::FullMatrixFuncs::EigenSolve(n, A, wr, wi, NullNekDouble1DArray);
    vector<NekDouble>EIG_R;
    for(int k = 0; k < sizemat-1; k++)
      {
	if(abs(wi[k])<1e-7 && abs(wr[k])<=1.0)
	  {
	    EIG_R.push_back(wr[k]);
	  }
      }
    
    return EIG_R;
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
    if(gamhold == 0)
      {
	exit(0);
      }
    boost::ignore_unused(testcoord2dlattice);
    int dim = 2;
    Array<OneD, NekDouble> g(dim);
    double inf = numeric_limits<double>::infinity();
    retarr = Array<OneD, Array<OneD, NekDouble> >(2);
    retarr[0] = Array<OneD,NekDouble>(1);
    retarr[1] = Array<OneD,NekDouble>(1);

    Array<OneD,NekDouble> xnew(dim) ,x(dim), denseeval(1e4+1), holdpq(1);;
    NekDouble gprev = inf;
    Array<OneD, NekDouble> nullarr(0), temp2(testcoord2d[0].size()), temp3, alltemp;
    vector<NekDouble> idxx, idxy;
    NekDouble gnew = inf;

    Array<OneD, Array<OneD, NekDouble> > tempdense (dim), coeffhold(dim);
    for(int k = 0; k < dim; k++)
      {
	tempdense[k] = Array<OneD, NekDouble>(1e4+1);
	coeffhold[k] = Array<OneD, NekDouble>(1);
	
      }
    int ct = 0;
    NekDouble allmin = inf;
    for(int y  = 0; y < 1e2; y++)
      {
	NekDouble valtmp = ( 1.0*y)/50 - 1.0;
	for(int u = 0; u < 1e2; u++)
	  {
	    tempdense[0][ct] = valtmp;
	    tempdense[1][ct] =(1.0*u)/50 - 1.0;
	    coeffhold[0][0] = valtmp;
	    coeffhold[1][0] =(1.0*u)/50 - 1.0;
	    Array<OneD, NekDouble > tmp  = E->PhysEvaluateBasis(coeffhold, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray  );
	    pq(uhats, coeffhold, tmp, nullarr, holdpq);

	    if(holdpq[0]<=allmin )
	      {
		if(abs(holdpq[0]-allmin)<1e-7)
		  {
		    idxx.push_back(tempdense[0][ct]);
		    idxy.push_back(tempdense[1][ct]);
		  }
		else
		  {
		    idxx.clear();
		    idxy.clear();
		    allmin = holdpq[0];
		    idxx.push_back(tempdense[0][ct]);
		    idxy.push_back(tempdense[1][ct]);

		  }
	      }
	    ct++;
	  }
      }
    tempdense[0][ct] = 1.0;
    tempdense[1][ct] = 1.0;

    coeffhold[0][0] = 1.0;
    coeffhold[1][0] = 1.0;

    Array<OneD, NekDouble> tmp  = E->PhysEvaluateBasis(coeffhold, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray  );

    pq(uhats, coeffhold, tmp, nullarr, holdpq);

    if(holdpq[0]<=allmin)
      {
	if(abs(holdpq[0]-allmin)<1e-7)
	  {
	    idxx.push_back(tempdense[0][ct]);
	    idxy.push_back(tempdense[1][ct]);
	  }
	else
	  {
	    idxx.clear();
	    idxy.clear();
	    allmin = holdpq[0];
	    idxx.push_back(tempdense[0][ct]);
	    idxy.push_back(tempdense[1][ct]);

	  }
      }
    cout<<"\n Dense lattice min val at ";
    for(int k = 0; k < idxx.size(); k++)
      {
	cout<<idxx[k] <<" ,"<<idxy[k]<<" val = "<<allmin<<"\n\n";
      }

    int sz = testcoord2d[0].size();
    Array<OneD, Array<OneD, NekDouble > > tempeval(4);
    for(int i = 0; i < 4; i++)
      {
	tempeval[i] = Array<OneD, NekDouble>(E->GetNcoeffs());
      }

    pq(uhats, testcoord2d, storage[0], nullarr, temp2 );

    gprev = Vmath::Vmin(sz, temp2, 1);
    pq(uhats, testcoord2dqmid, interioreval2dqmid, nullarr, temp3 );
    NekDouble gprevmid =  Vmath::Vmin(temp3.size(), temp3, 1);
    if(gprev > gprevmid)
      {
	gprev = gprevmid;
      }

    xnew[0] = 0.2;
    xnew[1] = -0.1;
    Array<OneD, NekDouble> xstart(2);
    xstart[0] = xnew[0];
    xstart[1] = xnew[1];
    Array<OneD, NekDouble> dereval(dim);
    fstream fio;
    fio.open("dumquad.txt", ios::app | ios::out | ios::in);
    fio<<"\n func = -2*x[i] + pow((x[i] + 0.6),3) + (pow(y[i],2) - 0.2); starting pt = "<<xnew[0]<<","<<xnew[1]<<" N = 5";

    if(gprev < 0 && abs(gprev)>1e-13)
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
	cout<<"\n begiin="<<t1[0];
	Array<OneD, NekDouble > holdxandval(3), saveholdxandval(3), savesavehold(3) ;//dim+1
	holdxandval[0] = xastaa[0][0];
	holdxandval[1] = xastaa[1][0];
	holdxandval[2] = t1[0];
	savesavehold[2] = inf;
	int truth_val = abs(g[0]+g[1]) > 1e-8;
	int ctr = 0, ct  = 0;

	NekDouble fac = 1,  gnew0 = gnew, gnew1 = gnew, gnew2 = gnew;
	if(truth_val)
	  {
	    NekDouble iter = secarg;
	    int counter_min_change = 0;

	    while(ctr < iter && (abs(g[0]) > 1e-9) && abs(g[1]) > 1e-9 && fac > 1e-6 && counter_min_change < 5 )
	      {
		if(ctr > 2 && abs(gnew0 - gnew2) < 1e-6)
		  break;
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
		cout<<"\n holdval = "<<holdval<<" saveholdxandval[2]="<<saveholdxandval[2]<<" fc = "<<fac<<" gsave[1] = "<<gsave[1]<<" gsave[0] = "<<gsave[0]<<" xastaa[0][0] = :"<<xastaa[0][0]		    <<" "<<xastaa[1][0]<<" saveholdxandval[0] =  "<<saveholdxandval[0]<<" saveholdxandval[1] = "<<saveholdxandval[1]<<"\n\n rhs="<<saveholdxandval[2] - c*fac*(gsave[0] + gsave[1]);

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
		    cout<<"\n holdval = "<<holdval<<"t1 = "<<t1[0]<<"\n";
		    if(abs(holdval -t1[0]) > 1e-5)
		      holdval = t1[0];
		    else
		      break;
		    cout<<"\n at ct = "<<ct<<" val="<<t1[0]<<" saveholdxandval[2] = "<<saveholdxandval[2]<<" fac = "<<fac<< " g = "<<g[0]<<"  "<<g[1]<<"\n";
		    cout<<"\n rhs="<<holdxandval[2] - c*fac*(g[0] + g[1]);
		  }
		// holdxandval[0] = xastaa[0][0];
		// holdxandval[1] = xastaa[1][0];
		// holdxandval[2] = t1[0];
		fio<<"\n ctr = "<< ctr<<" "<<xastaa[0][0]<<" "<<xastaa[1][0]<<" "<<t1[0]<<" ct \
= "<<ct;
		avgiterGD = ctr;
		if(saveholdval < holdval)
		  {
		    cout<<"\n 1... holdxandval="<<holdxandval[0]<<","<<holdxandval[1]<<","<< holdxandval[2]<< "\n before: min till now="<<savesavehold[2];
		    if((saveholdxandval[2]<savesavehold[2]) && abs(saveholdxandval[2] - savesavehold[2]) > 1e-5)
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
		    fio<<" "<<savesavehold[2];
		    cout<<"\n  after: min till now="<<savesavehold[2];
		    xastaa[0][0] = saveholdxandval[0];
		    xastaa[1][0] = saveholdxandval[1];
		    cout<<"\n saveholdxandval[2]="<<saveholdxandval[2]<< " t1[0] = "<<t1[0];

		    t1[0] = saveholdval;
		    cout<<"\n at ct = "<<ct<<" at ctr = "<<ctr<<"\n  returning:  "<<xastaa[0][0]
			<<","<< xastaa[1][0]<<" val = "<<saveholdval<<" \n g= "<<g[0]<<","<<g[1];

		  }
		else
		  {
 		    cout<<"\n 2... holdxandval="<<holdxandval[0]<<","<<holdxandval[1]<<","<< holdxandval[2];
		    cout<<"\n  before: min till now="<<savesavehold[2]<< "saveholdxandval[2] = "<<saveholdxandval[2]<<" ";;

		    if((saveholdxandval[2]<savesavehold[2]) && abs(saveholdxandval[2] - savesavehold[2]) > 1e-5)
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
		    cout<<"\n  after: min till now="<<savesavehold[2];
		    fio<<" "<<savesavehold[2];

		    cout<<"\n savesavehold[0] = "<<savesavehold[0]<<" "<<savesavehold[1]<<"  saveholdxandval[2]="<<saveholdxandval[2]<< " t1[0] = "<<t1[0];

		    t1[0] = holdval;
		    if(abs(xastaa[0][0])>1 || abs(xastaa[1][0]) > 1)
		      {
			xastaa[0][0] = savesavehold[0];
			xastaa[1][0] = savesavehold[1];
		      }
		    cout<<"\n 2.........at ct = "<<ct<<" at ctr = "<<ctr<<"  \n returning:  "<<xastaa[0][0]<<","<< xastaa[1][0]<<" val = "<<holdval<<" \n g= "<<g[0]<<","<<g[1];
		    fac = 1;

		  }
		int tval = ctr < iter && (abs(g[0]) > 1e-9) && abs(g[1]) > 1e-9 && fac > 1e-6;
		cout<<"\n tval ="<<tval<<"\n\n";
		cout<<"\n ctr  = "<<ctr<<" g[0] = "<<g[0]<<" g[1] = "<<g[1]<<" fac = "<<fac;
	      }
	  }
	cout<<"\n saveholdxandval[2]="<<saveholdxandval[2]<< " t1[0] = "<<t1[0];
	Array<OneD, NekDouble> t2;
	retarr = Array<OneD, Array<OneD, NekDouble> >(2);
	retarr[0] = Array<OneD, NekDouble>(1, xastaa[0][0]);
	retarr[1] = Array<OneD, NekDouble>(1, xastaa[1][0]);

	tempeval[0] = E->PhysEvaluateBasis(retarr, storage, tempeval[1], tempeval[2], tempeval[3]);
	pq(uhats, retarr, tempeval[0], nullarr, t1);
	cout<<"\n  val =at xastaa:"<<xastaa[0][0]<<" ,"<<xastaa[1][0]<<" == "<<t1[0];

	retarr[0][0] = savesavehold[0];
	retarr[1][0] = savesavehold[1];
	tempeval[0] = E->PhysEvaluateBasis(retarr, storage, tempeval[1], tempeval[2], tempeval[3]);
	pq(uhats, retarr, tempeval[0], nullarr, t2);
	cout<<"\n savesavehold="<<savesavehold[0]<<","<<savesavehold[1]<<"\n\nval at savesavehold="<<t2[0]<<" also, saved = "<<savesavehold[2];
	if(t1[0] < t2[0])
	  {
	    retarr[0][0] = xastaa[0][0];
	    retarr[1][0] = xastaa[1][0];
	  }
	tempeval[0] = E->PhysEvaluateBasis(retarr, storage, tempeval[1], tempeval[2], tempeval[3]);
	pq(uhats, retarr, tempeval[0], nullarr, t2);
	cout<<"\n returning:"<<retarr[0][0]<<","<<retarr[1][0]<<" and val = "<<t2[0];
	fio<<"\n*********************************************************\n";
	if(allmin > t2[0])
	  {
	    allmin = t2[0];
	  }
	if(abs(t2[0]-allmin) > 1e-4)
	  {
	    cout<<"\n \n fail!  gam ="<< gamhold<<", c = "<<chold<<"  "<<ctr <<" iters";
	    cout<<" \n expected min = "<<allmin<<" found = "<<t2[0]<<" error ="<< abs(t2[0]-allmin) ;
	    fio<<"\n \n fail!  gam ="<< gamhold<<", c = "<<chold<<"  "<<ctr <<" iters";
	    fio<<" \n expected min = "<<allmin<<" at "<< idxx[0]<<" "<<idxy[0]<<" found = "<<t2[0]<<" error ="<< abs(t2[0]-allmin) ;
	    

	  }
	else
	  {

	    cout<<"\n\n pass!  gam ="<< gamhold<<", c = "<<chold<<"  "<<ctr <<" iters";
	    cout<<" \n expected min = "<<allmin<<" found = "<<t2[0]<<" error ="<< abs(t2[0]-allmin) ;
	    fio<<"\n\n pass!  gam ="<< gamhold<<", c = "<<chold<<"  "<<ctr <<" iters";
	    fio<<" \n expected min = "<<allmin<<" at "<< idxx[0]<<" "<<idxy[0]<<" found = "<<t2[0]<<" error ="<< abs(t2[0]-allmin) ;

	  }
	fio.close();
	avgiterGD = ctr;

	return;
      }

  }

  void steepestgradient_descent2Dtri(Array<OneD, NekDouble> &uhats, StdExpansion *E, Array<OneD,Array<OneD, NekDouble> > &storage, Array<OneD, Array<OneD, NekDouble> >&retarr, NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> > &testcoord2d, Array<OneD, Array<OneD, NekDouble> > &testcoord2dqmid, Array<OneD, NekDouble> &interioreval2dqmid, Array<OneD, Array<OneD, NekDouble> > &testcoord2dlattice)
  {
    if(gamhold == 0)
      {
	exit(0);
      }
    boost::ignore_unused(testcoord2dlattice);
    int dim = 2;
    Array<OneD, NekDouble> g(dim);
    double inf = numeric_limits<double>::infinity();
    retarr = Array<OneD, Array<OneD, NekDouble> >(2);
    retarr[0] = Array<OneD,NekDouble>(1);
    retarr[1] = Array<OneD,NekDouble>(1);

    Array<OneD,NekDouble> xnew(dim) ,x(dim), denseeval(1e4+1), holdpq(1);;
    NekDouble gprev = inf;
    Array<OneD, NekDouble> nullarr(0), temp2(testcoord2d[0].size()), temp3, alltemp;
    vector<NekDouble> idxx, idxy;
    NekDouble gnew = inf;

    Array<OneD, Array<OneD, NekDouble> > tempdense (dim), coeffhold(dim);
    for(int k = 0; k < dim; k++)
      {
	tempdense[k] = Array<OneD, NekDouble>(1e4+1);
	coeffhold[k] = Array<OneD, NekDouble>(1);
      }
    int ct = 0;
    NekDouble allmin = inf;
    for(int y  = 0; y < 1e2; y++)
      {
	NekDouble valtmp = ( 1.0*y)/50 - 1.0;
	for(int u = 0; u < 1e2; u++)
	  {
	    tempdense[0][ct] = valtmp;
	    tempdense[1][ct] =(1.0*u)/50 - 1.0;
	    coeffhold[0][0] = valtmp;
	    coeffhold[1][0] =(1.0*u)/50 - 1.0;
	    if(coeffhold[0][0]+coeffhold[1][0] <= 1)
	      {
		Array<OneD, NekDouble > tmp  = E->PhysEvaluateBasis(coeffhold, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray  );
		pq(uhats, coeffhold, tmp, nullarr, holdpq);

		if(holdpq[0]<=allmin )
		  {
		    if(abs(holdpq[0]-allmin)<1e-7)
		      {
			idxx.push_back(tempdense[0][ct]);
			idxy.push_back(tempdense[1][ct]);
		      }
		    else{
		      idxx.clear();
		      idxy.clear();
		      allmin = holdpq[0];
		      idxx.push_back(tempdense[0][ct]);
		      idxy.push_back(tempdense[1][ct]);

		    }
		  }
	      }
	    ct++;
	  }
      }
    tempdense[0][ct] = 1.0;
    tempdense[1][ct] = 1.0;

    coeffhold[0][0] = 1.0;
    coeffhold[1][0] = 1.0;

    Array<OneD, NekDouble> tmp  = E->PhysEvaluateBasis(coeffhold, storage, NullNekDouble1DArray,NullNekDouble1DArray, NullNekDouble1DArray  );

    pq(uhats, coeffhold, tmp, nullarr, holdpq);

    if(holdpq[0]<=allmin)
      {
	if(abs(holdpq[0]-allmin)<1e-7)
	  {
	    idxx.push_back(1e2);
	    idxy.push_back(1e2);
	  }
	else{

	  idxx.clear();
	  idxy.clear();
	  allmin = holdpq[0];
	  idxx.push_back(1e2);
	  idxy.push_back(1e2);
	}
      }
    cout<<"\n Dense lattice min val at ";
    for(int k = 0; k < idxx.size(); k++)
      for(int k = 0; k < idxx.size(); k++)
	{
	  cout<<idxx[k] <<" ,"<<idxy[k]<<" val = "<<allmin<<"\n\n";
	}

    int sz = testcoord2d[0].size();
    Array<OneD, Array<OneD, NekDouble > > tempeval(4);
    for(int i = 0; i < 4; i++)
      {
	tempeval[i] = Array<OneD, NekDouble>(E->GetNcoeffs());
      }

    pq(uhats, testcoord2d, storage[0], nullarr, temp2 );

    gprev = Vmath::Vmin(sz, temp2, 1);
    pq(uhats, testcoord2dqmid, interioreval2dqmid, nullarr, temp3 );
    NekDouble gprevmid =  Vmath::Vmin(temp3.size(), temp3, 1);
    if(gprev > gprevmid)
      {
	gprev = gprevmid;
      }

    xnew[0] = 0.2;//-1.0;//0.1;//-0.9;
    xnew[1] = -0.2;//0.6;//0.1;//0.6;

    Array<OneD, NekDouble> xstart(2);
    xstart[0] = xnew[0];
    xstart[1] = xnew[1];
    //    NekDouble gnew = inf;
    cout<<"\n gprev = "<<gprev<<"\n\n";
    Array<OneD, NekDouble> dereval(dim);
    if(gprev < 0 && abs(gprev)>1e-13)
      {

	fstream fio;
	fio.open("dumtri.txt", ios::app | ios::out | ios::in);
	fio<<"\n func =-2*x[i] + pow((x[i] + 0.6),3) + (pow(y[i],2) - 0.2);;  = "<<xnew[0]<<","<<xnew[1]<<" N = 5";

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
	cout<<"\n begiin="<<t1[0];
	Array<OneD, NekDouble > holdxandval(3), saveholdxandval(3), savesavehold(3) ;

	holdxandval[0] = xastaa[0][0];
	holdxandval[1] = xastaa[1][0];
	holdxandval[2] = t1[0];
	savesavehold[2] = inf;
	int truth_val = abs(g[0]+g[1]) > 1e-8;
	cout<<"\n g = "<<g[0]<<" and "<<g[1];
	int ctr = 0, ct  = 0;
	NekDouble fac = 1, gnew0 = gnew, gnew1 = gnew, gnew2 = gnew;
	cout<<" dummy = "<<gnew0 + gnew1+gnew2<<"\n";;
	if(truth_val)
	  {
	    NekDouble iter = secarg;
	    int counter_min_change = 0;
	    cout<<"\n fac="<<fac<<" ter_min_change="<<counter_min_change<<" ctr="<<ctr<<"\n";
	    while(ctr < iter && (abs(g[0]) > 1e-9) && abs(g[1]) > 1e-9 && fac > 1e-6 && counter_min_change <5)
	      {
		cout<<"\n****counter_min_change="<<counter_min_change<<"****\n";
		if(ctr > 2 && abs(gnew0 - gnew2) < 1e-6)
		  {
		    cout<<"\n breaking!\n\n";
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
		cout<<"\n holdval = "<<holdval<<" saveholdxandval[2]="<<saveholdxandval[2]<<" fc = "<<fac<<" gsave[1] = "<<gsave[1]<<" gsave[0] = "<<gsave[0]<<" xastaa[0][0] = :"<<xastaa[0][0]<<" "<<xastaa[1][0]<<" saveholdxandval[0] =  "<<saveholdxandval[0]<<" saveholdxandval[1] = "<<saveholdxandval[1]<<"\n\n rhs="<<saveholdxandval[2] - c*fac*(gsave[0] + gsave[1]);
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
		    cout<<"\n holdval = "<<holdval<<"t1 = "<<t1[0]<<"\n";
		    if(abs(holdval -t1[0]) > 1e-5)
		      holdval = t1[0];
		    else
		      break;
		    cout<<"\n at ct = "<<ct<<" val="<<t1[0]<<" saveholdxandval[2] = "<<saveholdxandval[2]<<" fac = "<<fac<< " g = "<<g[0]<<"  "<<g[1]<<"\n";
		    cout<<"\n rhs="<<holdxandval[2] - c*fac*(g[0] + g[1]);
		  }
		fio<<"\n ctr = "<< ctr<<" "<<xastaa[0][0]<<" "<<xastaa[1][0]<<" "<<t1[0]<<" ct = "<<ct;
		avgiterGD = ctr;
		cout<<"\n faccccccc="<<fac;
		if(saveholdval < holdval )
		  {
		    cout<<"\n 1... holdxandval="<<holdxandval[0]<<","<<holdxandval[1]<<","<< holdxandval[2]<< "\n before: min till now="<<savesavehold[2];
		    if((saveholdxandval[2]<savesavehold[2]) && abs(saveholdxandval[2] - savesavehold[2]) > 1e-5)//saveholdxandval[2]< savesavehold[2])
		      {
			savesavehold[0] = saveholdxandval[0];
			savesavehold[1] = saveholdxandval[1];
			savesavehold[2] = saveholdxandval[2];
			counter_min_change = 0;
			cout<<"\n reset ctr min change\n";
		      }
		    else
		      {
			counter_min_change++;

			cout<<"\n ctr min change"<< counter_min_change<<"\n";
		      }
		    cout<<"\n  after: min till now="<<savesavehold[2];
		    fio<<" "<<savesavehold[2];
		    xastaa[0][0] = saveholdxandval[0];
		    xastaa[1][0] = saveholdxandval[1];
		    cout<<"\n saveholdxandval[2]="<<saveholdxandval[2]<< " t1[0] = "<<t1[0];

		    t1[0] = saveholdval;
		    cout<<"\n at ct = "<<ct<<" at ctr = "<<ctr<<"\n  returning:  "<<xastaa[0][0]<<","<<xastaa[1][0]<<" val = "<<saveholdval<<" \n g= "<<g[0]<<","<<g[1];

		  }
		else
		  {

		    cout<<"\n 2... holdxandval="<<holdxandval[0]<<","<<holdxandval[1]<<","<< holdxandval[2];
		    cout<<"\n  before: min till now="<<savesavehold[2]<< "saveholdxandval[2] = "<<saveholdxandval[2]<<" ";;

		    if((saveholdxandval[2]<savesavehold[2]) && abs(saveholdxandval[2] - savesavehold[2]) > 1e-5)
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
		    cout<<"\n counter_min_change = "<<counter_min_change<<" ct = "<<ct<<"\n";
		    cout<<"\n  after: min till now="<<savesavehold[2];
		    fio<<" "<<savesavehold[2];
		    cout<<"\n savesavehold[0] = "<<savesavehold[0]<<" "<<savesavehold[1]<<"  saveholdxandval[2]="<<saveholdxandval[2]<< " t1[0] = "<<t1[0];                                     
		    t1[0] = holdval;
		    if(abs(xastaa[0][0])>1 || abs(xastaa[1][0]) > 1 || (xastaa[0][0] + xastaa[1][0])>1)
		      {
			xastaa[0][0] = savesavehold[0];
			xastaa[1][0] = savesavehold[1];
		      }
		    cout<<"\n 2.........at ct = "<<ct<<" at ctr = "<<ctr<<"  \n returning:  "<<xastaa[0][0]<<","<< xastaa[1][0]<<" val = "<<holdval<<" \n g= "<<g[0]<<","<<g[1];
		    fac = 1;
		  }
		int tval = ctr < iter && (abs(g[0]) > 1e-9) && abs(g[1]) > 1e-9 && fac > 1e-6;
		cout<<"\n tval ="<<tval<<"\n\n";
		cout<<"\n ctr  = "<<ctr<<" g[0] = "<<g[0]<<" g[1] = "<<g[1]<<" fac = "<<fac<< " gneew="<<gnew;

	      }
	  }
	cout<<"\n saveholdxandval[2]="<<saveholdxandval[2]<< " t1[0] = "<<t1[0];
	Array<OneD, NekDouble> t2;
	retarr = Array<OneD, Array<OneD, NekDouble> >(2);
	retarr[0] = Array<OneD, NekDouble>(1, xastaa[0][0]);
	retarr[1] = Array<OneD, NekDouble>(1, xastaa[1][0]);

	tempeval[0] = E->PhysEvaluateBasis(retarr, storage, tempeval[1], tempeval[2], tempeval[3]);
	pq(uhats, retarr, tempeval[0], nullarr, t1);
	cout<<"\n  val =at xastaa:"<<xastaa[0][0]<<" ,"<<xastaa[1][0]<<" == "<<t1[0];

	retarr[0][0] = savesavehold[0];
	retarr[1][0] = savesavehold[1];
	tempeval[0] = E->PhysEvaluateBasis(retarr, storage, tempeval[1], tempeval[2], tempeval[3]);
	pq(uhats, retarr, tempeval[0], nullarr, t2);
	cout<<"\n savesavehold="<<savesavehold[0]<<","<<savesavehold[1]<<"\n\nval at savesavehold="<<t2[0]<<" also, saved = "<<savesavehold[2];

	if(t1[0] < t2[0])
	  {
	    retarr[0][0] = xastaa[0][0];
	    retarr[1][0] = xastaa[1][0];
	  }
	tempeval[0] = E->PhysEvaluateBasis(retarr, storage, tempeval[1], tempeval[2], tempeval[3]);
	pq(uhats, retarr, tempeval[0], nullarr, t2);
	cout<<"\n returning:"<<retarr[0][0]<<","<<retarr[1][0]<<" and val = "<<t2[0];


	fio<<"\n*********************************************************\n";
	if(allmin > t2[0])
	  {
	    allmin = t2[0];
	  }
	if(abs(t2[0]-allmin) > 1e-4)
	  {
	    cout<<"\n \n fail!  gam ="<< gamhold<<", c = "<<chold<<"  "<<ctr <<" iters";
	    cout<<" \n expected min = "<<allmin<<" found = "<<t2[0]<<" error ="<< abs(t2[0]-allmin) ;
	    fio<<"\n \n fail!  gam ="<< gamhold<<", c = "<<chold<<"  "<<ctr <<" iters";
	    fio<<" \n expected min = "<<allmin<<" at "<< idxx[0]<<" "<<idxy[0]<<" found = "<<t2[0]<<" error ="<< abs(t2[0]-allmin) ;

	  }
	else
	  {

	    cout<<"\n\n pass!  gam ="<< gamhold<<", c = "<<chold<<"  "<<ctr <<" iters";
	    cout<<" \n expected min = "<<allmin<<" found = "<<t2[0]<<" error ="<< abs(t2[0]-allmin) ;
	    fio<<"\n\n pass!  gam ="<< gamhold<<", c = "<<chold<<"  "<<ctr <<" iters";
	    fio<<" \n expected min = "<<allmin<<" at "<< idxx[0]<<" "<<idxy[0]<<" found = "<<t2[0]<<" error ="<< abs(t2[0]-allmin) ;

	  }
	fio.close();
	avgiterGD = ctr;

	return;
      }


  }


 void steepestgradientdescent3D(Array<OneD, NekDouble> &uhats, StdExpansion *E, Array<OneD, Array<OneD, NekDouble> > &storage, Array<OneD, Array<OneD, NekDouble> >&retarr, NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> > &testcoord3d, Array<OneD, Array<OneD, NekDouble> > &testcoord3dqmid, Array<OneD, NekDouble> &interioreval3dqmid, Array<OneD, Array<OneD, NekDouble> > &testcoord3dlattice)
    {
      int dim = 3;
      Array<OneD, NekDouble> g(dim);
      double inf = numeric_limits<double>::infinity();


      Array<OneD,NekDouble> xnew(dim) ,x(dim);
      NekDouble gprev = inf;
      int idxgprev;
      int sz = testcoord3d[0].size();
      Array<OneD, Array<OneD, NekDouble > > tempeval(4);

      for(int i = 0; i < 4; i++)
        {
          tempeval[i] = Array<OneD, NekDouble>(E->GetNcoeffs());
        }

      Array<OneD, NekDouble> nullarr(0), temp2, temp3, alltemp;
      pq(uhats, testcoord3d, storage[0], nullarr, temp2 );
      
      gprev = Vmath::Vmin(sz, temp2, 1);
      idxgprev = Vmath::Imin(temp2.size(), temp2, 1);

      pq(uhats, testcoord3dqmid, interioreval3dqmid, nullarr, temp3 );

      NekDouble gprevmid =  Vmath::Vmin(temp3.size(), temp3, 1);

      if(gprev > gprevmid)
        {
          gprev = gprevmid;
          idxgprev = sz + Vmath::Imin(temp3.size(), temp3, 1);
        }
      xnew[0] = testcoord3dlattice[0][idxgprev];
      xnew[1] = testcoord3dlattice[1][idxgprev];
      xnew[2] = testcoord3dlattice[2][idxgprev];

      Array<OneD, NekDouble> xstart(dim);
      xstart[0] = xnew[0];
      xstart[1] = xnew[1];
      xstart[2] = xnew[2];
      NekDouble   fnval, fnvalnew, dfval;
      NekDouble gnew = inf;
      Array<OneD, NekDouble> dereval(dim);
      Array<OneD, NekDouble> xsave(dim);
        
      if(gprev < 0 && abs(gprev)>1e-13)
        {
          Timer tt1;
          tt1.Start();

          Timer t12;
        

          NekDouble  iter = secarg;
          Array<OneD, Array<OneD, NekDouble> > xastaa(dim);
          Array<OneD, NekDouble> temp1 (uhats.size());
          Array<OneD, Array<OneD, NekDouble > > tempeval(4);
        

          for(int k = 0; k < tempeval.size(); k++)
            {
              tempeval[k] = Array<OneD, NekDouble>(uhats.size());
            }
          for(int p = 0; p < dim; p++)
            {
              xastaa[p] = Array<OneD, NekDouble>(1,xnew[p]);
            }

          xsave[0] = xnew[0];
   xsave[1] = xnew[1];
          xsave[2] = xnew[2];
        
          NekDouble c = chold;
           //find der at xnew:                                                                           
          xastaa[0][0] = xnew[0];
          xastaa[1][0] = xnew[1];
          xastaa[2][0] = xnew[2];
        
          tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2], tempeval[3]);
        

          for(int p = 0; p < dim; p++)
            {
              derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
              g[p] = dereval[p];
            }
          fnval = 0;
          Array<OneD, NekDouble> t1(1);
        
          pq(uhats, xastaa, tempeval[0], nullarr, t1);
           
          fnval = t1[0];
        
          gnew = pow((pow(g[0],2) + pow(g[1],2) + pow(g[2],2)) ,0.5);
          if(gnew <1e-6)
            {
              retarr[0][0] = (xastaa[0][0]);
              retarr[1][0] = (xastaa[1][0]);
              retarr[2][0] = (xastaa[2][0]);
              avgiterGD = 0;
              tt1.Stop();
              return;
            }
          dfval = c*(g[0] + g[1] + g[2]);

        
          xastaa[0][0] = xnew[0] - g[0]/gnew;
          xastaa[1][0] = xnew[1] - g[1]/gnew;
          xastaa[2][0] = xnew[2] - g[2]/gnew;
          NekDouble fac = gamhold;

          switch(E->DetShapeType())
            {
            case LibUtilities::eHexahedron:
              while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1)
                {
                  xastaa[0][0] = xnew[0] - fac*g[0]/gnew;
                  xastaa[1][0] = xnew[1] - fac*g[1]/gnew;
                  xastaa[2][0] = xnew[2] - fac*g[2]/gnew;
                  fac = fac*gamhold;

                }
              break;
            case LibUtilities::eTetrahedron:
              while(abs(xastaa[0][0])>1 || abs(xastaa[1][0])>1 || abs(xastaa[2][0] ) > 1 || xastaa[0][0] + xastaa[1][0] > 1e-9 || xastaa[0][0] + xastaa[2][0] > 1e-9 || xastaa[1][0] + xastaa[2][0] > 1e-9)
                {
        
                  xastaa[0][0] = xnew[0] - fac*g[0]/gnew;
                  xastaa[1][0] = xnew[1] - fac*g[1]/gnew;
                  xastaa[2][0] = xnew[2] - fac*g[2]/gnew;
                  fac = fac*gamhold;
                }
              break;
            default: cout<<"\n invalid element!\n";
              exit(0);
            }

          tempeval[0] = E->PhysEvaluateBasis(xastaa, storage, tempeval[1], tempeval[2], tempeval[3]);
          fnvalnew = 0;
          pq(uhats, xastaa, tempeval[0], nullarr, t1);
          fnvalnew = t1[0];
          NekDouble tls = 1, gam = fac;
          Array<OneD, NekDouble> gdnew(dim);
          Array<OneD, NekDouble> xold(dim);
          tls = tls*gam;
          int ctr = 0;

          while(ctr < iter && fnvalnew > fnval + tls*dfval/gnew)
            {
              ctr++;
              xastaa[0][0] = xnew[0] - tls*(g[0]/gnew);
              xastaa[1][0] = xnew[1] - tls*g[1]/gnew;
              xastaa[2][0] = xnew[2] - tls*g[2]/gnew;
              if(abs(xastaa[0][0] ) > 1 || abs(xastaa[1][0] ) > 1 || abs(xastaa[1][0] ) > 1 )
                {
                  tls = tls*gam;

                  continue;
                }

              tempeval[0] = E->PhysEvaluateBasis(xastaa,storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
              fnvalnew = 0;
              Array<OneD, NekDouble> t1(1);
              pq(uhats, xastaa, tempeval[0], nullarr, t1);
              fnvalnew = t1[0];

              tls= tls*gam;
            }
          xnew[0] = xastaa[0][0];//xnew[0] - tls*g[0]/gnew;
	  xnew[1] = xastaa[1][0];//xnew[1] - tls*g[1]/gnew;
	  xnew[2] = xastaa[2][0];//xnew[2] - tls*g[2]/gnew;
	   
          if(xnew[0] - xsave[0] < 1e-7 && xnew[1] - xsave[1] <1e-7  && xnew[2] - xsave[2] <1e-7 )
            {
              retarr[0][0] = (xnew[0]);
              retarr[1][0] = (xnew[1]);
              retarr[2][0] = (xnew[2]);
              iterGD = ctr;
              tt1.Stop();
        
              return;

            }
        

          avgiterGD = ctr;
          if(abs(xnew[0] ) > 1 || abs(xnew[1] ) > 1 || abs(xnew[2] ) > 1 )
            {

              retarr[0][0] = (xstart[0]);
              retarr[1][0] = (xstart[1]);
              retarr[2][0] = (xstart[2]);

            }
          else
            {

              retarr[0][0] = (xnew[0]);
              retarr[1][0] = (xnew[1]);
              retarr[2][0] = (xnew[2]);

            }
          tt1.Stop();
        

          return;
        }
      return;

    }



  
  
  // upon return, coords will have the only point with min value out of all vals
  Array<OneD, Array<OneD, NekDouble> >find_roots( Array<OneD, NekDouble> &uhats, StdExpansion *E , Array<OneD, Array<OneD, NekDouble> > &storage, NekDouble &avgiterGD, int d , int surfflag, int volflag, int surfid)
  {
    int dimension;
    if(surfflag == 0 && volflag == 0)
      dimension = 1;
    else if(volflag == 0)
      dimension = 2;
    else
      dimension = 3;

    Array<OneD, Array<OneD, NekDouble> > coords(dimension);

    //    if testing rootfinder only:
    if(1)
      {
	// for(NekDouble g =0.01; g < 0.95; g = g+0.04)
	//   {
	//     gamhold = g;
	//     for(NekDouble c = 0.01; c < 0.95; c = c+0.04)
	//       {
	//      chold = c;

	if(surfid == 0)
	  {
	    steepestgradient_descent2Dquad(uhats, E, storage, coords,  avgiterGD, testcoord2dqqpts, testcoord2dqqmidpts, interioreval2dqqmidpts, testcoord2dqlattice);
	    cout<<"\n avgitergd = "<<avgiterGD<<" ";

	  }
	else //tri
	  {

	    steepestgradient_descent2Dtri(uhats, E, storage, coords,  avgiterGD, testcoord2dtqpts, testcoord2dtqmidpts, interioreval2dtqmidpts, testcoord2dtlattice);
	    cout<<"\n avgitergd = "<<avgiterGD<<" ";

	  }
	//     }
	// }
	return coords;
      }
    else
      {
	boost::ignore_unused(surfid);
	//Confederate matrix approach
	if(surfflag == 0 && volflag == 0)
	  {

	    vector<NekDouble>ret;

	    while(true)
	      {

		int N = uhats.size();

		//vector<NekDouble> uhatsmon;
		while(abs(uhats[N-1])<1e-8 && N > 0)
		  {
		    N = N-1;

		    if(N == 0)
		      {

			ret.push_back(-1.0);
			ret.push_back(1.0);
			break;
		      }

		  }
		if(N == 0)
		  break;
		Array<OneD, NekDouble> uhatsmon;

		//          vector<NekDouble> temp(N);
		Array<OneD, NekDouble> temp2(N*N);
		int ct = 0;
		// convert uhats to monomial, find roots of uhats or der of uhats
		//          cout<<"\n uhatsmon=\n";

		for(int k = 0; k < N; k++)
		  {
		    for(int jj= 0 ; jj < N; jj++)
		      {
			//  temp[jj ] = C[jj][k];
			temp2[ct++] = this->C[jj][k];
		      }
		    //  Vmath::Vmul(N, &temp[0], 1,  &uhats[0], 1, &temp[0], 1);
		    //NekDouble temp2 =
		    //cout<<" "<<Vmath::Vsum(N, &temp[0], 1);
		    //  uhatsmon.push_back(temp2);
		  }

		//          cout<<"\n uhatsmon blas:\n";
		uhatsmon = blasmatvec(temp2, uhats, N, N);
		//for(int k = 0; k < uhatsmon.size(); k++)
		// cout<<" "<<uhatsmon[k]<<" ";
		if(abs(Vmath::Vmax(uhatsmon.size(), &uhatsmon[0], 1))<1e-10)
		  {
		    ret.push_back(-1.0);
		    ret.push_back(1.0);
		    break;
		  }

		// truncate trailing zeros
		while(abs(uhatsmon[N-1])<1e-8 && N > 0)
		  {
		    N = N-1;
		    if(N == 0)
		      {

			ret.push_back(-1.0);
			ret.push_back(1.0);
			break;

		      }

		  }
		if(N == 0)
		  break;
		//N = uhatsmon.size();

		// if(N == 0 || Vmath::Vmin(N,&uhatsmon[0],1)<-1e10)
		//   {
		//     ret[0].push_back(-1.0);
		//     ret[0].push_back(1.0);
		//     return ret;
		//   }
		// now size of uhatsmon = N;
		vector<NekDouble> uhatsdiff;
		// if d == 1,

		if(d == 1)
		  {
		    for(int k = 1; k < N; k++)
		      {
			uhatsdiff.push_back(k*uhatsmon[k]);
		      }
		    N = N-1;
		    if(N == 0)
		      {
			ret.push_back(-1.0);
			ret.push_back(1.0);

			break;
		      }

		  }
		else //d == 0
		  {
		    for(int k = 0; k<N; k++)
		      uhatsdiff.push_back(uhatsmon[k]);
		  }

		// if(N == 1)
		//   {
		//       // ret[0].push_back(-1.0);
		//       // ret[0].push_back(1.0);

		//     return ret;
		//   }

		Vmath::Smul(N, 1.0/uhatsdiff[N-1], &uhatsdiff[0], 1, &uhatsdiff[0], 1);

		vector<NekDouble> EIG_R = FindEigenval(uhatsdiff, N);


		for(int kk = 0; kk <EIG_R.size(); kk++)
		  {
		    ret.push_back( EIG_R[kk] );

		  }
		ret.push_back(-1.0);
		ret.push_back(1.0);

		break;

	      }
	    Array<OneD, Array<OneD, NekDouble> > retarr(1);
	    retarr[0] = Array<OneD, NekDouble> (ret.size(), ret.data()); // check if this is eq to ret

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
		cout<<"\n avgitergd = "<<avgiterGD<<" ";

	      }
	    else //tri
	      {

		steepestgradient_descent2Dtri(uhats, E, storage, coords,  avgiterGD, testcoord2dtqpts, testcoord2dtqmidpts, interioreval2dtqmidpts, testcoord2dtlattice);
		cout<<"\n avgitergd = "<<avgiterGD<<" ";

	      }
	    return coords;

	  }
	else if(volflag == 1)
	  {
	    for(int k = 0; k < dimension; k++)
	      {
		coords[k] = Array<OneD, NekDouble> (1);
	      }
	    //  coords = gradient_descent3D(uhats, E, storage,  avgiterGD, testcoord3dqpts, test
	    //coord3dqmidpts, interioreval3dqmidpts, testcoord3dlattice, sig);

	    steepestgradientdescent3D(uhats, E, storage, coords,  avgiterGD, testcoord3dqpts, testcoord3dqmidpts, interioreval3dqmidpts, testcoord3dlattice);
	    //cout<<"\n roots3d = "<<coords[0][0] <<" "<<coords[1][0]<<" "<<coords[2][0]<<"\n";
	  }

	return coords;
      }
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

  Array<OneD, Array<OneD, NekDouble>> GetQuadratureMidCoords(StdExpansion *E, Array<OneD, Array<\
							     OneD, NekDouble>> &coords)

  {

    int dimension = E->GetShapeDimension();
    //   const auto totPoints = (unsigned) E->GetTotPoints();
    int totPoints = (E->GetTotPoints()-1);//pow(E->GetBasis(0)->GetZ().size()-1, dimension);

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
    //    Array<OneD, Array<OneD, NekDouble>> coords = GetCoords( E );
    //Array<OneD, Array<OneD, NekDouble>> midcoords = GetQuadratureMidCoords( E );
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

protected:
  po::options_description m_desc;
  po::variables_map m_vm;
  vector<PointsKey> m_pkey;
  vector<BasisKey> m_bkey;
  NekDouble iterGD = 1e-2;
  NekDouble secarg = 1e3;
  NekDouble eps = -0.1;
  int avgnum = 1;
  NekDouble chold = 1.0;
  NekDouble gamhold = 0.1;
  std::string    m_shape;
  std::string    m_ntype;
  vector<string> m_basis{3, "NoBasisType"};
  vector<string> m_pointstype{3, "NoPointsType"};
  vector<int>    m_order;
  vector<int>    m_points;
  Array<OneD, Array<OneD, NekDouble> > C;
  LibUtilities::ShapeType stypeglo;
  
};

#endif
