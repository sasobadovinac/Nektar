///////////////////////////////////////////////////////////////////////////////
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

#ifndef DEMOS_STDREGIONS_STDDEMOSUPPORTNEW_HPP
#define DEMOS_STDREGIONS_STDDEMOSUPPORTNEW_HPP

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

#include <LibUtilities/Foundations/ManagerAccess.h>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>

using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;
namespace po = boost::program_options;
#include <LibUtilities/BasicUtils/Timer.h>

class DemoSupportNew
{
public:

  Array<OneD, Array<OneD, NekDouble> > testcoord3d; 
  Array<OneD, Array<OneD, NekDouble> > testcoord2dquad; 
  Array<OneD, Array<OneD, NekDouble> > testcoord2dtri; 
  Array<OneD, NekDouble > interioreval3d;
  Array<OneD, NekDouble > interioreval2dquad;
  Array<OneD, NekDouble > interioreval2dtri;
  Array<OneD, Array<OneD, NekDouble> > C;

  DemoSupportNew() : m_desc("Available options")
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
       "Max iters of gradient descent.")
      ("chold,c",
       po::value<NekDouble>(&chold)->multitoken(),
       "Max iters of gradient descent.")

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
    Array<OneD, Array<OneD, NekDouble> > C(J.size());
    for(i = 0; i < J.size(); i++)
      {
	C[i] = Array<OneD, NekDouble>(J[i].size(),J[i].data());
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

    if(m_pkey.size() == 0 && m_bkey.size() == 0)
      {
	for (int i = 0; i < dimension; ++i)
	  {
	    m_pkey.emplace_back(PointsKey(m_points[i], ptype[i]));
	    m_bkey.emplace_back(BasisKey(btype[i], m_order[i], m_pkey[i]));
	  }
      }
       

    switch (stype)
      {
      case ePoint:
	{
	  E = new StdPointExp(m_bkey[0]);
	  break;
	}
      case eSegment:
	{
	  E = new StdSegExp(m_bkey[0]);
	  break;
	}
      case eTriangle:
	{
	  E = nodaltype != eNoPointsType ? new StdNodalTriExp(m_bkey[0],
							      m_bkey[1],
							      nodaltype)
	    : new StdTriExp(m_bkey[0], m_bkey[1]);
	  break;
	}
      case eQuadrilateral:
	{
	  E = new StdQuadExp(m_bkey[0], m_bkey[1]);
	  break;
	}
      case eTetrahedron:
	{
	  E = nodaltype != eNoPointsType ? new StdNodalTetExp(m_bkey[0],
							      m_bkey[1],
							      m_bkey[2],
							      nodaltype)
	    : new StdTetExp(m_bkey[0], m_bkey[1],
			    m_bkey[2]);
	  break;
	}
      case ePyramid:
	{
	  E = new StdPyrExp(m_bkey[0], m_bkey[1], m_bkey[2]);
	  break;
	}
      case ePrism:
	{
	  E = nodaltype != eNoPointsType ? new StdNodalPrismExp(m_bkey[0],
								m_bkey[1],
								m_bkey[2],
								nodaltype)
	    : new StdPrismExp(m_bkey[0], m_bkey[1],
			      m_bkey[2]);
	  break;
	}
      case eHexahedron:
	{
	  E = new StdHexExp(m_bkey[0], m_bkey[1], m_bkey[2]);
	  break;
	}
      default:
	break;
      }
    C = formConf((dimension)*pow(3*m_order[0]+1,2)); 
    return E;
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

    
  vector< NekDouble> FindEigenval( vector<NekDouble> &uhatsdiff,
				   int N) 
  {    
    const unsigned int n = N-1;
    const int sz = (n)*(n);
    Array<OneD, NekDouble> CMdiff(sz,0.0);

    // for(int k  = 0 ; k < uhatsdiff.size(); k++)
    // {
    //     cout<<uhatsdiff[k]<<" ";

    // }
    // cout<<"\n";exit(0)
    int sizemat = uhatsdiff.size();
    for(int k = 0; k < sizemat-2; k++)
      {
	CMdiff[(sizemat-1)*k+k+1] = 1.0;
      }    
    // for(int k  = 0 ; k < CMdiff.size(); k++)
    // {
    //     cout<<CMdiff[k]<<" ";

    // }

    Vmath::Vcopy(uhatsdiff.size()-1, &uhatsdiff[0], 1, &CMdiff[sz-(sizemat-1)], 1);
    Vmath::Smul(uhatsdiff.size()-1, -1.0,  &CMdiff[sz-(sizemat-1)], 1, &CMdiff[sz-(sizemat-1)], 1);
        
     const Array<OneD, const double>& A = CMdiff;
    // cout<<"\n CMdiff:\n";
    // cout<<"\n";
    // exit(0);
    Array<OneD, NekDouble>  wr(n), wi(n);
    //    cout<<"\n 1 1\n\n";
    Nektar::FullMatrixFuncs::EigenSolve(n, A, wr, wi);
    //    cout<<"\n 1 1\n\n";
    vector<NekDouble>EIG_R;
    for(int k = 0; k < sizemat-1; k++)
      {
	if(abs(wi[k])<1e-7 && abs(wr[k])<=1.0)
	  {
	    EIG_R.push_back(wr[k]);
	  }
      }
    
    // EIG_R.push_back(-1.0);
    // EIG_R.push_back(1.0);


    return EIG_R;
  }



  Array<OneD, Array<OneD, NekDouble> > monomial_connection(int N)
  {

    double *a;
    double *b;

    Array<OneD, Array<OneD, NekDouble> > R(N);
    a = new double[N];
    b = new double[N];    


    // Initialize a and b to zero
    for(int i = 0; i < N; i++)
      {
	a[i] = 0.0;
	b[i] = 0.0;
      }
    
    // [a,b] = jacobi_recurrence(maxN+1, 0, 0);
    
    Polylib::RecCoeff(N, a, b, 0.0, 0.0);
    
    
    //     assert( (N >= numel(a)) && (N >= numel(b)) );

    //b = sqrt(b)
    transform(b, b+N, b, (double(*)(double)) sqrt);
    //b = 1./b
    Array<OneD, NekDouble> r(N, b);
    for(int i = 0; i<r.size(); i++)
      {
	r[i] = 1/r[i];
      }

    // local_leading_coeffs = cumprod(1./b);

    //r = cumprod(1./b) 
    partial_sum (r.begin(), r.end(), r.begin(), multiplies<double>());


    //Build the connection matrix for this dimension and then distribute it accordingly

    for(int i = 0; i<N; i++)
      {
	Array<OneD, NekDouble> row(N);
	for(int j = 0; j<N; j++)
	  {
	    if( i == j)
	      {
		row[i] = r[i];
	      }

	  }
	R[i] = row;

      }

    for(int i = 1; i<N; i++)
      {
	//"side" conditions (i.e. those w/o left/right boundary points)
	if(i<2)
	  {
	    R[i][0] = (1/b[i])*(-a[i-1]*R[i-1][0]);
	  }
	else
	  {
	    R[i][0] = (1/b[i])*(-a[i-1]*R[i-1][0] - b[i-1]*R[i-2][0]);
            
	    //All the rest are `vectorizable'
	    for(int j = 1; j<i-1; j++)
	      {
		R[i][j] = R[i-1][j-1] - a[i-1]*R[i-1][j] - b[i-1]*R[i-2][j];
		R[i][j] = R[i][j]/b[i];
	      }
            
	  }
      }

    return R;

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
  
    vector<NekDouble> Vsumsq;
    for(int i = 0; i < roots[0].size(); i++)
      {

	Vmath::Vmul(N, &V1[i], roots[0].size(), &uhats[0], 1, &w1[0], 1);
	fvals[i] = ( Vmath::Vsum(N, &w1[0], 1));
    
      }
    for( int i = 0; i < roots[0].size(); i++)
      {
	Vmath::Vmul(N, &V1[i], roots[0].size(), &V1[i], roots[0].size(), &w1[0], 1);

	Vsumsq.push_back(pow(Vmath::Vsum(N, &w1[0], 1),-0.5));
      }


    //    cts = vec_vec_dot(Vsumsq, ret);
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
    
    // //\sum(u*phi')
    // Vmath::Vmul(N, uhats, 1,Vxyd, 1, temp, 1);
    // NekDouble v1 = Vmath::Vsum(N, temp, 1);
  
    // v1 = v1*( Vmath::Vsum(N, temp, 1))*pow(tmp, -0.5);
    // //\sum(phi*phi')
    // Vmath::Vmul(N, Vxy, 1,Vxyd, 1, temp, 1);
    // NekDouble v2  = Vmath::Vsum(N, temp, 1);

    // //\sum(u*phi)
    // Vmath::Vmul(N, uhats, 1, Vxy, 1, temp, 1); 
    // v2 = v2*Vmath::Vsum(N, temp, 1);
    // //\sum(phi*phi)^(-3/2)
    // v2 = v2*( pow(tmp,-1.5));

    // ret = v1 + v2;
  } 

  
  //if flag = 0, call if not from optimizer
  //if flag = 1, call is from optimizer
  vector< vector< NekDouble> >  gradient_descent2Dquad(Array<OneD, NekDouble> &uhats, StdExpansion *E, NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> > &testcoord2dquad, Array<OneD, NekDouble> &interioreval2dquad, int flag)
  {

    // LibUtilities::Timer timer;
    // timer.Start();

    int dim = 2;
    Array<OneD, NekDouble> g(dim);
    double inf = numeric_limits<double>::infinity();
            

    Array<OneD,NekDouble> xnew(dim), x(dim) ;
    NekDouble gprev = inf; 
    int idxgprev;
    vector< vector< NekDouble> > ret(dim);
    int sz;
    Array<OneD, Array<OneD, NekDouble > > tempeval(4);

    sz = testcoord2dquad[0].size();
    
    for(int i = 0; i < 4; i++)
      {
	tempeval[i] = Array<OneD, NekDouble>(E->GetNcoeffs());  
      }
    Array<OneD, NekDouble> nullarr(0);
    Array<OneD, NekDouble> temp2(sz);

    if(flag==0)
      {
	Array<OneD, NekDouble> temp(uhats.size());
 
	for(int k = 0; k < sz; k++)
	  {

	    Vmath::Vmul(uhats.size(), &interioreval2dquad[k], sz, &uhats[0], 1, &temp[0], 1);          
	    temp2[k] = Vmath::Vsum(temp.size(), temp, 1);
	  }

	gprev = Vmath::Vmin(temp2.size(), temp2, 1);	
	idxgprev = Vmath::Imin(sz, temp2, 1);        
            
      }
    else
      {
	pq(uhats, testcoord2dquad, interioreval2dquad, nullarr, temp2 );
	gprev = Vmath::Vmin(sz, temp2, 1);
	idxgprev = Vmath::Imin(sz, temp2, 1);
      }
    xnew[0] = testcoord2dquad[0][idxgprev];
    xnew[1] = testcoord2dquad[1][idxgprev];
    //timer.Stop();
    //    cout<<"\n xnew="<<xnew[0]<<" "<<xnew[1]<<"  setup phase of GD 2d quad took "<<timer.TimePerTest(1)<<"s\n";

    int n = 0;
    //    NekDouble  prevepsl = 0;
    NekDouble saveminx = xnew[0], saveminy = xnew[1], savefnval = inf, fnval;;      
    if(gprev < 0 && abs(gprev)>1e-13)
      {
	NekDouble epsl = 1.0;

	NekDouble stepsize = iterGD, iter = secarg;
	Array<OneD, Array<OneD, NekDouble> > xastaa(dim);

	Array<OneD, NekDouble> dereval(dim);
	Array<OneD, NekDouble> xystart(xnew.size());
	for(int k = 0; k < dim; k++)
	  xystart[k] = xnew[k];
	
	// loop for gradient descent 2D:
        
	bool res2;
	Array<OneD, NekDouble> temp1 (uhats.size());

	//timer.Start();

	while(n < iter && epsl > 1e-8)
	  {
	    //prevepsl = epsl;
	    epsl = 0.0;
	    for(int p = 0; p < dim; p++)
	      {

		xastaa[p] = Array<OneD, NekDouble>(1,xnew[p]);
	      }   
           
	    E->PhysEvalBasisGrad(xastaa, tempeval[0], tempeval[1], tempeval[2], tempeval[3]);
	    // multiply tempeval by uhats to get g[p]
           
	    if(flag == 0) // from opt_needed
	      {
                
		    fnval = 0;
		    for(int p = 0; p < dim; p++)
		      {
			
			Vmath::Vmul(uhats.size(), &tempeval[p+1][0], 1, &uhats[0], 1, &temp1[0], 1 );
			g[p] = Vmath::Vsum(temp1.size(), temp1, 1);
			
			epsl += abs(g[p]);
			Vmath::Vmul(uhats.size(), &tempeval[p][0], 1, &uhats[0], 1, &temp1[0], 1 );
			fnval += Vmath::Vsum(temp1.size(), temp1, 1);
		      }
		  
	      }
	    else
	      {
		for(int p = 0; p < dim; p++)
		  {
		    fnval = 0;
		    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
                 
		    g[p] = dereval[p];
		    epsl += abs(g[p]);
		    
		  }
		Array<OneD, NekDouble> t1(1);
		pq(uhats, xastaa, tempeval[0], nullarr, t1);
		fnval = t1[0];

		
	      }
	    if(savefnval > fnval  && abs(xastaa[0][0] )<=1 && abs(xastaa[1][0])<=1)
	      {
		savefnval = fnval;
		saveminx = xastaa[0][0];
		saveminy = xastaa[1][0];
	      }
	    
	    x[0] = xnew[0] - stepsize*(g[0]);
	    x[1] = xnew[1] - stepsize*(g[1]);
	    res2 = true;
	    for(int k = 0; k < dim; k++)
	      {   
	    	if(abs(x[k]-xnew[k]) > 1e-4 || abs(x[k])>1)
	    	  {
	    	    res2 = res2 && false  ;
	    	  }
	      }
	    if(res2)
	   
	      {
		saveminx = x[0];
		saveminy = x[1];
		break;
	      }
	    if((abs(x[0]) >1.0) || (abs(x[1])>1.0) || (std::isnan(epsl)) )
	      {
		
		stepsize = stepsize/10;
		xnew[0] = xystart[0];
		xnew[1] = xystart[1];
		
                if(stepsize <1e-5)
		  {
		    xnew[0] = saveminx;
		    xnew[1] = saveminy;
		    break;
		  }
		continue;
	      }
	    for(int p = 0; p < dim; p++)
	      {
		xnew[p] = x[p];
	      }            
            
	    n = n+1;
            
	  }
   	// res2 = true;
	// for(int k = 0; k < dim; k++)
	//   {
	//     if(abs(x[k]) >1.0)
	//       {
	// 	res2 = false;
		
	//       }
	//   }
	// if((res2 && n<iter && (!std::isnan(epsl))))
	//   {
	// 	for(int p = 0; p < dim; p++)
	// 	  {
	// 	    ret[p].push_back( (xnew[p]));
	// 	  }
	//     timer.Stop();
	//     cout<<"\n 1.  GD iters took "<<timer.TimePerTest(1)<<"s and "<<n<<"  iters  ret = "<<xnew[0]<<" "<<xnew[1]<<"\n\n";
	//     avgiterGD = n;
	//     return ret;
	//   }

	    
	// ret[0].push_back(saveminx);
	// ret[1].push_back(saveminy);
	
	//	timer.Stop();
	ret[0].push_back(saveminx);
	ret[1].push_back(saveminy);
	avgiterGD = 0;
	//cout<<"\n quad  GD iters took "<<timer.TimePerTest(1)<<"s and "<<n<<"  iters. Returning="<<ret[0][0]<<" "<<ret[1][0]<<" x[0]="<<x[0]<<" x[1]="<<x[1]<<"  stepsize = "<<stepsize<<" ";//ret = "<<xnew[0	

      }
    return ret;

    
    
  }
  

  //if flag = 0, call if not from optimizer
  //if flag = 1, call is from optimizer
  vector< vector< NekDouble> >  gradient_descent2Dtri(Array<OneD, NekDouble> &uhats, StdExpansion *E, NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> > &testcoord2dtri, Array<OneD, NekDouble> &interioreval2dtri, int flag)
  {
    LibUtilities::Timer timer;
    timer.Start(); 
    int dim = 2;
    Array<OneD, NekDouble> g(dim);
    double inf = numeric_limits<double>::infinity();
            

    Array<OneD,NekDouble> xnew(dim), x(dim);
    NekDouble gprev = inf; 
    int idxgprev;
    vector< vector< NekDouble> > ret(dim);
    int sz;
    Array<OneD, Array<OneD, NekDouble > > tempeval(4);

    sz = testcoord2dtri[0].size();


    for(int i = 0; i < 4; i++)
      {
	tempeval[i] = Array<OneD, NekDouble>(E->GetNcoeffs());  
      }
    Array<OneD, NekDouble> nullarr(0);
    Array<OneD, NekDouble> temp2(sz);

    if(flag==0)
      {
	Array<OneD, NekDouble> temp(uhats.size());
	for(int k = 0; k < sz; k++)
	  {

	    Vmath::Vmul(uhats.size(), &interioreval2dtri[k], sz, &uhats[0], 1, &temp[0], 1);          
	    temp2[k] = Vmath::Vsum(temp.size(), temp, 1);
	  }

	gprev = Vmath::Vmin(temp2.size(), temp2, 1);	
	idxgprev = Vmath::Imin(sz, temp2, 1);        
            
      }
    else
      {
	pq(uhats, testcoord2dtri, interioreval2dtri, nullarr, temp2 );
	gprev = Vmath::Vmin(sz, temp2, 1);
	idxgprev = Vmath::Imin(sz, temp2, 1);
      }
    xnew[0] = testcoord2dtri[0][idxgprev];
    xnew[1] = testcoord2dtri[1][idxgprev];
    timer.Stop();
    //    cout<<"\n xnew="<<xnew[0]<<" "<<xnew[1]<<"  setup phase of GD 2d tri took "<<timer.TimePerTest(1)<<"s\n";
    //       cout<<"\n start: "<<xnew[0]<<" "<<xnew[1]<<" gprev ="<<gprev<<"\n";;
												       //	cout<<"\n start: ("<< xnew[0] <<" ,"<< xnew[1]<<") gprev="<<gprev<<" xnew[0]+xnew[1] = "<<xnew[0]+xnew[1]<<"\n\n";  
													if(xnew[0]+xnew[1] > 0 && abs(xnew[0]+xnew[1])>1e-10)
      {
	exit(0);
      }
    int n  = 0;
    //    NekDouble prevepsl = 0,
    NekDouble saveminx = xnew[0], saveminy = xnew[1], savefnval = inf, fnval;  

													 if(gprev < 0 && abs(gprev)>1e-13)
      {
	NekDouble epsl = 1.0;

	NekDouble stepsize = iterGD, iter = secarg;
	Array<OneD, Array<OneD, NekDouble> > xastaa(dim);

	Array<OneD, NekDouble> dereval(dim);
	Array<OneD, NekDouble> xystart(xnew.size());
	for(int k = 0; k < dim; k++)
	  xystart[k] = xnew[k];
	
	// loop for gradient descent 2D:
        
	bool res, res2;//,  res3 = true;
        
	Array<OneD, NekDouble> temp1 (uhats.size());
	timer.Start();    
	
	while(n < iter  && epsl > 1e-8)
	  {
	    
	    //prevepsl = epsl;
	    epsl = 0.0;
	    for(int p = 0; p < dim; p++)
	      {
		xastaa[p] = Array<OneD, NekDouble>(1,xnew[p]);
	      }   
           
	    E->PhysEvalBasisGrad(xastaa, tempeval[0], tempeval[1], tempeval[2], tempeval[3]);
	    // multiply tempeval by uhats to get g[p]
           
	    if(flag == 0) // from opt_needed
	      {
		fnval = 0;
                for(int p = 0; p < dim; p++)
		  {
                 
		    Vmath::Vmul(uhats.size(), &tempeval[p+1][0], 1, &uhats[0], 1, &temp1[0], 1 );
		    g[p] = Vmath::Vsum(temp1.size(), temp1, 1);
    
		    epsl += abs(g[p]);
		    Vmath::Vmul(uhats.size(), &tempeval[p][0], 1, &uhats[0], 1, &temp1[0], 1 );
		    fnval += Vmath::Vsum(temp1.size(), temp1, 1);
		  }
		
	      }
	    else
	      {
		fnval = 0;
		for(int p = 0; p < dim; p++)
		  {

		    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
                 
		    g[p] = dereval[p];
		    epsl += abs(g[p]);
		  }
		Array<OneD, NekDouble> t1(1);
		pq(uhats, xastaa, tempeval[0], nullarr, t1);
		fnval = t1[0];
	      }
	    if(savefnval > fnval && xastaa[0][0]+ xastaa[1][0] <= 0 && abs(xastaa[0][0] )<=1 && abs(xastaa[1][0])<=1)
	      {
		savefnval = fnval;
		saveminx = xastaa[0][0];
		saveminy = xastaa[1][0];
	      }
	    x[0] = xnew[0] - stepsize*(g[0]);
	    x[1] = xnew[1] - stepsize*(g[1]);

	    res = true;
	    res = (x[0]+x[1])<0;
	    res2 = true;
	    for(int k = 0; k < dim; k++)
	      {   
	    	if(abs(x[k]-xnew[k]) > 1e-5 || abs(x[k])>1)
	    	  {
	    	    res2 = res2 && false  ;
	    	  }
	      }
	    if(res2 && res)
	      {
		saveminx = x[0];
		saveminy = x[1];
		break;
	      }

	    // if( ((abs(epsl) <1e-8) || (abs(epsl)>1e2) || abs(prevepsl - epsl)<1e-5 || (abs(x[0]-xnew[0]) < 1e-5 && abs(x[1]-xnew[1]) < 1e-5)) && abs(x[1])<=1.0 && abs(x[0] )<=1.0)
	    //   {
	    // 	break;
	    //   }
	    if((abs(x[0]) >1.0) || (abs(x[1])>1.0) || (std::isnan(epsl)) )
	      {
		stepsize = stepsize/10;
		xnew[0] = xystart[0];
		xnew[1] = xystart[1];

                if(stepsize <1e-6)
		  {
		    xnew[0] = saveminx;
		    xnew[1] = saveminy;
		    break;
		  }
		continue;
	      }

	    for(int p = 0; p < dim; p++)
	      {
		xnew[p] = x[p];
	      }            
            
	    n = n+1;
 
	  }
	// //	cout<<"\n n = "<<n<<" ";	
	// res = true;
	// res = (x[0]+x[1])<=0;
        // //check
	// res2 = true;
	// for(int k = 0; k < dim; k++)
	//   {
	//     if(abs(x[k]) >1.0)
	//       {
	// 	res2 = false;

	//       }
	//   }
	// if((res && res2 && (n<iter) && (!std::isnan(epsl))))  
	// {
	//   // if( n <iter ) 
	//   //   {
	//       for(int p = 0; p < dim; p++)
	//       {
	// 	ret[p].push_back( (xnew[p]));
	//       }
	//   timer.Stop();   
	//   cout<<"\n 1. GD iters took "<<timer.TimePerTest(1)<<"s and "<<n<<"  iters  ret = "<<xnew[0]<<" "<<xnew[1]<<"\n";
	//   avgiterGD = n;	  

	// }
	// else
	//   {
	    
	//     ret[0].push_back(saveminx);
	//     ret[1].push_back(saveminy);  
	//     timer.Stop();
	//     cout<<"\n 2. GD iters took "<<timer.TimePerTest(1)<<"s and "<<n<<"  iters  ";//ret = "<<xnew[0]<<" "<<xnew[1]<< "  g[0] + g[1] ="<<g[0]+g[1] <<" gprev = "<<gprev<<" x[0] = "<<x[0]<<" x[1] = "<<x[1] <<" xnew[0] ="<<xnew[0]<<" xnew[1] = "<<xnew[1]<<"res ="<<((x[0]+x[1])<0)<<" saveminx ="<<saveminx<<" saveminy = "<<saveminy<<"\n"; 
	    
	//     avgiterGD = n;
	      
	//   }
	timer.Stop();
	ret[0].push_back(saveminx);
	ret[1].push_back(saveminy);
	avgiterGD = n;
	//	cout<<"\n tri  GD iters took "<<timer.TimePerTest(1)<<"s and "<<n<<"  iters  stepsize = "<<stepsize;//ret = "<<xnew[0]<<" "<<xnew[1]<< "  g[0] + g[1] ="<<g[0]+g[1] <<" gprev = "<<gprev<<" x[0] = "<<x[0]<<" x[1] = "<<x[1] <<" xnew[0] ="<<xnew[0]<<" xnew[1] = "<<xnew[1]<<"res ="<<((x[0]+x[1])<0)<<" saveminx ="<<saveminx<<" saveminy = "<<saveminy<<"\n"; 
	
      }
													 													 return ret;
  }

  vector< vector< NekDouble> >  steepestgradientdescent3D(Array<OneD, NekDouble> &uhats, StdExpansion *E, NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> > &testcoord3d, Array<OneD, NekDouble> &interioreval3d, int flag)
  {
    int dim = 3;

    Array<OneD, NekDouble> g(dim);
    double inf = numeric_limits<double>::infinity();
        
        
    Array<OneD,NekDouble>  x(dim);
    vector<NekDouble> xnew;
    NekDouble gprev = inf;
    int idxgprev;
    vector< vector< NekDouble> > ret(dim);
    int sz = testcoord3d[0].size();
    Array<OneD, Array<OneD, NekDouble > > tempeval(4);

    for(int i = 0; i < 4; i++)
      {
	tempeval[i] = Array<OneD, NekDouble>(E->GetNcoeffs());  
      }
    Array<OneD, NekDouble> nullarr(0);
    Array<OneD, NekDouble> temp2(sz);
    Timer t;
    t.Start();
    if(flag==0)
      {
            
	Array<OneD, NekDouble> temp(uhats.size());
            
	for(int k = 0; k < sz; k++)
	  {
                
	    Vmath::Vmul(uhats.size(), &interioreval3d[k], sz, &uhats[0], 1, &temp[0], 1);          
	    temp2[k] = Vmath::Vsum(temp.size(), temp, 1);
	  }
            	
	gprev = Vmath::Vmin(temp2.size(), temp2, 1);	
	idxgprev = Vmath::Imin(sz, temp2, 1);        
            
      }
    else
      {
	pq(uhats, testcoord3d, interioreval3d, nullarr, temp2 );
	gprev = Vmath::Vmin(sz, temp2, 1);
	idxgprev = Vmath::Imin(sz, temp2, 1);
      }
    t.Stop();
    cout<<"\n eval time:"<<t.TimePerTest(1)<<" s \n";;
        
    xnew.push_back(testcoord3d[0][idxgprev]);//[0] = ;
    xnew.push_back(testcoord3d[1][idxgprev]);
    xnew.push_back( testcoord3d[2][idxgprev]);
    Array<OneD, NekDouble> xstart(3);
    xstart[0] = xnew[0];
    xstart[1] = xnew[1];
    xstart[2] = xnew[2];
    cout<<"\n start: ("<< xnew[0] <<" ,"<< xnew[1]<<" ,"<< xnew[2]<<") gprev="<<gprev<<"\n";
    // if(xnew[0]+xnew[2] >0 || xnew[1]+xnew[2] >0)
    //   {
    //  
    //   }
    //timer.Stop();
    //    cout<<"\n xnew="<<xnew[0]<<" "<<xnew[1]<<" "<<xnew[2]<<"  setup phase of GD 3D took "<<timer.TimePerTest(1)<<"s gprev="<<gprev<<"\n";
    // prevepsl = 0,
    NekDouble   fnval, fnvalnew, dfval;   
    //NekDouble saveminx = xnew[0], saveminy = xnew[1], saveminz = xnew[2], savefnval = inf;
    // timer.Stop();
    // tmptime = timer.TimePerTest(1);
    // cout<<"\n setup for GD 3d: "<<tmptime<<"gprev = "<<gprev<<"\n";
    // tmptime = 0.0;
    // NekDouble pqdertime = 0.0;
    // timer.Start(); 
    NekDouble gnew = inf;
    Array<OneD, NekDouble> dereval(dim);
    Array<OneD, NekDouble> xsave(3);

    if(gprev < 0 && abs(gprev)>1e-13)
      {
	int n = 0;
	NekDouble  iter = secarg;
	Array<OneD, Array<OneD, NekDouble> > xastaa(dim);
	Array<OneD, NekDouble> temp1 (uhats.size());
	Array<OneD, Array<OneD, NekDouble > > tempeval(4);
	//NekDouble fnval, fnvalnew;

	for(int k = 0; k < tempeval.size(); k++)
	  {
	    tempeval[k] = Array<OneD, NekDouble>(uhats.size());
	  }
	for(int p = 0; p < dim; p++)
	  {
	    xastaa[p] = Array<OneD, NekDouble>(1,xnew[p]);
	  }   
	NekDouble c = chold;

	while (n<iter)
	  {
	    		
	    xsave[0] = xnew[0];
	    xsave[1] = xnew[1];
	    xsave[2] = xnew[2];

	    //find der at xnew:
	    xastaa[0][0] = xnew[0];
	    xastaa[1][0] = xnew[1];
	    xastaa[2][0] = xnew[2];
	    E->PhysEvalBasisGrad(xastaa, tempeval[0], tempeval[1], tempeval[2], tempeval[3]);

	    if(flag == 0) // from opt_needed
	      {
		for(int p = 0; p < dim; p++)
		  {
		    
		    Vmath::Vmul(uhats.size(), &tempeval[p+1][0], 1, &uhats[0], 1, &temp1[0], 1 );
		    g[p] = Vmath::Vsum(temp1.size(), temp1, 1);
		  }
		fnval = 0;
		Vmath::Vmul(uhats.size(), &tempeval[0][0], 1, &uhats[0], 1, &temp1[0], 1 );
		fnval = Vmath::Vsum(temp1.size(), temp1, 1);
		
		
	      }
	    
	    else // from optimizer
	      {
		for(int p = 0; p < dim; p++)
		  {
		    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
		    
		    g[p] = dereval[p];
		    
		  }
		fnval = 0;
		Array<OneD, NekDouble> t1(1);
		pq(uhats, xastaa, tempeval[0], nullarr, t1);
		fnval = t1[0];
		
	      }
	    //	cout<<"\n here: g[0]="<<g[0]<<" g[1]= "<<g[1]<<" g[2] = "<<g[2];	
	    //L2 norm
	    gnew = pow((pow(g[0],2) + pow(g[1],2) + pow(g[2],2)) ,0.5);
	    //cout<<"\n n = "<<n <<" gnew="<<gnew<<"\n";
	    if(gnew <1e-6)
	      {
	        ret[0].push_back(xastaa[0][0]);
	        ret[1].push_back(xastaa[1][0]);
	        ret[2].push_back(xastaa[2][0]);
	        return ret;
	      }
	    dfval = c*(g[0] + g[1] + g[2]);
	    
	    xastaa[0][0] = xastaa[0][0] - g[0]/gnew;
	    xastaa[1][0] = xastaa[1][0] - g[1]/gnew;
	    xastaa[2][0] = xastaa[2][0] - g[2]/gnew;
	    if(abs(xastaa[0][0] ) > 1 || abs(xastaa[1][0] ) > 1 || abs(xastaa[2][0] ) > 1)
	      {
		break;
		
	      }
	    E->PhysEvalBasisGrad(xastaa, tempeval[0], tempeval[1], tempeval[2], tempeval[3]);
	    if(flag == 0) // from opt_needed
	      {
		fnvalnew = 0;
		Vmath::Vmul(uhats.size(), &tempeval[0][0], 1, &uhats[0], 1, &temp1[0], 1 );
		fnvalnew = Vmath::Vsum(temp1.size(), temp1, 1);
		
		
	      }
	    
	    else // from optimizer
	      {
		fnvalnew = 0;
		Array<OneD, NekDouble> t1(1);
		pq(uhats, xastaa, tempeval[0], nullarr, t1);
		fnvalnew = t1[0];
		
	      }
	    NekDouble t = 1, gam = gamhold;	
	    
	    
	    Array<OneD, NekDouble> gdnew(3);
	    Array<OneD, NekDouble> xold(3);
	    t = t*gam;
	    //	    cout<<" t = "<<t<<" fnvalnew = "<<fnvalnew <<" fnval + t*dfval = "<<fnval - t*dfval<<"\n";
	    while(fnvalnew > fnval + t*dfval/gnew)
	      {
		xastaa[0][0] = xnew[0] - t*(g[0]/gnew);
		xastaa[1][0] = xnew[1] - t*g[1]/gnew;
		xastaa[2][0] = xnew[2] - t*g[2]/gnew;
		if(abs(xastaa[0][0] ) > 1 || abs(xastaa[1][0] ) > 1 || abs(xastaa[2][0] ) > 1)
		  {
		    t = t*gam;

		    continue;
		  }
		E->PhysEvalBasisGrad(xastaa, tempeval[0], NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
		if(flag == 0) // from opt_needed  
		  {
		    fnvalnew = 0;
		    Vmath::Vmul(uhats.size(), &tempeval[0][0], 1, &uhats[0], 1, &temp1[0], 1 );
		    fnvalnew = Vmath::Vsum(temp1.size(), temp1, 1);
		    
		  }
		else
		  {
		    
		    fnvalnew = 0;
		    Array<OneD, NekDouble> t1(1);
		    pq(uhats, xastaa, tempeval[0], nullarr, t1);
		    fnvalnew = t1[0];
		    
		  }
		
		t = t*gam;
	      }
	    xnew[0] = xnew[0] - t*g[0]/gnew;
	    xnew[1] = xnew[1] - t*g[1]/gnew;
	    xnew[2] = xnew[2] - t*g[2]/gnew;
	    //	    cout<<"\n here: t="<<t<<" xnew = "<<xnew[0]<<" "<<xnew[1]<<" "<<xnew[2]<<" n= "<<n<<" ";
	   
	    if(xnew[0] - xsave[0] < 1e-6 && xnew[1] - xsave[1] <1e-6 && xnew[2] - xsave[2] < 1e-6)
	      {
		break;
		
	      }
	    n = n + 1;
	  }
	avgiterGD = n;
	if(abs(xnew[0] ) > 1 || abs(xnew[1] ) > 1 || abs(xnew[2] ) > 1)
	  {
	    ret[0].push_back(xstart[0]);
	    ret[1].push_back(xstart[1]);
	    ret[2].push_back(xstart[2]);
	  }
	else
	  {
	    ret[0].push_back(xnew[0]);
	    ret[1].push_back(xnew[1]);
	    ret[2].push_back(xnew[2]);
	  }
	return ret;
	
      }
    return ret;
  }
  
  //   }
	    
    // 		cout<<"\n here fvalnew = "<<fnvalnew<<" fnval + c*gam*(g[0] + g[1] + g[2])) = "<<fnval+c*gam*(g[0] + g[1] + g[2])<<"\n\n";
		
    // 		if(fnvalnew < fnval + c*gam*(g[0] + g[1] + g[2]))// || std::isnan(fnvalnew)|| abs(fnvalnew) <1e-14)
    // 		  {
    // 		    break;
    // 		  }
	      
    // 		a = a*gam;
    // 		//fnval = fnvalnew;
    // 		// if(a <1e-12)
    // 		//   break;
    // 	      }
    // 	    cout<<" xastaa="<<xastaa[0][0]<<" "<<xastaa[1][0]<<" "<<xastaa[2][0]<<"\n";
    // 	    // if(abs(xastaa[0][0]) > 1 || abs(xastaa[1][0]) > 1 || abs(xastaa[2][0]) > 1)
    // 	    //   {
    // 	    // 	// xastaa[0][0] = xstart[0];
    // 	    // 	// xastaa[1][0] = xstart[1];
    // 	    // 	// xastaa[2][0] = xstart[2];
    // 	    // 	break;
    // 	    //   }
    // 	    // else
    // 	    //   {
    // 	    // 	xnew[0] = xastaa[0][0];	    
    // 	    // 	xnew[1] = xastaa[1][0];	    
    // 	    // 	xnew[2] = xastaa[2][0];	    
    // 	    //   }

    // 	    //	    cout<<"\n out of while\n\n";
    // 	    //replace fnval with fnvalnew?
    // 	    E->PhysEvalBasisGrad(xastaa, tempeval[0], tempeval[1], tempeval[2], tempeval[3]);
    // 	    if(flag == 0) // from opt_needed  
    // 	      {
    // 		for(int p = 0; p < dim; p++)
    // 		  {
		    
    // 		    Vmath::Vmul(uhats.size(), &tempeval[p+1][0], 1, &uhats[0], 1, &temp1[0], 1 );
    // 		    gdnew[p] = Vmath::Vsum(temp1.size(), temp1, 1);
    // 		  }
		
    // 	      }
    // 	    else
    // 	      {
		
    // 		for(int p = 0; p < dim; p++)
    // 		  {
    // 		    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
		    
    // 		    gdnew[p] = dereval[p];
		    
    // 		  }
		
        	
    // 	      }
	
    // 	    g = gdnew;
	    
    // 	    //L2 norm
    // 	    gnew = pow((pow(g[0],2) + pow(g[1],2) + pow(g[2],2)) ,0.5);
    // 	    cout<<"\n gnew = "<<gnew<<" n = "<<n<<" xsave[0]="<<xsave[0]<<" xsave[1]="<<xsave[1]<<" xsave[2]"<<xsave[2]<<"\n";
    // 	    if(gnew <1e-6||(xsave[0] - xastaa[0][0] <1e-6 && xsave[1] - xastaa[1][0] <1e-6 && xsave[2] - xastaa[2][0] <1e-6))
    // 	      {
    // 		ret[0].push_back(xastaa[0][0]);
    // 		ret[1].push_back(xastaa[1][0]);
    // 		ret[2].push_back(xastaa[2][0]);
    // 		cout<<"\n ret = "<<ret[0][0]<<" "<<ret[1][0]<<" "<<ret[2][0];
    // 		return ret;
    // 	      }
    // 	    //	    xsave = xnew;
    // 	    xnew[0] = xastaa[0][0];
    // 	    xnew[1] = xastaa[1][0];
    // 	    xnew[2] = xastaa[2][0];
    // 	    fnval = fnvalnew;
    // 	    n = n+1;
    // 	  }

    // 	cout<<"\n ret = "<<xstart[0]<<" "<<xstart[1]<<" "<<xstart[2]<<"\n";
    // 	ret[0].push_back(xstart[0]);
    // 	ret[2].push_back(xstart[2]);
    // 	ret[1].push_back(xstart[1]);
    // 	avgiterGD = n;
    // 	cout<<" n = "<<n<<" ";

    //   }
    // return ret;

	// 	//stepsize = gam;
	




    // // // find g = grad(f)
    // // cout<<"\n calc g inside:"<<xastaa[0][0]<<" "<<xastaa[1][0]<<" "<<xastaa[2][0];
    // 	// E->PhysEvalBasisGrad(xastaa, tempeval[0], tempeval[1], tempeval[2], tempeval[3]);
    // 	// if(flag == 0) // from opt_needed
    // 	//   {
    // 	//     fnval = 0;
    // 	//     for(int p = 0; p < dim; p++)
    // 	//       {
                        
    // 	// 	Vmath::Vmul(uhats.size(), &tempeval[p+1][0], 1, &uhats[0], 1, &temp1[0], 1 );
    // 	// 	g[p] = Vmath::Vsum(temp1.size(), temp1, 1);
    // 	//       }
    // 	//     Vmath::Vmul(uhats.size(), &tempeval[0][0], 1, &uhats[0], 1, &temp1[0], 1 );
    // 	//     fnval = Vmath::Vsum(temp1.size(), temp1, 1);
		  
    // 	//   }
    // 	// else // from optimizer
    // 	//   {
    // 	//     fnval = 0;     
    // 	//     for(int p = 0; p < dim; p++)
    // 	//       {
    // 	// 	derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
                        
    // 	// 	g[p] = dereval[p];
    // 	// 	//		    epsl += abs(g[p]);
		    
    // 	//       }
    // 	//     //tm.Stop();
    // 	//     //pqdertime += tm.TimePerTest(1);
    // 	//     //		timer.Stop();
    // 	//     Array<OneD, NekDouble> t1(1);
    // 	//     pq(uhats, xastaa, tempeval[0], nullarr, t1);
    // 	//     fnval = t1[0];
    // 	//   }

    // 	gnew = pow((pow(g[0],2) + pow(g[1],2) + pow(g[2],2)) ,0.5); 
    // 	cout<<"\n beg gnew = "<<gnew<<" ";
    // 	NekDouble a = 1;
    // 	while(gnew > 1e-6 && n < iter)
    // 	  {

    // 	    int flg = 1;
    // 	    for(int k = 0; k < dim; k ++)
    // 	      {
    // 		NekDouble val = xnew[k] - a*(g[k]);
    // 		if(abs(val) > 1.0 )
    // 		  {
    // 		    flg = 0;
    // 		    //		    cout<<"\n2 n = "<<n<<"\n\n";
    // 		        break;
    // 		  }
    // 	      }
    // 	    if(flg == 1)
    // 	      {
    // 		// find x - ag
    // 		xastaa[0][0] = xnew[0] - a*(g[0]);
    // 		xastaa[1][0] = xnew[1] - a*(g[1]);
    // 		xastaa[2][0] = xnew[2] - a*(g[2]);
    // 		a = 1;
    // 	      }
    // 	    else
    // 	      {
    // 		//a = a*beta;
    // 		xastaa[0][0] = xnew[0];
    // 		xastaa[1][0] = xnew[1];
    // 		xastaa[2][0] = xnew[2];
    // 		a = a*beta;
    // 		continue;
    // 		//break;
    // 	      }

    // 	    // find func(x)
    // 	    E->PhysEvalBasisGrad(xastaa, tempeval[0], NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
    // 	    	    cout<<"\n	    xastaa[0][0] = "<<xastaa[0][0]<<"   xastaa[1][0] = "<<xastaa[1][0]<<"   xastaa[2][0] = "<<xastaa[2][0] << " xnew[0] = "<<xnew[0]<<" xnew[1] = "<<xnew[1]<<" xnew[2] = "<<xnew[2]<<" g[0] = "<<g[0]<<" g[1] = "<<g[1]<<" g[2] = "<<g[2]<<" ";
    // 	    if(flag == 0) // from opt_needed  
    // 	      {
    // 		fnvalnew = 0;
    // 		Vmath::Vmul(uhats.size(), &tempeval[0][0], 1, &uhats[0], 1, &temp1[0], 1 );
    // 		fnvalnew = Vmath::Vsum(temp1.size(), temp1, 1);
			
    // 	      }
    // 	    else
    // 	      {
    // 		fnvalnew = 0;
    // 		Array<OneD, NekDouble> t1(1);
    // 		pq(uhats, xastaa, tempeval[0], nullarr, t1);
    // 		fnvalnew = t1[0];

    // 	      }
    // 	    // if(flg == 0)
    // 	    //   cout<<"\n here\n\n";
	
    // 	    //	    cout<<"\nfnvalnew ="<<fnvalnew<<"fval="<<fnval<<"( fnvalnew - fnval)/a="<<( fnvalnew - fnval)/a<<" -sigma*(pow(g[0],2) + pow(g[1],2) + pow(g[2],2))="<<-sigma*(pow(g[0],2) + pow(g[1],2) + pow(g[2],2))<<"\n";;//)/a >= -sigma*(pow(g[0],2) + pow(g[1],2) + pow(g[2],2));

    // 	    cout<<"\nxastaa[0][0] = "<<xastaa[0][0]<<" xastaa[1][0]="<<xastaa[1][0]<<" xastaa[2][0] ="<< xastaa[2][0]<<" a="<<a<<" n = "<<n<< "gnew = "<<gnew<<" ";
    // 	   //	   int flg1 = 1;
    // 	    a =1;
    // 	   while(( fnvalnew - fnval)/a >= -sigma*(pow(g[0],2) + pow(g[1],2) + pow(g[2],2)))
    // 	      {

    // 		a = a*beta;
	        
    // 		//		if(fnvalnew - fnval <1e-10)
    // 		//  break;
    // 		// find x - ag after changing a
    // 		xastaa[0][0] = xnew[0] - a*(g[0]);
    // 		xastaa[1][0] = xnew[1] - a*(g[1]);
    // 		xastaa[2][0] = xnew[2] - a*(g[2]);

    // 		if(abs( xastaa[0][0]) > 1.0 || abs( xastaa[1][0]) > 1.0 || abs( xastaa[2][0])>1.0)
    // 		  {
    // 		    continue;
    // 		    // xastaa[0][0] = xstart[0];
    // 		    // xastaa[1][0] = xstart[1];
    // 		    // xastaa[2][0] = xstart[2];
		   
    // 		    // a = a*beta;
    // 		    // xnew[0] = xstart[0];
    // 		    // xnew[2] = xstart[2];
    // 		    // xnew[1] = xstart[1];
    // 		      //		    flg1 = 0;
    // 		    //		    break;    

    // 		  }
    // 		// find func(x)
    // 		E->PhysEvalBasisGrad(xastaa, tempeval[0], NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
    // 		if(flag == 0) // from opt_needed  
    // 		  {
    // 		    fnvalnew = 0;
    // 		    Vmath::Vmul(uhats.size(), &tempeval[0][0], 1, &uhats[0], 1, &temp1[0], 1 );
    // 		    fnvalnew = Vmath::Vsum(temp1.size(), temp1, 1);
		    
    // 		  }
    // 		else
    // 		  {
    // 		    fnvalnew = 0;
    // 		    Array<OneD, NekDouble> t1(1);
    // 		    pq(uhats, xastaa, tempeval[0], nullarr, t1);
    // 		    fnvalnew = t1[0];
		    
    // 		  }
    // 		cout<<"\n (fnvalnew - fnval)/a="<<(fnvalnew - fnval)/a<<" -sigma*(pow(g[0],2) + pow(g[1],2) + pow(g[2],2)) = "<<-sigma*(pow(g[0],2) + pow(g[1],2) + pow(g[2],2))<<" a = "<<a<< " fnval = "<<fnval<<" fnvalnew = "<<fnvalnew<<" n = "<<n<<" fnvalnew="<<fnvalnew<< "xastaa[0][0]="<<xastaa[0][0]<<" xastaa[1][0]="<<xastaa[1][0]<<" xastaa[2][0] = "<<xastaa[2][0];;
		
		
		
    // 	      }
    // 	   // if(flg1 ==0)
    // 	   //   {
    // 	   //     cout<<"\n here\n\n";
    // 	   //   }
    // 	    xnew[0] = xnew[0] - a*(g[0]);
    // 	    xnew[1] = xnew[1] - a*(g[1]);
    // 	    xnew[2] = xnew[2] - a*(g[2]);


    // 	      fnval = fnvalnew;

    // 	      // find grad(f at new x)
    // 	      xastaa[0][0] = xnew[0];
    // 	      xastaa[1][0] = xnew[1];
    // 	      xastaa[2][0] = xnew[2];
	      	        
    // 	      E->PhysEvalBasisGrad(xastaa, tempeval[0], tempeval[1], tempeval[2], tempeval[3]);
    // 	      if(flag == 0) // from opt_needed
    // 		{
    // 		  for(int p = 0; p < dim; p++)
    // 		    {
		    
    // 		      Vmath::Vmul(uhats.size(), &tempeval[p+1][0], 1, &uhats[0], 1, &temp1[0], 1 );
    // 		      g[p] = Vmath::Vsum(temp1.size(), temp1, 1);
    // 		    }
		
    // 		}

    // 	      else // from optimizer
    // 		{
    // 		  for(int p = 0; p < dim; p++)
    // 		    {
    // 		      derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
                        
    // 		      g[p] = dereval[p];
		    
    // 		    }
    // 		}

    // 	      //L2 norm
    // 	      gnew = pow((pow(g[0],2) + pow(g[1],2) + pow(g[2],2)) ,0.5);
    // 	      cout<<"\n gnew = "<<gnew<<" n = "<<n<<"\n";
    // 	      n  = n + 1;
	      
	      
	      
    // 	  }
  //  }
	      
  
  //if flag = 0, call if not from optimizer
  //if flag = 1, call is from optimizer
    vector< vector< NekDouble> >  gradient_descent3D(Array<OneD, NekDouble> &uhats, StdExpansion *E, NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> > &testcoord3d, Array<OneD, NekDouble> &interioreval3d, int flag)
  {
    // LibUtilities::Timer timer;
    int dim = 3;
    //    NekDouble tmptime =0.0;
    //timer.Start();
    NekDouble stepsize = iterGD;
    Array<OneD, NekDouble> g(dim);
    double inf = numeric_limits<double>::infinity();
        
        
    Array<OneD,NekDouble> xnew(dim) ,x(dim);
    NekDouble gprev = inf;
    int idxgprev;
    vector< vector< NekDouble> > ret(dim);
    int sz = testcoord3d[0].size();
    Array<OneD, Array<OneD, NekDouble > > tempeval(4);

    for(int i = 0; i < 4; i++)
      {
	tempeval[i] = Array<OneD, NekDouble>(E->GetNcoeffs());  
      }
        
    Array<OneD, NekDouble> nullarr(0);
    Array<OneD, NekDouble> temp2(sz);
        
    if(flag==0)
      {
            
	Array<OneD, NekDouble> temp(uhats.size());
            
	for(int k = 0; k < sz; k++)
	  {
                
	    Vmath::Vmul(uhats.size(), &interioreval3d[k], sz, &uhats[0], 1, &temp[0], 1);          
	    temp2[k] = Vmath::Vsum(temp.size(), temp, 1);
	  }
            	
gprev = Vmath::Vmin(temp2.size(), temp2, 1);	
	idxgprev = Vmath::Imin(sz, temp2, 1);        
            
      }
    else
      {
	pq(uhats, testcoord3d, interioreval3d, nullarr, temp2 );
	gprev = Vmath::Vmin(sz, temp2, 1);
	idxgprev = Vmath::Imin(sz, temp2, 1);
      }
        
    xnew[0] = testcoord3d[0][idxgprev];
    xnew[1] = testcoord3d[1][idxgprev];
    xnew[2] = testcoord3d[2][idxgprev];
    // if(xnew[0]+xnew[2] >0 || xnew[1]+xnew[2] >0)
    //   {
    // 	cout<<"\n start: ("<< xnew[0] <<" ,"<< xnew[1]<<" ,"<< xnew[2]<<") gprev="<<gprev<<"\n";  
    //   }
    //timer.Stop();
    //    cout<<"\n xnew="<<xnew[0]<<" "<<xnew[1]<<" "<<xnew[2]<<"  setup phase of GD 3D took "<<timer.TimePerTest(1)<<"s gprev="<<gprev<<"\n";
    int n  = 0;
    // prevepsl = 0,
    NekDouble saveminx = xnew[0], saveminy = xnew[1], saveminz = xnew[2], savefnval = inf, fnval;   
    // timer.Stop();
    // tmptime = timer.TimePerTest(1);
    // cout<<"\n setup for GD 3d: "<<tmptime<<"gprev = "<<gprev<<"\n";
    // tmptime = 0.0;
    // NekDouble pqdertime = 0.0;
    // timer.Start(); 
    
    if(gprev < 0 && abs(gprev)>1e-13)
      {
	n = 0;
	NekDouble epsl = 1.0;

	NekDouble  iter = secarg;
	Array<OneD, Array<OneD, NekDouble> > xastaa(dim);

	Array<OneD, NekDouble> dereval(dim);
	Array<OneD, NekDouble> xyzstart(xnew.size());
	for(int k = 0; k < dim; k++)
	  xyzstart[k] = xnew[k];
	bool res2, res3 = true;
            
	Array<OneD, NekDouble> temp1 (uhats.size());
	
	while(n < iter)
	  {
	    //	    prevepsl = epsl;  
	    epsl = 0.0;
	    for(int p = 0; p < dim; p++)
	      {
		xastaa[p] = Array<OneD, NekDouble>(1,xnew[p]);
	      }   
            //Timer tm;
	    //tm.Start();
	        
	    E->PhysEvalBasisGrad(xastaa, tempeval[0], tempeval[1], tempeval[2], tempeval[3]);
	    //tm.Stop();  
	    //	    NekDouble disp = tm.TimePerTest(1);
	    //	    cout<<"\n at n = "<<n<<" eval time = "<<disp;
	    //	    pqdertime += disp;

	    // multiply tempeval by uhats to get g[p]
	                

	    if(flag == 0) // from opt_needed
	      {
		fnval = 0;
		for(int p = 0; p < dim; p++)
		  {
                        
		    Vmath::Vmul(uhats.size(), &tempeval[p+1][0], 1, &uhats[0], 1, &temp1[0], 1 );
		    g[p] = Vmath::Vsum(temp1.size(), temp1, 1);
                    
		    epsl += abs(g[p]);
		  }
		Vmath::Vmul(uhats.size(), &tempeval[0][0], 1, &uhats[0], 1, &temp1[0], 1 );
		fnval = Vmath::Vsum(temp1.size(), temp1, 1);
		  
	      }
	    else
	      {
		//		Timer tm;
		fnval = 0;     
		//timer.Start();
		//tm.Start();
		for(int p = 0; p < dim; p++)
		  {
		    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
                        
		    g[p] = dereval[p];
		    epsl += abs(g[p]);
		    
		  }
		//tm.Stop();
		//pqdertime += tm.TimePerTest(1);
		//		timer.Stop();
		Array<OneD, NekDouble> t1(1);
		pq(uhats, xastaa, tempeval[0], nullarr, t1);
		fnval = t1[0];
	      }
	    
	    bool r1 = true;
	    // if(E->DetShapeType() == LibUtilities::ePrism)
	    //   r1 = xastaa[0][0]+ xastaa[2][0] <= 0;

	    // if(E->DetShapeType() == LibUtilities::eTetrahedron)
	    //   r1 = xastaa[0][0]+ xastaa[2][0] <= 0 &&   xastaa[0][0]+ xastaa[1][0] <= 0 &&  xastaa[2][0]+ xastaa[1][0] ;
	    if(r1  && savefnval > fnval && abs(xastaa[2][0]) <= 1 && abs(xastaa[0][0] )<=1 && abs(xastaa[1][0])<=1)
	      {
		savefnval = fnval;
		saveminx = xastaa[0][0];
		saveminy = xastaa[1][0];
		saveminz = xastaa[2][0];
	      }
	    
	    x[0] = xnew[0] - stepsize*(g[0]);
	    x[1] = xnew[1] - stepsize*(g[1]);
	    x[2] = xnew[2] - stepsize*(g[2]);
	    res2 = true; res3 = true;
	    // if(E->DetShapeType() == LibUtilities::ePrism)
	    //   {
	    // 	res3 = (x[0]+x[2])<=0;
	    //   }
	    // if(E->DetShapeType() == LibUtilities::eTetrahedron)
	    //   {
	    // 	res3 = (x[0]+x[2])<=0 && (x[0]+x[1])<=0 && (x[1]+x[2])<=0;
	    //   }
	    // else if(E->DetShapeType() == LibUtilities::ePyramid)
	    //   {
            //       res3 = (x[0]+x[2])<=0 &&  (x[1]+x[2])<=0;;

	    //   }
	    
	    for(int k = 0; k < dim; k++)
	      {   
                    
		if(abs(x[k]-xnew[k]) > 1e-4 || abs(x[k])>1)
		  {
		    res2 = res2 && false  ;
		  }
	      }

            if(res2 && (res3)  )    
	      {
	    //   saveminx = x[0];
	    //   saveminy = x[1];
	    //   saveminz = x[2];
	       break;
             }
	    //if(prism)
	    //{
	    //check if x+z <= 0
	    //  res3 = (x[0]+x[2])<=0 ;   

	    //if(pyr)
	    //{
	    //check if x+z <= 0 and y+z<=0
	    //  res3 = (x[0]+x[2])<=0 ;   
	    // if( ((std::isnan(epsl)) || (epsl <1e-10) || (epsl>1e2) || abs(prevepsl - epsl)<1e-5 || (abs(x[0]-xnew[0]) < 1e-5 && abs(x[1]-xnew[1]) < 1e-5 && abs(x[2]-xnew[2]) < 1e-5)  ) && abs(x[0]) <=1.0 && abs(x[1])<=1.0 &&  abs(x[2])<=1.0  )
	    //   {
	    // 	//		cout<<"\n break at n = "<<n<<"\n";
	    // 	break;
	    //   }

	    if((abs(x[0]) >1.0) || (abs(x[1])>1.0) || (abs(x[2])>1.0) || (std::isnan(epsl)) )
	      {
		stepsize = stepsize/10;
		xnew[0] = xyzstart[0];
		xnew[1] = xyzstart[1];
		xnew[2] = xyzstart[2];

                if(stepsize <1e-4)
		  {
		    // xnew[0] = saveminx;
		    // xnew[1] = saveminy;
		    // xnew[2] = saveminz;

		    break;
		  }
		continue;
	      }

	    
	    xnew[0] = x[0];
	    xnew[1] = x[1];
	    xnew[2] = x[2];
            
	    n = n+1;
	  }
	//	timer.Stop();
	//tmptime =  timer.TimePerTest(1);
		
	//	cout<<"\n time taken by GD iters = "<<tmptime<<" time taken by derpq() = "<< pqdertime<<" n = "<<n<<"\n";;

	//	timer.Stop();
	ret[0].push_back(saveminx);
	ret[1].push_back(saveminy);
	ret[2].push_back(saveminz);
	avgiterGD = n;
      }
    // if(ret[0].size() >0 && ret[1].size()>0 && ret[2].size() >0)
    //   cout<<"\n  3D GD iters took "<<timer.TimePerTest(1)<<"s and "<<n<<"  iters  << returning = "<<ret[0][0]<<" "<<ret[1][0]<<" "<<ret[2][0]<<" savefnval= "<<savefnval<<" fnval = "<<fnval<<" stepsize = "<<stepsize<<" ";

    return ret;
  }
    
    
  //rename to findrootsedges
  vector<vector<NekDouble> > find_roots(Array<OneD, NekDouble> &uhats, StdExpansion *E, NekDouble &avgiterGD, int d ,int sig, int surfflag, int volflag, int surfid)
  {
    int dimension; 
    if(surfflag == 0 && volflag == 0)
      dimension = 1;
    else if(volflag == 0)
      dimension = 2;
    else
      dimension = 3;
    vector<vector< NekDouble> > ret(dimension);
    //Confederate matrix approach
    if(surfflag == 0 && volflag == 0)
      {
        int N = uhats.size();

        vector<NekDouble> uhatsmon;
        while(abs(uhats[N-1])<1e-8)
	  {
            N = N-1;
            
            if(N == 0)
	      {
        
                ret[0].push_back(-1.0);
                ret[0].push_back(1.0);
                
                return ret;
	      }

	  }
	
        vector<NekDouble> temp(N);
        // convert uhats to monomial, find roots of uhats or der of uhats
        //        cout<<"\n uhatsmon=\n";
        for(int k = 0; k < N; k++)
	  {
            for(int jj= 0 ; jj < N; jj++)
	      {
                temp[jj ] = C[jj][k];
	      }
            Vmath::Vmul(N, &temp[0], 1,  &uhats[0], 1, &temp[0], 1);
            NekDouble temp2 = Vmath::Vsum(N, &temp[0], 1);
            uhatsmon.push_back(temp2);
            //  cout<< temp2<< " ";
	  }

        if(abs(Vmath::Vmax(uhatsmon.size(), &uhatsmon[0], 1))<1e-10)
	  {
            ret[0].push_back(-1.0);
            ret[0].push_back(1.0);
            
            return ret;
         
	  }
	
        // truncate trailing zeros
        while(abs(uhatsmon[N-1])<1e-8)
	  {
            N = N-1;
            if(N == 0)
	      {
       
                  ret[0].push_back(-1.0);
                  ret[0].push_back(1.0);
                
                return ret;
	      }

	  }
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
                  ret[0].push_back(-1.0);
                  ret[0].push_back(1.0);
                
                  return ret;
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
            ret[0].push_back( EIG_R[kk] );
            
	  }


      }
    else if(dimension > 1 && surfflag == 1 && volflag == 0 )
      {
	
	if(surfid == 0)
	  ret = gradient_descent2Dquad(uhats, E, avgiterGD, testcoord2dquad, interioreval2dquad, sig);
	else
	  ret = gradient_descent2Dtri(uhats, E, avgiterGD, testcoord2dtri, interioreval2dtri, sig);
	
      }
    else if(volflag == 1)
      {
	//	ret = gradient_descent3D(uhats, E, avgiterGD, testcoord3d, interioreval3d, sig);

	ret = steepestgradientdescent3D(uhats, E, avgiterGD, testcoord3d, interioreval3d, sig);

      }
    
    return ret;    
    
  }
 

  Array<OneD, Array<OneD, NekDouble>> GetCoords(StdExpansion *E)
  {
    int dimension = E->GetShapeDimension();
    const auto totPoints = (unsigned) E->GetTotPoints();
    Array<OneD, NekDouble> x(totPoints), y(totPoints), z(totPoints);
    Array<OneD, Array<OneD, NekDouble> > coords(dimension);

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
  NekDouble eps = 1e-3;
  int avgnum = 1;
  NekDouble chold = 0.3;
  NekDouble gamhold = 0.5;
  std::string    m_shape;
  std::string    m_ntype;
  vector<string> m_basis{3, "NoBasisType"};
  vector<string> m_pointstype{3, "NoPointsType"};
  vector<int>    m_order;
  vector<int>    m_points;

};

#endif
