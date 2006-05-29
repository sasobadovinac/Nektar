///////////////////////////////////////////////////////////////////////////////
//
// File $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/LocalRegions/TriExp.h,v $ 
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef TRIEXP_H
#define TRIEXP_H

#include <StdRegions/StdBasis.h>
#include <StdRegions/StdTriExp.h>

#include <SpatialDomains/TriGeom.h>

#include <LocalRegions/MetricRelatedInfo.h>

#include <StdRegions/StdRegions.hpp>

namespace Nektar
{
    namespace LocalRegions 
    {

	class TriExp: public StdRegions::StdTriExp
	{

	public:
	    
	    /** \brief Constructor using BasisKey class for quadrature
		points and order definition */
	    TriExp(const StdRegions::BasisKey &Ba, const StdRegions::BasisKey &Bb,
		   SpatialDomains::SharedTriGeomPtr geom);

	    
	    /** \brief Constructor using BasisKey class for quadrature
		points and order definition where _coeffs and _phys are all
		set. */
	    TriExp(const StdRegions::BasisKey &Ba, const StdRegions::BasisKey &Bb,
		   double *coeffs, double *phys, SpatialDomains::SharedTriGeomPtr geom);
	    
	    /// Copy Constructor
	    TriExp(const TriExp &T);

	    /// Destructor
	    ~TriExp();
    
	    /// Return Shape of region, using  ShapeType enum list. i.e. Triangle
	    StdRegions::ShapeType DetShapeType() 
	    {
		return StdRegions::eTriangle;
	    }
	    
	    SharedMetricRelatedInfoPtr TriExp::GenGeoFac();
	    
	    
	    inline void SetGeoFac(SharedMetricRelatedInfoPtr minfo)
	     {
		 m_minfo = minfo;
	     }


	    void GetCoords(double **coords);
      

	    void GetCoord(const double *Lcoords, double *coords);
    

	    virtual StdRegions::GeomType V_GeoFacType()
	    {
		return m_minfo->GetGtype();
	    }
	    
	    void WriteToFile(FILE *outfile);
	    
	    //----------------------------
	    // Integration Methods
	    //----------------------------
	    
	    /// \brief Integrate the physical point list \a inarray over region
	    double Integral(const double *inarray);

	    /** \brief  Inner product of \a inarray over region with respect to the 
		expansion basis (this)->_Base[0] and return in \a outarray */
	    void IProductWRTBase(const double * inarray, double * outarray);

	    //----------------------------------
	    // Local Matrix Routines 
	    //----------------------------------
	    
	    /** \brief Get the mass matrix attached to this expansion by using
		the StdMatrix manager _ElmtMats and return the standard Matrix
		container */
	    StdRegions::StdMatContainer *GetMassMatrix(); 
      
	    /** \brief Get the weak Laplacian matrix attached to this
		expansion by using the StdMatrix manager _ElmtMats and return
		the standard Matrix container */
	    StdRegions::StdMatContainer * GetLapMatrix();
	    
	    //-----------------------------
	    // Differentiation Methods
	    //-----------------------------
	    
	    void Deriv(double * outarray_d1, double *outarray_d2)
	    {
		double *out[2];
		out [0] = outarray_d1;  
		out [1] = outarray_d2;
		Deriv(2, this->m_phys, out);
	    }
	    
	    void Deriv(const double *inarray, double * outarray_d1, 
		       double *outarray_d2)
	    {
		double *out[2];
		out [0] = outarray_d1;  
		out [1] = outarray_d2;
		Deriv(2, inarray, out);
	    }
	    
	    void Deriv(const int n, double **outarray);
	    
	    void Deriv(const int n, const double *inarray, double ** outarray);
	    
	    //----------------------------
	    // Evaluations Methods
	    //---------------------------
	    
	    /** \brief Forward transform from physical quadrature space
		stored in \a inarray and evaluate the expansion coefficients and
		store in \a (this)->_coeffs  */
	    void FwdTrans(const double * inarray);

	    double Evaluate(const double *coord);
	    
	protected:
	    int m_id;
	    int m_field;
	    
	    SpatialDomains::SharedTriGeomPtr m_geom;
	    SharedMetricRelatedInfoPtr       m_minfo;
      
	    /** \brief  Inner product of \a inarray over region with respect to
		the expansion basis \a base and return in \a outarray */
	    inline void IProductWRTBase(const double *base0, const double *base1, 
					const double *inarray, double *outarray);

	private:

	    virtual StdRegions::ShapeType v_DetShapeType() 
	    {
		DetShapeType();
	    }

	    virtual SharedMetricRelatedInfoPtr v_GenGeoFac()
	    {
		return GenGeoFac();
	    }

	    virtual void v_SetGeoFac(SharedMetricRelatedInfoPtr minfo)
	    {
		SetGeoFac(minfo);
	    }

	    virtual void v_GetCoords(double **coords)
	    {
		GetCoords(coords);
	    }

	    virtual void v_GetCoord(const double *Lcoords, double *coords)
	    {
		GetCoord(Lcoords, coords);
	    }

	    virtual void v_WriteToFile(FILE *outfile)
	    {
		WriteToFile(outfile);
	    }

	    /** \brief Virtual call to integrate the physical point list \a inarray
		over region (see SegExp::Integral) */
	    virtual double v_Integral(const double *inarray)
	    {
		return Integral(inarray);
	    }

	    /** \brief Virtual call to TriExp::IProduct_WRT_B */
	    virtual void v_IProductWRTBase(const double * inarray, double * outarray)
	    {
		IProductWRTBase(inarray,outarray);
	    }
	    
	    /// virtual call to GetMassMatrix
	    virtual StdRegions::StdMatContainer *v_GetMassMatrix() 
	    {
		return GetMassMatrix();
	    }
	    
	    /// virtual call to GetLapatrix
	    virtual StdRegions::StdMatContainer *v_GetLapMatrix() 
	    {
		return GetLapMatrix();
	    }
      
	    virtual void v_Deriv(double * outarray_d1, double *outarray_d2)
	    {
		Deriv(this->m_phys, outarray_d1, outarray_d2);
	    }

	    virtual void v_StdDeriv(double * outarray_d1, double *outarray_d2)
	    {
		StdTriExp::Deriv(this->m_phys, outarray_d1, outarray_d2);
	    }
    
	    virtual void v_Deriv(const double *inarray, double * outarray_d1, 
				 double *outarray_d2)
	    {
		Deriv(inarray, outarray_d1, outarray_d2);
	    }

	    virtual void v_StdDeriv(const double *inarray, double * outarray_d1, 
				    double *outarray_d2)
	    {
		StdTriExp::Deriv(inarray, outarray_d1, outarray_d2);
	    }
  
	    virtual void v_Deriv(const int n,  double ** outarray)
	    {
		Deriv(n, outarray);
	    }

	    virtual void v_Deriv(const int n, const double *inarray,
				 double ** outarray)
	    {
		Deriv(n, inarray, outarray);
	    }
    
	    /// Virtual call to SegExp::FwdTrans
	    virtual void v_FwdTrans(const double * inarray)
	    {
		FwdTrans(inarray);
	    }
	    
	    /// Virtual call to TriExp::Evaluate
	    virtual double v_Evaluate(const double * coords)
	    {
		return Evaluate(coords);
	    }
	};
	
	
    } //end of namespace
} //end of namespace

#endif // TRIEXP_H

/** 
 *    $Log: TriExp.h,v $
 *    Revision 1.1  2006/05/04 18:58:47  kirby
 *    *** empty log message ***
 *
 *    Revision 1.13  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.12  2006/03/12 07:43:33  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
