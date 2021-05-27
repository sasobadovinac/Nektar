///////////////////////////////////////////////////////////////////////////////
//
// File FilterEnergy.h
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
// Description: Outputs solution fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////   
#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERSTRUCTUREPRES_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERSTRUCTUREPRES_H

#include <SolverUtils/Filters/Filter.h>
#include <SolverUtils/Filters/StdDemoSupport.hpp>
/* #include <LibUtilities/BasicUtils/NekFactory.hpp> */
/* #include <LibUtilities/BasicUtils/SessionReader.h> */

namespace Nektar
{
  
  namespace SolverUtils
  {
    class FilterStructurePres : public Filter
    {
    public:
      friend class MemoryManager<FilterStructurePres>;
      /// Creates an instance of this class
      static FilterSharedPtr create(
				    const LibUtilities::SessionReaderSharedPtr         &pSession,
				    const std::weak_ptr<EquationSystem> &pEquation,
				    const std::map<std::string, std::string>   &pParams)
      {
	FilterSharedPtr p = MemoryManager<FilterStructurePres>
	  ::AllocateSharedPtr(pSession, pEquation, pParams);
	return p;
      }

      ///Name of the class
      static std::string className;

      SOLVER_UTILS_EXPORT FilterStructurePres(
					      const LibUtilities::SessionReaderSharedPtr &pSession,
					      const std::weak_ptr<EquationSystem>      &pEquation,
					      const ParamMap                             &pParams);
      SOLVER_UTILS_EXPORT  ~FilterStructurePres();
      

    protected:

      Array<OneD, Array<OneD, NekDouble> > FormConf(NekDouble N);
      
      SOLVER_UTILS_EXPORT virtual void v_Initialise(
						    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pField,
						    const NekDouble &time);
      SOLVER_UTILS_EXPORT virtual void v_Update(
						const Array<OneD, const MultiRegions::ExpListSharedPtr> &pField,
						const NekDouble &time);
      SOLVER_UTILS_EXPORT virtual void v_Finalise(
						  const Array<OneD, const MultiRegions::ExpListSharedPtr> &pField, const NekDouble &time);
      SOLVER_UTILS_EXPORT virtual bool v_IsTimeDependent();

      Array<OneD, Array<OneD, NekDouble> > FindIds(
						   const MultiRegions::ExpListSharedPtr &pFields);

      void  Optimize(LocalRegions::ExpansionSharedPtr exp, Array<OneD, NekDouble> &coefflocal, Array<OneD, NekDouble> &physlocal, Array<OneD, NekDouble> tmpcoord);

      Array<OneD, Array<OneD, NekDouble> > formConf(NekDouble N); 

      void project_edges( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret, StdExpansion *E);
      void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret, int modes,  Array <OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2);
	
      Array<OneD, NekDouble>  call_find_roots(Array<OneD,  NekDouble> &uhats , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats, NekDouble &minv, StdExpansion *E, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime);


      void surfuhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz, StdExpansion *Etemp);

      void project_surfaces( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret, StdExpansion *E);

      void call_setup_seg(LocalRegions::ExpansionSharedPtr exp, int N);
      void call_setup_tri(LocalRegions::ExpansionSharedPtr exp);
      void call_setup_quad(LocalRegions::ExpansionSharedPtr exp);
      void call_setup_hex(LocalRegions::ExpansionSharedPtr exp);
      void call_setup_tet(LocalRegions::ExpansionSharedPtr exp);
      void call_setup_pyr(LocalRegions::ExpansionSharedPtr exp);
      void DoOrthonormalize(const MultiRegions::ExpListSharedPtr &pFields,
			    int flag = 0);
    private:
      int dimension;
      unsigned int                m_index;
      unsigned int                m_outputFrequency;
      std::ofstream               m_outFile;
      bool                        m_homogeneous;
      NekDouble                   m_homogeneousLength;
      NekDouble                   m_area;
      LibUtilities::CommSharedPtr m_comm;
      Array<OneD, unsigned int>   m_planes;

      
      // Structure pres filter
      //Orthogonal ver of elements 
      StdExpansion *E3seg = nullptr; 
      StdExpansion *Etet = nullptr;
      StdExpansion *Ehex = nullptr;
      StdExpansion *Epyr = nullptr;
      StdExpansion *Equad = nullptr;
      StdExpansion *Etri = nullptr;

      Array<OneD, Array<OneD, NekDouble> > C;
      
      Array<OneD, const MultiRegions::ExpListSharedPtr> saveFields;
      /* Array<OneD, Array<OneD, NekDouble > > storage3dhex, storage3dtet, storage3dpyr, storage2dt, storage2dq, coordmidhex, coordmidtet, coordmidquad,  coordmidtri, coordmidpyr, coordhex, coordtet, coordpyr, coordtri, coordquad; */
      
      /* Array<OneD, NekDouble> midptevalhex, midptevaltet, midptevalpyr, midptevalquad, midptevaltri;  */
      Array<OneD, Array<OneD, NekDouble > > x, y, z;// coords of glo domain (all elements)
      
      DemoSupport demo;
      int nelmt;
      int nfields;

      // for edges root finding quad:
      // left (x = -1)
      Array<OneD, NekDouble> Vxm1q;
      Array<OneD, NekDouble> Vdyxm1q;// =(totszedges);
      Array<OneD, NekDouble> Vdxxm1q;// = Array<OneD, NekDouble>(totszedges);

      // bot (y = -1)
      Array<OneD, NekDouble> Vym1q;//  = Array<OneD, NekDouble>(totszedges);
      Array<OneD, NekDouble> Vdxym1q;//  = Array<OneD, NekDouble>(totszedges);
      Array<OneD, NekDouble> Vdyym1q;//  = Array<OneD, NekDouble>(totszedges);
      
      // right quad (x = 1)
      Array<OneD, NekDouble> Vx1q;// = Array<OneD, NekDouble>(totszedges);
      Array<OneD, NekDouble> Vdxx1q;// = Array<OneD, NekDouble>(totszedges);
      Array<OneD, NekDouble> Vdyx1q;// = Array<OneD, NekDouble>(totszedges);

      // top quad (y = 1)
      Array<OneD, NekDouble> Vdxy1q;// = Array<OneD, NekDouble>(totszedges);
      Array<OneD, NekDouble> Vdyy1q;// = Array<OneD, NekDouble>(totszedges);
      Array<OneD, NekDouble> Vy1q;// = Array<OneD, NekDouble>(totszedges);

      // for edges root finding tri:
      // left (x = -1)
      Array<OneD, NekDouble> Vxm1t;//(totszedges);
      Array<OneD, NekDouble> Vdyxm1t;// = Array<OneD, NekDouble>(totszedges);
      Array<OneD, NekDouble> Vdxxm1t;// = Array<OneD, NekDouble>(totszedges);

      // bot (y = -1)
      Array<OneD, NekDouble> Vym1t;//  = (totszedges);
      Array<OneD, NekDouble> Vdxym1t;//  = Array<OneD, NekDouble>(totszedges);
      Array<OneD, NekDouble> Vdyym1t;//  = Array<OneD, NekDouble>(totszedges);

      // hypt tri (y = -x)
      Array<OneD, NekDouble> Vxyhypt;// = (totszedges);
      Array<OneD, NekDouble> Vdxxyhypt;// = Array<OneD, NekDouble>(totszedges);
      Array<OneD, NekDouble> Vdyxyhypt;// = Array<OneD, NekDouble>(totszedges);
      

      // for edges root finding hex:
      // edge front left (x = -1) (y = -1)
      Array<OneD, NekDouble> Vxm1ym1z   ;
      Array<OneD, NekDouble> Vdxxm1ym1z   ;
      Array<OneD, NekDouble> Vdyxm1ym1z   ;
      Array<OneD, NekDouble> Vdzxm1ym1z   ;
                
      //edge front right (x = 1) (y = -1)
      Array<OneD, NekDouble> Vx1ym1z   ;
      Array<OneD, NekDouble> Vdxx1ym1z   ;
      Array<OneD, NekDouble> Vdyx1ym1z   ;
      Array<OneD, NekDouble> Vdzx1ym1z   ;
               
      //edge front top (y = -1) (z = 1)
      Array<OneD, NekDouble> Vym1xz1   ;
      Array<OneD, NekDouble> Vdxym1xz1   ;
      Array<OneD, NekDouble> Vdyym1xz1   ;
      Array<OneD, NekDouble> Vdzym1xz1   ;
               
      //edge front bot (y = -1) (z = -1)
      Array<OneD, NekDouble> Vym1xzm1   ;
      Array<OneD, NekDouble> Vdxym1xzm1   ;
      Array<OneD, NekDouble> Vdyym1xzm1   ;
      Array<OneD, NekDouble> Vdzym1xzm1   ;

      // edge back left (y = 1), (x = -1)
      Array<OneD, NekDouble> Vxm1y1z   ;
      Array<OneD, NekDouble> Vdxxm1y1z   ;
      Array<OneD, NekDouble> Vdyxm1y1z   ;
      Array<OneD, NekDouble> Vdzxm1y1z   ;
                
      //edge back right (x = 1), (y = 1))
      Array<OneD, NekDouble> Vx1y1z   ;
      Array<OneD, NekDouble> Vdxx1y1z   ;
      Array<OneD, NekDouble> Vdyx1y1z   ;
      Array<OneD, NekDouble> Vdzx1y1z   ;
               
      //edge back top ( y = 1) (z = 1)
      Array<OneD, NekDouble> Vy1xz1   ;
      Array<OneD, NekDouble> Vdxy1xz1   ;
      Array<OneD, NekDouble> Vdyy1xz1   ;
      Array<OneD, NekDouble> Vdzy1xz1   ;
               
      //edge back bot (y = 1) (z = -1))
      Array<OneD, NekDouble> Vy1xzm1   ;
      Array<OneD, NekDouble> Vdxy1xzm1   ;
      Array<OneD, NekDouble> Vdyy1xzm1   ;
      Array<OneD, NekDouble> Vdzy1xzm1   ;

      // edge left bot (z = -1), (x = -1)
      Array<OneD, NekDouble> Vxm1yzm1   ;
      Array<OneD, NekDouble> Vdxxm1yzm1   ;
      Array<OneD, NekDouble> Vdyxm1yzm1   ;
      Array<OneD, NekDouble> Vdzxm1yzm1   ;
                
      //edge left top (x = -1), (z = 1))
      Array<OneD, NekDouble> Vxm1yz1   ;
      Array<OneD, NekDouble> Vdxxm1yz1   ;
      Array<OneD, NekDouble> Vdyxm1yz1   ;
      Array<OneD, NekDouble> Vdzxm1yz1   ;
               
      //edge right bot ( z = -1) (x = 1)
      Array<OneD, NekDouble> Vx1yzm1   ;
      Array<OneD, NekDouble> Vdxx1yzm1   ;
      Array<OneD, NekDouble> Vdyx1yzm1   ;
      Array<OneD, NekDouble> Vdzx1yzm1   ;
               
      //edge right top (z  1) (x  1))
      Array<OneD, NekDouble>Vx1yz1 ;
      Array<OneD, NekDouble>Vdxx1yz1 ;
      Array<OneD, NekDouble>Vdyx1yz1 ;
      Array<OneD, NekDouble> Vdzx1yz1;

      //hex edges:

      //surface bot z = -1
      Array<OneD, NekDouble>Vxyzm1 ;

      //surface right x = 1
      Array<OneD, NekDouble>Vx1yz ;

      //surface top z = 1
      Array<OneD, NekDouble>Vxyz1 ; 

      //surface left x = -1
      Array<OneD, NekDouble>Vxm1yz ;

      //surface front y = -1
      Array<OneD, NekDouble>Vxym1z ;

      //surface back y = 1
      Array<OneD, NekDouble>Vxy1z ;

      // only tets edges:

      // edge front left (x = -1) (y = -1)
      Array<OneD, NekDouble> Vxm1ym1ztet   ;
      Array<OneD, NekDouble> Vdxxm1ym1ztet   ;
      Array<OneD, NekDouble> Vdyxm1ym1ztet   ;
      Array<OneD, NekDouble> Vdzxm1ym1ztet   ;


      //edge front bot (y = -1) (z = -1)
      Array<OneD, NekDouble> Vym1xzm1tet   ;
      Array<OneD, NekDouble> Vdxym1xzm1tet   ;
      Array<OneD, NekDouble> Vdyym1xzm1tet   ;
      Array<OneD, NekDouble> Vdzym1xzm1tet   ;


      // edge left bot (z = -1), (x = -1)
      Array<OneD, NekDouble> Vxm1yzm1tet   ;
      Array<OneD, NekDouble> Vdxxm1yzm1tet   ;
      Array<OneD, NekDouble> Vdyxm1yzm1tet   ;
      Array<OneD, NekDouble> Vdzxm1yzm1tet   ;


      //edge front hypt (DB) (y = -1) (z = -x)                                                       
      Array<OneD, NekDouble> Vym1xmztet   ;
      Array<OneD, NekDouble> Vdxym1xmztet   ;
      Array<OneD, NekDouble> Vdyym1xmztet   ;
      Array<OneD, NekDouble> Vdzym1xmztet   ;

      //edge left hypt (DC) ( x = -1) (z = -y)                                                       
      Array<OneD, NekDouble> Vxm1ymztet   ;
      Array<OneD, NekDouble> Vdxxm1ymztet   ;
      Array<OneD, NekDouble> Vdyxm1ymztet   ;
      Array<OneD, NekDouble> Vdzxm1ymztet   ;

      // edge bot diag (BC) (z = -1) (y = -x)                                                        
      Array<OneD, NekDouble> Vxmyzm1tet   ;
      Array<OneD, NekDouble> Vdxxmyzm1tet   ;
      Array<OneD, NekDouble> Vdyxmyzm1tet   ;
      Array<OneD, NekDouble> Vdzxmyzm1tet   ;

      //tet surfaces:
      //surface bot z = -1 restriction on GD: x+y = 0 (ABC)
      Array<OneD, NekDouble>  Vxyzm1tet ;

      //surface left x = -1 restriction on GD: y + z = 0 (DAC)
      Array<OneD, NekDouble>Vxm1ypz0tet;

      //surf front y = -1   restriction on GD: x + z = 0 (DAB)
      Array<OneD, NekDouble> Vxpz0ym1tet;

      //surf DCB restriction  on GD: (x + y + z = -1)  
      Array<OneD, NekDouble> Vxpypzm1tet;

      //pyr dges:
      // edge front left (x = -1) (y = -1) EA
      Array<OneD, NekDouble> Vxm1ym1zpyr;
      Array<OneD, NekDouble> Vdyxm1ym1zpyr;
      Array<OneD, NekDouble> Vdxxm1ym1zpyr;
      Array<OneD, NekDouble> Vdzxm1ym1zpyr;

      //edge front hypt (y = -1) (z = -x) EB
      Array<OneD, NekDouble> Vym1xmzpyr   ;
      Array<OneD, NekDouble> Vdxym1xmzpyr   ;
      Array<OneD, NekDouble> Vdyym1xmzpyr   ;
      Array<OneD, NekDouble> Vdzym1xmzpyr   ;


      //edge back hypt (z = -x  or x = -z) (x = y) EC
      Array<OneD, NekDouble> Vxeyxmzpyr   ;
      Array<OneD, NekDouble> Vdxxeyxmzpyr   ;
      Array<OneD, NekDouble> Vdyxeyxmzpyr   ;
      Array<OneD, NekDouble> Vdzxeyxmzpyr   ;

      //edge back left (y = -z) (x = -1) ED
      Array<OneD, NekDouble> Vx1ymzpyr   ;
      Array<OneD, NekDouble> Vdxx1ymzpyr   ;
      Array<OneD, NekDouble> Vdyx1ymzpyr   ;
      Array<OneD, NekDouble> Vdzx1ymzpyr   ;


      //edge front bot (y = -1) (z = -1) AB
      Array<OneD, NekDouble> Vym1xzm1pyr   ;
      Array<OneD, NekDouble> Vdyym1xzm1pyr   ;
      Array<OneD, NekDouble> Vdxym1xzm1pyr   ;
      Array<OneD, NekDouble> Vdzym1xzm1pyr   ;


      //edge back bot (y = 1) (z = -1)) DC
      Array<OneD, NekDouble> Vy1xzm1pyr   ;
      Array<OneD, NekDouble> Vdyy1xzm1pyr   ;
      Array<OneD, NekDouble> Vdxy1xzm1pyr   ;
      Array<OneD, NekDouble> Vdzy1xzm1pyr   ;

      // edge left bot (z = -1), (x = -1) AD
      Array<OneD, NekDouble> Vxm1yzm1pyr   ;
      Array<OneD, NekDouble> Vdyxm1yzm1pyr   ;
      Array<OneD, NekDouble> Vdxxm1yzm1pyr   ;
      Array<OneD, NekDouble> Vdzxm1yzm1pyr   ;

      //edge right bot ( z = -1) (x = 1) BC
      Array<OneD, NekDouble> Vx1yzm1pyr   ;
      Array<OneD, NekDouble> Vdyx1yzm1pyr   ;
      Array<OneD, NekDouble> Vdxx1yzm1pyr   ;
      Array<OneD, NekDouble> Vdzx1yzm1pyr   ;

      //pyr surfaces
      //surface bot z = -1
      Array<OneD, NekDouble> Vxyzm1pyr ;

      //surface hypt x+z = 0
      Array<OneD, NekDouble> Vxmzypyr;

      //surface left x = -1
      Array<OneD, NekDouble> Vxm1yzpyr;


      //surface front y = -1
      Array<OneD, NekDouble> Vxym1zpyr;

      //surface back y +z = 0
      Array<OneD, NekDouble> Vxymzpyr;


    };

  }
}
#endif
