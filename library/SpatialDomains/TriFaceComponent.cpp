////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/SpatialDomains/TriFaceComponent.cpp,v $
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific 
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:  
//
//
////////////////////////////////////////////////////////////////////////////////
#include "pchSpatialDomains.h"

#include <SpatialDomains/TriFaceComponent.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        TriFaceComponent::TriFaceComponent()
        {
        }

        TriFaceComponent::TriFaceComponent(const int coordim):
        Geometry2D(coordim)
        {
            const StdRegions::BasisKey B0(StdRegions::eModified_A, 2,
                StdRegions::eLobatto, 3,0,0);
            const StdRegions::BasisKey B1(StdRegions::eModified_B, 2,
                StdRegions::eRadauM, 3,1,0);
            m_xmap = new StdRegions::StdTriExp * [m_coordim];

            for(int i = 0; i < m_coordim; ++i)
            {
#pragma message("Fix this! " __FILE__)
                m_xmap[i] = NULL; //new StdRegions::StdTriExp(B0,B1);
            }

            m_state = eNotFilled;
        }


        TriFaceComponent::~TriFaceComponent()
        {
            if(m_xmap)
            {
                for(int i = 0; i < m_coordim; ++i)
                {
#pragma message("Fix this! " __FILE__)
                    //delete m_xmap[i]; 
                }

                delete[] m_xmap;
                m_xmap = (StdRegions::StdTriExp **) NULL;
            }
        }


        TriFaceComponent::TriFaceComponent(const TriFaceComponent &T)
        {
            m_fid = T.m_fid;

            std::list<CompToElmt>::const_iterator def;
            for(def = T.m_elmtmap.begin(); def != T.m_elmtmap.end(); def++)
            {
                m_elmtmap.push_back(*def);    
            }
        }


        void TriFaceComponent::AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtmap.push_back(ee);
        }


        int TriFaceComponent::NumElmtConnected() const
        {
            return int(m_elmtmap.size());
        }


        bool TriFaceComponent::IsElmtConnected(int gvo_id, int locid) const
        {
            std::list<CompToElmt>::const_iterator def;
            CompToElmt ee(gvo_id,locid);

            def = find(m_elmtmap.begin(),m_elmtmap.end(),ee);

            // Found the element connectivity object in the list
            return (def != m_elmtmap.end());
        }

        /** given local collapsed coordinate Lcoord return the value of
        physical coordinate in direction i **/


        double TriFaceComponent::GetCoord(const int i, const double *Lcoord)
        {
            ASSERTL1(m_state == ePtsFilled,
                "Goemetry is not in physical space");

            return m_xmap[i]->Evaluate(Lcoord);
        }
    }; //end of namespace
}; //end of namespace

//
// $Log: TriFaceComponent.cpp,v $
// Revision 1.1  2006/05/04 18:59:04  kirby
// *** empty log message ***
//
// Revision 1.9  2006/04/09 02:08:36  jfrazier
// Added precompiled header.
//
// Revision 1.8  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.7  2006/03/12 14:20:44  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.6  2006/03/12 11:06:40  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.5  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
