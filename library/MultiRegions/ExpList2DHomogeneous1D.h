///////////////////////////////////////////////////////////////////////////////
//
// File ExpList2DHomogeneous1D.h
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
// Description: A 2D field which is homogeneous in 1 direction and so
// uses much of the functionality from a ExpList1D and its daughters
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPLIST2DHOMO1D_H
#define EXPLIST2DHOMO1D_H

#include <vector>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpListHomogeneous1D.h>
#include <MultiRegions/ExpList1D.h>
#include <SpatialDomains/MeshGraph1D.h>

namespace Nektar
{
    namespace MultiRegions
    {

        // Forward declaration for typedefs
        class ExpList2DHomogeneous1D;

        /// Shared pointer to an ExpList2DHomogeneous1D object.
        typedef boost::shared_ptr<ExpList2DHomogeneous1D>      ExpList2DHomogeneous1DSharedPtr;
        /// Vector of pointers to ExpList2DHomogeneous1D objects.
        typedef std::vector< ExpList2DHomogeneous1DSharedPtr > ExpList2DHomogeneous1DVector;
        /// Iterator for the vector of ExpList2DHomogeneous1D pointers.
        typedef std::vector< ExpList2DHomogeneous1DSharedPtr >::iterator ExpList2DHomogeneous1DVectorIter;

        /// Abstraction of a two-dimensional multi-elemental expansion which
        /// is merely a collection of local expansions.
        class ExpList2DHomogeneous1D: public ExpListHomogeneous1D
        {
        public:
            /// Default constructor.
            ExpList2DHomogeneous1D();

            /// Sets up a list of local expansions based on an input mesh.
            ExpList2DHomogeneous1D(const LibUtilities::BasisKey &HomoBasis,
                                   const NekDouble lz, 
                                   SpatialDomains::MeshGraph1D &graph1D);
            
            ExpList2DHomogeneous1D(const LibUtilities::BasisKey &HomoBasis,
                                   const NekDouble lhom,
                                   boost::shared_ptr<StdRegions::StdExpansionVector> &exp,
                                   Array<OneD, ExpList1DSharedPtr> &planes);

            /// Copy constructor.
            ExpList2DHomogeneous1D(const ExpList2DHomogeneous1D &In);

            /// Destructor.
            ~ExpList2DHomogeneous1D();

            /// This function calculates the coordinates of all the elemental
            /// quadrature points \f$\boldsymbol{x}_i\f$.
            inline void GetCoords(Array<OneD, NekDouble> &coord_0,
                                  Array<OneD, NekDouble> &coord_1 = NullNekDouble1DArray,
                                  Array<OneD, NekDouble> &coord_2 = NullNekDouble1DArray);

            void GetCoords(const int eid,
                           Array<OneD, NekDouble> &xc0,
                           Array<OneD, NekDouble> &xc1,
                           Array<OneD, NekDouble> &xc2);
        protected:
            
            /// Definition of the total number of degrees of freedom and
            /// quadrature points. Sets up the storage for \a m_coeff and \a
            ///  m_phys.
            void             SetCoeffPhys(void);

            //  virtual functions

            virtual void v_GetCoords(Array<OneD, NekDouble> &coord_0,
                                     Array<OneD, NekDouble> &coord_1,
                                     Array<OneD, NekDouble> &coord_2);

            virtual void v_WriteTecplotZone(std::ofstream &outfile, 
                                            int expansion);

        private:
        };
        
        inline void ExpList2DHomogeneous1D::GetCoords(Array<OneD, NekDouble> &coord_0,
                                                      Array<OneD, NekDouble> &coord_1,
                                                      Array<OneD, NekDouble> &coord_2)
            
        {
            v_GetCoords(coord_0,coord_1,coord_2);
        }
            
    } //end of namespace
} //end of namespace

#endif//EXPLIST3DHOMO1D_H

/**
* $Log: v $
*
**/

