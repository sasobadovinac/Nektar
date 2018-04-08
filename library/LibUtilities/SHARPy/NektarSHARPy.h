///////////////////////////////////////////////////////////////////////////////
//
// File NektarSHARPy.h
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
// Description: Header file for the SHARPy class in Nektar++
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILIITIES_SHARPY_NEKTARSHARPY_H
#define NEKTAR_LIB_UTILIITIES_SHARPY_NEKTARSHARPY_H

#include <LibUtilities/BasicUtils/SessionReader.h>
			
namespace Nektar
{
    template <typename Dim, typename DataType>
    class Array;
    
	namespace LibUtilities
	{
		/**
		 * The NektarSHARPy class manages the use of the third-party code SHARPy to solve nonlinear structural dynamics.
		 * The function here defined will link to a proper implementation of the SHARPy routine.
		 * Depending on the user definition the functions can link to a class which is a wrapper around the SHARPy
		 * library or to a specific SHARPy implementation.
		 */
		class NektarSHARPy;
		
		// A shared pointer to the NektarSHARPy object
		typedef std::shared_ptr<NektarSHARPy>  NektarSHARPySharedPtr;
		
		/// Datatype of the NekFactory used to instantiate classes derived from
		/// the NektarFFT class.
		typedef LibUtilities::NekFactory< std::string, NektarSHARPy> NektarSHARPyFactory;
		
		LIB_UTILITIES_EXPORT NektarSHARPyFactory& GetNektarSHARPyFactory();

		class NektarSHARPy
		{
		public:
			/// Initialises NektarFFT class members.
			LIB_UTILITIES_EXPORT NektarSHARPy();
			
			// Distructor
			LIB_UTILITIES_EXPORT  ~NektarSHARPy();
			
		
			/**
			 *
			 **/
			LIB_UTILITIES_EXPORT void InitialiseStatic(const LibUtilities::SessionReaderSharedPtr &pSession);

			/**
			 *
			 **/
			LIB_UTILITIES_EXPORT void InitialiseDynamic(const LibUtilities::SessionReaderSharedPtr& pSession, const Array<OneD, Array<OneD, NekDouble> > &CdCl);
		

			/**
			 *
			 **/
			LIB_UTILITIES_EXPORT void SolvenlnStatic(const LibUtilities::SessionReaderSharedPtr& pSession, const Array<OneD, Array<OneD, NekDouble> >& CdCl);
			

			/**
             *
			 **/
			LIB_UTILITIES_EXPORT void SolvenlnDynamic(const LibUtilities::SessionReaderSharedPtr& pSession, const Array<OneD, Array<OneD, NekDouble> > &CdCl, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &Vmotions);
	
		protected:
			
			virtual void v_InitialiseStatic(const LibUtilities::SessionReaderSharedPtr &pSession);	

			virtual void v_InitialiseDynamic(const LibUtilities::SessionReaderSharedPtr& pSession, const Array<OneD, Array<OneD, NekDouble> > &CdCl);

			virtual void v_SolvenlnStatic(const LibUtilities::SessionReaderSharedPtr& pSession, const Array<OneD, Array<OneD, NekDouble> >& CdCl);

			virtual void v_SolvenlnDynamic(const LibUtilities::SessionReaderSharedPtr& pSession, const Array<OneD, Array<OneD, NekDouble> > &CdCl, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &Vmotions);

		private:
			
		};
	}//end namespace LibUtilities
}//end of namespace Nektar
#endif //NEKTAR_LIB_UTILIITIES_SHARPY_NEKTARSHARPY_H
