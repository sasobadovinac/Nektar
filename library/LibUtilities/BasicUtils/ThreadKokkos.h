///////////////////////////////////////////////////////////////////////////////
//
// File ThreadKokkos.h
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

#ifndef NEKTAR_LIBUTILITIES_THREADKOKKOS_H_ 
#define NEKTAR_LIBUTILITIES_THREADKOKKOS_H_ 

#include "LibUtilities/BasicUtils/Thread.h"
#include <vector>

#include <Kokkos_Core.hpp>

#include "LibUtilities/Memory/NekMemoryManager.hpp"

namespace Nektar
{
namespace Thread
{

/**
 * @brief Implementation of ThreadManager using Kokkos threads.
 */
class ThreadManagerKokkos: public ThreadManager
{
    public:
        /// Constructs a ThreadManagerBoost.
        ThreadManagerKokkos(unsigned int numWorkers);
        /// Shuts down threading.
        virtual ~ThreadManagerKokkos();

        virtual void QueueJobs(std::vector<ThreadJob*> &joblist);
        virtual void QueueJob(ThreadJob* job);
        virtual unsigned int GetNumWorkers();
        virtual unsigned int GetWorkerNum();
        virtual void SetNumWorkers(unsigned int num);
        virtual void SetNumWorkers();
        virtual unsigned int GetMaxNumWorkers();
        virtual void Wait();
        virtual void SetChunkSize(unsigned int chnk);
        virtual void SetSchedType(SchedType s);
        virtual bool InThread();
        virtual void Hold();
        virtual const std::string& GetType() const;

        /// Called by the factory method.
        static ThreadManagerSharedPtr Create(unsigned int numT)
        {
            return boost::shared_ptr<ThreadManager>(
                new ThreadManagerKokkos(numT));
        }

    private:
        ThreadManagerKokkos();
        ThreadManagerKokkos(const ThreadManagerKokkos&);
        
        // Member variables
        const unsigned int                         m_numThreads;
        unsigned int                               m_numWorkers;
        std::vector<ThreadJob*> *                  m_masterJob;        

        static std::string                         className;
        std::string                                m_type;
        bool                                       m_Initialised;
};

} // Thread
} /* namespace Nektar */


#endif /* THREADKOKKOS_H_ */
