///////////////////////////////////////////////////////////////////////////////
//
// File ThreadKokkoscpp
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

#define NOMINMAX

#include "LibUtilities/BasicUtils/ThreadKokkos.h"
#include <iostream>

#include <Kokkos_Core.hpp>


namespace Nektar
{
namespace Thread
{

std::string ThreadManagerKokkos::className =
    GetThreadManagerFactory().RegisterCreatorFunction("ThreadManagerKokkos",
            ThreadManagerKokkos::Create, "Threading using Kokkos.");


/**
 * @param numWorkers The number of threads to start (including master thread).
 * @note Do not use, use factory instead.
 */
ThreadManagerKokkos::ThreadManagerKokkos(unsigned int numT) :
        m_numThreads(numT), m_numWorkers(numT)
{        
    m_type = "Threading with Kokkos";
    m_Initialised = false;
    
    if ( !Kokkos::HostSpace::execution_space::is_initialized())
    {
        // ** Initialise Kokkos with number of threads **
        Kokkos::InitArguments args;
        args.num_threads = m_numThreads;
        Kokkos::DefaultHostExecutionSpace::initialize(args.num_threads);
        std::cout << "ThreadManager::Initialise" << std::endl;
        std::cout << "Number of threads: = " << args.num_threads << std::endl;
        m_Initialised = true;
    }
}


/**
 * Terminates all running threads (they will finish their current job),
 * releases resources and destructs.
 */
ThreadManagerKokkos::~ThreadManagerKokkos()
{
    // ** finalize kokkos **
    if ( m_Initialised)
    {
        Kokkos::DefaultHostExecutionSpace::finalize();
        std::cout << "ThreadManager::Finalise" << std::endl;
    }
}


/**
 *
 */
void ThreadManagerKokkos::QueueJobs(std::vector<ThreadJob*> &joblist)
{
    m_masterJob = &joblist;    
}


/*
 *
 */
void ThreadManagerKokkos::QueueJob(ThreadJob *job)
{
    //empty
}


/**
 *
 */
void ThreadManagerKokkos::SetChunkSize(unsigned int chnk)
{
    //empty
}


/**
 *
 */
void ThreadManagerKokkos::SetSchedType(SchedType s)
{
    //empty
}


/**
 *
 */
bool ThreadManagerKokkos::InThread()
{
    return false;
}


/**
 *
 */
void ThreadManagerKokkos::Wait()
{
    // ** retrieve queued jobs **
    std::vector<ThreadJob*> *jobs = m_masterJob;

    // ** establish size of parallel_for loop **
    int jobSize = jobs->size();
    //std::cout << "Number of threads in Kokkos::parallel_for loop: " << N << std::endl;

    // ** call Kokkos parallel_for loop **
    typedef Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace> range_policy;
    Kokkos::parallel_for ( range_policy(0,jobSize), KOKKOS_LAMBDA (const int& i) {
        jobs->at(i)->Run();
        //printf("Hello Kokkos %i\n", i);
    });
}


/**
 *
 */
unsigned int ThreadManagerKokkos::GetNumWorkers()
{
    return m_numWorkers;
}


/**
 *
 */
unsigned int ThreadManagerKokkos::GetWorkerNum()
{
    return 0;
}


/**
 *
 */
void ThreadManagerKokkos::SetNumWorkers(unsigned int num)
{
    m_numWorkers = num;
}


/**
 *
 */
void ThreadManagerKokkos::SetNumWorkers()
{
    SetNumWorkers(m_numThreads);
}


/**
 *
 */
unsigned int ThreadManagerKokkos::GetMaxNumWorkers()
{
    return m_numThreads;
}


/**
 *
 */
void ThreadManagerKokkos::Hold()
{
    //empty
}


/**
 *
 */
const std::string& ThreadManagerKokkos::GetType() const
{
    return m_type;
}




} // Thread
} /* namespace Nektar */
