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
}


/**
 * Terminates all running threads (they will finish their current job),
 * releases resources and destructs.
 */
ThreadManagerKokkos::~ThreadManagerKokkos()
{
    //empty
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
    Kokkos::InitArguments args;
    args.num_threads = GetNumWorkers();
    Kokkos::initialize(args);
    //std::cout << "Number of threads: = " << args.num_threads << std::endl;

    std::vector<ThreadJob*> *jobs = m_masterJob;
    int jobSize = jobs->size();
    //std::cout << "Number of threads in Kokkos::parallel_for loop: " << N << std::endl;

    Kokkos::parallel_for ( jobSize, KOKKOS_LAMBDA (const int& i) {
        jobs->at(i)->Run();
        //printf("Hello Kokkos %i\n", i);
    });

    Kokkos::finalize();
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
