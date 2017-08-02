////////////////////////////////////////////////////////////////////////////////
//
//  File: PowerMonitor.h
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
//  Description: Class for reading PM counter files
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_POWERMONITOR_H
#define NEKTAR_LIB_UTILITIES_POWERMONITOR_H


namespace Nektar
{
namespace LibUtilities
{

// Supported power management counters
typedef enum pm_counters
{
    PM_COUNTER_FRESHNESS = 0,  // The freshness counter MUST be first in the list
    PM_COUNTER_POWER,
    PM_COUNTER_ENERGY,
    PM_COUNTER_ACCEL_POWER,
    PM_COUNTER_ACCEL_ENERGY,
    PM_COUNTER_STARTUP,
    PM_COUNTER_POWER_CAP,
    PM_COUNTER_ACCEL_POWER_CAP,
    PM_NCOUNTERS
} pm_counter_e ;

/// Base class for unsteady solvers.
class PowerMonitor
{
public:
    static void Initialise(const char *out_fn);
    static unsigned int Record(const int nstep, const int sstep, const bool initial_sync, const bool initial_rec);
    static void Finalise(void);

private:
    PowerMonitor(void) {};
    ~PowerMonitor(void) {};

    static void OpenCounterFiles(void);
    static void CloseCounterFiles(void);
    static bool IsAcceleratorCounter(const unsigned int i);
    static int GetFirstLine(const unsigned int i, char* line, const unsigned int len);
    static long int GetCounterValue(const unsigned int i);
    static int GetNodeNumber(void);
    
    static bool IsInitialised(void);

    static unsigned int RecordCounterValues(const int nstep, const int sstep);


    static const char ver[];

    static const unsigned int PM_RECORD_OK;
    static const unsigned int PM_RECORD_UNINITIALISED;
    static const unsigned int PM_RECORD_BLADE_RESTART;
    static const unsigned int PM_RECORD_COUNTER_FILE_ERROR;
  
    static const unsigned int MAX_FPATH_LEN;
    static const unsigned int MAX_FLINE_LEN;

    static const char sys_pm_cnt_dir[];
    static const char* cnt_fname[];

    static int rank, min_node_rank;
    static int monitor_cnt, non_monitor_cnt;
    static bool first_record;
    static int mpi_comm_monitor;

    static FILE* cnt_fp[PM_NCOUNTERS];
    static FILE* log_fp;
    static double tm0;
    static long int entot0;
    static int last_nstep;
    static long int init_startup;

    static bool all_initialised;

    static int system_error;
};

}
}

#endif
