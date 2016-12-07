/*
 Copyright (c) 2014 Harvey Richardson, Michael Bareford
 All rights reserved.

 See the LICENSE file elsewhere in this distribution for the
 terms under which this software is made available.

*/

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
    static void Record(const int nstep, const int sstep);
    static void Finalise(void);

private:
    PowerMonitor(void) {};
    ~PowerMonitor(void) {};

    static void OpenCounterFiles(void);
    static void CloseCounterFiles(void);
    static bool IsAcceleratorCounter(const unsigned int i);
    static void GetFirstLine(const unsigned int i, char* line, const unsigned int len);
    static unsigned int GetCounterValue(const unsigned int i);
    static unsigned long long int GetLongCounterValue(const unsigned int i);

    static bool IsInitialised(void);


    static const char ver[];
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
    static double tm0, entot0;
    static int last_nstep;

    static bool all_initialised;
};

}
}

#endif
