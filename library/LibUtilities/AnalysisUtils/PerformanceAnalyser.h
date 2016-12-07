/*
 Copyright (c) 2016 Michael Bareford
 All rights reserved.

 See the LICENSE file elsewhere in this distribution for the
 terms under which this software is made available.

*/

#ifndef NEKTAR_LIB_UTILITIES_PERFORMANCEANALYSER_H
#define NEKTAR_LIB_UTILITIES_PERFORMANCEANALYSER_H

namespace Nektar
{
namespace LibUtilities
{
        
class PerformanceAnalyser
{
public:
    static void Initialise(const char *out_fn);
    static void Record(const int nstep, const int sstep);
    static void Finalise(void);

private:
    PerformanceAnalyser(void) {};
    ~PerformanceAnalyser(void) {};
            
    static bool IsInitialised(void);
    static int GetCategoryValue(const char* str);
            
    static const char ver[];
    static const unsigned int MAX_NAME_LEN;
    static const char PAT_RT_SEPARATOR[];
    static const unsigned int PAT_REGION_OPEN;
    static const unsigned int PAT_REGION_MONITOR;

    static int rank;
    static int root_rank;
    static FILE* log_fp;
    static bool first_monitor;
    static double tm0;
    static double tm;
    static int last_nstep;
  
    static int ncat;
    static int* cat_id;
    static int* cat_ncntr;
    static int* cat_comm;
  
    static int ncntr;
    static char*** cat_cntr_name;
    static unsigned long** cat_cntr_value;
    static unsigned long long** cat_cntr_value_tot;
		
    static bool open;
    static bool debug; 
};

}
}

#endif
