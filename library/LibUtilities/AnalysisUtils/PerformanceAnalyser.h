////////////////////////////////////////////////////////////////////////////////
//
//  File: PerformanceAnalyser.h
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
//  Description: Class for reading hardware counters via CrayPAT API
//
////////////////////////////////////////////////////////////////////////////////

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
    static unsigned int Record(const int nstep, const int sstep, const bool initial_sync, const bool initial_rec);
    static void Finalise(void);

private:
    PerformanceAnalyser(void) {};
    ~PerformanceAnalyser(void) {};
            
    static bool IsInitialised(void);
    static int GetCategoryValue(const char* str);
    static unsigned int RecordCounterValues(const int nstep, const int sstep);
            
    static const char ver[];

    static const unsigned int PAT_RECORD_OK;
    static const unsigned int PAT_RECORD_UNINITIALISED;
    static const unsigned int PAT_RECORD_ERROR;
    
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
