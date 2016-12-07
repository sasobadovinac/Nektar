/*

 Copyright (c) 2016 Michael Bareford
 All rights reserved.

 See the LICENSE file elsewhere in this distribution for the
 terms under which this software is made available.

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <mpi.h>
#include "pat_api.h"
#include <LibUtilities/AnalysisUtils/PerformanceAnalyser.h> 

namespace Nektar
{	
namespace LibUtilities
{

// initialise static constant attributes of PerformanceAnalyser class
const char PerformanceAnalyser::ver[] = "2.0.0";
const unsigned int PerformanceAnalyser::MAX_NAME_LEN(128);
const char PerformanceAnalyser::PAT_RT_SEPARATOR[] = ",";
const unsigned int PerformanceAnalyser::PAT_REGION_OPEN(1);
const unsigned int PerformanceAnalyser::PAT_REGION_MONITOR(2);

// initialise static attributes of PerformanceAnalyser class
int PerformanceAnalyser::rank(-1);
int PerformanceAnalyser::root_rank(-1);
FILE* PerformanceAnalyser::log_fp(NULL);
bool PerformanceAnalyser::first_monitor(true);
double PerformanceAnalyser::tm0(0.0);
double PerformanceAnalyser::tm(0.0);
int PerformanceAnalyser::last_nstep(0);
  
int PerformanceAnalyser::ncat(0);
int* PerformanceAnalyser::cat_id(NULL);
int* PerformanceAnalyser::cat_ncntr(NULL);
MPI_Comm* PerformanceAnalyser::cat_comm(NULL);
  
int PerformanceAnalyser::ncntr(0);
char*** PerformanceAnalyser::cat_cntr_name(NULL);
unsigned long** PerformanceAnalyser::cat_cntr_value(NULL);
unsigned long long** PerformanceAnalyser::cat_cntr_value_tot(NULL);
		
bool PerformanceAnalyser::open(false);
bool PerformanceAnalyser::debug(false);  

        
// private methods
//////////////////////////////////////////////////////////////////////////////////

// return true if PerformanceAnalyser::Initialise has been called successfully
bool PerformanceAnalyser::IsInitialised(void)
{
    bool ok(false);

    if (-1 != rank)
    {
        ok = true;
        if (ncat > 0)
	{
            ok = (ok && cat_id && cat_ncntr && cat_comm);
            if (!ok)
	    {
                if (cat_id) free(cat_id);
                if (cat_ncntr) free(cat_ncntr);
                if (cat_comm) free(cat_comm);
                cat_id = NULL;
                cat_ncntr = NULL;
                cat_comm = NULL;
            }
            else
	    {
                if (ncntr > 0)
		{
                    ok = (ok && cat_cntr_name && cat_cntr_value && cat_cntr_value_tot);
                    if (!ok)
		    {
                        if (cat_cntr_name)
			{
                            for (int i = 0; i < ncat; i++)
			    {
                                if (cat_cntr_name[i])
				{
                                    for (int j = 0; j < cat_ncntr[i]; j++)
				    {
                                        if (cat_cntr_name[i][j]) free(cat_cntr_name[i][j]);
                                    }
                                    free(cat_cntr_name[i]);
                                }
                            }
                            free(cat_cntr_name);
                            cat_cntr_name = NULL;
                        }

                        if (cat_cntr_value)
			{
                            for (int i = 0; i < ncat; i++)
			    {
                                if (cat_cntr_value[i]) free(cat_cntr_value[i]);
                            }
                            free(cat_cntr_value);
                            cat_cntr_value = NULL;
                        }

                        if (cat_cntr_value_tot)
			{
                            for (int i = 0; i < ncat; i++)
			    {
                                if (cat_cntr_value_tot[i]) free(cat_cntr_value_tot[i]);
                            }
                            free(cat_cntr_value_tot);
                            cat_cntr_value_tot = NULL;
                        }
                    }

                } // end of <if (ncntr > 0)> clause

            }

        } // end of <if (ncat > 0)> clause

        if (rank == root_rank)
        {
	    ok = (ok && (NULL != log_fp));
        }

    } // end of <if (-1 != rank)> clause

    return ok;
}
  

// convert string to counter category number
int PerformanceAnalyser::GetCategoryValue(const char* str)
{
    int cat(0);

    if (NULL == str)
    {
        cat = 0;
    }
    else if (0 == strcmp("PAT_CTRS_CPU", str))
    {
        // HWPCs on processor
        cat = PAT_CTRS_CPU;
    }
    else if (0 == strcmp("PAT_CTRS_NETWORK", str))
    {
        // NWPCs on network router
        cat = PAT_CTRS_NETWORK;
    }
    else if (0 == strcmp("PAT_CTRS_ACCEL", str))
    {
        // HWPCs on attached GPUs
        cat = PAT_CTRS_ACCEL;
    }
    else if (0 == strcmp("PAT_CTRS_NB", str))
    {
        // AMD NorthBridge on nodes
        cat = PAT_CTRS_NB;
    }
    else if (0 == strcmp("PAT_CTRS_RAPL", str))
    {
        // Running Avr Power Level on package
        cat = PAT_CTRS_RAPL;
    }
    else if (0 == strcmp("PAT_CTRS_PM", str))
    {
        // Cray Power Management
        cat = PAT_CTRS_PM;
    }
    else if (0 == strcmp("PAT_CTRS_UNCORE", str))
    {
        // Intel Uncore on socket
        cat = PAT_CTRS_UNCORE;
    }
    else
    {
        // next counter category?
        cat = PAT_CTRS_UNCORE+1;
    }

    return cat;
}
	

// public methods
//////////////////////////////////////////////////////////////////////////////////

// determine which processes are reading which counters
// determine the root process
// call pat_mpi_monitor(-1,1)
void PerformanceAnalyser::Initialise(const char* log_fpath)
{
    int str_len(0), nrank(0);
    int ncntr_max(0);
    int pat_res(0);
    char* cat_list_str(NULL);
    char* cat_name(NULL);

    if (open)
    {
        return;
    }

    ncat = 0;
    ncntr = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    // get the number of counter categories
    ///////////////////////////////////////////////
    cat_list_str = getenv("MY_RT_CTRCAT");
    str_len = strlen(cat_list_str);
    ncat = str_len > 0 ? 1 : 0;
    for (int i = 0; i < str_len; i++)
    {
        if (PAT_RT_SEPARATOR[0] == cat_list_str[i])
	{
            ncat++;
        }
    }

    if (debug && 0 == rank)
    {
        printf("%d: %d categories listed.\n", rank, ncat);
    }

    if (ncat <= 0)
    {
        return;
    }
    ///////////////////////////////////////////////


    // for each category, get the id, number of counters and setup an mpi communicator
    // then determine the root rank
    // first allocate the necessary arrays
    ///////////////////////////////////////////////////////////////////////////////////
    cat_id = (int*) calloc(ncat, sizeof(int));
    cat_ncntr = (int*) calloc(ncat, sizeof(int));
    cat_comm = (MPI_Comm*) calloc(ncat, sizeof(MPI_Comm));
    if (!PerformanceAnalyser::IsInitialised())
    {
        ncat = 0;
        ncntr = 0;
    }
    else
    {
        pat_res = PAT_record(PAT_STATE_ON);
        if (PAT_API_OK != pat_res)
	{
            printf("%d: PAT_record(PAT_STATE_ON) failed with error %d.\n", rank, pat_res);
        }

        cat_name = strtok(cat_list_str, PAT_RT_SEPARATOR);
        ncntr = 0;
        pat_res = PAT_region_begin(PAT_REGION_OPEN, "PerformanceAnalyser::Initialise");
        if (PAT_API_OK != pat_res)
	{
            printf("%d: PAT_region_begin(PAT_REGION_OPEN) failed with error %d.\n", rank, pat_res);
        }
        for (int i = 0; i < ncat; i++)
	{
            cat_id[i] = PerformanceAnalyser::GetCategoryValue(cat_name);
            cat_name = strtok(NULL, PAT_RT_SEPARATOR);

            pat_res = PAT_counters(cat_id[i], 0, 0, &cat_ncntr[i]);
            if (1 == debug)
	    {
                if (PAT_API_OK != pat_res)
		{
                    printf("%d: PAT_counters failed with error %d within pat_mpi_open.\n",
                        rank, pat_res);
                }
                else
		{
                    printf("%d: counter category %d has %d counter(s).\n",
                        rank, cat_id[i], cat_ncntr[i]);
                }
            }
            MPI_Comm_split(MPI_COMM_WORLD, cat_ncntr[i], rank, &cat_comm[i]);

            ncntr += cat_ncntr[i];
        }
        pat_res = PAT_region_end(PAT_REGION_OPEN);
        if (PAT_API_OK != pat_res)
	{
            printf("%d: PAT_region_end(PAT_REGION_OPEN) failed with error %d.\n", rank, pat_res);
        }
    }

    ncntr_max = 0;
    MPI_Allreduce(&ncntr, &ncntr_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);
    root_rank = (ncntr == ncntr_max) ? rank : nrank;
    MPI_Allreduce(MPI_IN_PLACE, &root_rank, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD);
    ///////////////////////////////////////////////////////////////////////////////////


    // root rank opens file for counter data
    ////////////////////////////////////////
    if (rank == root_rank)
    {
        if (NULL != log_fp)
	{
            fclose(log_fp);
            log_fp = NULL;
        }
        // open performance counter data file
        if (NULL != log_fpath)
	{
            log_fp = fopen(log_fpath, "w");
        }
        if (NULL == log_fp)
	{
            log_fp = fopen("./pat_log.out", "w");
        }
    }
    ////////////////////////////////////////


    if (ncat > 0 && ncntr > 0)
    {
        // for each category, allocate the counter name, value and total arrays
        ////////////////////////////////////////////////////////////////////////////////////////
        cat_cntr_name = (char***) calloc(ncat, sizeof(char**));
        cat_cntr_value = (unsigned long**) calloc(ncat, sizeof(unsigned long*));
        cat_cntr_value_tot = (unsigned long long**) calloc(ncat, sizeof(unsigned long long*));
        if (cat_cntr_name && cat_cntr_value && cat_cntr_value_tot)
	{
            for (int i = 0; i < ncat; i++)
	    {
                cat_cntr_name[i] = (char**) calloc(cat_ncntr[i], sizeof(char*));
                if (cat_cntr_name[i])
		{
                    for (int j = 0; j < cat_ncntr[i]; j++)
		    {
                        cat_cntr_name[i][j] = (char*) calloc(MAX_NAME_LEN, sizeof(char));
                    }
                }

                cat_cntr_value[i] = (unsigned long*) calloc(cat_ncntr[i], sizeof(unsigned long));
                cat_cntr_value_tot[i] = (unsigned long long*) calloc(cat_ncntr[i], sizeof(unsigned long long));
            }
        }

        if (!PerformanceAnalyser::IsInitialised())
	{
            ncntr = 0;
        }
        /////////////////////////////////////////////////////////////////////////////////////////
    }


    bool all_ok = PerformanceAnalyser::IsInitialised();
    if (debug)
    {
        printf("%d: pat_ok() returned %d.\n", rank, all_ok);
    }
    MPI_Allreduce(MPI_IN_PLACE, &all_ok, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD);
    if (!all_ok)
    {
        open = false;
        PerformanceAnalyser::Finalise();
    }
    else
    {
        // do initial monitoring call, which ends with MPI_Barrier
        open = true;
        first_monitor = true;
	PerformanceAnalyser::Record(-1, 1);
    }
  
} // end of PerformanceAnalyser::Initialise method



// read counter values if first rank on node,
// and output those values if rank zero
void PerformanceAnalyser::Record(const int nstep, const int sstep) {

    int ncntr_test(0), pat_res(0);
    unsigned long long ncntr_val(0);

    if (!open)
    {
        return;
    }

    if (rank == root_rank)
    {
        tm = MPI_Wtime();
        if (first_monitor)
	{
            tm0 = tm;
            first_monitor = false;
        }
    }


    if (ncat > 0 && ncntr > 0)
    {
        // get counter values
        ///////////////////////////////////////////////////////////////////////////////////////////
        pat_res = PAT_region_begin(PAT_REGION_MONITOR, "pat_mpi_monitor");
        if (PAT_API_OK != pat_res)
	{
            printf("%d: PAT_region_begin(PAT_REGION_MONITOR) failed with error %d.\n", rank, pat_res);
        }
        for (int i = 0; i < ncat; i++)
	{
            pat_res = PAT_counters(cat_id[i], (const char**) cat_cntr_name[i], cat_cntr_value[i], &ncntr_test);
            if (1 == debug)
	    {
                if (PAT_API_OK != pat_res)
		{
                    printf("%d: PAT_counters failed with error %d within pat_mpi_monitor.\n", rank, pat_res);
                }
                else if (cat_ncntr[i] != ncntr_test)
		{
                    printf("%d: counter category %d has %d counter(s) when %d expected.\n",
                        rank, cat_id[i], ncntr_test, cat_ncntr[i]);
                }
            }
        }
        pat_res = PAT_region_end(PAT_REGION_MONITOR);
        if (PAT_API_OK != pat_res)
	{
            printf("%d: PAT_region_end(PAT_REGION_MONITOR) failed with error %d.\n", rank, pat_res);
        }
        ///////////////////////////////////////////////////////////////////////////////////////////

        // get counter value totals
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for (int i = 0; i < ncat; i++)
	{
            for (int j = 0; j < cat_ncntr[i]; j++)
	    {
                ncntr_val = (unsigned long long) cat_cntr_value[i][j];
                MPI_Reduce(&ncntr_val, &cat_cntr_value_tot[i][j], 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, root_rank, cat_comm[i]);
            }
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }


    if (rank == root_rank)
    {
        if (NULL != log_fp)
	{
            // output data
            if (tm0 == tm)
	    {
                fprintf(log_fp, "PerformanceAnalyser v%s: time (s), step, substep", ver);
                for (int i = 0; i < ncat; i++)
		{
                    for (int j = 0; j < cat_ncntr[i]; j++)
		    {
                        fprintf(log_fp, ", %s", cat_cntr_name[i][j]);
                    }
                }
                fprintf(log_fp, "\n");
            }

            // update counter data file
            fprintf(log_fp, "%f %d %d", tm-tm0, nstep, sstep);
            for (int i = 0; i < ncat; i++)
	    {
                for (int j = 0; j < cat_ncntr[i]; j++)
		{
                    fprintf(log_fp, " %llu", cat_cntr_value_tot[i][j]);
                }
            }
            fprintf(log_fp, "\n");
        }
    }

    last_nstep = nstep;

    MPI_Barrier(MPI_COMM_WORLD);
    			
} // end of PerformanceAnalyser::Record method


// close the files used to read and record counter data
void PerformanceAnalyser::Finalise(void)
{
    int pat_res(0);

    if (open)
    {
        // do the last monitoring call
        PerformanceAnalyser::Record(last_nstep+1, 1);
    }

    // turn off recording and free memory   
    pat_res = PAT_record(PAT_STATE_OFF);
    if (PAT_API_OK != pat_res)
    {
        printf("%d: PAT_record(PAT_STATE_OFF) failed with error %d.\n", rank, pat_res);
    }

    if (ncat > 0)
    {
        if (ncntr > 0)
	{
            for (int i = 0; i < ncat; i++)
	    {
                for (int j = 0; j < cat_ncntr[i]; j++)
		{
                    free(cat_cntr_name[i][j]);
                }
                free(cat_cntr_name[i]);
                free(cat_cntr_value[i]);
                free(cat_cntr_value_tot[i]);
            }

            free(cat_cntr_name);
            free(cat_cntr_value);
            free(cat_cntr_value_tot);

            cat_cntr_name = NULL;
            cat_cntr_value = NULL;
            cat_cntr_value_tot = NULL;

            ncntr = 0;
        }

        free(cat_id);
        free(cat_ncntr);
        free(cat_comm);

        cat_id = NULL;
        cat_ncntr = NULL;
        cat_comm = NULL;

        ncat = 0;
    }


    if (rank == root_rank)
    {
        if (NULL != log_fp)
	{
            // close performance counter data file
            fclose(log_fp);
            log_fp = NULL;
        }
    }

    root_rank = -1;
    rank = -1;

    open = false;
    MPI_Barrier(MPI_COMM_WORLD);  
  		  
} // end of PerformanceAnalyser::Finalise method

}
}
