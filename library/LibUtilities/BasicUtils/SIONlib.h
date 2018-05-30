///////////////////////////////////////////////////////////////////////////////
//
// File SIONlib.h
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
// Description: Simple OO wrapper around SIONlib
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_SIONLIB_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_SIONLIB_H

#include "mpi.h"
#include "sion.h"

#include <exception>
#include <string>
#include <string.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <fstream>

#include <LibUtilities/Foundations/BasisType.h>

namespace Nektar
{
namespace LibUtilities
{
namespace SIONlib
{
  
class SION_Base
{

public:
    static const unsigned int LUSTRE_STRIPE_BASE_SIZE;
    
    SION_Base() : _sid(-9999) {}
    ~SION_Base() {}

    // Common functions
    char *getSionFileName() const;

    void setDebug(bool debug);
    
    void setMode(std::string mode);
    std::string getMode() const;
    bool isModeCollective() const;
    bool isModeCollectiveMerge() const;

    void setNumberOfFiles(int num_files);
    int getNumberOfFiles() const;

    void setNumberOfTasks(int num_tasks);
    int getNumberOfTasks() const;

    void setRank(int rank);
    int getRank() const;

    void setChunkSize(sion_int64 chunk_size);
    sion_int64 getChunkSize() const;

    void setChunkSizes(sion_int64 *chunk_sizes);
    sion_int64 *getChunkSizes() const;

    void setGlobalRanks(int *global_ranks);
    int *getGlobalRanks() const;

    void setFileSystemBlockSize(sion_int32 fs_blk_size);
    sion_int32 getFileSystemBlockSize() const;

    int getNumberOfSuccessfulReadElements() const;
   
    int getSid() const;

    int getReturnCode() const;

    void seek();
    
    /* get information (with sion datatypes) */
    int getFileEndianness() const;
    sion_int64 getBytesWritten() const;
    sion_int64 getBytesRead() const;
    sion_int64 getBytesAvailInBlock() const;
    sion_int64 getBytesAvailInChunk() const;
    sion_int64 getPosition() const;

protected:
    // Common attributes
    char *_sion_file_name;
    std::string _mode;
    bool _is_mode_collective;
    bool _is_mode_collectivemerge;
    int _num_files;
    int _num_tasks;
    int _rank;
    sion_int64 *_chunk_sizes;
    sion_int64 _chunk_size;
    sion_int32 _fs_blk_size;
    int *_global_ranks;
    FILE *_file_ptr;
    int _number_of_elements_sucessfully_read;
    int _return_code;
    int _sid;
    bool _debug;

    //  get information (with sion datatypes) 
    int _file_endianness;
    sion_int64 _bytes_written;
    sion_int64 _bytes_read;
    sion_int64 _bytes_avail_in_block;
    sion_int64 _bytes_avail_in_chunk;
    sion_int64 _position;
};


class SIONFile : public SION_Base
{

public:
      
    SIONFile();
    SIONFile(std::string sion_file_name, std::string mode = "bw",
        int num_files = 1, sion_int64 chunk_size=LUSTRE_STRIPE_BASE_SIZE, sion_int32 block_size=-1,
        int global_rank = 0, MPI_Comm gComm=MPI_COMM_WORLD, MPI_Comm lComm=MPI_COMM_WORLD);
    virtual ~SIONFile();

    char *getNewSionFileName() const;
    void setLocalCommunicator(MPI_Comm lComm);
    MPI_Comm getLocalCommunicator() const;

    void setGlobalCommunicator(MPI_Comm gComm);
    MPI_Comm getGlobalCommunicator() const;

    void setGlobalRank(int global_rank);
    int getGlobalRank() const;

    void open();
    void close();

    void ensureFreeSpace(long numbytes);
    void endOfFile();


    template<class T>
    void write(const T rhs);

    friend SIONFile &operator<<(SIONFile &sf, const LibUtilities::BasisType &rhs);
    friend SIONFile &operator<<(SIONFile &sf, const unsigned int &rhs);
    friend SIONFile &operator<<(SIONFile &sf, const size_t &rhs);
    friend SIONFile &operator<<(SIONFile &sf, const double &rhs);

    void write_string(const std::string &rhs);
    
    friend SIONFile &operator<<(SIONFile &sf, const std::string &rhs);

    template<class T>
    void write_vector(const std::vector<T> &rhs);

    friend SIONFile &operator<<(SIONFile &sf, const std::vector<LibUtilities::BasisType> &rhs);
    friend SIONFile &operator<<(SIONFile &sf, const std::vector<unsigned int> &rhs);
    friend SIONFile &operator<<(SIONFile &sf, const std::vector<double> &rhs);
    friend SIONFile &operator<<(SIONFile &sf, const std::vector<unsigned char> &rhs);

    
    template<class T>
    void read(T *rhs);

    friend SIONFile &operator>>(SIONFile &sf, LibUtilities::BasisType &rhs);
    friend SIONFile &operator>>(SIONFile &sf, unsigned int &rhs);
    friend SIONFile &operator>>(SIONFile &sf, size_t &rhs);
    friend SIONFile &operator>>(SIONFile &sf, double &rhs);

    void read_string(std::string &rhs);
    
    friend SIONFile &operator>>(SIONFile &sf, std::string &rhs);

    template<class T>
    void read_vector(std::vector<T> &rhs);

    friend SIONFile &operator>>(SIONFile &sf, std::vector<LibUtilities::BasisType> &rhs);
    friend SIONFile &operator>>(SIONFile &sf, std::vector<unsigned int> &rhs);
    friend SIONFile &operator>>(SIONFile &sf, std::vector<double> &rhs);

    
private:
    char *_new_sion_file_name;
    MPI_Comm _g_comm;
    MPI_Comm _l_comm;
    int _global_rank;
};

 
}
}
}

#endif
