////////////////////////////////////////////////////////////////////////////////
//
//  File: FieldIOBenchmarker.cpp
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
//  Description: Measure the performance of FieldIO XML, HDF5 and SIONlib classes.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/bimap.hpp>
#include <LibUtilities/BasicUtils/H5.h>
#include <LibUtilities/BasicUtils/SIONlib.h>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/FieldIOXml.h>
#include <LibUtilities/BasicUtils/FieldIOHdf5.h>
#include <LibUtilities/BasicUtils/FieldIOSIONlib.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/Communication/CommMpi.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

// Below, we'd like to use an unordered set for its faster lookup performance
// However this is only available if C++11 is.
//
#if __cplusplus >= 201103L
#include <unordered_set>
typedef std::unordered_set<int> IntSet;
#else
#include <set>
typedef std::set<int> IntSet;
#endif

using namespace Nektar;
using namespace LibUtilities;

namespace po = boost::program_options;
namespace berrc = boost::system::errc;

typedef std::vector<FieldDefinitionsSharedPtr> DefVec;
typedef std::vector<std::vector<NekDouble> > DatVec;

typedef boost::shared_ptr<FieldIOXml> FieldIOXmlSharedPtr;
typedef boost::shared_ptr<FieldIOHdf5> FieldIOHdf5SharedPtr;
typedef boost::shared_ptr<FieldIOSIONlib> FieldIOSIONlibSharedPtr;

struct Experiment
{
    /// Test read (false) or write (true)
    bool write;
    /// Test the HDF5 (true) or XML/SIONlib (false) reader/writer.
    bool hdf;
    /// Test the SIONlib (true) or XML/HDF5 (false) reader/writer.
    bool sionlib;
    /// If true, print additional debugging information.
    bool verbose;
    /// Number of writes to perform.
    int n;
    /// Input file source name
    std::string dataSource;
    /// Output filename
    std::string dataDest;
    /// Communicator for the experiment.
    CommSharedPtr comm;
};

IOSettingsSharedPtr hdf5_settings(new IOSettings);
IOSettingsSharedPtr sionlib_settings(new IOSettings);

typedef std::vector<double> Results;

const double MB = 1000000.0;
const int AORTIC_ARCH = 0;
const int RACING_CAR = 1;

std::string SIONlib_GetIOBlocksPerChunk(int nprocs, int testid);
Results TestRead(Experiment &exp);
Results TestWrite(Experiment &exp);


int main(int argc, char *argv[])
{
    Experiment exp;
    exp.write   = false;
    exp.hdf     = false;
    exp.sionlib = false;
    exp.verbose = false;
    exp.n       = 3;
    exp.comm    = GetCommFactory().CreateInstance("ParallelMPI", argc, argv);

    hdf5_settings->insert( IOSettings::value_type("Reformatting", "1") );
    
    sionlib_settings->insert( IOSettings::value_type("IOBlockSize", "65536") );
    sionlib_settings->insert( IOSettings::value_type("IOBlocksPerChunk", SIONlib_GetIOBlocksPerChunk(exp.comm->GetSize(), AORTIC_ARCH)) );
    sionlib_settings->insert( IOSettings::value_type("IOReadMode", "br,collective,collsize=-1") );
    sionlib_settings->insert( IOSettings::value_type("IOWriteMode", "bw,collective,collectivemerge,collsize=-1") );
        
    po::options_description desc("Available options");
    desc.add_options()("help,h", "Produce this help message.")(
        "mode,m", po::value<char>(),
        "Choose r[ead] (default), x[ml write] or h[df5 write] or s[sionlib write]")(
        "number,n", po::value<unsigned>(),
        "Number of iterations to perform, default 3")("verbose,v",
                                                      "Enable verbose mode.");

    po::options_description hidden("Hidden options");
    hidden.add_options()("input-file", po::value<std::string>(),
                         "Input filename")("output-file",
                                           po::value<std::string>());

    po::options_description cmdline_options;
    cmdline_options.add(hidden).add(desc);

    po::options_description visible("Allowed options");
    visible.add(desc);

    po::positional_options_description p;
    p.add("input-file", 1).add("output-file", 1);

    po::variables_map vm;

    try
    {
        po::store(po::command_line_parser(argc, argv)
                      .options(cmdline_options)
                      .positional(p)
                      .run(),
                  vm);
        po::notify(vm);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << desc;
        return 1;
    }

    if (vm.count("help") || vm.count("input-file") != 1)
    {
        std::cerr
            << "Usage: FieldIOBenchmarker [options] inputfile [outputfile]"
            << endl;
        std::cout << desc;
        std::cout << endl;
        return 1;
    }

    ASSERTL0(vm.count("input-file"), "Must specify input.");

    exp.dataSource = vm["input-file"].as<std::string>();

    if (vm.count("verbose") && exp.comm->GetRank() == 0)
    {
        exp.verbose = true;
    }

    if (vm.count("number"))
    {
        exp.n = vm["number"].as<unsigned>();
    }

    if (vm.count("mode"))
    {
        char mode = vm["mode"].as<char>();
        switch (mode)
        {
            case 'r':
                exp.write   = false;
                break;
            case 'x':
                exp.write   = true;
                exp.hdf     = false;
                exp.sionlib = false;
                break;
            case 'h':
                exp.write   = true;
                exp.hdf     = true;
                exp.sionlib = false;
                break;
            case 's':
                exp.write   = true;
                exp.hdf     = false;
                exp.sionlib = true;
                break;
            default:
                std::cout << "Unrecognised mode: " << mode << std::endl;
                std::cout << desc << endl;
                return 1;
                break;
        }
    }

    Results res;
    if (exp.write)
    {
        if (vm.count("output-file"))
        {
            exp.dataDest = vm["output-file"].as<std::string>();
        }
        else
        {
            exp.dataDest = exp.dataSource + ".tmp";
        }

        res = TestWrite(exp);
    }
    else
    {
        res = TestRead(exp);
    }

    exp.comm->Finalise();
}


/**
 * @brief Return the number of file system blocks per SIONlib chunk.
 *
 * This result depends on the number of processes and on the size of the checkpoint
 * file associated with the test case.
 *
 * @param nprocs  The number of MPI processes.
 * @param testid  The test case --- this determines the size of the checkpoint file.
 *
 * @return string.
 */
std::string SIONlib_GetIOBlocksPerChunk(int nprocs, int testid)
{
    std::string bpc = "1.0";
  
    if (AORTIC_ARCH == testid)
    {
        switch (nprocs)
        {
            case 48:
                bpc = "8.0";
                break;
            case 96:
                bpc = "4.0";
                break;
            case 192:
                bpc = "2.0";
                break;
        }
    }
    else if (RACING_CAR == testid)
    {
        switch (nprocs)
        {
            case 768:
                bpc = "110.0";
                break;   
            case 1536:
                bpc = "55.0";
                break;
            case 3072:
                bpc = "29.0";
                break;
            case 6144:
	        bpc = "16.0";
                break;
        }
    }

    return bpc;
}
    
/**
 * @brief Read elemental IDs from XML field file format for this rank.
 *
 * Here we read the Info.xml to figure out in which file the elements are
 * stored.  The elements are then divided amongst the ranks (in a trivial
 * decomposition based on order in the Info.xml).
 *
 * @param exp  Experiment to be run.
 * @param fio  FieldIO object to perform test on.
 *
 * @return Array containing element IDs for this rank.
 */
Array<OneD, int> ReadIDsForThisRank(Experiment &exp, FieldIOSharedPtr fio)
{
    std::vector<std::string> fileNames;
    std::vector<std::vector<unsigned int> > elementList;
    FieldMetaDataMap fieldmetadatamap;

    std::string infoFile = exp.dataSource + "/Info.xml";

    std::shared_ptr<FieldIOXml> fioXml =
        std::dynamic_pointer_cast<FieldIOXml>(fio);
    fioXml->ImportMultiFldFileIDs(infoFile, fileNames, elementList,
                                  fieldmetadatamap);

    unsigned totalEls = 0;
    std::vector<unsigned> elStartFile(elementList.size(), 0);
    std::vector<unsigned> elStopFile(elementList.size(), 0);

    for (unsigned i = 0; i < elementList.size(); ++i)
    {
        elStartFile[i] = totalEls;
        totalEls += elementList[i].size();
        elStopFile[i] = totalEls;
    }
    NekDouble elemPerNode = totalEls / (NekDouble)exp.comm->GetSize();
    unsigned elStart      = elemPerNode * exp.comm->GetRank();
    unsigned elStop       = elemPerNode * (exp.comm->GetRank() + 1);
    unsigned nEls         = elStop - elStart;

    Array<OneD, int> ElementIDs(nEls);

    for (unsigned iFile = 0, iEl = elStart; iEl < elStop;)
    {
        // Find the index of the file that contains the element index we want
        while (!(elStartFile[iFile] <= iEl && iEl < elStopFile[iFile]))
            iFile++;

        unsigned startInFile = iEl - elStartFile[iFile];
        unsigned stopInFile;

        // Determine how much of the file we want
        if (elStop > elStopFile[iFile])
        {
            // Need some of the next one too
            // Copy to the end
            stopInFile = elStopFile[iFile] - elStartFile[iFile];
        }
        else
        {
            // This is the last file we need
            // Copy up to elStop
            stopInFile = elStop - elStartFile[iFile];
        }

        // Copy the chunk
        std::memcpy(&ElementIDs[iEl - elStart],
                    &elementList[iFile][startInFile],
                    (stopInFile - startInFile) * sizeof(int));

        iEl += stopInFile - startInFile;
    }
    return ElementIDs;
}

/**
 * @brief  Get the IDs of the elements managed by this rank.
 *
 * @param exp  Experimental setup.
 *
 * @return IDs of the elements managed by this rank.
 */
Array<OneD, int> ReadIDsForThisRank(Experiment& exp)
{
    std::string ft = FieldIO::GetFileType(exp.dataSource, exp.comm);
    FieldIOSharedPtr fio =
        GetFieldIOFactory().CreateInstance(ft, exp.comm, true);

    if (fs::is_directory(exp.dataSource))
    {
        return ReadIDsForThisRank(exp, fio);
    }
    else if (exp.hdf)
    {
        // ensure auto_selective
        Array<OneD, int> ElementIDs;
        return ElementIDs;
    }

    return NullInt1DArray;
}

/**
 * @brief Test read speed for this experimental setup.
 *
 * This routine performs Experiment::n reads, timing the read times and
 * returning a Results struct containing the results.
 *
 * @param exp  Experimental setup.
 *
 * @return Resulting timings.
 */
Results TestRead(Experiment &exp)
{
    const std::string ft = FieldIO::GetFileType(exp.dataSource, exp.comm);
    if (exp.verbose)
    {
        std::cout << "Beginning read experiment..."
                  << exp.n << " iterations." << std::endl;
        std::cout << "Determining file type... ";
        std::cout << ft << std::endl;
    }

    if (0 == ft.compare("Hdf5"))
    {
        exp.hdf = true;
    }
    else if (0 == ft.compare("SIONlib"))
    {
        exp.sionlib = true;
    }

    std::vector<FieldDefinitionsSharedPtr> fielddefs;
    std::vector < std::vector<NekDouble> > fielddata;
    FieldMetaDataMap fieldmetadatamap;
    FieldIOSharedPtr fio = GetFieldIOFactory().CreateInstance(
            ft, exp.comm, true);

    if (exp.sionlib)
    {
        fio->InitFromBenchmarker(sionlib_settings);
    }
    else if (exp.hdf)
    {
        fio->InitFromBenchmarker(hdf5_settings);
    }
    
    Array<OneD, int> ElementIDs = ReadIDsForThisRank(exp);
        
    Results res(exp.n, 0.0);
    for (unsigned i = 0; i < exp.n; ++i)
    {
        if (exp.verbose)
        {
            std::cout << "Test " << i+1 << " of " << exp.n;
        }
                                
        double t0 = MPI_Wtime();
        unsigned long nRead = fio->Import(exp.dataSource,
            fielddefs, fielddata,
            fieldmetadatamap, ElementIDs);        
        double t1 = MPI_Wtime() - t0;

        exp.comm->AllReduce(nRead, LibUtilities::ReduceSum);

        if (exp.verbose)
        {
            std::cout << ": read " << nRead/MB << " MB in "
                      << t1 << " s." << std::endl;
        }
        
        res[i] = (nRead/MB) / t1;
    }
    return res;
}

/**
 * @brief Test write speed for this experimental setup.
 *
 * This routine performs Experiment::n writes, timing the write times and
 * returning a Results struct containing the results.
 *
 * @param exp  Experimental setup.
 *
 * @return Resulting timings.
 */
Results TestWrite(Experiment &exp)
{
    const std::string ft = FieldIO::GetFileType(exp.dataSource, exp.comm);
    if (exp.verbose)
    {
        std::cout << "Beginning write experiment..."
                  << exp.n << " iterations." << std::endl;
        std::cout << "Determining file type... ";
        std::cout << ft << std::endl;
    }
    
    std::vector<FieldDefinitionsSharedPtr> fielddefs;
    std::vector < std::vector<NekDouble> > fielddata;
    FieldMetaDataMap fieldmetadatamap;
    FieldIOSharedPtr fio = GetFieldIOFactory().CreateInstance(
            ft, exp.comm, true);

    if (exp.sionlib)
    {
        fio->InitFromBenchmarker(sionlib_settings);
    }
    else if (exp.hdf)
    {
        fio->InitFromBenchmarker(hdf5_settings);
    }
    
    Array<OneD, int> ElementIDs = ReadIDsForThisRank(exp);
    fio->Import(exp.dataSource, fielddefs, fielddata,
        fieldmetadatamap, ElementIDs);

    exp.comm->Block();

    Results res(exp.n, 0);
    fs::path specPath(exp.dataDest);

    for (unsigned i = 0; i < exp.n; ++i)
    {
        if (exp.verbose)
        {
            std::cout << "Test " << i+1 << " of " << exp.n;
        }
        
        double t0 = MPI_Wtime();
        unsigned long nWritten = fio->Write(exp.dataDest, fielddefs, fielddata);
        double t1 = MPI_Wtime() - t0;

        exp.comm->AllReduce(nWritten, LibUtilities::ReduceSum);
         
        if (exp.verbose)
        {
            std::cout << ": written " << nWritten/MB << " MB in "
                      << t1 << " s."  << std::endl;
        }
        
        res[i] = (nWritten/MB) / t1;
        /*
        if (exp.comm->GetRank() == 0)
        {
            try
            {
                // remove any files that have just been written
                // as part of test
                fs::remove_all(specPath);  
            }
            catch (fs::filesystem_error& e)
            {
                ASSERTL0(
                    e.code().value()
                    == berrc::no_such_file_or_directory,
                    "Filesystem error: " + string(e.what()));
            }
        }
        */
        exp.comm->Block();
    }
    return res;
}


