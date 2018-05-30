///////////////////////////////////////////////////////////////////////////////
//
// File FieldIOSIONlib.h
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
// Description: Field IO to/from SIONlib files.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_FIELDIOSIONLIB_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_FIELDIOSIONLIB_H

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/FileSystem.h>

namespace Nektar
{
namespace LibUtilities
{

/**
 * @class Class encapsulating simple SIONlib data source.
 */
class SIONlibDataSource : public DataSource
{
public:
    /// Constructor based on filename.
    SIONlibDataSource(const std::string &fn)
    {
    }

    std::string Get()
    {
        return fn;
    }

    const std::string Get() const
    {
        return fn;
    }

    /// Static constructor for this data source.
    static DataSourceSharedPtr create(const std::string &fn)
    {
        return DataSourceSharedPtr(new SIONlibDataSource(fn));
    }

private:
    /// SION document.
    std::string fn;
};
typedef std::shared_ptr<SIONlibDataSource> SIONlibDataSourceSharedPtr;

 
class FieldIOSIONlib : public FieldIO
{
public:
    static const unsigned int FORMAT_VERSION;
    
    /// Creates an instance of this class
    LIB_UTILITIES_EXPORT static FieldIOSharedPtr create(
        LibUtilities::CommSharedPtr pComm, bool sharedFilesystem)
    {
        return MemoryManager<FieldIOSIONlib>::AllocateSharedPtr(pComm,
                                                             sharedFilesystem);
    }

    /// Name of class
    LIB_UTILITIES_EXPORT static std::string className;

    LIB_UTILITIES_EXPORT FieldIOSIONlib(
        LibUtilities::CommSharedPtr pComm,
        bool sharedFilesystem);

    LIB_UTILITIES_EXPORT virtual ~FieldIOSIONlib()
    {
    }
    
    /// Get class name
    inline virtual const std::string &GetClassName() const
    {
        return className;
    }

private:
    static const std::string ATTRNAME_FIELDS;
    static const std::string ATTRNAME_BASIS;
    static const std::string ATTRNAME_SHAPE;
    static const std::string ATTRNAME_HOMOLENS;
    static const std::string ATTRNAME_NUMMODES;
    static const std::string ATTRNAME_POINTSTYPE;
    static const std::string ATTRNAME_NUMPOINTS;
    
    static const std::string ATTRVALUE_MIXORDER;
    static const std::string ATTRVALUE_UNIORDER;

    static sion_int32 block_size;
    static sion_int64 chunk_size;
    static std::string read_mode;
    static std::string write_mode;
    
    SIONlib::SIONFile *OpenFileForWriting(const std::string &outFilename);

    LIB_UTILITIES_EXPORT virtual void v_Init(
        const LibUtilities::SessionReaderSharedPtr session);

    LIB_UTILITIES_EXPORT virtual void v_InitFromBenchmarker(
        const LibUtilities::IOSettingsSharedPtr iosettings);


    void WriteRawDataToBinaryVector(std::vector<NekByte> &v, const NekByte* dptr, const size_t dsize);
    void WriteStringToBinaryVector(std::vector<NekByte> &v, const std::string &s);

    template<class T>   
    void WriteDatumToBinaryVector(std::vector<NekByte> &v, const T d);

    template<class T>   
    void WriteDataToBinaryVector(std::vector<NekByte> &v, const std::vector<T> &d);
 
    LIB_UTILITIES_EXPORT virtual unsigned long v_Write(
        const std::string &outFileSuffix,
        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        std::vector<std::vector<NekDouble> > &fielddata,
        const FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap,
        const bool backup = false);

    
    SIONlib::SIONFile *OpenFileForReading(const std::string &inFilename);

    unsigned long DirectImport(SIONlib::SIONFile &f,
        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        std::vector<std::vector<NekDouble> > &fielddata,
        FieldMetaDataMap &fieldinfomap,
        const Array<OneD, int> &ElementiDs);
        
    LIB_UTILITIES_EXPORT virtual unsigned long v_Import(const std::string &infilename,
        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        std::vector<std::vector<NekDouble> > &fielddata = NullVectorNekDoubleVector,
        FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap,
        const Array<OneD, int> &ElementiDs = NullInt1DArray);

    
    LIB_UTILITIES_EXPORT virtual DataSourceSharedPtr v_ImportFieldMetaData(
        const std::string &filename, FieldMetaDataMap &fieldmetadatamap);
};

}
}

#endif
