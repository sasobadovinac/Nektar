///////////////////////////////////////////////////////////////////////////////
//
// File FieldIOHdf5.h
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
// Description: Field IO to/from HDF5
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_FIELDIOHDF5_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_FIELDIOHDF5_H

#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicUtils/H5.h>

namespace Nektar
{
namespace LibUtilities
{

namespace H5
{
class Group;
typedef std::shared_ptr<Group> GroupSharedPtr;
}

/**
 * @class Class encapsulating simple HDF5 data source using H5 reader utilities.
 */
class H5DataSource : public DataSource
{
public:
    /// Constructor based on filename.
    H5DataSource(const std::string &fn, H5::PListSharedPtr parallelProps)
        : doc(H5::File::Open(fn, H5F_ACC_RDONLY, parallelProps))
    {
    }

    /// Get H5::FileSharedPtr reference to file.
    H5::FileSharedPtr Get()
    {
        return doc;
    }

    /// Get H5::FileSharedPtr reference to file.
    const H5::FileSharedPtr Get() const
    {
        return doc;
    }

    /// Static constructor for this data source.
    static DataSourceSharedPtr create(
        const std::string &fn, H5::PListSharedPtr parallelProps)
    {
        return DataSourceSharedPtr(new H5DataSource(fn, parallelProps));
    }

private:
    /// HDF5 document.
    H5::FileSharedPtr doc;
};
typedef std::shared_ptr<H5DataSource> H5DataSourceSharedPtr;

/**
 * @class Simple class for writing hierarchical data using HDF5.
 */
class H5TagWriter : public TagWriter
{
public:
    /// Default constructor.
    H5TagWriter(H5::GroupSharedPtr grp) : m_Group(grp) {}

    /// Add a child node.
    TagWriterSharedPtr AddChild(const std::string &name)
    {
        H5::GroupSharedPtr child = m_Group->CreateGroup(name);
        return TagWriterSharedPtr(new H5TagWriter(child));
    }

    /// Set an attribute key/value pair on this tag.
    void SetAttr(const std::string &key, const std::string &val)
    {
        m_Group->SetAttribute(key, val);
    }

private:
    /// HDF5 group for this tag.
    H5::GroupSharedPtr m_Group;
};
typedef std::shared_ptr<H5TagWriter> H5TagWriterSharedPtr;

/**
 * @class Class for operating on HDF5-based FLD files.
 *
 * This class implements a HDF5 reader/writer based on MPI/O that is designed to
 * operate on a single file across all processors of a simulation. The
 * definition follows vaguely similar lines to XML output but is stored somewhat
 * differently to accommodate parallel reading and writing. At a basic level
 * metadata is organised as follows:
 *
 *   - Nektar++ data lies in the root `/NEKTAR` group.
 *   - The contents of a FieldDefinitions object is hashed to construct a unique
 *     identifier for each object, which becomes the name of a group within the
 *     root group. We then use the H5TagWriter to assign the field definitions
 *     to each group.
 *   - In a similar fashion, we create a `Metadata` group to contain field
 *     metadata that is written.
 *
 * We then define five data sets to contain field data:
 *
 *   - The `DATA` dataset contains the double-precision modal coefficient data.
 *   - The `IDS` dataset contains the element IDs of the elements that are
 *     written out.
 *   - The `POLYORDERS` dataset is written if the field data contains variable
 *     polynomial order, and contains the (possibly hetergeneous) mode orders in
 *     each direction for each of the elements.
 *   - The `HOMOGENEOUSZIDS` dataset contains the IDs of z-planes for
 *     homogeneous simulations, if the data are homogeneous.
 *   - The `HOMOGENEOUSYIDS` dataset contains the IDs of y-planes for
 *     homogeneous simulations, if the data are homogeneous.
 *   - The `HOMOGENEOUSSIDS` dataset contains the strip IDs for
 *     homogeneous simulations, if the data are homogeneous and use strips.
 *
 * The ordering is defined according to the `DECOMPOSITION` dataset. A
 * `decomposition' in this class is essentially a single field definition with
 * its accompanying data. Data are written into each dataset by the order of
 * each decomposition. Each decomposition contains the following seven integers
 * that define it per field definition per processor:
 *
 *   - Number of elements in this field definition (index #ELEM_DCMP_IDX).
 *   - Number of entries in the `DATA` array for this field definition
 *     (index #VAL_DCMP_IDX)
 *   - Number of entries in the `POLYORDERS` array for this field definition
 *     (index #ORDER_DCMP_IDX)
 *   - Number of entries in the `HOMOGENEOUSZIDS` array (index #HOMZ_DCMP_IDX).
 *   - Number of entries in the `HOMOGENEOUSYIDS` array (index #HOMY_DCMP_IDX).
 *   - Number of entries in the `HOMOGENEOUSSIDS` array (index #HOMS_DCMP_IDX).
 *   - Hash of the field definition, represented as a 32-bit integer, which
 *     describes the name of the attribute that contains the rest of the field
 *     definition information (e.g. field names, basis type, etc).
 *
 * The number of decompositions is therefore calculated as the field size
 * divided by #MAX_DCMPS which allows us to calculate the offsets of the data
 * for each field definition within the arrays.
 */
class FieldIOHdf5 : public FieldIO
{
public:
    static const unsigned int FORMAT_VERSION;

    static const unsigned int FIELD_COUNT_IDS;
    static const unsigned int FIELD_COUNT_DATA;
    static const unsigned int FIELD_COUNT_HOMY;
    static const unsigned int FIELD_COUNT_HOMZ;
    static const unsigned int FIELD_COUNT_HOMS;
    static const unsigned int FIELD_COUNT_ORDER;
    static const unsigned int FIELD_COUNT_SIZE;

    static const unsigned int FIELD_DECOMP_OFF;
    static const unsigned int FIELD_DECOMP_CNT;
    static const unsigned int FIELD_DECOMP_SIZE;

    static const unsigned int DATA_DECOMP_FIELD_HASH;
    static const unsigned int DATA_DECOMP_IDS_OFF;
    static const unsigned int DATA_DECOMP_IDS_CNT;
    static const unsigned int DATA_DECOMP_DATA_OFF;
    static const unsigned int DATA_DECOMP_DATA_CNT;
    static const unsigned int DATA_DECOMP_HOMY_OFF;
    static const unsigned int DATA_DECOMP_HOMY_CNT;
    static const unsigned int DATA_DECOMP_HOMZ_OFF;
    static const unsigned int DATA_DECOMP_HOMZ_CNT;
    static const unsigned int DATA_DECOMP_HOMS_OFF;
    static const unsigned int DATA_DECOMP_HOMS_CNT;
    static const unsigned int DATA_DECOMP_ORDER_OFF;
    static const unsigned int DATA_DECOMP_ORDER_CNT;
    static const unsigned int DATA_DECOMP_SIZE;
      
    
    /// Creates an instance of this class
    LIB_UTILITIES_EXPORT static FieldIOSharedPtr create(
        LibUtilities::CommSharedPtr pComm, bool sharedFilesystem)
    {
        return MemoryManager<FieldIOHdf5>::AllocateSharedPtr(pComm,
                                                             sharedFilesystem);
    }

    /// Name of class
    LIB_UTILITIES_EXPORT static std::string className;

    LIB_UTILITIES_EXPORT FieldIOHdf5(
        LibUtilities::CommSharedPtr pComm,
        bool sharedFilesystem);

    LIB_UTILITIES_EXPORT virtual ~FieldIOHdf5()
    {
    }

    /// Get class name
    inline virtual const std::string &GetClassName() const
    {
        return className;
    }

private:

    bool reformatting;
    
    LIB_UTILITIES_EXPORT virtual void v_Init(
        const LibUtilities::SessionReaderSharedPtr session);

    LIB_UTILITIES_EXPORT virtual void v_InitFromBenchmarker(
        const LibUtilities::IOSettingsSharedPtr iosettings);
    
    LIB_UTILITIES_EXPORT virtual uint64_t v_Write(
        const std::string &outFilename,
        std::vector<FieldDefinitionsSharedPtr> &fieldDefs,
        std::vector<std::vector<NekDouble> > &fieldData,
        const FieldMetaDataMap &fieldMetaDataMap = NullFieldMetaDataMap,
        const bool backup = false);

    LIB_UTILITIES_EXPORT uint64_t CreateDataSets(const std::string &outFilename,
        const int rkFormatter,
        const std::size_t nTotFields,
        std::vector<uint64_t> &allFieldHashes,
        std::vector<uint64_t> &allFieldCounts,
        std::vector<uint64_t> &allFieldDecomps,
	std::vector<uint64_t> &allFirstDataDecomps,
        std::vector<uint64_t> &allDataDecomps,
        std::vector<FieldDefinitionsSharedPtr> &fieldDefs,
        std::vector<std::vector<unsigned int> > &homoYIDs,
        std::vector<std::vector<unsigned int> > &homoZIDs,
        std::vector<std::vector<unsigned int> > &homoSIDs,
        std::vector<std::vector<unsigned int> > &numModesPerDirVar,
	std::vector<std::vector<NekDouble> > &fieldData,
        const FieldMetaDataMap &fieldMetaDataMap);
    LIB_UTILITIES_EXPORT uint64_t WriteDecompositionData(const std::string &outFilename,
        const int rkFormatter,
        std::vector<uint64_t> &allFieldDecomps,
	std::vector<uint64_t> &allDataDecomps);
    LIB_UTILITIES_EXPORT uint64_t WriteFieldAttributes(const std::string &outFilename,
        const uint64_t nFields,
	const int varOrder,
        std::vector<uint64_t> &writeFieldHashes,
        std::vector<FieldDefinitionsSharedPtr> &fieldDefs,
	std::vector<std::string> &fieldNames,
	std::vector<std::string> &shapeStrings,
	std::vector<std::vector<NekDouble> > &homoLengths,
	std::vector<std::string> &numModesPerDirUni);
    LIB_UTILITIES_EXPORT uint64_t WriteData(
        const std::string &outFilename,
	const std::size_t nFields,
	std::vector<uint64_t> &totalFieldCounts,
	std::vector<uint64_t> &firstDataDecomps,
        std::vector<FieldDefinitionsSharedPtr> &fieldDefs,
	std::vector<std::vector<unsigned int> > &homoYIDs,
        std::vector<std::vector<unsigned int> > &homoZIDs,
        std::vector<std::vector<unsigned int> > &homoSIDs,
        std::vector<std::vector<unsigned int> > &numModesPerDirVar,
        std::vector<std::vector<NekDouble> > &fieldData);
      
    template <class T>
    LIB_UTILITIES_EXPORT uint64_t WriteFieldDataInd(std::size_t nFields,
        H5::DataSpaceSharedPtr &space, H5::DataSetSharedPtr &dset,
        uint64_t data_i, std::vector<std::vector<T> > &data);
    
    template <class T>
    LIB_UTILITIES_EXPORT uint64_t WriteFieldData(std::size_t nMinFields, std::size_t nFields,
        H5::DataSpaceSharedPtr &space, H5::DataSetSharedPtr &dset,
        uint64_t data_i, std::vector<std::vector<T> > &data);

    LIB_UTILITIES_EXPORT virtual uint64_t v_Import(const std::string &inFilename,
        std::vector<FieldDefinitionsSharedPtr> &fieldDefs,
        std::vector<std::vector<NekDouble> > &fieldData = NullVectorNekDoubleVector,
        FieldMetaDataMap &fieldMetaInfoMap = NullFieldMetaDataMap,
        const Array<OneD, int> &elementIDs = NullInt1DArray);
    
    template <class T>
    LIB_UTILITIES_EXPORT uint64_t ImportFieldDataInd(H5::GroupSharedPtr root,
        std::string dsetName, std::string dataTag, uint64_t nDataItems,
        uint64_t offset, std::vector<T> &data);

    LIB_UTILITIES_EXPORT uint64_t ImportFieldDef(H5::GroupSharedPtr root,
        std::string               group,
        FieldDefinitionsSharedPtr def);
    
    LIB_UTILITIES_EXPORT virtual DataSourceSharedPtr v_ImportFieldMetaData(
        const std::string &filename, FieldMetaDataMap &fieldmetadatamap);

    LIB_UTILITIES_EXPORT void ImportHDF5FieldMetaData(
        DataSourceSharedPtr dataSource, FieldMetaDataMap &fieldmetadatamap);

    
};
}
}
#endif
