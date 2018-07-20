////////////////////////////////////////////////////////////////////////////////
//
//  File: FieldIOHdf5.cpp
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
//  Description: I/O routines relating to Fields into HDF
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/FieldIOHdf5.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>

#include <unordered_set>
#include <functional>

namespace berrc = boost::system::errc;

namespace Nektar
{
namespace LibUtilities
{
namespace H5
{

template <> inline DataTypeSharedPtr DataTypeTraits<BasisType>::GetType()
{
    return PredefinedDataType::Native<int>();
}

}

std::string FieldIOHdf5::className =
    GetFieldIOFactory().RegisterCreatorFunction(
        "Hdf5", FieldIOHdf5::create, "HDF5-based output of field data.");

/// Version of the Nektar++ HDF5 format, which is embedded into the main NEKTAR
/// group as an attribute.
const unsigned int FieldIOHdf5::FORMAT_VERSION = 1;

// The following definitions allow us to consistently refer to indexes pulled
// out of the various datasets.

const unsigned int FieldIOHdf5::FIELD_COUNT_IDS   = 0;
const unsigned int FieldIOHdf5::FIELD_COUNT_DATA  = 1;
const unsigned int FieldIOHdf5::FIELD_COUNT_HOMY  = 2;
const unsigned int FieldIOHdf5::FIELD_COUNT_HOMZ  = 3;
const unsigned int FieldIOHdf5::FIELD_COUNT_HOMS  = 4;
const unsigned int FieldIOHdf5::FIELD_COUNT_ORDER = 5;
const unsigned int FieldIOHdf5::FIELD_COUNT_SIZE  = 6;

const unsigned int FieldIOHdf5::FIELD_DECOMP_OFF  = 0;
const unsigned int FieldIOHdf5::FIELD_DECOMP_CNT  = 1;
const unsigned int FieldIOHdf5::FIELD_DECOMP_SIZE = 2;

const unsigned int FieldIOHdf5::DATA_DECOMP_FIELD_HASH = 0;
const unsigned int FieldIOHdf5::DATA_DECOMP_IDS_OFF    = 1;
const unsigned int FieldIOHdf5::DATA_DECOMP_IDS_CNT    = 2;
const unsigned int FieldIOHdf5::DATA_DECOMP_DATA_OFF   = 3;
const unsigned int FieldIOHdf5::DATA_DECOMP_DATA_CNT   = 4;
const unsigned int FieldIOHdf5::DATA_DECOMP_HOMY_OFF   = 5;
const unsigned int FieldIOHdf5::DATA_DECOMP_HOMY_CNT   = 6;
const unsigned int FieldIOHdf5::DATA_DECOMP_HOMZ_OFF   = 7;
const unsigned int FieldIOHdf5::DATA_DECOMP_HOMZ_CNT   = 8;
const unsigned int FieldIOHdf5::DATA_DECOMP_HOMS_OFF   = 9;
const unsigned int FieldIOHdf5::DATA_DECOMP_HOMS_CNT   = 10;
const unsigned int FieldIOHdf5::DATA_DECOMP_ORDER_OFF  = 11;
const unsigned int FieldIOHdf5::DATA_DECOMP_ORDER_CNT  = 12;
const unsigned int FieldIOHdf5::DATA_DECOMP_SIZE       = 13;

  
/**
 * @brief Construct the FieldIO object for HDF5 output.
 *
 * @param pComm              Communicator object.
 * @param sharedFilesystem   True if this system has a shared filesystem.
 */
FieldIOHdf5::FieldIOHdf5(LibUtilities::CommSharedPtr pComm,
                         bool sharedFilesystem)
  : reformatting(true), FieldIO(pComm, sharedFilesystem)
{
}


void FieldIOHdf5::v_Init(const LibUtilities::SessionReaderSharedPtr session)
{
    reformatting = true;
}

void FieldIOHdf5::v_InitFromBenchmarker(const LibUtilities::IOSettingsSharedPtr iosettings)
{
    LibUtilities::IOSettings::iterator it;
    
    reformatting = true;
    it = iosettings->find("Reformatting");
    if (iosettings->end() != it)
    {
        std::istringstream(it->second) >> reformatting;
    }
}


/**
 * @brief Write a HDF5 file to @p outFile given the field definitions @p
 * fielddefs, field data @p fielddata and metadata @p fieldmetadatamap.
 *
 * The writing strategy for HDF5 output is as follows:
 *
 *   - Each rank determines the amount of data needed to be written into each
 *     dataset.
 *   - Each rank communicates its decomposition information to the formatting process.
 *   - The formatting processor initialises the output structure, writes the
 *     decomposition datasets and all the field definition information.
 *   - Other ranks may have field definitions that do not belong to the formatting
 *     process, in which case they open the file and append this (since
 *     attributes cannot be written in parallel).
 *   - Each of the other ranks writes their data contributions to the rest of
 *     the set.
 *
 * @param outFilename       Output filename.
 * @param fieldDefs         Input field definitions.
 * @param fieldData         Input field data.
 * @param fieldMetaDataMap  Field metadata.
 * @return The number of bytes written.
 */
uint64_t FieldIOHdf5::v_Write(const std::string &outFilename,
                              std::vector<FieldDefinitionsSharedPtr> &fieldDefs,
                              std::vector<std::vector<NekDouble> > &fieldData,
                              const FieldMetaDataMap &fieldMetaDataMap,
                              const bool backup)
{
    int nranks = m_comm->GetSize();
    int rk = m_comm->GetRank();
        
    std::stringstream prfx;
    prfx << rk << ": FieldIOHdf5::v_Write(): ";

    /*
    double tm0 = 0.0;
    if (0 == rk)
    {
        tm0 = m_comm->Wtime();
    }
    */

    if (reformatting)
    {
        SetUpOutput(outFilename, false, backup);
    }

    // We make a number of assumptions in this code.
    //   1. All element ids have the same type: unsigned int.
    //   2. All elements within a given field have the same number of values.
    //   3. All element values have the same type, NekDouble.

    ASSERTL1(fieldDefs.size() == fieldData.size(),
        prfx.str() + "fielddefs and fielddata have incompatible lengths.");

    uint64_t nWritten = 0;
    std::size_t nFields = fieldDefs.size();
    int rkFormatter = 0;

    std::vector<uint64_t> fieldHashes(nFields);
    std::vector<uint64_t> fieldCounts(FIELD_COUNT_SIZE*nFields, 0);
    std::vector<uint64_t> totalFieldCounts(FIELD_COUNT_SIZE, 0);

    std::vector<std::string> fieldNames(nFields);
    std::vector<std::string> shapeStrings(nFields);
    std::vector<std::vector<NekDouble> > homoLengths(nFields);
    std::vector<std::vector<unsigned int> > homoSIDs(nFields),
        homoYIDs(nFields), homoZIDs(nFields);
    std::vector<std::vector<unsigned int> > numModesPerDirVar(nFields);
    std::vector<std::string> numModesPerDirUni(nFields);

    int homDim = -1;
    int varOrder = 0;
    for (std::size_t f = 0; f < nFields; ++f)
    {
        if (!fieldDefs[f]->m_uniOrder)
        {
            varOrder = 1;
            break;
        }
    }
    m_comm->AllReduce(varOrder, LibUtilities::ReduceMax);

    // Calculate the total number of elements handled by this MPI process and
    // the total number of bytes required to store the elements. Base the name
    // of each field as the hash of the field definition.
    for (std::size_t f = 0; f < nFields; ++f)
    {
        std::size_t fco = FIELD_COUNT_SIZE*f;
        ASSERTL1(fieldData[f].size() > 0,
            prfx.str() + "fieldData vector must contain at least one value.");
        ASSERTL1(fieldData[f].size() ==
            fieldDefs[f]->m_fields.size() * CheckFieldDefinition(fieldDefs[f]),
            prfx.str() + "fieldData vector has invalid size.");

        std::size_t nElemIds = fieldDefs[f]->m_elementIDs.size();
        fieldCounts[fco + FIELD_COUNT_IDS]  = nElemIds;
        fieldCounts[fco + FIELD_COUNT_DATA] = fieldData[f].size();
	totalFieldCounts[FIELD_COUNT_IDS]  += nElemIds;
	totalFieldCounts[FIELD_COUNT_DATA] += fieldData[f].size();

        // Hash the field specification
        std::stringstream hashStream;
        std::size_t nSubFields = fieldDefs[f]->m_fields.size();
        for (std::size_t sf = 0; sf < nSubFields; ++sf)
        {
            hashStream << fieldDefs[f]->m_fields[sf];
        }

        nSubFields = fieldDefs[f]->m_basis.size();
        for (std::size_t sf = 0; sf < nSubFields; ++sf)
        {
            hashStream << fieldDefs[f]->m_basis[sf];
        }

        // Determine SHAPE attribute
        std::stringstream shapeStringStream;
        shapeStringStream << ShapeTypeMap[fieldDefs[f]->m_shapeType];

        if (fieldDefs[f]->m_numHomogeneousDir > 0)
        {
            if (homDim == -1)
            {
                homDim = fieldDefs[f]->m_numHomogeneousDir;
            }

            ASSERTL1(homDim == fieldDefs[f]->m_numHomogeneousDir,
                     "HDF5 does not support variable homogeneous directions in "
                     "the same file.");

            shapeStringStream << "-HomogenousExp"
                              << fieldDefs[f]->m_numHomogeneousDir << "D";
        }

        if (fieldDefs[f]->m_homoStrips)
        {
            shapeStringStream << "-Strips";
        }

        shapeStrings[f] = shapeStringStream.str();
        hashStream << shapeStringStream.str();

        // Determine HOMOGENEOUS attributes
        if (fieldDefs[f]->m_numHomogeneousDir)
        {
            nSubFields = fieldDefs[f]->m_homogeneousLengths.size();
            homoLengths[f].resize(nSubFields);
            for (std::size_t sf = 0; sf < nSubFields; ++sf)
            {
                NekDouble len = fieldDefs[f]->m_homogeneousLengths[sf];
                hashStream << len;
                homoLengths[f][sf] = len;
            }

            nSubFields = fieldDefs[f]->m_homogeneousYIDs.size();
            if (nSubFields > 0)
            {
                homoYIDs[f].resize(nSubFields);
                fieldCounts[fco + FIELD_COUNT_HOMY] = nSubFields;
		totalFieldCounts[FIELD_COUNT_HOMY] += nSubFields;
                for (std::size_t sf = 0; sf < nSubFields; ++sf)
                {
                    homoYIDs[f][sf] = fieldDefs[f]->m_homogeneousYIDs[sf];
                }
            }

            nSubFields = fieldDefs[f]->m_homogeneousZIDs.size();
            if (nSubFields > 0)
            {
                homoZIDs[f].resize(nSubFields);
                fieldCounts[fco + FIELD_COUNT_HOMZ] = nSubFields;
		totalFieldCounts[FIELD_COUNT_HOMZ] += nSubFields;
                for (std::size_t sf = 0; sf < nSubFields; ++sf)
                {
                    homoZIDs[f][sf] = fieldDefs[f]->m_homogeneousZIDs[sf];
                }
            }

            nSubFields = fieldDefs[f]->m_homogeneousSIDs.size();
            if (nSubFields > 0)
            {
                homoSIDs[f].resize(nSubFields);
                fieldCounts[fco + FIELD_COUNT_HOMS] = nSubFields;
		totalFieldCounts[FIELD_COUNT_HOMS] += nSubFields;
                for (std::size_t sf = 0; sf < nSubFields; ++sf)
                {
                    homoSIDs[f][sf] = fieldDefs[f]->m_homogeneousSIDs[sf];
                }
            }
        }

        if (fieldDefs[f]->m_uniOrder)
        {
            std::vector<unsigned int> elemModes(fieldDefs[f]->m_basis.size());

            for (std::vector<int>::size_type i = 0;
                 i < fieldDefs[f]->m_basis.size(); ++i)
            {
                elemModes[i] = fieldDefs[f]->m_numModes[i];
            }

            if (varOrder)
            {
                for (std::vector<int>::size_type i = 0; i < nElemIds; ++i)
                {
                    std::copy(elemModes.begin(), elemModes.end(),
                              std::back_inserter(numModesPerDirVar[f]));
                }
                fieldCounts[fco + FIELD_COUNT_ORDER] = nElemIds * elemModes.size();
		totalFieldCounts[FIELD_COUNT_ORDER] += fieldCounts[fco + FIELD_COUNT_ORDER];
            }
            else
            {
                stringstream numModesStringStream;
                numModesStringStream << "UNIORDER:";
                for (std::vector<int>::size_type i = 0;
                     i < elemModes.size(); i++)
                {
                    if (i > 0)
                    {
                        numModesStringStream << ",";
                    }
                    numModesStringStream << elemModes[i];
                }

                numModesPerDirUni[f] = numModesStringStream.str();
                hashStream << numModesPerDirUni[f];
            }
        }
        else
        {
            numModesPerDirVar[f] = fieldDefs[f]->m_numModes;
            fieldCounts[fco + FIELD_COUNT_ORDER] = fieldDefs[f]->m_numModes.size();
	    totalFieldCounts[FIELD_COUNT_ORDER] += fieldCounts[fco + FIELD_COUNT_ORDER];
        }

        std::hash<std::string> stringHasher;
        std::stringstream fieldNameStream;
        uint64_t fieldDefHash = stringHasher(hashStream.str());

	if (reformatting)
	{
            fieldHashes[f] = fieldDefHash;
	}

        fieldNameStream << fieldDefHash;
        fieldNames[f] = fieldNameStream.str();
	
    }  // end of <for (std::size_t f = 0; f < nFields; ++f)> loop

    std::vector<uint64_t> firstDataDecomps;
    
    if (reformatting)
    {
        std::size_t nTotFields = 0;
        std::vector<int> sizeMap, offsetMap;
	std::vector<int> sizeMap2, offsetMap2;
	std::vector<int> vnFields(1, nFields);
        std::vector<int> allnFields = m_comm->Gather(rkFormatter, vnFields);

	if (rkFormatter == rk)
	{
	    sizeMap.resize(nranks, 0);
	    offsetMap.resize(nranks, 0);
	    sizeMap[0] = allnFields[0];
	    for (int r = 1; r < nranks; ++r)
	    {
	        sizeMap[r]   = allnFields[r];
	        offsetMap[r] = offsetMap[r-1] + sizeMap[r-1];
	    }
	    nTotFields = offsetMap[nranks-1] + sizeMap[nranks-1];
	}
	std::vector<uint64_t> allFieldHashes = m_comm->Gatherv(rkFormatter, fieldHashes, sizeMap, offsetMap);

	std::vector<uint64_t> allFieldDecomps;
	if (rkFormatter == rk)
	{
	    allFieldDecomps.resize(FIELD_DECOMP_SIZE*nranks, 0);
	    sizeMap2.resize(nranks, 0);
	    offsetMap2.resize(nranks, 0);
	    for (int r = 0; r < nranks; ++r)
	    {
	        std::size_t fdo = FIELD_DECOMP_SIZE*r;
		allFieldDecomps[fdo + FIELD_DECOMP_CNT] = sizeMap[r];
                allFieldDecomps[fdo + FIELD_DECOMP_OFF] = offsetMap[r];
	        sizeMap2[r]   = FIELD_COUNT_SIZE*sizeMap[r];
	        offsetMap2[r] = FIELD_COUNT_SIZE*offsetMap[r];
	    }
	}
        std::vector<uint64_t> allFieldCounts = m_comm->Gatherv(rkFormatter, fieldCounts, sizeMap2, offsetMap2);

	std::vector<uint64_t> allFirstDataDecomps;
	
	if (rkFormatter == rk)
        {
            std::vector<uint64_t> allDataDecomps(DATA_DECOMP_SIZE*nTotFields);
	    
	    nWritten += CreateDataSets(outFilename, rkFormatter, nTotFields,
	        allFieldHashes, allFieldCounts, allFieldDecomps,
	        allFirstDataDecomps, allDataDecomps,
	        fieldDefs, homoYIDs, homoZIDs, homoSIDs, numModesPerDirVar,
		fieldData, fieldMetaDataMap);

	    nWritten += WriteDecompositionData(outFilename, rkFormatter,
	        allFieldDecomps, allDataDecomps);

	    std::set<uint64_t> hashesAssigned;
	    for (int r = 0; r < nranks; ++r)
            {
	        std::size_t fdo = FIELD_DECOMP_SIZE*r;
                std::size_t fdc = allFieldDecomps[fdo + FIELD_DECOMP_CNT];
		std::size_t fho = allFieldDecomps[fdo + FIELD_DECOMP_OFF];
				
	        for (std::size_t f = 0; f < fdc; ++f)
                {
		    uint64_t hash = allFieldHashes[fho + f];
		    
                    if (hashesAssigned.find(hash) == hashesAssigned.end())
		    {
		        hashesAssigned.insert(hash);
		    }
		    else
                    {
		        // This field hash has been assigned to be written by another process.
		        allFieldHashes[fho + f] = 0;
                    }
                }
	    }

        } // end of <if (rkFormatter == rk)> clause

	// Scatter those field hashes that indicate which the field definitions that each process writes to file.
        std::vector<uint64_t> writeFieldHashes = m_comm->Scatterv(rkFormatter, allFieldHashes, sizeMap, offsetMap, nFields);

	// Write field definitions to checkpoint file.
	nWritten += WriteFieldAttributes(outFilename, nFields, varOrder, writeFieldHashes, fieldDefs,
            fieldNames, shapeStrings, homoLengths, numModesPerDirUni);

	firstDataDecomps = m_comm->Scatter(rkFormatter, allFirstDataDecomps, DATA_DECOMP_SIZE);
        
    } // end of <if (reformatting)> clause

    nWritten = WriteData(outFilename, nFields, totalFieldCounts,
	firstDataDecomps, fieldDefs, homoYIDs, homoZIDs, homoSIDs,
	numModesPerDirVar, fieldData);           

    m_comm->Block();

    /*
    if (0 == rk)
    {
        cout << " (" << m_comm->Wtime() - tm0 << "s, HDF5)" << endl;
    }
    */

    return nWritten;
}


/**
 * Create the HDF5 datasets for the checkpoint file.
 *
 * @param outFilename          Output filename.
 * @param rkFormatter          The rank of the process responsible for formatting the checkpoint file.
 * @param nTotFields           The total number of fields handled by all processes.
 * @param allFieldHashes       Hashes of all field definitions.
 * @param allFieldCounts       The counts associated with all the fields handled by all the processes. 
 * @param allFieldDecomps      Field decomposition per process.
 * @param allFirstDataDecomps  Data decomposition of first field per process.
 * @param allDataDecomps       Data decomposition per field per process.
 * @param fieldDefs            Input field definitions.
 * @param homoYIDs             Homogeneous YIDs handled by this process.
 * @param homoZIDs             Homogeneous ZIDs handled by this process.
 * @param homoSIDs             Homogeneous SIDs handled by this process.
 * @param numModesPerDirVar    The number of (universal) modes associated with the fields handled by this process.
 * @param fieldData            Input field data.
 * @param fieldMetaDataMap     Field metadata.
 * @return The number of bytes written.
 */
uint64_t FieldIOHdf5::CreateDataSets(const std::string &outFilename,
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
				     const FieldMetaDataMap &fieldMetaDataMap)
{
    int nranks = m_comm->GetSize();
    int rk = m_comm->GetRank();
    uint64_t nWritten = 0;

    std::stringstream prfx;
    prfx << rk << ": FieldIOHdf5:: CreateDataSets(): ";

    if (rkFormatter != rk)
    {
        return nWritten;
    }
    
    // The root rank creates the file layout from scratch.
    H5::FileSharedPtr outFile = H5::File::Create(outFilename, H5F_ACC_TRUNC);
    ASSERTL1(outFile, prfx.str() + "cannot create HDF5 file.");
    H5::GroupSharedPtr root = outFile->CreateGroup("NEKTAR");
    ASSERTL1(root, prfx.str() + "cannot create root group.");
    TagWriterSharedPtr info_writer(new H5TagWriter(root));
    AddInfoTag(info_writer, fieldMetaDataMap);

    // Record file format version as attribute in main group.
    root->SetAttribute("FORMAT_VERSION", FORMAT_VERSION);
    nWritten += sizeof(FORMAT_VERSION);

    // Calculate the indexes to be used by each MPI process when reading the
    // IDS and DATA datasets.
    std::size_t runningElemIdCnt = 0, runningElemDataCnt = 0, runningOrderCnt = 0;
    std::size_t runningHomyCnt = 0, runningHomzCnt = 0, runningHomsCnt = 0;
    std::size_t ddo = 0;
    
    for (int r = 0; r < nranks; ++r)
    {
        std::size_t fdo = FIELD_DECOMP_SIZE*r;
        std::size_t rkFieldCnt = allFieldDecomps[fdo + FIELD_DECOMP_CNT];
	std::size_t rkFieldOff = allFieldDecomps[fdo + FIELD_DECOMP_OFF];
	std::size_t fco = FIELD_COUNT_SIZE*rkFieldOff;
	
        for (int f = 0; f < rkFieldCnt; ++f)
        {
	    std::size_t elemIdCnt, elemDataCnt, orderCnt;
            std::size_t homyCnt, homzCnt, homsCnt;

            allDataDecomps[ddo + DATA_DECOMP_FIELD_HASH] = allFieldHashes[rkFieldOff + f];

            elemIdCnt = allFieldCounts[fco + FIELD_COUNT_IDS];
            allDataDecomps[ddo + DATA_DECOMP_IDS_CNT] = elemIdCnt;
            allDataDecomps[ddo + DATA_DECOMP_IDS_OFF] = runningElemIdCnt;
            runningElemIdCnt += elemIdCnt;

            elemDataCnt = allFieldCounts[fco + FIELD_COUNT_DATA];
            allDataDecomps[ddo + DATA_DECOMP_DATA_CNT] = elemDataCnt;
            allDataDecomps[ddo + DATA_DECOMP_DATA_OFF] = runningElemDataCnt;
            runningElemDataCnt += elemDataCnt;

            homyCnt = allFieldCounts[fco + FIELD_COUNT_HOMY];
            allDataDecomps[ddo + DATA_DECOMP_HOMY_CNT] = homyCnt;
            allDataDecomps[ddo + DATA_DECOMP_HOMY_OFF] = runningHomyCnt;
            runningHomyCnt += homyCnt;

            homzCnt = allFieldCounts[fco + FIELD_COUNT_HOMZ];
            allDataDecomps[ddo + DATA_DECOMP_HOMZ_CNT] = homzCnt;
            allDataDecomps[ddo + DATA_DECOMP_HOMZ_OFF] = runningHomzCnt;
            runningHomzCnt += homzCnt;

            homsCnt = allFieldCounts[fco + FIELD_COUNT_HOMS];
            allDataDecomps[ddo + DATA_DECOMP_HOMS_CNT] = homsCnt;
            allDataDecomps[ddo + DATA_DECOMP_HOMS_OFF] = runningHomsCnt;
            runningHomsCnt += homsCnt;

	    orderCnt = allFieldCounts[fco + FIELD_COUNT_ORDER];
            allDataDecomps[ddo + DATA_DECOMP_ORDER_CNT] = orderCnt;
            allDataDecomps[ddo + DATA_DECOMP_ORDER_OFF] = runningOrderCnt;
            runningOrderCnt += orderCnt;

	    if (0 == f)
	    {
	        allFirstDataDecomps.insert(allFirstDataDecomps.end(),
		    allDataDecomps.begin()+ddo,
		    allDataDecomps.begin()+ddo+DATA_DECOMP_SIZE);
	    }
	            
            fco += FIELD_COUNT_SIZE;
	    ddo += DATA_DECOMP_SIZE;
        } // end of <for (int f = 0; f < rkFieldCnt; ++f)> loop
	
    } // end of <for (int r = 0; r < nranks; ++r)> loop

    // Create FIELD_DECOMPOSITION dataset: basic field info for each MPI process.
    H5::DataTypeSharedPtr fieldDecompsType = H5::DataType::OfObject(allFieldDecomps[0]);
    H5::DataSpaceSharedPtr fieldDecompsSpace = H5::DataSpace::OneD(FIELD_DECOMP_SIZE*nranks);
    H5::DataSetSharedPtr fieldDecompsDset =
        root->CreateDataSet("FIELD_DECOMPOSITION", fieldDecompsType, fieldDecompsSpace);
    ASSERTL1(fieldDecompsDset, prfx.str() + "cannot create FIELD_DECOMPOSITION dataset.");
	
    // Create DATA_DECOMPOSITION dataset: basic info for each field handled by each MPI process.
    H5::DataTypeSharedPtr dataDecompsType = H5::DataType::OfObject(allDataDecomps[0]);
    H5::DataSpaceSharedPtr dataDecompsSpace = H5::DataSpace::OneD(DATA_DECOMP_SIZE*nTotFields);
    H5::DataSetSharedPtr dataDecompsDset =
        root->CreateDataSet("DATA_DECOMPOSITION", dataDecompsType, dataDecompsSpace);
    ASSERTL1(dataDecompsDset, prfx.str() + "cannot create DATA_DECOMPOSITION dataset.");
	
    // Create IDS dataset: element ids.
    H5::DataTypeSharedPtr elemIdsType = H5::DataType::OfObject(fieldDefs[0]->m_elementIDs[0]);
    H5::DataSpaceSharedPtr elemIdsSpace = H5::DataSpace::OneD(runningElemIdCnt);
    H5::DataSetSharedPtr elemIdsDset =
        root->CreateDataSet("ELEM_IDS", elemIdsType, elemIdsSpace);
    ASSERTL1(elemIdsDset, prfx.str() + "cannot create ELEM_IDS dataset.");

    // Create DATA dataset: element data.
    H5::DataTypeSharedPtr elemDataType = H5::DataType::OfObject(fieldData[0][0]);
    H5::DataSpaceSharedPtr elemDataSpace = H5::DataSpace::OneD(runningElemDataCnt);
    H5::DataSetSharedPtr elemDataDset =
        root->CreateDataSet("ELEM_DATA", elemDataType, elemDataSpace);
    ASSERTL1(elemDataDset, prfx.str() + "cannot create ELEM_DATA dataset.");

    // Create HOMOGENEOUSYIDS dataset: homogeneous y-plane IDs.
    if (runningHomyCnt > 0)
    {
        H5::DataTypeSharedPtr homyType = H5::DataType::OfObject(homoYIDs[0][0]);
        H5::DataSpaceSharedPtr homySpace = H5::DataSpace::OneD(runningHomyCnt);
        H5::DataSetSharedPtr homyDset =
            root->CreateDataSet("HOMOGENEOUSY_IDS", homyType, homySpace);
        ASSERTL1(homyDset, prfx.str() + "cannot create HOMOGENEOUSYIDS dataset.");
    }

    // Create HOMOGENEOUSYIDS dataset: homogeneous z-plane IDs.
    if (runningHomzCnt > 0)
    {
        H5::DataTypeSharedPtr homzType = H5::DataType::OfObject(homoZIDs[0][0]);
        H5::DataSpaceSharedPtr homzSpace = H5::DataSpace::OneD(runningHomzCnt);
        H5::DataSetSharedPtr homzDset =
            root->CreateDataSet("HOMOGENEOUSZ_IDS", homzType, homzSpace);
        ASSERTL1(homzDset, prfx.str() + "cannot create HOMOGENEOUSZIDS dataset.");
    }

    // Create HOMOGENEOUSSIDS dataset: homogeneous strip IDs.
    if (runningHomsCnt > 0)
    {
        H5::DataTypeSharedPtr homsType = H5::DataType::OfObject(homoSIDs[0][0]);
        H5::DataSpaceSharedPtr homsSpace = H5::DataSpace::OneD(runningHomsCnt);
        H5::DataSetSharedPtr homsDset =
            root->CreateDataSet("HOMOGENEOUSS_IDS", homsType, homsSpace);
        ASSERTL1(homsDset, prfx.str() + "cannot create HOMOGENEOUSSIDS dataset.");
    }

    // Create POLYORDERS dataset: elemental polynomial orders.
    if (runningOrderCnt)
    {
        H5::DataTypeSharedPtr orderType = H5::DataType::OfObject(numModesPerDirVar[0][0]);
        H5::DataSpaceSharedPtr orderSpace = H5::DataSpace::OneD(runningOrderCnt);
        H5::DataSetSharedPtr orderDset =
            root->CreateDataSet("POLYORDERS", orderType, orderSpace);
        ASSERTL1(orderDset, prfx.str() + "cannot create POLYORDERS dataset.");
    }
    
    return nWritten;
}

/**
 * Write the FIELD_DECOMPOSITION and DATA_DECOMPOSITION datasets.
 *
 * @param outFilename      Output filename.
 * @param rkFormatter      The rank of the process responsible for formatting the checkpoint file.
 * @param allFieldDecomps  Field decomposition per process.
 * @param allDataDecomps   Data decomposition per field per process.
 * @return The number of bytes written.
 */
uint64_t FieldIOHdf5::WriteDecompositionData(const std::string &outFilename,
                                             const int rkFormatter,
					     std::vector<uint64_t> &allFieldDecomps,
					     std::vector<uint64_t> &allDataDecomps)
{
    int rk = m_comm->GetRank();
    uint64_t nWritten = 0;

    std::stringstream prfx;
    prfx << rk << ": FieldIOHdf5:: WriteDecompositionData(): ";

    if (rkFormatter != rk)
    {
        return nWritten;
    }
    
    H5::PListSharedPtr serialProps = H5::PList::Default();
    H5::PListSharedPtr writeSR     = H5::PList::Default();

    // Reopen the file.
    H5::FileSharedPtr outFile = H5::File::Open(outFilename, H5F_ACC_RDWR, serialProps);
    ASSERTL1(outFile, prfx.str() + "cannot open HDF5 file.");
    H5::GroupSharedPtr root = outFile->OpenGroup("NEKTAR");
    ASSERTL1(root, prfx.str() + "cannot open root group.");

    // Write the FIELD_DECOMPOSITION dataset.
    H5::DataSetSharedPtr fieldDecompsDset = root->OpenDataSet("FIELD_DECOMPOSITION");
    ASSERTL1(fieldDecompsDset, prfx.str() + "cannot open FIELD_DECOMPOSITION dataset.");

    H5::DataSpaceSharedPtr fieldDecompsSpace = fieldDecompsDset->GetSpace();
    ASSERTL1(fieldDecompsSpace, prfx.str() + "cannot open FIELD_DECOMPOSITION filespace.");

    fieldDecompsSpace->SelectRange(0, allFieldDecomps.size());
    fieldDecompsDset->Write(allFieldDecomps, fieldDecompsSpace, writeSR);
    nWritten += allFieldDecomps.size();

    // Write the DATA_DECOMPOSITION dataset.
    H5::DataSetSharedPtr dataDecompsDset = root->OpenDataSet("DATA_DECOMPOSITION");
    ASSERTL1(dataDecompsDset, prfx.str() + "cannot open DATA_DECOMPOSITION dataset.");

    H5::DataSpaceSharedPtr dataDecompsSpace = dataDecompsDset->GetSpace();
    ASSERTL1(dataDecompsSpace, prfx.str() + "cannot open DATA_DECOMPOSITION filespace.");

    dataDecompsSpace->SelectRange(0, allDataDecomps.size());
    dataDecompsDset->Write(allDataDecomps, dataDecompsSpace, writeSR);
    nWritten += allDataDecomps.size();

    return nWritten;
}

/**
 * Determine which process will write the group representing the field description
 * in the HDF5 file. This is done by finding all unique hashes and then determining
 * one process that will create (possibly more than one) group for that hash.
 * An alternative would be to communicate the field information to the formatting process,
 * but that method is somewhat convoluted.
 *
 * @param outFilename       Output filename.
 * @param nFields           The number of fields handled by process.
 * @param varOrder          Is 1 if at least one of the field definitions is not universal order.
 * @param writeFieldHashes  Data structure containing the hashes of the fields to be written by
 *                          this process. A hash is a hash of the field definition.
 * @param fieldDefs         Input field definitions.
 * @param fieldNames        The hashes of the fields handled by this process
 * @param shapeString       The element shapes associated with the fields handled by this process.
 * @param homoLengths       The homogeneous lengths associated with the fields handled by this process.
 * @param numModesPerDirUni The number of (universal) modes associated with the fields handled by this process.
 * @return The number of bytes written.
 */
uint64_t FieldIOHdf5::WriteFieldAttributes(const std::string &outFilename,
					   const uint64_t nFields,
					   const int varOrder,
					   std::vector<uint64_t> &writeFieldHashes,
					   std::vector<FieldDefinitionsSharedPtr> &fieldDefs,
					   std::vector<std::string> &fieldNames,
                                           std::vector<std::string> &shapeStrings,
                                           std::vector<std::vector<NekDouble> > &homoLengths,
                                           std::vector<std::string> &numModesPerDirUni)
{
    int nranks = m_comm->GetSize();
    int rk = m_comm->GetRank();
    uint64_t nWritten = 0;
    bool allZeros = std::all_of(writeFieldHashes.begin(), writeFieldHashes.end(), [](uint64_t h) { return h==0; });

    std::stringstream prfx;
    prfx << rk << ": FieldIOHdf5:: WriteFieldAttributes(): ";
    
    for (int r = 0; r < nranks; ++r)
    {
        // Write out this rank's field definitions.
        if (rk == r && !allZeros)
        {
            H5::PListSharedPtr serialProps = H5::PList::Default();
            
            // Open the file.
            H5::FileSharedPtr outFile =
                H5::File::Open(outFilename, H5F_ACC_RDWR, serialProps);
            ASSERTL1(outFile, prfx.str() + "cannot open HDF5 file.");
            H5::GroupSharedPtr root = outFile->OpenGroup("NEKTAR");
            ASSERTL1(root, prfx.str() + "cannot open root group.");

            for (std::size_t f = 0; f < nFields; ++f)
	    {
	        if (0 == writeFieldHashes[f]) continue;
	        
	        H5::GroupSharedPtr fieldGroup = root->CreateGroup(fieldNames[f]);
                ASSERTL1(fieldGroup, prfx.str() + "cannot create field group.");

                fieldGroup->SetAttribute("FIELDS", fieldDefs[f]->m_fields);
                fieldGroup->SetAttribute("BASIS", fieldDefs[f]->m_basis);
                fieldGroup->SetAttribute("SHAPE", shapeStrings[f]);

                for (std::size_t j = 0; j < fieldDefs[f]->m_fields.size(); ++j)
                {
                    nWritten += fieldDefs[f]->m_fields[j].size();
                }
                nWritten += fieldDefs[f]->m_basis.size()*sizeof(LibUtilities::BasisType);
                nWritten += shapeStrings[f].size();

                if (homoLengths[f].size() > 0)
                {
                    fieldGroup->SetAttribute("HOMOGENEOUSLENGTHS", homoLengths[f]);
                    nWritten += homoLengths[f].size()*sizeof(NekDouble);
                }

                // If the field has only uniform order, we write the order
                // into the NUMMODESPERDIR attribute. Otherwise, we'll go
                // ahead and assume everything is mixed and fix this in the
                // read later if required.
                if (!varOrder)
                {
                    fieldGroup->SetAttribute("NUMMODESPERDIR", numModesPerDirUni[f]);
                    nWritten += numModesPerDirUni[f].size();
                }
                else
                {
                    std::string numModesPerDir = "MIXORDER";
                    fieldGroup->SetAttribute("NUMMODESPERDIR", numModesPerDir);
                    nWritten += numModesPerDir.size();
                }
		    
            } // end of <for (std::size_t f = 0; f < nFields; ++f)> loop
		
        } // end of <if (rk == r)> clause

        // We block to avoid more than one processor opening the file at a time.
        m_comm->Block();
	
    } // end of <for (r = 0; r < nranks; ++r)> loop

    return nWritten;
}

/**
 * Write the data to the checkpoint file.
 *
 * @param outFilename         Output filename.
 * @param nFields             The number of fields handled by this process.
 * @param totalFieldCounts    The cumulative counts over all the fields handled by this process.
 * @param firstDataDecomps    Data decomposition of first field for this process.
 * @param fieldDefs           Input field definitions.
 * @param homoYIDs            Homogeneous YIDs handled by this process.
 * @param homoZIDs            Homogeneous ZIDs handled by this process.
 * @param homoSIDs            Homogeneous SIDs handled by this process.
 * @param numModesPerDirUni   The number of (universal) modes associated with the fields handled by this process.
 * @param fieldData           Input field data.
 * @return The number of bytes written.
 */
uint64_t FieldIOHdf5::WriteData(const std::string &outFilename,
				const std::size_t nFields,
	                        std::vector<uint64_t> &totalFieldCounts,
	                        std::vector<uint64_t> &firstDataDecomps,
                                std::vector<FieldDefinitionsSharedPtr> &fieldDefs,
	                        std::vector<std::vector<unsigned int> > &homoYIDs,
                                std::vector<std::vector<unsigned int> > &homoZIDs,
                                std::vector<std::vector<unsigned int> > &homoSIDs,
                                std::vector<std::vector<unsigned int> > &numModesPerDirVar,
                                std::vector<std::vector<NekDouble> > &fieldData)
{
    int nranks = m_comm->GetSize();
    int rk = m_comm->GetRank();
    uint64_t nWritten = 0;

    std::stringstream prfx;
    prfx << rk << ": FieldIOHdf5::WriteData(): ";


    // Set properties for parallel file access (if we're in parallel).
    H5::PListSharedPtr parallelProps = H5::PList::Default();
    if (nranks > 1)
    {
        // Use MPI/O to access the file
        parallelProps = H5::PList::FileAccess();
        parallelProps->SetMpio(m_comm);
    }

    // Reopen the file.
    H5::FileSharedPtr outFile = H5::File::Open(outFilename, H5F_ACC_RDWR, parallelProps);
    ASSERTL1(outFile, prfx.str() + "cannot open HDF5 file.");
    H5::GroupSharedPtr root = outFile->OpenGroup("NEKTAR");
    ASSERTL1(root, prfx.str() + "cannot open root group.");

    if (!reformatting)
    {
        H5::PListSharedPtr readPL = H5::PList::Default();
        if (nranks > 1)
        {
	    // Use collective IO
            readPL = H5::PList::DatasetXfer();
            readPL->SetDxMpioCollective();
	}
	
        H5::DataSetSharedPtr fieldDecompsDset = root->OpenDataSet("FIELD_DECOMPOSITION");
        ASSERTL1(fieldDecompsDset, prfx.str() + "cannot open FIELD_DECOMPOSITION dataset.");
        H5::DataSpaceSharedPtr fieldDecompsSpace = fieldDecompsDset->GetSpace();
        ASSERTL1(fieldDecompsSpace, prfx.str() + "cannot open FIELD_DECOMPOSITION filespace.");

        std::vector<uint64_t> fieldDecomps(FIELD_DECOMP_SIZE, 0);
        fieldDecompsSpace->SelectRange(FIELD_DECOMP_SIZE*rk, FIELD_DECOMP_SIZE);
        fieldDecompsDset->Read(fieldDecomps, fieldDecompsSpace, readPL);
    
        H5::DataSetSharedPtr dataDecompsDset = root->OpenDataSet("DATA_DECOMPOSITION");
        ASSERTL1(dataDecompsDset, prfx.str() + "cannot open DATA_DECOMPOSITION dataset.");
        H5::DataSpaceSharedPtr dataDecompsSpace = dataDecompsDset->GetSpace();
        ASSERTL1(dataDecompsSpace, prfx.str() + "cannot open DATA_DECOMPOSITION filespace.");

	dataDecompsSpace->SelectRange(DATA_DECOMP_SIZE*fieldDecomps[FIELD_DECOMP_OFF], DATA_DECOMP_SIZE);
        dataDecompsDset->Read(firstDataDecomps, dataDecompsSpace, readPL);
    }  

    
    // Open the optional datasets and data spaces.
    // We use independent (non-collective) writes for all datasets that are actually part of the field definition,
    // i.e., {order_dset,homz_dset,homy_dset,homs_dset}; otherwise, we would have to assume that all fields
    // were defined such that they all used the same selection of field definition datasets.

    if (totalFieldCounts[FIELD_COUNT_HOMY] > 0)
    {
        H5::DataSetSharedPtr homyDset = root->OpenDataSet("HOMOGENEOUSYIDS");
        ASSERTL1(homyDset, prfx.str() + "cannot open HOMOGENEOUSYIDS dataset.");
        H5::DataSpaceSharedPtr homySpace = homyDset->GetSpace();
        ASSERTL1(homySpace, prfx.str() + "cannot open HOMOGENEOUSYIDS filespace.");
	nWritten += WriteFieldDataInd(nFields, homySpace, homyDset, firstDataDecomps[DATA_DECOMP_HOMY_OFF], homoYIDs);
    }

    if (totalFieldCounts[FIELD_COUNT_HOMZ] > 0)
    {
        H5::DataSetSharedPtr homzDset = root->OpenDataSet("HOMOGENEOUSZIDS");
        ASSERTL1(homzDset, prfx.str() + "cannot open HOMOGENEOUSZIDS dataset.");
        H5::DataSpaceSharedPtr homzSpace = homzDset->GetSpace();
        ASSERTL1(homz_fspace, prfx.str() + "cannot open HOMOGENEOUSZIDS filespace.");
	nWritten += WriteFieldDataInd(nFields, homzSpace, homzDset, firstDataDecomps[DATA_DECOMP_HOMZ_OFF], homoZIDs);
    }

    if (totalFieldCounts[FIELD_COUNT_HOMS] > 0)
    {
        H5::DataSetSharedPtr homsDset = root->OpenDataSet("HOMOGENEOUSSIDS");
        ASSERTL1(homsDset, prfx.str() + "cannot open HOMOGENEOUSSIDS dataset.");
        H5::DataSpaceSharedPtr homsSpace = homsDset->GetSpace();
        ASSERTL1(homsSpace, prfx.str() + "cannot open HOMOGENEOUSSIDS filespace.");
	nWritten += WriteFieldDataInd(nFields, homsSpace, homsDset, firstDataDecomps[DATA_DECOMP_HOMS_OFF], homoSIDs);
    }

    if (totalFieldCounts[FIELD_COUNT_ORDER] > 0)
    {
        H5::DataSetSharedPtr orderDset = root->OpenDataSet("POLYORDERS");
        ASSERTL1(orderDset, prfx.str() + "cannot open POLYORDERS dataset.");
        H5::DataSpaceSharedPtr orderSpace = orderDset->GetSpace();
        ASSERTL1(orderSpace, prfx.str() + "cannot open POLYORDERS filespace.");
	nWritten += WriteFieldDataInd(nFields, orderSpace, orderDset, firstDataDecomps[DATA_DECOMP_ORDER_OFF], numModesPerDirVar);
    }


    if (reformatting)
    {
        // Open the element ids dataset.
        H5::DataSetSharedPtr elemIdsDset = root->OpenDataSet("ELEM_IDS");
        ASSERTL1(elemIdsDset, prfx.str() + "cannot open ELEM_IDS dataset.");
        H5::DataSpaceSharedPtr elemIdsSpace = elemIdsDset->GetSpace();
        ASSERTL1(elemIdsSpace, prfx.str() + "cannot open ELEM_IDS filespace.");

        // Write all the element IDs and element data collectively, taking into account the fact that
        // different processes maybe handling different numbers of fields.
        std::vector<std::vector<unsigned int> > elemIDs(nFields);
        for (std::size_t f = 0; f < nFields; ++f)
        {
            elemIDs[f] = fieldDefs[f]->m_elementIDs;
        }

        nWritten += WriteFieldData(nFields, elemIdsSpace, elemIdsDset, firstDataDecomps[DATA_DECOMP_IDS_OFF], elemIDs);
    }

    
    // Open the element data dataset and associated data space.
    H5::DataSetSharedPtr elemDataDset = root->OpenDataSet("ELEM_DATA");
    ASSERTL1(elemDataDset, prfx.str() + "cannot open ELEM_DATA dataset.");
    H5::DataSpaceSharedPtr elemDataSpace = elemDataDset->GetSpace();
    ASSERTL1(elemDataSpace, prfx.str() + "cannot open ELEM_DATA filespace.");
    
    nWritten += WriteFieldData(nFields, elemDataSpace, elemDataDset, firstDataDecomps[DATA_DECOMP_DATA_OFF], fieldData);

    return nWritten;
}

/**
 * @brief Write the contributions from a number of fields to a specific dataset.
 *
 * The data are written independently.
 *
 * @param nFields  the number of fields handled by this process
 * @param fspace   hdf5 file space
 * @param dset     hdf5 data set
 * @param data_i   starting offset into hdf5 data set for this process
 * @param data     vector of data vectors (one for each field) to be written to data set
 * @return The number of bytes written.
 */
template <class T>
uint64_t FieldIOHdf5::WriteFieldDataInd(std::size_t nFields,
                                        H5::DataSpaceSharedPtr &space, H5::DataSetSharedPtr &dset,
                                        uint64_t data_i, std::vector<std::vector<T> > &data)
{
    if (!space || !dset) return 0;

    std::size_t nDataItems = 0;
    uint64_t nWritten = 0;

    H5::PListSharedPtr writePL = H5::PList::DatasetXfer();
    writePL->SetDxMpioIndependent();

    for (std::size_t f = 0; f < nFields; ++f)
    {
        nDataItems = data[f].size();
        if (nDataItems > 0)
        {
            space->SelectRange(data_i, nDataItems);
            dset->Write(data[f], space, writePL);
            data_i += nDataItems;
            nWritten += nDataItems*sizeof(T);
        }
    }

    return nWritten;
}

/**
 * @brief Write the contributions from a number of fields to a specific dataset.
 *
 * The data are written collectively; hence, this method allows for the fact that
 * different processes might be handling different numbers of fields.
 *
 * @param nFields    the number of fields handled by this process
 * @param space      hdf5 file space
 * @param dset       hdf5 data set
 * @param data_i     starting offset into hdf5 data set for this process
 * @param data       vector of data vectors (one for each field) to be written to data set
 * @return The number of bytes written.
 */
template <class T>
uint64_t FieldIOHdf5::WriteFieldData(std::size_t nFields,
                                     H5::DataSpaceSharedPtr &space, H5::DataSetSharedPtr &dset,
                                     uint64_t data_i, std::vector<std::vector<T> > &data)
{
    if (!space || !dset) return 0;

    H5::PListSharedPtr writePL = H5::PList::DatasetXfer();
    writePL->SetDxMpioCollective();

    std::vector<T> last_data;
    H5S_seloper_t h5_select = H5S_SELECT_SET;
    
    for (std::size_t f = 0; f < nFields; ++f)
    {
        std::size_t nDataItems = data[f].size();
        space->SelectRange(h5_select, data_i, nDataItems);
	data_i += nDataItems;
	
	last_data.insert(last_data.end(), data[f].begin(), data[f].end());

	h5_select = H5S_SELECT_OR;
    }

    dset->Write(last_data, space, writePL);
    return last_data.size()*sizeof(T);
}


/**
 * @brief Import a HDF5 format file.
 *
 * @param inFilename       Input filename
 * @param fieldDefs         Field definitions of resulting field
 * @param fieldData         Field data of resulting field
 * @param fieldMetaInfoMap      Field metadata of resulting field
 * @param elementIDs        If specified, contains the list of element IDs on
 *                          this rank. The resulting field definitions will only
 *                          contain data for the element IDs specified in this
 *                          array.
 * @return The number of bytes read.
 */
uint64_t FieldIOHdf5::v_Import(const std::string &inFilename,
                               std::vector<FieldDefinitionsSharedPtr> &fieldDefs,
                               std::vector<std::vector<NekDouble> > &fieldData,
                               FieldMetaDataMap &fieldMetaInfoMap,
                               const Array<OneD, int> &elementIDs)
{
    int nranks = m_comm->GetSize();
    int rk = m_comm->GetRank();
        
    std::stringstream prfx;
    prfx << rk << ": FieldIOHdf5::v_Import(): ";
    
    /*
    double tm0 = 0.0;
    if (0 == rk)
    {
        tm0 = m_comm->Wtime();
    }
    */

    uint64_t nRead = 0;
    
    // Set properties for parallel file access (if we're in parallel)
    H5::PListSharedPtr parallelProps = H5::PList::Default();
    H5::PListSharedPtr readPL = H5::PList::Default();
    
    if (nranks > 1)
    {
        // Use MPI/O to access the file
        parallelProps = H5::PList::FileAccess();
        parallelProps->SetMpio(m_comm);
        // Use collective IO
        readPL = H5::PList::DatasetXfer();
        readPL->SetDxMpioCollective();
    }

    DataSourceSharedPtr dataSource = H5DataSource::create(
        inFilename, parallelProps);

    // Open the root group of the hdf5 file
    H5DataSourceSharedPtr h5 =
        std::static_pointer_cast<H5DataSource>(dataSource);
    ASSERTL1(h5, prfx.str() + "cannot open HDF5 file.");
    H5::GroupSharedPtr root = h5->Get()->OpenGroup("NEKTAR");
    ASSERTL1(root, prfx.str() + "cannot open root group.");

    // Check format version
    unsigned int formatVersion;
    H5::Group::AttrIterator attrIt  = root->attr_begin();
    H5::Group::AttrIterator attrEnd = root->attr_end();
    for (; attrIt != attrEnd; ++attrIt)
    {
        if (*attrIt == "FORMAT_VERSION")
        {
            break;
        }
    }

    ASSERTL0(attrIt != attrEnd,
             "Unable to determine Nektar++ HDF5 file version");
    root->GetAttribute("FORMAT_VERSION", formatVersion);

    nRead += sizeof(formatVersion);
    
    ASSERTL0(formatVersion <= FORMAT_VERSION,
             "File format if " + inFilename + " is higher than supported in "
             "this version of Nektar++"); 

    // Obtain the decompositional data for this process.
    H5::DataSetSharedPtr fieldDecompsDset = root->OpenDataSet("FIELD_DECOMPOSITION");
    ASSERTL1(fieldDecompsDset, prfx.str() + "cannot open FIELD_DECOMPOSITION dataset.");
    H5::DataSpaceSharedPtr fieldDecompsSpace = fieldDecompsDset->GetSpace();
    ASSERTL1(fieldDecompsSpace, prfx.str() + "cannot open FIELD_DECOMPOSITION filespace.");

    std::vector<uint64_t> fieldDecomps(FIELD_DECOMP_SIZE, 0);
    fieldDecompsSpace->SelectRange(FIELD_DECOMP_SIZE*rk, FIELD_DECOMP_SIZE);
    fieldDecompsDset->Read(fieldDecomps, fieldDecompsSpace, readPL);
    std::size_t nFields = fieldDecomps[FIELD_DECOMP_CNT];
    
    H5::DataSetSharedPtr dataDecompsDset = root->OpenDataSet("DATA_DECOMPOSITION");
    ASSERTL1(dataDecompsDset, prfx.str() + "cannot open DATA_DECOMPOSITION dataset.");
    H5::DataSpaceSharedPtr dataDecompsSpace = dataDecompsDset->GetSpace();
    ASSERTL1(dataDecompsSpace, prfx.str() + "cannot open DATA_DECOMPOSITION filespace.");

    std::vector<uint64_t> dataDecomps;
    dataDecompsSpace->SelectRange(DATA_DECOMP_SIZE*fieldDecomps[FIELD_DECOMP_OFF],
        DATA_DECOMP_SIZE*fieldDecomps[FIELD_DECOMP_CNT]);
    dataDecompsDset->Read(dataDecomps, dataDecompsSpace, readPL);

    
    // Open the element ids dataset.
    H5::DataSetSharedPtr elemIdsDset = root->OpenDataSet("ELEM_IDS");
    ASSERTL1(elemIdsDset, "Cannot open ELEM_IDS dataset.");
    H5::DataSpaceSharedPtr elemIdsSpace = elemIdsDset->GetSpace();
    ASSERTL1(ids_fspace, "Cannot open ELEM_IDS filespace.");

    // Open the element data dataset.
    H5::DataSetSharedPtr elemDataDset = root->OpenDataSet("ELEM_DATA");
    ASSERTL1(elemDataDset, prfx.str() + "cannot open ELEM_DATA dataset.");
    H5::DataSpaceSharedPtr elemDataSpace = elemDataDset->GetSpace();
    ASSERTL1(elemDataSpace, prfx.str() + "cannot open ELEM_DATA filespace.");

    
    // Chain together all the range selections for the ELEM_IDS and ELEM_DATA datasets.
    uint64_t nIdItems = 0, nDataItems = 0;
    uint64_t nRealIdItems = 0, nRealDataItems = 0;
    H5S_seloper_t h5_select = H5S_SELECT_SET;
    for (std::size_t f = 0; f < nFields; ++f)
    {
        std::size_t ddo = DATA_DECOMP_SIZE*f;

	nIdItems = dataDecomps[ddo + DATA_DECOMP_IDS_CNT];
        elemIdsSpace->SelectRange(h5_select,
	    dataDecomps[ddo + DATA_DECOMP_IDS_OFF],
	    nIdItems);
	nRealIdItems += nIdItems;
	
	nDataItems = dataDecomps[ddo + DATA_DECOMP_DATA_CNT];
        elemDataSpace->SelectRange(h5_select,
	    dataDecomps[ddo + DATA_DECOMP_DATA_OFF],
	    nDataItems);
	nRealDataItems += nDataItems;
	
	h5_select = H5S_SELECT_OR; 
    }

    // Now read the element ids and data.
    std::vector<unsigned int> decompRealIds;
    if (nRealIdItems > 0)
    {
        elemIdsDset->Read(decompRealIds, elemIdsSpace, readPL);
        ASSERTL0(decompRealIds.size() == nRealIdItems,
            prfx.str() + "unexpected number of id items.");
        nRead += nRealIdItems;
    }
    std::vector<NekDouble> decompRealData;
    if (nRealDataItems > 0)
    {
        elemDataDset->Read(decompRealData, elemDataSpace, readPL);
        ASSERTL0(decompRealData.size() == nRealDataItems,
            prfx.str() + "unexpected number of data items.");
        nRead += nRealDataItems;
    }

    
    uint64_t elemIdsOff = 0;
    uint64_t elemDataOff = 0;
    for (std::size_t f = 0; f < nFields; ++f)
    {
        std::size_t ddo = DATA_DECOMP_SIZE*f;

	std::stringstream fieldNameStream;
        fieldNameStream << dataDecomps[ddo + DATA_DECOMP_FIELD_HASH];
        FieldDefinitionsSharedPtr fieldDef =
                MemoryManager<FieldDefinitions>::AllocateSharedPtr();
        nRead += ImportFieldDef(root, fieldNameStream.str(), fieldDef);
	    
	if (fieldDef->m_numHomogeneousDir >= 1)
        {
            nRead += ImportFieldDataInd(root, "HOMOGENEOUSZIDS", "homogeneousZIDs",
                         dataDecomps[ddo + DATA_DECOMP_HOMZ_CNT],
			 dataDecomps[ddo + DATA_DECOMP_HOMZ_OFF],
                         fieldDef->m_homogeneousZIDs);
        }

        if (fieldDef->m_numHomogeneousDir >= 2)
        {
            nRead += ImportFieldDataInd(root, "HOMOGENEOUSYIDS", "homogeneousYIDs",
                         dataDecomps[ddo + DATA_DECOMP_HOMY_CNT],
			 dataDecomps[ddo + DATA_DECOMP_HOMY_OFF],
                         fieldDef->m_homogeneousYIDs);
        }

        if (fieldDef->m_homoStrips)
        {
            nRead += ImportFieldDataInd(root, "HOMOGENEOUSSIDS", "homogeneousSIDs",
                         dataDecomps[ddo + DATA_DECOMP_HOMS_CNT],
			 dataDecomps[ddo + DATA_DECOMP_HOMS_OFF],
                         fieldDef->m_homogeneousSIDs);
        }

        if (!fieldDef->m_uniOrder)
        {
            nRead += ImportFieldDataInd(root, "POLYORDERS", "numModes",
                         dataDecomps[ddo + DATA_DECOMP_ORDER_CNT],
			 dataDecomps[ddo + DATA_DECOMP_ORDER_OFF],
                         fieldDef->m_numModes);
        }

	// Copy element ids to the field definition.
	uint64_t elemIdsCnt = dataDecomps[ddo + DATA_DECOMP_IDS_CNT];
	if (elemIdsCnt > 0)
        {
            vector<unsigned int>::const_iterator idsSelStart = decompRealIds.begin() + elemIdsOff; 
            vector<unsigned int>::const_iterator idsSelEnd = idsSelStart + elemIdsCnt;
            vector<unsigned int> idsSel(idsSelStart, idsSelEnd);

	    fieldDef->m_elementIDs.insert(fieldDef->m_elementIDs.begin(), idsSel.begin(), idsSel.end());
            fieldDefs.push_back(fieldDef);
	    
	    elemIdsOff += elemIdsCnt;
	}

	// Copy the element data to the fielddata structure.
	uint64_t elemDataCnt = dataDecomps[ddo + DATA_DECOMP_DATA_CNT];
        if (fieldData != NullVectorNekDoubleVector && elemDataCnt > 0)
        {
            vector<NekDouble>::const_iterator dataSelStart = decompRealData.begin() + elemDataOff; 
            vector<NekDouble>::const_iterator dataSelEnd = dataSelStart + elemDataCnt;
            vector<NekDouble> dataSel(dataSelStart, dataSelEnd);

            uint64_t dataSize = CheckFieldDefinition(fieldDef);
            ASSERTL0(dataSel.size() == dataSize * fieldDef->m_fields.size(),
                prfx.str() + "input data is not the same length as header information.");
               
            fieldData.push_back(dataSel);
	    
	    elemDataOff += elemDataCnt;
        }
	
    } // end of <for (std::size_t f = 0; f < nFields; ++f)> loop
    
          
    m_comm->Block();

    /*
    if (0 == rk)
    {
        cout << " (" << m_comm->Wtime() - tm0 << "s, HDF5)" << endl;
    }
    */
    
    return nRead;
    
}


/**
 * @brief Import a number of data items from a specified position within a hdf5 dataset.
 *
 * The data are imported independently.
 *
 * @param root        root group containing field definitions
 * @param dsetName    the name of the hdf5 dataset
 * @param dataTag     data item descriptive tag
 * @param nDataItems  the number of items to be imported from the dataset
 * @param offset      the position within dataset from where reading should commence
 * @param data        the data vector for storing data read
 * @return The number of bytes read.
 */
template <class T>
uint64_t FieldIOHdf5::ImportFieldDataInd(
    H5::GroupSharedPtr root,
    std::string dsetName,
    std::string dataTag,
    uint64_t nDataItems,
    uint64_t offset,
    std::vector<T> &data)
{
    if (nDataItems == 0) return 0;

    std::stringstream prfx;
    prfx << m_comm->GetRank() << ": FieldIOHdf5::ImportFieldDataInd(): ";

    H5::PListSharedPtr readPL = H5::PList::Default();
    readPL = H5::PList::DatasetXfer();
    readPL->SetDxMpioIndependent();

    H5::DataSetSharedPtr dset = root->OpenDataSet(dsetName);
    ASSERTL1(dset, prfx.str() + "cannot open " + dsetName + " dataset.");
    H5::DataSpaceSharedPtr fspace = dset->GetSpace();
    ASSERTL1(fspace, prfx.str() + "cannot open " + dsetName + " filespace.");

    fspace->SelectRange(offset, nDataItems);

    dset->Read(data, fspace, readPL);
    ASSERTL0(data.size() == nDataItems,
        prfx.str() + "unexpected number of " + dataTag + " items.");

    return nDataItems*sizeof(T);
}

/**
 * @brief Import field definitions from a HDF5 document.
 *
 * @param readPL       Reading parameter list.
 * @param root         Root group containing field definitions.
 * @param group        Group name to process.
 * @param def          On output contains field definitions.
 * @return The number of bytes read.
 */
uint64_t FieldIOHdf5::ImportFieldDef(
    H5::GroupSharedPtr        root,
    std::string               group,
    FieldDefinitionsSharedPtr def)
{
    std::stringstream prfx;
    prfx << m_comm->GetRank() << ": FieldIOHdf5::ImportFieldDef(): ";
    
    uint64_t nRead = 0;
    H5::GroupSharedPtr field = root->OpenGroup(group);
    ASSERTL1(field, prfx.str() + "cannot open field group, " + group + '.');

    def->m_uniOrder = false;

    H5::Group::AttrIterator attrIt  = field->attr_begin();
    H5::Group::AttrIterator attrEnd = field->attr_end();
    for (; attrIt != attrEnd; ++attrIt)
    {
        const std::string &attrName = *attrIt;
        if (attrName == "FIELDS")
        {
            field->GetAttribute(attrName, def->m_fields);
            for (int i = 0; i < def->m_fields.size(); i++)
            {  
                nRead += def->m_fields[i].size();
            }
        }
        else if (attrName == "SHAPE")
        {
            std::string shapeString;
            field->GetAttribute(attrName, shapeString);
            nRead += shapeString.size();

            // check to see if homogeneous expansion and if so
            // strip down the shapeString definition
            int loc;
            //---> this finds the first location of 'n'!
            if (shapeString.find("Strips") != string::npos)
            {
                def->m_homoStrips = true;
            }

            if ((loc = shapeString.find_first_of("-")) != string::npos)
            {
                if (shapeString.find("Exp1D") != string::npos)
                {
                    def->m_numHomogeneousDir = 1;
                }
                else // HomogeneousExp1D
                {
                    def->m_numHomogeneousDir = 2;
                }

                shapeString.erase(loc, shapeString.length());
            }

            // get the geometrical shape
            bool valid = false;
            for (unsigned int j = 0; j < SIZE_ShapeType; j++)
            {
                if (ShapeTypeMap[j] == shapeString)
                {
                    def->m_shapeType = (ShapeType)j;
                    valid = true;
                    break;
                }
            }

            ASSERTL0(valid, prfx.str() + std::string(
                         "unable to correctly parse the shape type: ")
                     .append(shapeString).c_str());
        }
        else if (attrName == "BASIS")
        {
            field->GetAttribute(attrName, def->m_basis);
            nRead += def->m_basis.size()*sizeof(LibUtilities::BasisType);
            
            // check the basis is in range
            std::vector<BasisType>::const_iterator bIt  = def->m_basis.begin();
            std::vector<BasisType>::const_iterator bEnd = def->m_basis.end();
            for (; bIt != bEnd; ++bIt)
            {
                BasisType bt = *bIt;
                ASSERTL0(bt >= 0 && bt < SIZE_BasisType,
                         prfx.str() +
                         "unable to correctly parse the basis types.");
            }
        }
        else if (attrName == "HOMOGENEOUSLENGTHS")
        {
            field->GetAttribute(attrName, def->m_homogeneousLengths);
            nRead += def->m_homogeneousLengths.size()*sizeof(NekDouble);
        }
        else if (attrName == "NUMMODESPERDIR")
        {
            std::string numModesPerDir;
            field->GetAttribute(attrName, numModesPerDir);
            nRead += numModesPerDir.size();

            if (strstr(numModesPerDir.c_str(), "UNIORDER:"))
            {
                def->m_uniOrder = true;
                bool valid = ParseUtils::GenerateVector(
                    numModesPerDir.substr(9), def->m_numModes);
                ASSERTL0(valid,
                         prfx.str() +
                         "unable to correctly parse the number of modes.");
            }
        }
        else if (attrName == "POINTSTYPE")
        {
            std::string pointsString;
            field->GetAttribute(attrName, pointsString);
            nRead += pointsString.size();
            def->m_pointsDef = true;

            std::vector<std::string> pointsStrings;
            bool valid = ParseUtils::GenerateVector(
                pointsString, pointsStrings);
            ASSERTL0(valid,
                     prfx.str() +
                     "unable to correctly parse the points types.");
            for (std::vector<std::string>::size_type i = 0;
                 i < pointsStrings.size();
                 i++)
            {
                valid = false;
                for (unsigned int j = 0; j < SIZE_PointsType; j++)
                {
                    if (kPointsTypeStr[j] == pointsStrings[i])
                    {
                        def->m_points.push_back((PointsType)j);
                        valid = true;
                        break;
                    }
                }

                ASSERTL0(
                    valid,
                    prfx.str() +
                    std::string(
                        "unable to correctly parse the points type: ")
                    .append(pointsStrings[i])
                    .c_str());
            }
        }
        else if (attrName == "NUMPOINTSPERDIR")
        {
            std::string numPointsPerDir;
            field->GetAttribute(attrName, numPointsPerDir);
            nRead += numPointsPerDir.size();
            def->m_numPointsDef = true;

            bool valid = ParseUtils::GenerateVector(
                numPointsPerDir, def->m_numPoints);
            ASSERTL0(valid,
                     prfx.str() +
                     "unable to correctly parse the number of points.");
        }
        else
        {
            std::string errstr("unknown attribute: ");
            errstr += attrName;
            ASSERTL1(false, prfx.str() + errstr.c_str());
        }
    }

    return nRead;
}


/**
 * @brief Import field metadata from @p filename and return the data source
 * which wraps @p filename.
 *
 * @param filename          Input filename.
 * @param fieldmetadatamap  Resulting field metadata from @p dataSource.
 */
DataSourceSharedPtr FieldIOHdf5::v_ImportFieldMetaData(
    const std::string &filename,
    FieldMetaDataMap  &fieldmetadatamap)
{
    H5::PListSharedPtr parallelProps = H5::PList::Default();
    DataSourceSharedPtr ans = H5DataSource::create(filename, parallelProps);
    ImportHDF5FieldMetaData(ans, fieldmetadatamap);

    return ans;
}

/**
 * @brief Import field metadata from @p dataSource.
 *
 * @param dataSource        Input datasource, which should be a H5DataSource.
 * @param fieldmetadatamap  Resulting field metadata from @p dataSource.
 */
void FieldIOHdf5::ImportHDF5FieldMetaData(DataSourceSharedPtr dataSource,
                                          FieldMetaDataMap &fieldmetadatamap)
{
    H5DataSourceSharedPtr hdf =
        std::static_pointer_cast<H5DataSource>(dataSource);

    H5::GroupSharedPtr master = hdf->Get()->OpenGroup("NEKTAR");
    // New metadata format only in HDF
    H5::GroupSharedPtr metadata = master->OpenGroup("Metadata");

    if (metadata)
    {
        H5::Group::AttrIterator param = metadata->attr_begin(),
                                pEnd = metadata->attr_end();
        for (; param != pEnd; ++param)
        {
            std::string paramString = *param;
            if (paramString != "Provenance")
            {
                std::string paramBodyStr;
                metadata->GetAttribute(paramString, paramBodyStr);
                fieldmetadatamap[paramString] = paramBodyStr;
            }
        }
    }
}

}
}
