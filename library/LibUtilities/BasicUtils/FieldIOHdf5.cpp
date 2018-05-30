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

/// A helper for FieldIOHdf5::v_Write and FieldIOHdf5::v_Import. Describes the
/// position of the number of elements in decomposition (i.e. field definition).
const unsigned int FieldIOHdf5::ELEM_DCMP_IDX  = 0;
/// A helper for FieldIOHdf5::v_Write and FieldIOHdf5::v_Import. Describes the
/// position of the number of data points in decomposition (i.e. field
/// definition).
const unsigned int FieldIOHdf5::VAL_DCMP_IDX   = 1;
/// A helper for FieldIOHdf5::v_Write and FieldIOHdf5::v_Import. Describes the
/// position of the number of elements multiplied by the dimension of the
/// element, giving number of modes when variable polynomial order is defined.
const unsigned int FieldIOHdf5::ORDER_DCMP_IDX = 2;
/// A helper for FieldIOHdf5::v_Write and FieldIOHdf5::v_Import. Describes the
/// position of the number of the number of y-planes for homogeneous
/// simulations.
const unsigned int FieldIOHdf5::HOMY_DCMP_IDX  = 3;
/// A helper for FieldIOHdf5::v_Write and FieldIOHdf5::v_Import. Describes the
/// position of the number of the number of z-planes for homogeneous
/// simulations.
const unsigned int FieldIOHdf5::HOMZ_DCMP_IDX  = 4;
/// A helper for FieldIOHdf5::v_Write and FieldIOHdf5::v_Import. Describes the
/// position of the number of the number of strips for homogeneous simulations.
const unsigned int FieldIOHdf5::HOMS_DCMP_IDX  = 5;
/// The hash of the field definition information, which defines the name of the
/// attribute containing the field definition itself.
const unsigned int FieldIOHdf5::HASH_DCMP_IDX  = 6;
/// A helper for FieldIOHdf5::v_Write. Describes the maximum number of items in
/// the decomposition per field definition.
const unsigned int FieldIOHdf5::MAX_DCMPS      = FieldIOHdf5::HASH_DCMP_IDX + 1;

/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// elements in the cnt array.
const unsigned int FieldIOHdf5::ELEM_CNT_IDX   = 0;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// data points in the cnt array.
const unsigned int FieldIOHdf5::VAL_CNT_IDX    = 1;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// order points in the cnt array.
const unsigned int FieldIOHdf5::ORDER_CNT_IDX  = 2;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// homogeneous y-planes in the cnt array.
const unsigned int FieldIOHdf5::HOMY_CNT_IDX   = 3;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// homogeneous z-planes in the cnt array.
const unsigned int FieldIOHdf5::HOMZ_CNT_IDX   = 4;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// homogeneous strips in the cnt array.
const unsigned int FieldIOHdf5::HOMS_CNT_IDX   = 5;
/// A helper for FieldIOHdf5::v_Write. Describes the maximum number of items in
/// the cnt array per field definition.
const unsigned int FieldIOHdf5::MAX_CNTS       = FieldIOHdf5::HOMS_CNT_IDX + 1;

/// A helper for FieldIOHdf5::v_Write. Describes the position of the element IDs
/// within the indexing set.
const unsigned int FieldIOHdf5::IDS_IDX_IDX    = 0;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the data size
/// within the indexing set.
const unsigned int FieldIOHdf5::DATA_IDX_IDX   = 1;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the element
/// order within the indexing set.
const unsigned int FieldIOHdf5::ORDER_IDX_IDX  = 2;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// y-planes within the indexing set.
const unsigned int FieldIOHdf5::HOMY_IDX_IDX   = 3;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// z-planes within the indexing set.
const unsigned int FieldIOHdf5::HOMZ_IDX_IDX   = 4;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// homogeneous strips within the indexing set.
const unsigned int FieldIOHdf5::HOMS_IDX_IDX   = 5;
/// A helper for FieldIOHdf5::v_Write. Describes the maximum number of items in
/// the indexing set.
const unsigned int FieldIOHdf5::MAX_IDXS       = FieldIOHdf5::HOMS_IDX_IDX + 1;

/**
 * @brief Construct the FieldIO object for HDF5 output.
 *
 * @param pComm              Communicator object.
 * @param sharedFilesystem   True if this system has a shared filesystem.
 */
FieldIOHdf5::FieldIOHdf5(LibUtilities::CommSharedPtr pComm,
                         bool sharedFilesystem)
    : FieldIO(pComm, sharedFilesystem)
{
}

/**
 * @brief Write a HDF5 file to @p outFile given the field definitions @p
 * fielddefs, field data @p fielddata and metadata @p fieldmetadatamap.
 *
 * The writing strategy for HDF5 output is as follows:
 *
 *   - Each rank determines the amount of data needed to be written into each
 *     dataset.
 *   - Each rank communicates its decomposition information to the root process.
 *   - The root processor initialises the output structure, writes the
 *     decomposition dataset and all the field definition information.
 *   - Other ranks may have field definitions that do not belong to the root
 *     process, in which case they open the file and append this (since
 *     attributes cannot be written in parallel).
 *   - Each of the other ranks writes their data contributions to the rest of
 *     the set.
 *
 * @param outFile           Output filename.
 * @param fielddefs         Input field definitions.
 * @param fielddata         Input field data.
 * @param fieldmetadatamap  Field metadata.
 * @return The number of bytes written.
 */
uint64_t FieldIOHdf5::v_Write(const std::string &outFile,
                              std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                              std::vector<std::vector<NekDouble> > &fielddata,
                              const FieldMetaDataMap &fieldmetadatamap,
                              const bool backup)
{
    std::stringstream prfx;
    prfx << m_comm->GetRank() << ": FieldIOHdf5::v_Write(): ";

    uint64_t nWritten = 0;

    double tm0 = 0.0, tm1 = 0.0;

    if (m_comm->GetRank() == 0)
    {
        tm0 = m_comm->Wtime();
    }
    
    SetUpOutput(outFile, false, backup);

    // We make a number of assumptions in this code:
    //   1. All element ids have the same type: unsigned int
    //   2. All elements within a given field have the same number of values
    //   3. All element values have the same type, NekDouble

    // Determine the root MPI process, i.e., the lowest ranked process handling
    // nMaxFields fields, that will be responsible for writing our file.
    ASSERTL1(fielddefs.size() == fielddata.size(),
             prfx.str() + "fielddefs and fielddata have incompatible lengths.");

    std::size_t nFields    = fielddefs.size();
    std::size_t nMaxFields = nFields;
    std::size_t nMinFields = nFields;
    m_comm->AllReduce(nMaxFields, LibUtilities::ReduceMax);
    m_comm->AllReduce(nMinFields, LibUtilities::ReduceMin);
    //cout << prfx.str() << "nFields " << nFields << endl;
    
    int root_rank = -1;
    bool amRoot = false;
    LibUtilities::CommSharedPtr max_fields_comm;

    if (m_comm->GetSize() > 1)
    {
        max_fields_comm = m_comm->CommCreateIf((nFields == nMaxFields) ? 1 : 0);
    }
    else
    {
        max_fields_comm = m_comm;
    }

    if (max_fields_comm)
    {
        int rank  = m_comm->GetRank();
        root_rank = rank;
        max_fields_comm->AllReduce(root_rank, LibUtilities::ReduceMin);
        amRoot = (rank == root_rank);
        if (!amRoot)
        {
            root_rank = -1;
        }
    }

    m_comm->AllReduce(root_rank, LibUtilities::ReduceMax);
    ASSERTL1(root_rank >= 0 && root_rank < m_comm->GetSize(),
             prfx.str() + "invalid root rank.");

    std::vector<uint64_t> decomps(nMaxFields * MAX_DCMPS, 0);
    std::vector<uint64_t> all_hashes(nMaxFields * m_comm->GetSize(), 0);
    std::vector<uint64_t> cnts(MAX_CNTS, 0);
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
        if (!fielddefs[f]->m_uniOrder)
        {
            varOrder = 1;
            break;
        }
    }

    m_comm->AllReduce(varOrder, LibUtilities::ReduceMax);

    // Calculate the total number of elements handled by this MPI process and
    // the total number of bytes required to store the elements. Base the name
    // of each field on the hash of the field definition.
    for (std::size_t f = 0; f < nFields; ++f)
    {
        ASSERTL1(fielddata[f].size() > 0,
                 prfx.str() +
                     "fielddata vector must contain at least one value.");
        ASSERTL1(fielddata[f].size() ==
                     fielddefs[f]->m_fields.size() *
                         CheckFieldDefinition(fielddefs[f]),
                 prfx.str() + "fielddata vector has invalid size.");

        std::size_t nFieldElems = fielddefs[f]->m_elementIDs.size();
        std::size_t nElemVals   = fielddata[f].size();

        decomps[f * MAX_DCMPS + ELEM_DCMP_IDX] = nFieldElems;
        decomps[f * MAX_DCMPS + VAL_DCMP_IDX]  = nElemVals;

        cnts[ELEM_CNT_IDX] += nFieldElems;
        cnts[VAL_CNT_IDX]  += nElemVals;

        // Hash the field specification
        std::stringstream hashStream;
        std::size_t nSubFields = fielddefs[f]->m_fields.size();
        for (std::size_t sf = 0; sf < nSubFields; ++sf)
        {
            hashStream << fielddefs[f]->m_fields[sf];
        }

        nSubFields = fielddefs[f]->m_basis.size();
        for (std::size_t sf = 0; sf < nSubFields; ++sf)
        {
            hashStream << fielddefs[f]->m_basis[sf];
        }

        // Determine SHAPE attribute
        std::stringstream shapeStringStream;
        shapeStringStream << ShapeTypeMap[fielddefs[f]->m_shapeType];

        if (fielddefs[f]->m_numHomogeneousDir > 0)
        {
            if (homDim == -1)
            {
                homDim = fielddefs[f]->m_numHomogeneousDir;
            }

            ASSERTL1(homDim == fielddefs[f]->m_numHomogeneousDir,
                     "HDF5 does not support variable homogeneous directions in "
                     "the same file.");

            shapeStringStream << "-HomogenousExp"
                              << fielddefs[f]->m_numHomogeneousDir << "D";
        }

        if (fielddefs[f]->m_homoStrips)
        {
            shapeStringStream << "-Strips";
        }

        shapeStrings[f] = shapeStringStream.str();
        hashStream << shapeStringStream.str();

        // Determine HOMOGENEOUS attributes
        if (fielddefs[f]->m_numHomogeneousDir)
        {
            nSubFields = fielddefs[f]->m_homogeneousLengths.size();
            homoLengths[f].resize(nSubFields);
            for (std::size_t sf = 0; sf < nSubFields; ++sf)
            {
                NekDouble len = fielddefs[f]->m_homogeneousLengths[sf];
                hashStream << len;
                homoLengths[f][sf] = len;
            }

            nSubFields = fielddefs[f]->m_homogeneousYIDs.size();
            if (nSubFields > 0)
            {
                homoYIDs[f].resize(nSubFields);
                decomps[f * MAX_DCMPS + HOMY_DCMP_IDX] = nSubFields;
                cnts[HOMY_CNT_IDX] += nSubFields;
                for (std::size_t sf = 0; sf < nSubFields; ++sf)
                {
                    homoYIDs[f][sf] = fielddefs[f]->m_homogeneousYIDs[sf];
                }
            }

            nSubFields = fielddefs[f]->m_homogeneousZIDs.size();
            if (nSubFields > 0)
            {
                homoZIDs[f].resize(nSubFields);
                decomps[f * MAX_DCMPS + HOMZ_DCMP_IDX] = nSubFields;
                cnts[HOMZ_CNT_IDX] += nSubFields;
                for (std::size_t sf = 0; sf < nSubFields; ++sf)
                {
                    homoZIDs[f][sf] = fielddefs[f]->m_homogeneousZIDs[sf];
                }
            }

            nSubFields = fielddefs[f]->m_homogeneousSIDs.size();
            if (nSubFields > 0)
            {
                homoSIDs[f].resize(nSubFields);
                decomps[f * MAX_DCMPS + HOMS_DCMP_IDX] = nSubFields;
                cnts[HOMS_CNT_IDX] += nSubFields;
                for (std::size_t sf = 0; sf < nSubFields; ++sf)
                {
                    homoSIDs[f][sf] = fielddefs[f]->m_homogeneousSIDs[sf];
                }
            }
        }

        if (fielddefs[f]->m_uniOrder)
        {
            std::vector<unsigned int> elemModes(fielddefs[f]->m_basis.size());

            for (std::vector<int>::size_type i = 0;
                 i < fielddefs[f]->m_basis.size(); ++i)
            {
                elemModes[i] = fielddefs[f]->m_numModes[i];
            }

            if (varOrder)
            {
                for (std::vector<int>::size_type i = 0; i < nFieldElems; ++i)
                {
                    std::copy(elemModes.begin(), elemModes.end(),
                              std::back_inserter(numModesPerDirVar[f]));
                }
                decomps[f * MAX_DCMPS + ORDER_DCMP_IDX] =
                    nFieldElems * elemModes.size();
                cnts[ORDER_CNT_IDX] += nFieldElems * elemModes.size();
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
            numModesPerDirVar[f] = fielddefs[f]->m_numModes;
            decomps[f * MAX_DCMPS + ORDER_DCMP_IDX] =
                fielddefs[f]->m_numModes.size();
            cnts[ORDER_CNT_IDX] += fielddefs[f]->m_numModes.size();
        }

        std::hash<std::string> string_hasher;
        std::stringstream fieldNameStream;
        uint64_t fieldDefHash = string_hasher(hashStream.str());

        decomps[f * MAX_DCMPS + HASH_DCMP_IDX]         = fieldDefHash;
        all_hashes[m_comm->GetRank() * nMaxFields + f] = fieldDefHash;

        fieldNameStream << fieldDefHash;
        fieldNames[f] = fieldNameStream.str();
    }

    // Gather information from all MPI processes
    std::vector<uint64_t> all_cnts = m_comm->Gather(root_rank, cnts);
    std::vector<uint64_t> all_idxs(m_comm->GetSize() * MAX_IDXS, 0);
    std::vector<uint64_t> all_decomps = m_comm->Gather(root_rank, decomps);
    std::vector<uint64_t> all_dsetsize(MAX_CNTS, 0);

    // The root rank creates the file layout from scratch
    if (amRoot)
    {
        H5::FileSharedPtr outfile = H5::File::Create(outFile, H5F_ACC_TRUNC);
        ASSERTL1(outfile, prfx.str() + "cannot create HDF5 file.");
        H5::GroupSharedPtr root = outfile->CreateGroup("NEKTAR");
        ASSERTL1(root, prfx.str() + "cannot create root group.");
        TagWriterSharedPtr info_writer(new H5TagWriter(root));
        AddInfoTag(info_writer, fieldmetadatamap);

        // Record file format version as attribute in main group.
        root->SetAttribute("FORMAT_VERSION", FORMAT_VERSION);
        nWritten += sizeof(FORMAT_VERSION);

        // Calculate the indexes to be used by each MPI process when reading the
        // IDS and DATA datasets
        std::size_t nTotElems = 0, nTotVals = 0, nTotOrder = 0;
        std::size_t nTotHomY = 0, nTotHomZ = 0, nTotHomS = 0;
        int nRanks = m_comm->GetSize();
        for (int r = 0; r < nRanks; ++r)
        {
            all_idxs[r * MAX_IDXS + IDS_IDX_IDX]   = nTotElems;
            all_idxs[r * MAX_IDXS + DATA_IDX_IDX]  = nTotVals;
            all_idxs[r * MAX_IDXS + ORDER_IDX_IDX] = nTotOrder;
            all_idxs[r * MAX_IDXS + HOMY_IDX_IDX]  = nTotHomY;
            all_idxs[r * MAX_IDXS + HOMZ_IDX_IDX]  = nTotHomZ;
            all_idxs[r * MAX_IDXS + HOMS_IDX_IDX]  = nTotHomS;

            nTotElems += all_cnts[r * MAX_CNTS + ELEM_CNT_IDX];
            nTotVals  += all_cnts[r * MAX_CNTS + VAL_CNT_IDX];
            nTotOrder += all_cnts[r * MAX_CNTS + ORDER_CNT_IDX];
            nTotHomY  += all_cnts[r * MAX_CNTS + HOMY_CNT_IDX];
            nTotHomZ  += all_cnts[r * MAX_CNTS + HOMZ_CNT_IDX];
            nTotHomS  += all_cnts[r * MAX_CNTS + HOMS_CNT_IDX];
        }

        all_dsetsize[ELEM_CNT_IDX ] = nTotElems;
        all_dsetsize[VAL_CNT_IDX  ] = nTotVals;
        all_dsetsize[ORDER_CNT_IDX] = nTotOrder;
        all_dsetsize[HOMY_CNT_IDX ] = nTotHomY;
        all_dsetsize[HOMZ_CNT_IDX ] = nTotHomZ;
        all_dsetsize[HOMS_CNT_IDX ] = nTotHomS;

        // Create DECOMPOSITION dataset: basic field info for each MPI process
        H5::DataTypeSharedPtr decomps_type =
            H5::DataType::OfObject(all_decomps[0]);
        H5::DataSpaceSharedPtr decomps_space =
            H5::DataSpace::OneD(all_decomps.size());
        H5::DataSetSharedPtr decomps_dset =
            root->CreateDataSet("DECOMPOSITION", decomps_type, decomps_space);
        ASSERTL1(decomps_dset,
                 prfx.str() + "cannot create DECOMPOSITION dataset.");

        // Create IDS dataset: element ids
        H5::DataTypeSharedPtr ids_type =
            H5::DataType::OfObject(fielddefs[0]->m_elementIDs[0]);
        H5::DataSpaceSharedPtr ids_space = H5::DataSpace::OneD(nTotElems);
        H5::DataSetSharedPtr ids_dset =
            root->CreateDataSet("ELEMENTIDS", ids_type, ids_space);
        ASSERTL1(ids_dset, prfx.str() + "cannot create ELEMENTIDS dataset.");

        // Create DATA dataset: element data
        H5::DataTypeSharedPtr data_type =
            H5::DataType::OfObject(fielddata[0][0]);
        H5::DataSpaceSharedPtr data_space = H5::DataSpace::OneD(nTotVals);
        H5::DataSetSharedPtr data_dset =
            root->CreateDataSet("DATA", data_type, data_space);
        ASSERTL1(data_dset, prfx.str() + "cannot create DATA dataset.");

        // Create HOMOGENEOUSYIDS dataset: homogeneous y-plane IDs
        if (nTotHomY > 0)
        {
            H5::DataTypeSharedPtr homy_type =
                H5::DataType::OfObject(homoYIDs[0][0]);
            H5::DataSpaceSharedPtr homy_space = H5::DataSpace::OneD(nTotHomY);
            H5::DataSetSharedPtr homy_dset =
                root->CreateDataSet("HOMOGENEOUSYIDS", homy_type, homy_space);
            ASSERTL1(homy_dset,
                     prfx.str() + "cannot create HOMOGENEOUSYIDS dataset.");
        }

        // Create HOMOGENEOUSYIDS dataset: homogeneous z-plane IDs
        if (nTotHomZ > 0)
        {
            H5::DataTypeSharedPtr homz_type =
                H5::DataType::OfObject(homoZIDs[0][0]);
            H5::DataSpaceSharedPtr homz_space = H5::DataSpace::OneD(nTotHomZ);
            H5::DataSetSharedPtr homz_dset =
                root->CreateDataSet("HOMOGENEOUSZIDS", homz_type, homz_space);
            ASSERTL1(homz_dset,
                     prfx.str() + "cannot create HOMOGENEOUSZIDS dataset.");
        }

        // Create HOMOGENEOUSSIDS dataset: homogeneous strip IDs
        if (nTotHomS > 0)
        {
            H5::DataTypeSharedPtr homs_type =
                H5::DataType::OfObject(homoSIDs[0][0]);
            H5::DataSpaceSharedPtr homs_space = H5::DataSpace::OneD(nTotHomS);
            H5::DataSetSharedPtr homs_dset =
                root->CreateDataSet("HOMOGENEOUSSIDS", homs_type, homs_space);
            ASSERTL1(homs_dset,
                     prfx.str() + "cannot create HOMOGENEOUSSIDS dataset.");
        }

        // Create POLYORDERS dataset: elemental polynomial orders
        if (varOrder)
        {
            H5::DataTypeSharedPtr order_type =
                H5::DataType::OfObject(numModesPerDirVar[0][0]);
            H5::DataSpaceSharedPtr order_space = H5::DataSpace::OneD(nTotOrder);
            H5::DataSetSharedPtr order_dset =
                root->CreateDataSet("POLYORDERS", order_type, order_space);
            ASSERTL1(order_dset,
                     prfx.str() + "cannot create POLYORDERS dataset.");
        }
    }

    m_comm->Bcast(all_dsetsize, root_rank);

    // Datasets, root group and HDF5 file are all closed automatically since
    // they are now out of scope. Now we need to determine which process will
    // write the group representing the field description in the HDF5 file. This
    // next block of code performs this by finding all unique hashes and then
    // determining one process that will create (possibly more than one) group
    // for that hash. An alternative would be to communicate the field
    // information to the root processor, but this is a bit convoluted.

    // This set stores the unique hashes.
    std::set<uint64_t> hashToProc;
    // This map takes ranks to hashes this process will write.
    std::map<int, std::vector<uint64_t> > writingProcs;

    // Gather all field hashes to every processor.
    m_comm->AllReduce(all_hashes, LibUtilities::ReduceMax);

    for (int n = 0; n < m_comm->GetSize(); ++n)
    {
        for (std::size_t i = 0; i < nMaxFields; ++i)
        {
            uint64_t hash = all_hashes[n*nMaxFields + i];

            // Note hash can be zero if, on this process, nFields < nMaxFields.
            if (hashToProc.find(hash) != hashToProc.end() || hash == 0)
            {
                continue;
            }
            hashToProc.insert(hash);
            writingProcs[n].push_back(hash);
        }
    }

    // Having constructed the map, go ahead and write the attributes out.
    map<int, std::vector<uint64_t> >::iterator sIt;
    for (sIt = writingProcs.begin(); sIt != writingProcs.end(); sIt++)
    {
        int rank = sIt->first;

        // Write out this rank's groups.
        if (m_comm->GetRank() == rank)
        {
            H5::PListSharedPtr serialProps = H5::PList::Default();
            H5::PListSharedPtr writeSR     = H5::PList::Default();

            // Reopen the file
            H5::FileSharedPtr outfile =
                H5::File::Open(outFile, H5F_ACC_RDWR, serialProps);
            ASSERTL1(outfile, prfx.str() + "cannot open HDF5 file.");
            H5::GroupSharedPtr root = outfile->OpenGroup("NEKTAR");
            ASSERTL1(root, prfx.str() + "cannot open root group.");

            // Write a HDF5 group for each field
            hashToProc.clear();
            for (std::size_t i = 0; i < sIt->second.size(); ++i)
            {
                for (std::size_t f = 0; f < nFields; ++f)
                {
                    if (sIt->second[i] !=
                        all_hashes[m_comm->GetRank() * nMaxFields + f] ||
                        hashToProc.find(sIt->second[i]) != hashToProc.end())
                    {
                        continue;
                    }

                    hashToProc.insert(sIt->second[i]);

                    // Just in case we've already written this
                    H5::GroupSharedPtr field_group =
                        root->CreateGroup(fieldNames[f]);
                    ASSERTL1(field_group,
                             prfx.str() + "cannot create field group.");

                    field_group->SetAttribute("FIELDS", fielddefs[f]->m_fields);
                    field_group->SetAttribute("BASIS", fielddefs[f]->m_basis);
                    field_group->SetAttribute("SHAPE", shapeStrings[f]);

                    for (std::size_t j = 0; j < fielddefs[f]->m_fields.size(); ++j)
                    {
                        nWritten += fielddefs[f]->m_fields[j].size();
                    }
                    nWritten += fielddefs[f]->m_basis.size()*sizeof(LibUtilities::BasisType);
                    nWritten += shapeStrings[f].size();

                    if (homoLengths[f].size() > 0)
                    {
                        field_group->SetAttribute("HOMOGENEOUSLENGTHS",
                                                  homoLengths[f]);
                        nWritten += homoLengths[f].size()*sizeof(NekDouble);
                    }

                    // If the field has only uniform order, we write the order
                    // into the NUMMODESPERDIR attribute. Otherwise, we'll go
                    // ahead and assume everything is mixed and fix this in the
                    // read later if required.
                    if (!varOrder)
                    {
                        field_group->SetAttribute("NUMMODESPERDIR",
                                                  numModesPerDirUni[f]);
                        nWritten += numModesPerDirUni[f].size();
                    }
                    else
                    {
                        std::string numModesPerDir = "MIXORDER";
                        field_group->SetAttribute("NUMMODESPERDIR",
                                                  numModesPerDir);
                        nWritten += numModesPerDir.size();
                    }
                }
            }
        }

        // We block to avoid more than one processor opening the file at a time.
        m_comm->Block();
    }

    // Write the DECOMPOSITION dataset
    if (amRoot)
    {
        H5::PListSharedPtr serialProps = H5::PList::Default();
        H5::PListSharedPtr writeSR     = H5::PList::Default();

        // Reopen the file
        H5::FileSharedPtr outfile =
            H5::File::Open(outFile, H5F_ACC_RDWR, serialProps);
        ASSERTL1(outfile, prfx.str() + "cannot open HDF5 file.");
        H5::GroupSharedPtr root = outfile->OpenGroup("NEKTAR");
        ASSERTL1(root, prfx.str() + "cannot open root group.");

        // Write the DECOMPOSITION dataset
        H5::DataSetSharedPtr decomps_dset = root->OpenDataSet("DECOMPOSITION");
        ASSERTL1(decomps_dset,
                 prfx.str() + "cannot open DECOMPOSITION dataset.");

        H5::DataSpaceSharedPtr decomps_fspace = decomps_dset->GetSpace();
        ASSERTL1(decomps_fspace,
                 prfx.str() + "cannot open DECOMPOSITION filespace.");

        decomps_fspace->SelectRange(0, all_decomps.size());
        decomps_dset->Write(all_decomps, decomps_fspace, writeSR);
        nWritten += all_decomps.size()*sizeof(uint64_t);
    }

    // Initialise the dataset indexes for all MPI processes
    std::vector<uint64_t> idx = m_comm->Scatter(root_rank, all_idxs);
    uint64_t ids_i            = idx[IDS_IDX_IDX];
    uint64_t data_i           = idx[DATA_IDX_IDX];
    uint64_t order_i          = idx[ORDER_IDX_IDX];
    uint64_t homy_i           = idx[HOMY_IDX_IDX];
    uint64_t homz_i           = idx[HOMZ_IDX_IDX];
    uint64_t homs_i           = idx[HOMS_IDX_IDX];

    // Set properties for parallel file access (if we're in parallel)
    H5::PListSharedPtr parallelProps = H5::PList::Default();
    if (m_comm->GetSize() > 1)
    {
        // Use MPI/O to access the file
        parallelProps = H5::PList::FileAccess();
        parallelProps->SetMpio(m_comm);
    }

    // Reopen the file
    H5::FileSharedPtr outfile =
        H5::File::Open(outFile, H5F_ACC_RDWR, parallelProps);
    ASSERTL1(outfile, prfx.str() + "cannot open HDF5 file.");
    H5::GroupSharedPtr root = outfile->OpenGroup("NEKTAR");
    ASSERTL1(root, prfx.str() + "cannot open root group.");

    m_comm->Block();

    // all HDF5 groups have now been created. Open the IDS dataset and
    // associated data space
    H5::DataSetSharedPtr ids_dset = root->OpenDataSet("ELEMENTIDS");
    ASSERTL1(ids_dset, prfx.str() + "cannot open ELEMENTIDS dataset.");
    H5::DataSpaceSharedPtr ids_fspace = ids_dset->GetSpace();
    ASSERTL1(ids_fspace, prfx.str() + "cannot open ELEMENTIDS filespace.");

    // Open the DATA dataset and associated data space
    H5::DataSetSharedPtr data_dset = root->OpenDataSet("DATA");
    ASSERTL1(data_dset, prfx.str() + "cannot open DATA dataset.");
    H5::DataSpaceSharedPtr data_fspace = data_dset->GetSpace();
    ASSERTL1(data_fspace, prfx.str() + "cannot open DATA filespace.");

    // Open the optional datasets and data spaces.
    H5::DataSetSharedPtr order_dset, homy_dset, homz_dset, homs_dset;
    H5::DataSpaceSharedPtr order_fspace , homy_fspace, homz_fspace, homs_fspace;

    if (all_dsetsize[ORDER_CNT_IDX])
    {
        order_dset = root->OpenDataSet("POLYORDERS");
        ASSERTL1(order_dset, prfx.str() + "cannot open POLYORDERS dataset.");
        order_fspace = order_dset->GetSpace();
        ASSERTL1(order_fspace, prfx.str() + "cannot open POLYORDERS filespace.");
    }

    if (all_dsetsize[HOMY_CNT_IDX])
    {
        homy_dset = root->OpenDataSet("HOMOGENEOUSYIDS");
        ASSERTL1(homy_dset, prfx.str() + "cannot open HOMOGENEOUSYIDS dataset.");
        homy_fspace = homy_dset->GetSpace();
        ASSERTL1(homy_fspace, prfx.str() + "cannot open HOMOGENEOUSYIDS filespace.");
    }

    if (all_dsetsize[HOMZ_CNT_IDX])
    {
        homz_dset   = root->OpenDataSet("HOMOGENEOUSZIDS");
        ASSERTL1(homz_dset, prfx.str() + "cannot open HOMOGENEOUSZIDS dataset.");
        homz_fspace = homz_dset->GetSpace();
        ASSERTL1(homz_fspace, prfx.str() + "cannot open HOMOGENEOUSZIDS filespace.");
    }

    if (all_dsetsize[HOMS_CNT_IDX])
    {
        homs_dset   = root->OpenDataSet("HOMOGENEOUSSIDS");
        ASSERTL1(homs_dset, prfx.str() + "cannot open HOMOGENEOUSSIDS dataset.");
        homs_fspace = homs_dset->GetSpace();
        ASSERTL1(homs_fspace, prfx.str() + "cannot open HOMOGENEOUSSIDS filespace.");
    }

    
    // we use independent (non-collective) writes for all datasets that are actually part of the field definition,
    // i.e., {order_dset,homz_dset,homy_dset,homs_dset}; otherwise, we would have to assume that all fields
    // were defined such that they all used the same selection of field definition datasets
    nWritten += WriteFieldDataInd(nFields, order_fspace, order_dset, order_i, numModesPerDirVar);
    nWritten += WriteFieldDataInd(nFields, homy_fspace, homy_dset, homy_i, homoYIDs);
    nWritten += WriteFieldDataInd(nFields, homz_fspace, homz_dset, homz_i, homoZIDs);
    nWritten += WriteFieldDataInd(nFields, homs_fspace, homs_dset, homs_i, homoSIDs);
    
    // write all the element IDs and element data collectively, taking into account the fact that
    // different processes maybe handling different numbers of fields
    std::vector<std::vector<unsigned int> > elemIDs(nFields);
    for (std::size_t f = 0; f < nFields; ++f)
    {
        elemIDs[f] = fielddefs[f]->m_elementIDs;
    }
    nWritten += WriteFieldData(nMinFields, nFields, ids_fspace, ids_dset, ids_i, elemIDs);
    nWritten += WriteFieldData(nMinFields, nFields, data_fspace, data_dset, data_i, fielddata);
            
    
    m_comm->Block();

    // all data has been written
    //if (m_comm->GetRank() == 0)
    //{
        //tm1 = m_comm->Wtime();
        //cout << " (" << tm1 - tm0 << "s, HDF5)" << endl;
    //}

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
                                        H5::DataSpaceSharedPtr &fspace, H5::DataSetSharedPtr &dset,
                                        uint64_t data_i, std::vector<std::vector<T> > &data)
{
    if (!fspace || !dset) return 0;

    std::size_t nDataItems = 0;
    uint64_t nWritten = 0;

    H5::PListSharedPtr writePL = H5::PList::DatasetXfer();
    writePL->SetDxMpioIndependent();

    for (std::size_t f = 0; f < nFields; ++f)
    {
        nDataItems = data[f].size();
        if (nDataItems > 0)
        {
            fspace->SelectRange(data_i, nDataItems);
            dset->Write(data[f], fspace, writePL);
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
 * @param nMinFields the lowest number of fields handled by any process
 * @param nFields    the number of fields handled by this process
 * @param fspace     hdf5 file space
 * @param dset       hdf5 data set
 * @param data_i     starting offset into hdf5 data set for this process
 * @param data       vector of data vectors (one for each field) to be written to data set
 * @return The number of bytes written.
 */
template <class T>
uint64_t FieldIOHdf5::WriteFieldData(std::size_t nMinFields, std::size_t nFields,
                                     H5::DataSpaceSharedPtr &fspace, H5::DataSetSharedPtr &dset,
                                     uint64_t data_i, std::vector<std::vector<T> > &data)
{
    if (!fspace || !dset) return 0;

    bool concatenate_last_fields = nMinFields < nFields;
    std::size_t nFirstFields = concatenate_last_fields ? nMinFields-1 : nFields;
    std::size_t nDataItems = 0;
    std::size_t f = 0;
    uint64_t nWritten = 0;

    H5::PListSharedPtr writePL = H5::PList::DatasetXfer();
    writePL->SetDxMpioCollective();

    for (; f < nFirstFields; ++f)
    {
        nDataItems = data[f].size();
        fspace->SelectRange(data_i, nDataItems);
        dset->Write(data[f], fspace, writePL);
        data_i += nDataItems;
        nWritten += nDataItems*sizeof(T);
    }

    if (!concatenate_last_fields) return nWritten;

    std::vector<T> last_data;

    nDataItems = data[f].size();
    fspace->SelectRange(H5S_SELECT_SET, data_i, nDataItems);
    last_data.insert(last_data.end(), data[f].begin(), data[f].end());
    data_i += nDataItems;
    f++;
        
    for (; f < nFields; ++f)
    {
        nDataItems = data[f].size();
        fspace->SelectRange(H5S_SELECT_OR, data_i, nDataItems);
        last_data.insert(last_data.end(), data[f].begin(), data[f].end());
        data_i += nDataItems;
    }

    dset->Write(last_data, fspace, writePL);
    return nWritten + last_data.size()*sizeof(T);
}


/**
 * @brief Import a HDF5 format file.
 *
 * @param finfilename       Input filename
 * @param fielddefs         Field definitions of resulting field
 * @param fielddata         Field data of resulting field
 * @param fieldinfomap      Field metadata of resulting field
 * @param ElementIDs        If specified, contains the list of element IDs on
 *                          this rank. The resulting field definitions will only
 *                          contain data for the element IDs specified in this
 *                          array.
 * @return The number of bytes read.
 */
uint64_t FieldIOHdf5::v_Import(const std::string &infilename,
                               std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                               std::vector<std::vector<NekDouble> > &fielddata,
                               FieldMetaDataMap &fieldinfomap,
                               const Array<OneD, int> &ElementIDs)
{
    std::stringstream prfx;
    prfx << m_comm->GetRank() << ": FieldIOHdf5::v_Import(): ";
    int nRanks = m_comm->GetSize();

    uint64_t nRead = 0;
    
    // Set properties for parallel file access (if we're in parallel)
    H5::PListSharedPtr parallelProps = H5::PList::Default();
    H5::PListSharedPtr readPL = H5::PList::Default();
    
    if (nRanks > 1)
    {
        // Use MPI/O to access the file
        parallelProps = H5::PList::FileAccess();
        parallelProps->SetMpio(m_comm);
        // Use collective IO
        readPL = H5::PList::DatasetXfer();
        readPL->SetDxMpioCollective();
    }

    DataSourceSharedPtr dataSource = H5DataSource::create(
        infilename, parallelProps);

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
             "File format if " + infilename + " is higher than supported in "
             "this version of Nektar++");

    std::unordered_set<uint64_t> selectedIDs;
    bool selective = false;
    bool auto_selective = false;
    if (&ElementIDs != &NullInt1DArray)
    {
        if (ElementIDs.num_elements() > 0)
        {
            selective = true;
            for (uint64_t i = 0; i < ElementIDs.num_elements(); ++i)
            {
                selectedIDs.insert(ElementIDs[i]);
            }
        }
        else
        {
            auto_selective = true;
        }
    }

    // Open and read in the entire the decomps dataset.
    H5::DataSetSharedPtr decomps_dset = root->OpenDataSet("DECOMPOSITION");
    ASSERTL1(decomps_dset, "Cannot open DECOMPOSITION dataset.");
    H5::DataSpaceSharedPtr decomps_fspace = decomps_dset->GetSpace();
    ASSERTL1(decomps_fspace, "Cannot open DECOMPOSITION filespace.");
    
    std::vector<uint64_t> decomps;
    decomps_dset->Read(decomps, decomps_fspace, readPL);
    nRead += decomps.size()*sizeof(uint64_t);
    
    uint64_t nDecomps = decomps.size() / MAX_DCMPS;
    uint64_t i_autosel_min = 0, i_autosel_lim = 0;
    if (auto_selective)
    {
        uint64_t nDecompsPerRank = nDecomps / nRanks;
        i_autosel_min = m_comm->GetRank()*nDecompsPerRank;
        i_autosel_lim = i_autosel_min + nDecompsPerRank; 
    }


    // Open the element ids dataset.
    H5::DataSetSharedPtr ids_dset = root->OpenDataSet("ELEMENTIDS");
    ASSERTL1(ids_dset, "Cannot open ELEMENTIDS dataset.");
    H5::DataSpaceSharedPtr ids_fspace = ids_dset->GetSpace();
    ASSERTL1(ids_fspace, "Cannot open ELEMENTIDS filespace.");


    // Mapping from each decomposition to offsets in the data array.
    vector<OffsetHelper> decompsToOffsets(nDecomps);
    
    // Mapping from each group's hash to a vector of element IDs. Note this has
    // to be unsigned int, since that's what we use in FieldDefinitions.
    map<uint64_t, vector<unsigned int> > groupsToElems;
    
    // Mapping from each group's hash to each of the decompositions.
    map<uint64_t, set<uint64_t> > groupsToDecomps;
    
    // Counters for data offsets
    OffsetHelper running;

    H5S_seloper_t h5_select = H5S_SELECT_SET;

    uint64_t nAutoSelectedElems = 0;
    uint64_t cnt = 0;
    for (uint64_t i = 0; i < nDecomps; ++i, cnt += MAX_DCMPS)
    {
        uint64_t nElems = decomps[cnt + ELEM_DCMP_IDX];
	uint64_t groupHash = decomps[cnt + HASH_DCMP_IDX];

	if (auto_selective)
        {
            // auto_selective means the number of ranks is equal to
            // that used to write the checkpoint data.
            if (i < i_autosel_min)
            {
                // Have not yet reached this rank's
                // decomposition data.
                nElems = 0;
            }
            else if (i >= i_autosel_lim)
            {
                // Have reached the end of this rank's
                // decomposition data.
                break;
            }
        }

        if (nElems > 0)
        {
	    if (auto_selective)
	    {
	        nAutoSelectedElems += nElems;
	        ids_fspace->SelectRange(h5_select, running.ids, nElems);
                h5_select = H5S_SELECT_OR;

                groupsToDecomps[groupHash].insert(i);
                decompsToOffsets[i] = running;
	    }
            else
	    {
	        vector<unsigned int> ids(nElems);
		ids_fspace->SelectRange(running.ids, nElems);
                ids_dset->Read(ids, ids_fspace, readPL);
                nRead += ids.size()*sizeof(unsigned int);

                bool found = !selective;
                if (selective)
                {
		    for (auto &it : ids)
                    {
		        if (selectedIDs.find(it) != selectedIDs.end())
                        {
			    found = true; 
			    break;
                        }
		    }    
                }

                if (found)
                {  
                    groupsToDecomps[groupHash].insert(i);
                    groupsToElems[i] = ids;
                    decompsToOffsets[i] = running;

                    if (selective && selectedIDs.size() == 0)
                    {
                        break;
                    }
                }
	    } // end of <if (auto_selective)> else clause
        } // end of <if (nElems > 0)> clause

	running.ids   += decomps[cnt + ELEM_DCMP_IDX];
        running.data  += decomps[cnt + VAL_DCMP_IDX];
        running.order += decomps[cnt + ORDER_DCMP_IDX];
        running.homy  += decomps[cnt + HOMY_DCMP_IDX];
        running.homz  += decomps[cnt + HOMZ_DCMP_IDX];
        running.homs  += decomps[cnt + HOMS_DCMP_IDX];
	
    } // end of <for (uint64_t i = 0; i < nDecomps; ++i, cnt += MAX_DCMPS)> loop

    if (auto_selective)
    {
        vector<unsigned int> ids(nAutoSelectedElems);
        ids_dset->Read(ids, ids_fspace, readPL);

	uint64_t elemOffset = 0;
	for (auto &gIt : groupsToDecomps)
        {
            for (auto &di : gIt.second)
            {
                uint64_t nElems = decomps[di*MAX_DCMPS + ELEM_DCMP_IDX];
                vector<unsigned int> sids(&ids[elemOffset],&ids[elemOffset+nElems]);
                groupsToElems[di] = sids;
		elemOffset += nElems;
            }  
        }
    }
    
 
    H5::DataSetSharedPtr data_dset;
    H5::DataSpaceSharedPtr data_fspace;
    data_dset = root->OpenDataSet("DATA");
    ASSERTL1(data_dset, prfx.str() + "cannot open DATA dataset.");
    data_fspace = data_dset->GetSpace();
    ASSERTL1(data_fspace, prfx.str() + "cannot open DATA filespace.");
    

    // chain together all the necessary HDF5 range selections
    uint64_t nDataItems = 0;
    uint64_t nRealDataItems = 0;
    h5_select = H5S_SELECT_SET;
    map<uint64_t, uint64_t> dataCountMap;
    
    for (auto &gIt : groupsToDecomps)
    {
        for (auto &di : gIt.second)
        {
            OffsetHelper offset = decompsToOffsets[di];

            nDataItems = decomps[di*MAX_DCMPS + VAL_DCMP_IDX];
            data_fspace->SelectRange(h5_select, offset.data, nDataItems);
            dataCountMap[offset.data] = nDataItems;
            nRealDataItems += nDataItems;
            h5_select = H5S_SELECT_OR;            
        }  
    }

    // now read the actual HDF5 data
    std::vector<NekDouble> decompRealData;
    if (nRealDataItems > 0)
    {
        data_dset->Read(decompRealData, data_fspace, readPL);
        ASSERTL0(decompRealData.size() == nRealDataItems,
            prfx.str() + "unexpected number of data items.");
        nRead += nRealDataItems;
    }


    // calculate the adjusted offsets for the recently-read data
    map<uint64_t, uint64_t> dataOffsetMap;
    uint64_t offVal;

    offVal = 0;
    for (auto &offIt : dataCountMap)
    {
        dataOffsetMap[offIt.first] = offVal;
        offVal += offIt.second;
    }

    // transfer the recently read data to the fielddata data structures
    // and read in the rest of the field definition data
    for (auto &gIt : groupsToDecomps)
    {
        std::stringstream fieldNameStream;
        fieldNameStream << gIt.first;

        for (auto &di : gIt.second)
        {
            FieldDefinitionsSharedPtr fielddef =
                MemoryManager<FieldDefinitions>::AllocateSharedPtr();
            nRead += ImportFieldDef(root, fieldNameStream.str(), fielddef);

            OffsetHelper offset = decompsToOffsets[di];

            if (fielddef->m_numHomogeneousDir >= 1)
            {
                nRead += ImportFieldDataInd(root, "HOMOGENEOUSZIDS", "homogeneousZIDs",
                             decomps[di*MAX_DCMPS + HOMZ_DCMP_IDX], offset.homz,
                             fielddef->m_homogeneousZIDs);
            }

            if (fielddef->m_numHomogeneousDir >= 2)
            {
                nRead += ImportFieldDataInd(root, "HOMOGENEOUSYIDS", "homogeneousYIDs",
                             decomps[di*MAX_DCMPS + HOMY_DCMP_IDX], offset.homy,
                             fielddef->m_homogeneousYIDs);
            }

            if (fielddef->m_homoStrips)
            {
                nRead += ImportFieldDataInd(root, "HOMOGENEOUSSIDS", "homogeneousSIDs",
                             decomps[di*MAX_DCMPS + HOMS_DCMP_IDX], offset.homs,
                             fielddef->m_homogeneousSIDs);
            }

            if (!fielddef->m_uniOrder)
            {
                nRead += ImportFieldDataInd(root, "POLYORDERS", "numModes",
                             decomps[di*MAX_DCMPS + ORDER_DCMP_IDX], offset.order,
                             fielddef->m_numModes);
            }
            
            fielddef->m_elementIDs = groupsToElems[di];
            fielddefs.push_back(fielddef);

            if (fielddata != NullVectorNekDoubleVector && nRealDataItems > 0)
            {
                vector<NekDouble>::const_iterator dataSelStart = decompRealData.begin() + dataOffsetMap[offset.data]; 
                vector<NekDouble>::const_iterator dataSelEnd = dataSelStart + dataCountMap[offset.data];
                vector<NekDouble> dataSel(dataSelStart, dataSelEnd);

                uint64_t datasize = CheckFieldDefinition(fielddef);
                ASSERTL0(dataSel.size() == datasize * fielddef->m_fields.size(),
                    prfx.str() + "input data is not the same length as header information.");
                
                fielddata.push_back(dataSel);
            }
            
        }  // end of <for (auto &di : gIt.second)> loop
        
    } // end of <for (auto &gIt : groupsToDecomps)> loop

    m_comm->Block();

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
