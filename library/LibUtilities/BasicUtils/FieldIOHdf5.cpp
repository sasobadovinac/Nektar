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
#include <boost/unordered_set.hpp>

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
H5TagWriter::H5TagWriter(H5::GroupSharedPtr grp) : m_Group(grp)
{
}
TagWriterSharedPtr H5TagWriter::AddChild(const std::string &name)
{
    H5::GroupSharedPtr child = m_Group->CreateGroup(name);
    return TagWriterSharedPtr(new H5TagWriter(child));
};

void H5TagWriter::SetAttr(const std::string &key, const std::string &val)
{
    m_Group->SetAttribute(key, val);
}
class H5DataSource : public DataSource
{
    H5::FileSharedPtr doc;

public:
    H5DataSource(const std::string &fn)
        : doc(H5::File::Open(fn, H5F_ACC_RDONLY))
    {
    }

    H5::FileSharedPtr Get()
    {
        return doc;
    }
    const H5::FileSharedPtr Get() const
    {
        return doc;
    }

    static DataSourceSharedPtr create(const std::string &fn)
    {
        return DataSourceSharedPtr(new H5DataSource(fn));
    }
};

typedef boost::shared_ptr<H5DataSource> H5DataSourceSharedPtr;

std::string FieldIOHdf5::className =
    GetFieldIOFactory().RegisterCreatorFunction(
        "Hdf5", FieldIOHdf5::create, "HDF5-based output of field data.");

const unsigned int FieldIOHdf5::ELEM_DCMP_IDX = 0;
const unsigned int FieldIOHdf5::VAL_DCMP_IDX  = 1;
const unsigned int FieldIOHdf5::HASH_DCMP_IDX = 2;
const unsigned int FieldIOHdf5::MAX_DCMPS     = FieldIOHdf5::HASH_DCMP_IDX + 1;
const unsigned int FieldIOHdf5::ELEM_CNT_IDX  = 0;
const unsigned int FieldIOHdf5::VAL_CNT_IDX   = 1;
const unsigned int FieldIOHdf5::MAX_CNTS      = FieldIOHdf5::VAL_CNT_IDX + 1;
const unsigned int FieldIOHdf5::IDS_IDX_IDX   = 0;
const unsigned int FieldIOHdf5::DATA_IDX_IDX  = 1;
const unsigned int FieldIOHdf5::MAX_IDXS      = FieldIOHdf5::DATA_IDX_IDX + 1;

FieldIOHdf5::FieldIOHdf5(LibUtilities::CommSharedPtr pComm,
                         bool sharedFilesystem)
    : FieldIO(pComm, sharedFilesystem)
{
}

void FieldIOHdf5::v_Write(const std::string &outFile,
                          std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                          std::vector<std::vector<NekDouble> > &fielddata,
                          const FieldMetaDataMap &fieldmetadatamap)
{
    std::stringstream prfx;
    prfx << m_comm->GetRank() << ": FieldIOHdf5::v_Write(): ";
    double tm0 = 0.0, tm1 = 0.0;
    if (0 == m_comm->GetRank())
    {
        cout << prfx.str() << "entering..." << endl;
        tm0 = m_comm->Wtime();
    }

    ///////////////////////////////////////////////////////////////////////////////////
    // ASSUMPTION 1: all element ids have the same type, unsigned int
    ///////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////
    // ASSUMPTION 2: all elements within a given field have the same number of
    // values
    ///////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////
    // ASSUMPTION 3: all element values have the same type, NekDouble (i.e.,
    // double)
    ///////////////////////////////////////////////////////////////////////////////////

    // determine the root MPI process, i.e., the lowest ranked process handling
    // nMaxFields fields
    /////////////////////////////////////////////////////////////////////////////////////////////
    ASSERTL1(fielddefs.size() == fielddata.size(),
             prfx.str() + "fielddefs and fielddata have incompatible lengths.");

    size_t nFields    = fielddefs.size();
    size_t nMaxFields = nFields;
    m_comm->AllReduce(nMaxFields, LibUtilities::ReduceMax);

    int root_rank = -1;
    bool amRoot = false;
    LibUtilities::CommSharedPtr max_fields_comm =
        m_comm->CommCreateIf((nFields == nMaxFields) ? 1 : 0);
    if (max_fields_comm)
    {
        int rank  = m_comm->GetRank();
        root_rank = rank;
        max_fields_comm->AllReduce(root_rank, LibUtilities::ReduceMin);
        amRoot = (rank == root_rank);
        if (!amRoot)
            root_rank = -1;
    }
    m_comm->AllReduce(root_rank, LibUtilities::ReduceMax);
    ASSERTL1(root_rank >= 0 && root_rank < m_comm->GetSize(),
             prfx.str() + "invalid root rank.");
    /////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<std::size_t> decomps(nMaxFields * MAX_DCMPS, 0);
    std::vector<std::size_t> all_hashes(nMaxFields * m_comm->GetSize(), 0);
    std::vector<std::size_t> cnts(MAX_CNTS, 0);
    std::vector<std::string> fieldNames(nFields);
    std::vector<std::string> shapeStrings(nFields);
    std::vector<std::vector<NekDouble> > homoLengths(nFields);
    std::vector<std::vector<unsigned int> > homoSIDs(nFields),
        homoYIDs(nFields), homoZIDs(nFields);
    std::vector<std::string> numModesPerDirs(nFields);

    // calculate the total number of elements handled by this MPI process
    // and the total number of bytes required to store the elements,
    // base the name of each field on the hash of the field definition
    /////////////////////////////////////////////////////////////////////////////////////////
    for (int f = 0; f < nFields; ++f)
    {
        ASSERTL1(fielddata[f].size() > 0,
                 prfx.str() +
                     "fielddata vector must contain at least one value.");
        ASSERTL1(fielddata[f].size() ==
                     fielddefs[f]->m_fields.size() *
                         CheckFieldDefinition(fielddefs[f]),
                 prfx.str() + "fielddata vector has invalid size.");

        std::size_t nFieldElems = fielddefs[f]->m_elementIDs.size();
        std::size_t nElemVals   = fielddata[f].size() / nFieldElems;

        decomps[f * MAX_DCMPS + ELEM_DCMP_IDX] = nFieldElems;
        decomps[f * MAX_DCMPS + VAL_DCMP_IDX]  = nElemVals;

        cnts[ELEM_CNT_IDX] += nFieldElems;
        cnts[VAL_CNT_IDX] += nElemVals * nFieldElems;

        // hash the field specification
        /////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////
        std::stringstream hashStream;
        std::size_t nSubFields = fielddefs[f]->m_fields.size();
        for (int sf = 0; sf < nSubFields; ++sf)
            hashStream << fielddefs[f]->m_fields[sf];

        nSubFields = fielddefs[f]->m_basis.size();
        for (int sf = 0; sf < nSubFields; ++sf)
            hashStream << fielddefs[f]->m_basis[sf];

        // determine SHAPE attribute
        /////////////////////////////////////////////////////////////////////////////////////////
        std::stringstream shapeStringStream;
        shapeStringStream << ShapeTypeMap[fielddefs[f]->m_shapeType];
        if (fielddefs[f]->m_numHomogeneousDir == 1)
        {
            shapeStringStream << "-HomogenousExp1D";
        }
        else if (fielddefs[f]->m_numHomogeneousDir == 2)
        {
            shapeStringStream << "-HomogenousExp2D";
        }
        if (fielddefs[f]->m_homoStrips)
        {
            shapeStringStream << "-Strips";
        }
        shapeStrings[f] = shapeStringStream.str();
        hashStream << shapeStringStream.str();
        /////////////////////////////////////////////////////////////////////////////////////////

        // determine HOMOGENEOUS attributes
        /////////////////////////////////////////////////////////////////////////////////////////
        if (fielddefs[f]->m_numHomogeneousDir)
        {
            nSubFields = fielddefs[f]->m_homogeneousLengths.size();
            homoLengths[f].resize(nSubFields);
            for (int sf = 0; sf < nSubFields; ++sf)
            {
                std::size_t len = fielddefs[f]->m_homogeneousLengths[sf];
                hashStream << len;
                homoLengths[f][sf] = len;
            }

            nSubFields = fielddefs[f]->m_homogeneousYIDs.size();
            if (nSubFields > 0)
            {
                homoYIDs[f].resize(nSubFields);
                for (int sf = 0; sf < nSubFields; ++sf)
                {
                    unsigned int yid = fielddefs[f]->m_homogeneousYIDs[sf];
                    hashStream << yid;
                    homoYIDs[f][sf] = yid;
                }
            }

            nSubFields = fielddefs[f]->m_homogeneousZIDs.size();
            if (nSubFields > 0)
            {
                homoZIDs[f].resize(nSubFields);
                for (int sf = 0; sf < nSubFields; ++sf)
                {
                    unsigned int zid = fielddefs[f]->m_homogeneousZIDs[sf];
                    hashStream << zid;
                    homoZIDs[f][sf] = zid;
                }
            }

            nSubFields = fielddefs[f]->m_homogeneousSIDs.size();
            if (nSubFields > 0)
            {
                homoSIDs[f].resize(nSubFields);
                for (int sf = 0; sf < nSubFields; ++sf)
                {
                    unsigned int sid = fielddefs[f]->m_homogeneousSIDs[sf];
                    hashStream << sid;
                    homoSIDs[f][sf] = sid;
                }
            }
        }
        /////////////////////////////////////////////////////////////////////////////////////////

        // determine NUMMODESPERDIR attribute
        /////////////////////////////////////////////////////////////////////////////////////////
        std::stringstream numModesStringStream;
        if (fielddefs[f]->m_uniOrder)
        {
            numModesStringStream << "UNIORDER:";
            for (std::vector<int>::size_type i = 0;
                 i < fielddefs[f]->m_basis.size();
                 i++)
            {
                if (i > 0)
                    numModesStringStream << ",";
                numModesStringStream << fielddefs[f]->m_numModes[i];
            }
        }
        else
        {
            numModesStringStream << "MIXORDER:";
            for (std::vector<int>::size_type i = 0;
                 i < fielddefs[f]->m_numModes.size();
                 i++)
            {
                if (i > 0)
                    numModesStringStream << ",";
                numModesStringStream << fielddefs[f]->m_numModes[i];
            }
        }
        numModesPerDirs[f] = numModesStringStream.str();
        hashStream << numModesStringStream.str();
        /////////////////////////////////////////////////////////////////////////////////////////

        boost::hash<std::string> string_hasher;
        std::stringstream fieldNameStream;
        std::size_t fieldDefHash               = string_hasher(hashStream.str());
        decomps[f * MAX_DCMPS + HASH_DCMP_IDX] = fieldDefHash;
        all_hashes[m_comm->GetRank() * nMaxFields + f] = fieldDefHash;
        fieldNameStream << fieldDefHash;
        fieldNames[f] = fieldNameStream.str();
        /////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////
    } // end of <for (int f = 0; f < nFields; ++f)>
    /////////////////////////////////////////////////////////////////////////////////////////

    // gather information from all MPI processes
    /////////////////////////////////////////////////////////////////////////////////////////
    std::vector<std::size_t> all_cnts = m_comm->Gather(root_rank, cnts);
    std::vector<std::size_t> all_idxs(m_comm->GetSize() * MAX_IDXS, 0);
    std::vector<std::size_t> all_decomps = m_comm->Gather(root_rank, decomps);
    /////////////////////////////////////////////////////////////////////////////////////////

    // the root rank creates the file layout from scratch
    /////////////////////////////////////////////////////////////////////////////////////////
    if (amRoot)
    {
        H5::FileSharedPtr outfile = H5::File::Create(outFile, H5F_ACC_TRUNC);
        ASSERTL1(outfile, prfx.str() + "cannot create HDF5 file.");
        H5::GroupSharedPtr root = outfile->CreateGroup("NEKTAR");
        ASSERTL1(root, prfx.str() + "cannot create root group.");
        TagWriterSharedPtr info_writer(new H5TagWriter(root));
        AddInfoTag(info_writer, fieldmetadatamap);

        // calculate the indexes to be used by each MPI process when reading the
        // IDS and DATA datasets
        /////////////////////////////////////////////////////////////////////////////////////////
        std::size_t nTotElems = 0, nTotVals = 0;
        int nRanks = m_comm->GetSize();
        for (int r = 0; r < nRanks; ++r)
        {
            all_idxs[r * MAX_IDXS + IDS_IDX_IDX]  = nTotElems;
            all_idxs[r * MAX_IDXS + DATA_IDX_IDX] = nTotVals;

            nTotElems += all_cnts[r * MAX_CNTS + ELEM_CNT_IDX];
            nTotVals += all_cnts[r * MAX_CNTS + VAL_CNT_IDX];
        }
        /////////////////////////////////////////////////////////////////////////////////////////

        // create DECOMPOSITION dataset: basic field info for each MPI process
        /////////////////////////////////////////////////////////////////////////////////////////
        H5::DataTypeSharedPtr decomps_type =
            H5::DataType::OfObject(all_decomps[0]);
        H5::DataSpaceSharedPtr decomps_space =
            H5::DataSpace::OneD(all_decomps.size());
        H5::DataSetSharedPtr decomps_dset =
            root->CreateDataSet("DECOMPOSITION", decomps_type, decomps_space);
        ASSERTL1(decomps_dset,
                 prfx.str() + "cannot create DECOMPOSITION dataset.");
        /////////////////////////////////////////////////////////////////////////////////////////

        // create INDEXES dataset: MPI process indexes into IDS and DATA
        // datasets
        /////////////////////////////////////////////////////////////////////////////////////////
        H5::DataTypeSharedPtr idxs_type = H5::DataType::OfObject(all_idxs[0]);
        H5::DataSpaceSharedPtr idxs_space =
            H5::DataSpace::OneD(all_idxs.size());
        H5::DataSetSharedPtr idxs_dset =
            root->CreateDataSet("INDEXES", idxs_type, idxs_space);
        ASSERTL1(idxs_dset, prfx.str() + "cannot create INDEXES dataset.");
        /////////////////////////////////////////////////////////////////////////////////////////

        // create IDS dataset: element ids
        /////////////////////////////////////////////////////////////////////////////////////////
        H5::DataTypeSharedPtr ids_type =
            H5::DataType::OfObject(fielddefs[0]->m_elementIDs[0]);
        H5::DataSpaceSharedPtr ids_space = H5::DataSpace::OneD(nTotElems);
        H5::DataSetSharedPtr ids_dset =
            root->CreateDataSet("IDS", ids_type, ids_space);
        ASSERTL1(ids_dset, prfx.str() + "cannot create IDS dataset.");
        /////////////////////////////////////////////////////////////////////////////////////////

        // create DATA dataset: element data
        /////////////////////////////////////////////////////////////////////////////////////////
        H5::DataTypeSharedPtr data_type =
            H5::DataType::OfObject(fielddata[0][0]);
        H5::DataSpaceSharedPtr data_space = H5::DataSpace::OneD(nTotVals);
        H5::DataSetSharedPtr data_dset =
            root->CreateDataSet("DATA", data_type, data_space);
        ASSERTL1(data_dset, prfx.str() + "cannot create DATA dataset.");
        /////////////////////////////////////////////////////////////////////////////////////////

        // RAII => field group is closed automatically

    } // if (amRoot)
    // RAII => datasets (ids_dset, data_dset), root group and HDF5 file are all
    // closed automatically
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Determine which process will write the group representing the field
    // description in the HDF5 file.
    std::map<size_t, int> hashToProc;
    std::map<int, std::vector<size_t> > writingProcs;

    // Gather all field hashes to every processor.
    m_comm->AllReduce(all_hashes, LibUtilities::ReduceMax);

    std::cout << "RANK " << m_comm->GetRank() << " HAS: ";
    for (int i = 0; i < all_hashes.size(); ++i)
    {
        std::cout << all_hashes[i] << " ";
    }
    std::cout << std::endl;

    for (int n = 0; n < m_comm->GetSize(); ++n)
    {
        for (int i = 0; i < nMaxFields; ++i)
        {
            size_t hash = all_hashes[n*nMaxFields + i];
            if (hashToProc.find(hash) != hashToProc.end() || hash == 0)
            {
                continue;
            }
            hashToProc[hash] = m_comm->GetRank();
            writingProcs[n].push_back(hash);
        }
    }

    map<int, std::vector<size_t> >::iterator sIt;
    for (sIt = writingProcs.begin(); sIt != writingProcs.end(); sIt++)
    {
        int rank = sIt->first;
        if (m_comm->GetRank() == rank)
        {
            H5::PListSharedPtr serialProps = H5::PList::Default();
            H5::PListSharedPtr writeSR     = H5::PList::Default();

            // reopen the file
            H5::FileSharedPtr outfile =
                H5::File::Open(outFile, H5F_ACC_RDWR, serialProps);
            ASSERTL1(outfile, prfx.str() + "cannot open HDF5 file.");
            H5::GroupSharedPtr root = outfile->OpenGroup("NEKTAR");
            ASSERTL1(root, prfx.str() + "cannot open root group.");

            // write a hdf5 group for each field
            for (int i = 0; i < sIt->second.size(); ++i)
            {
                for (int f = 0; f < nFields; ++f)
                {
                    if (sIt->second[i] != all_hashes[m_comm->GetRank() * nMaxFields + f])
                    {
                        continue;
                    }

                    std::cout << "RANK " << rank << " WRITING HASH " << sIt->second[i] << std::endl;
                    H5::GroupSharedPtr field_group = root->CreateGroup(fieldNames[f]);
                    ASSERTL1(field_group, prfx.str() + "cannot create field group.");
                    field_group->SetAttribute("FIELDS", fielddefs[f]->m_fields);
                    field_group->SetAttribute("BASIS", fielddefs[f]->m_basis);
                    field_group->SetAttribute("SHAPE", shapeStrings[f]);
                    if (homoLengths[f].size() > 0)
                        field_group->SetAttribute("HOMOGENEOUSLENGTHS", homoLengths[f]);
                    if (homoYIDs[f].size() > 0)
                        field_group->SetAttribute("HOMOGENEOUSYIDS", homoYIDs[f]);
                    if (homoZIDs[f].size() > 0)
                        field_group->SetAttribute("HOMOGENEOUSZIDS", homoZIDs[f]);
                    if (homoSIDs[f].size() > 0)
                        field_group->SetAttribute("HOMOGENEOUSSIDS", homoSIDs[f]);
                    field_group->SetAttribute("NUMMODESPERDIR", numModesPerDirs[f]);
                }
            }
        }
        m_comm->Block();
    }

    // write the DECOMPOSITION and INDEXES datasets
    /////////////////////////////////////////////////////////////////////////////////////////
    if (amRoot)
    {
        H5::PListSharedPtr serialProps = H5::PList::Default();
        H5::PListSharedPtr writeSR     = H5::PList::Default();

        // reopen the file
        H5::FileSharedPtr outfile =
            H5::File::Open(outFile, H5F_ACC_RDWR, serialProps);
        ASSERTL1(outfile, prfx.str() + "cannot open HDF5 file.");
        H5::GroupSharedPtr root = outfile->OpenGroup("NEKTAR");
        ASSERTL1(root, prfx.str() + "cannot open root group.");

        // write the DECOMPOSITION dataset
        /////////////////////////////////////////////////////////////////////////////////////////
        H5::DataSetSharedPtr decomps_dset = root->OpenDataSet("DECOMPOSITION");
        ASSERTL1(decomps_dset,
                 prfx.str() + "cannot open DECOMPOSITION dataset.");

        H5::DataSpaceSharedPtr decomps_fspace = decomps_dset->GetSpace();
        ASSERTL1(decomps_fspace,
                 prfx.str() + "cannot open DECOMPOSITION filespace.");

        decomps_fspace->SelectRange(0, all_decomps.size());
        decomps_dset->Write(all_decomps, decomps_fspace, writeSR);
        /////////////////////////////////////////////////////////////////////////////////////////

        // write the INDEXES dataset
        /////////////////////////////////////////////////////////////////////////////////////////
        H5::DataSetSharedPtr idxs_dset = root->OpenDataSet("INDEXES");
        ASSERTL1(idxs_dset, prfx.str() + "cannot open INDEXES dataset.");

        H5::DataSpaceSharedPtr idxs_fspace = idxs_dset->GetSpace();
        ASSERTL1(idxs_fspace, prfx.str() + "cannot open INDEXES filespace.");

        idxs_fspace->SelectRange(0, all_idxs.size());
        idxs_dset->Write(all_idxs, idxs_fspace, writeSR);
        /////////////////////////////////////////////////////////////////////////////////////////

    } // if (amRoot)
    // RAII => datasets (idxs_dset, decomps_dset), root group and HDF5 file are
    // all closed automatically
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // initialise the dataset indexes for all MPI processes
    std::vector<std::size_t> idx = m_comm->Scatter(root_rank, all_idxs);
    std::size_t ids_i            = idx[IDS_IDX_IDX];
    std::size_t data_i           = idx[DATA_IDX_IDX];

    // set properties for parallel file access (if we're in parallel)
    H5::PListSharedPtr parallelProps = H5::PList::Default();
    H5::PListSharedPtr writePL = H5::PList::Default();
    if (m_comm->GetSize() > 1)
    {
        // use MPI/O to access the file
        parallelProps = H5::PList::FileAccess();
        parallelProps->SetMpio(m_comm);
        // use collective IO
        writePL = H5::PList::DatasetXfer();
        writePL->SetDxMpioCollective();
    } // end of <if (m_comm->GetSize() > 1)>

    // reopen the file
    H5::FileSharedPtr outfile =
        H5::File::Open(outFile, H5F_ACC_RDWR, parallelProps);
    ASSERTL1(outfile, prfx.str() + "cannot open HDF5 file.");
    H5::GroupSharedPtr root = outfile->OpenGroup("NEKTAR");
    ASSERTL1(root, prfx.str() + "cannot open root group.");

    // RAII => field group is closed automatically
    m_comm->Block();
    // all hdf5 groups have now been created
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // open the IDS dataset and associated data space
    H5::DataSetSharedPtr ids_dset = root->OpenDataSet("IDS");
    ASSERTL1(ids_dset, prfx.str() + "cannot open IDS dataset.");
    H5::DataSpaceSharedPtr ids_fspace = ids_dset->GetSpace();
    ASSERTL1(ids_fspace, prfx.str() + "cannot open IDS filespace.");

    // open the DATA dataset and associated data space
    H5::DataSetSharedPtr data_dset = root->OpenDataSet("DATA");
    ASSERTL1(data_dset, prfx.str() + "cannot open DATA dataset.");
    H5::DataSpaceSharedPtr data_fspace = data_dset->GetSpace();
    ASSERTL1(data_fspace, prfx.str() + "cannot open DATA filespace.");

    // write the data
    ////////////////////////////////////////////////////////////////////////////////////////////////
    for (int f = 0; f < nFields; ++f)
    {
        // write the element ids
        std::size_t nFieldElems = fielddefs[f]->m_elementIDs.size();
        ids_fspace->SelectRange(ids_i, nFieldElems);
        ids_dset->Write(fielddefs[f]->m_elementIDs, ids_fspace, writePL);
        ids_i += nFieldElems;

        // write the element values
        std::size_t nElemVals  = fielddata[f].size() / nFieldElems;
        std::size_t nFieldVals = nElemVals * nFieldElems;
        data_fspace->SelectRange(data_i, nFieldVals);
        data_dset->Write(fielddata[f], data_fspace, writePL);
        data_i += nFieldVals;
    }

    for (int f = nFields; f < nMaxFields; ++f)
    {
        // this MPI process is handling fewer than nMaxFields fields
        // so, since this is a collective operation
        // just rewrite the element ids and values of the last field
        ids_dset->Write(
            fielddefs[nFields - 1]->m_elementIDs, ids_fspace, writePL);
        data_dset->Write(fielddata[nFields - 1], data_fspace, writePL);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////

    m_comm->Block();
    // all data has been written

    if (0 == m_comm->GetRank())
    {
        tm1 = m_comm->Wtime();
        cout << prfx.str() << "leaving after " << tm1 - tm0 << " secs..."
             << endl;
    }
} // end of <FieldIOHdf5::v_Write> method
// RAII => filespaces (ids_fspace, data_fspace), datasets (ids_dset, data_dset),
//         root group and HDF5 file are all closed automatically

void FieldIOHdf5::v_Import(const std::string &infilename,
                           std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                           std::vector<std::vector<NekDouble> > &fielddata,
                           FieldMetaDataMap &fieldinfomap,
                           const Array<OneD, int> ElementIDs)
{
    std::stringstream prfx;
    prfx << m_comm->GetRank() << ": FieldIOHdf5::v_ImportFile(): ";
    double tm0 = 0.0, tm1 = 0.0;
    if (0 == m_comm->GetRank())
    {
        cout << prfx.str() << "entering..." << endl;
        tm0 = m_comm->Wtime();
    }

    int nRanks = m_comm->GetSize();

    DataSourceSharedPtr dataSource = H5DataSource::create(infilename);

    // set properties for parallel file access (if we're in parallel)
    ////////////////////////////////////////////////////////////////////////////////////
    H5::PListSharedPtr parallelProps = H5::PList::Default();
    H5::PListSharedPtr readPL = H5::PList::Default();
    if (nRanks > 1)
    {
        // use MPI/O to access the file
        parallelProps = H5::PList::FileAccess();
        parallelProps->SetMpio(m_comm);
        // use collective IO
        readPL = H5::PList::DatasetXfer();
        readPL->SetDxMpioCollective();
    }
    ////////////////////////////////////////////////////////////////////////////////////

    // open the root group of the hdf5 file
    ////////////////////////////////////////////////////////////////////////////////////
    H5DataSourceSharedPtr h5 =
        boost::static_pointer_cast<H5DataSource>(dataSource);
    ASSERTL1(h5, prfx.str() + "cannot open HDF5 file.");
    H5::GroupSharedPtr root = h5->Get()->OpenGroup("NEKTAR");
    ASSERTL1(root, prfx.str() + "cannot open root group.");
    ////////////////////////////////////////////////////////////////////////////////////

    // open the datasets
    ////////////////////////////////////////////////////////////////////////////////////
    H5::DataSetSharedPtr decomps_dset = root->OpenDataSet("DECOMPOSITION");
    ASSERTL1(decomps_dset, prfx.str() + "cannot open DECOMPOSITION dataset.");

    H5::DataSetSharedPtr idxs_dset = root->OpenDataSet("INDEXES");
    ASSERTL1(idxs_dset, prfx.str() + "cannot open INDEXES dataset.");

    H5::DataSetSharedPtr ids_dset = root->OpenDataSet("IDS");
    ASSERTL1(ids_dset, prfx.str() + "cannot open IDS dataset.");

    H5::DataSetSharedPtr data_dset = root->OpenDataSet("DATA");
    ASSERTL1(data_dset, prfx.str() + "cannot open DATA dataset.");
    ////////////////////////////////////////////////////////////////////////////////////

    // open the dataset file spaces
    ////////////////////////////////////////////////////////////////////////////////////
    H5::DataSpaceSharedPtr decomps_fspace = decomps_dset->GetSpace();
    ASSERTL1(decomps_fspace,
             prfx.str() + "cannot open DECOMPOSITION filespace.");

    H5::DataSpaceSharedPtr idxs_fspace = idxs_dset->GetSpace();
    ASSERTL1(idxs_fspace, prfx.str() + "cannot open INDEXES filespace.");

    H5::DataSpaceSharedPtr ids_fspace = ids_dset->GetSpace();
    ASSERTL1(ids_fspace, prfx.str() + "cannot open IDS filespace.");

    H5::DataSpaceSharedPtr data_fspace = data_dset->GetSpace();
    ASSERTL1(data_fspace, prfx.str() + "cannot open DATA filespace.");
    ////////////////////////////////////////////////////////////////////////////////////

    // Read entire IDS data set to get list of global IDs.
    std::vector<std::size_t> ids;

    ids_dset->Read(ids, ids_fspace, readPL);

    boost::unordered_set<size_t> toread;
    if (ElementIDs != NullInt1DArray)
    {
        for (size_t i = 0; i < ElementIDs.num_elements(); ++i)
        {
            toread.insert(ElementIDs[i]);
        }
    }

    std::vector<std::size_t> decomps;
    decomps_dset->Read(decomps, decomps_fspace, readPL);

    size_t nDecomps = decomps.size() / MAX_DCMPS;
    size_t cnt = 0, cnt2 = 0;

    // Mapping from each decomposition to offsets in the data array.
    vector<size_t> decompsToDataOffsets(nDecomps);
    // Mapping from each group's hash to a vector of element IDs. Note this has
    // to be unsigned int, since that's what we use in FieldDefinitions.
    map<size_t, vector<unsigned int> > groupsToElmts;
    // Mapping from each group's hash to each of the decompositions.
    map<size_t, set<size_t> > groupsToDecomps;

    // True if we are pulling element IDs from ElementIDs.
    bool selective = toread.size() > 0;

    // Counter for data offsets
    size_t dataOffset = 0;

    for (size_t i = 0; i < nDecomps; ++i, cnt += MAX_DCMPS)
    {
        size_t nElmt     = decomps[cnt + ELEM_DCMP_IDX];
        size_t dataSize  = decomps[cnt + VAL_DCMP_IDX];
        size_t groupHash = decomps[cnt + HASH_DCMP_IDX];

        vector<size_t> tmp;

        if (selective)
        {
            for (size_t j = 0; j < nElmt; ++j)
            {
                size_t elmtId = ids[cnt2 + j];
                if (toread.find(elmtId) != toread.end())
                {
                    tmp.push_back(elmtId);
                }
            }
        }
        else
        {
            tmp.insert(
                tmp.begin(), ids.begin() + cnt2, ids.begin() + cnt2 + nElmt);
        }

        vector<unsigned int> tmp2(nElmt);
        for (size_t j = 0; j < nElmt; ++j)
        {
            tmp2[j] = ids[cnt2+j];
        }

        cnt2 += nElmt;

        if (tmp.size() > 0)
        {
            groupsToDecomps[groupHash].insert(i);
        }

        groupsToElmts[i] = tmp2;
        decompsToDataOffsets[i] = dataOffset;
        dataOffset += nElmt * dataSize;
    }

    map<size_t, set<size_t> >::iterator gIt;
    for (gIt = groupsToDecomps.begin(); gIt != groupsToDecomps.end(); ++gIt)
    {
        // Select region from dataset for this decomposition.
        set<size_t>::iterator sIt;
        for (sIt = gIt->second.begin(); sIt != gIt->second.end(); ++sIt)
        {
            std::stringstream fieldNameStream;
            fieldNameStream << gIt->first;

            FieldDefinitionsSharedPtr fielddef =
                MemoryManager<FieldDefinitions>::AllocateSharedPtr();
            ImportFieldDef(readPL, root, fieldNameStream.str(), fielddef);

            fielddef->m_elementIDs = groupsToElmts[*sIt];
            fielddefs.push_back(fielddef);

            if (fielddata != NullVectorNekDoubleVector)
            {
                std::vector<NekDouble> decompFieldData;
                ImportFieldData(
                    readPL, data_dset, data_fspace, decompsToDataOffsets[*sIt],
                    decomps, *sIt, fielddef, decompFieldData);
                fielddata.push_back(decompFieldData);
            }
        }
    }

    m_comm->Block();
    // all data has been read

    if (0 == m_comm->GetRank())
    {
        tm1 = m_comm->Wtime();
    }
}
// RAII => filespaces (ids_fspace, data_fspace, idxs_fspace, decomps_fspace),
//         datasets (ids_dset, data_dset, idxs_dset, decomps_dset),
//         root group and HDF5 file are all closed automatically

void FieldIOHdf5::ImportFieldDef(
    H5::PListSharedPtr readPL,
    H5::GroupSharedPtr root,
    std::string group,
    FieldDefinitionsSharedPtr def)
{
    std::stringstream prfx;
    prfx << m_comm->GetRank() << ": FieldIOHdf5::ImportFieldDefsHdf5(): ";

    H5::GroupSharedPtr field = root->OpenGroup(group);
    ASSERTL1(field, prfx.str() + "cannot open field group, " + group + '.');

    H5::Group::AttrIterator attrIt  = field->attr_begin();
    H5::Group::AttrIterator attrEnd = field->attr_end();
    for (; attrIt != attrEnd; ++attrIt)
    {
        const std::string &attrName = *attrIt;
        if (attrName == "FIELDS")
        {
            field->GetAttribute(attrName, def->m_fields);
        }
        else if (attrName == "SHAPE")
        {
            std::string shapeString;
            field->GetAttribute(attrName, shapeString);

            // check to see if homogeneous expansion and if so
            // strip down the shapeString definition
            size_t loc;
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

            ASSERTL0(valid,
                     prfx.str() +
                     std::string(
                         "unable to correctly parse the shape type: ")
                     .append(shapeString)
                     .c_str());
        }
        else if (attrName == "BASIS")
        {
            field->GetAttribute(attrName, def->m_basis);
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
        }
        else if (attrName == "HOMOGENEOUSSIDS")
        {
            field->GetAttribute(attrName, def->m_homogeneousSIDs);
        }
        else if (attrName == "HOMOGENEOUSZIDS")
        {
            field->GetAttribute(attrName, def->m_homogeneousZIDs);
        }
        else if (attrName == "HOMOGENEOUSYIDS")
        {
            field->GetAttribute(attrName, def->m_homogeneousYIDs);
        }
        else if (attrName == "NUMMODESPERDIR")
        {
            std::string numModesPerDir;
            field->GetAttribute(attrName, numModesPerDir);

            if (strstr(numModesPerDir.c_str(), "UNIORDER:"))
            {
                def->m_uniOrder = true;
            }

            bool valid = ParseUtils::GenerateOrderedVector(
                numModesPerDir.c_str() + 9, def->m_numModes);

            ASSERTL0(valid,
                     prfx.str() +
                     "unable to correctly parse the number of modes.");
        }
        else if (attrName == "POINTSTYPE")
        {
            std::string pointsString;
            field->GetAttribute(attrName, pointsString);
            def->m_pointsDef = true;

            std::vector<std::string> pointsStrings;
            bool valid = ParseUtils::GenerateOrderedStringVector(
                pointsString.c_str(), pointsStrings);
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
            def->m_numPointsDef = true;

            bool valid = ParseUtils::GenerateOrderedVector(
                numPointsPerDir.c_str(), def->m_numPoints);
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
}

void FieldIOHdf5::ImportFieldData(
    H5::PListSharedPtr readPL,
    H5::DataSetSharedPtr data_dset,
    H5::DataSpaceSharedPtr data_fspace,
    size_t data_i,
    std::vector<std::size_t> &decomps,
    size_t decomp,
    const FieldDefinitionsSharedPtr fielddef,
    std::vector<NekDouble> &fielddata)
{
    std::stringstream prfx;
    prfx << m_comm->GetRank() << ": FieldIOHdf5::ImportFieldData(): ";

    size_t nFieldElems = decomps[decomp * MAX_DCMPS + ELEM_DCMP_IDX];
    size_t nElemVals   = decomps[decomp * MAX_DCMPS + VAL_DCMP_IDX];
    size_t nFieldVals  = nFieldElems * nElemVals;

    data_fspace->SelectRange(data_i, nFieldVals);
    data_dset->Read(fielddata, data_fspace, readPL);

    int datasize = CheckFieldDefinition(fielddef);
    ASSERTL0(
        fielddata.size() == datasize * fielddef->m_fields.size(),
        prfx.str() +
        "input data is not the same length as header information.");
}

DataSourceSharedPtr FieldIOHdf5::v_ImportFieldMetaData(
    std::string filename, FieldMetaDataMap &fieldmetadatamap)
{
    DataSourceSharedPtr ans = H5DataSource::create(filename);
    v_ImportFieldMetaData(ans, fieldmetadatamap);
    return ans;
}

void FieldIOHdf5::v_ImportFieldMetaData(DataSourceSharedPtr dataSource,
                                        FieldMetaDataMap &fieldmetadatamap)
{
    H5DataSourceSharedPtr hdf =
        boost::static_pointer_cast<H5DataSource>(dataSource);

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
