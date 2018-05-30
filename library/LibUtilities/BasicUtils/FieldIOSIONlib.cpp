////////////////////////////////////////////////////////////////////////////////
//
//  File: FieldIOSIONlib.cpp
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
//  Description: I/O routines relating to Fields into SIONlib files.
//
////////////////////////////////////////////////////////////////////////////////

#include "sion.h"

#include <LibUtilities/BasicUtils/SIONlib.h>
#include <LibUtilities/BasicUtils/FieldIOSIONlib.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <boost/unordered_set.hpp>
#include <LibUtilities/Communication/CommMpi.h>


namespace berrc = boost::system::errc;

namespace Nektar
{
namespace LibUtilities
{ 

std::string FieldIOSIONlib::className =
    GetFieldIOFactory().RegisterCreatorFunction(
        "SIONlib", FieldIOSIONlib::create, "SIONlib-based output of field data.");

  
const unsigned int FieldIOSIONlib::FORMAT_VERSION = 1;
  
const std::string FieldIOSIONlib::ATTRNAME_FIELDS = "FIELDS";
const std::string FieldIOSIONlib::ATTRNAME_BASIS = "BASIS";
const std::string FieldIOSIONlib::ATTRNAME_SHAPE = "SHAPE";
const std::string FieldIOSIONlib::ATTRNAME_HOMOLENS = "HOMOGENEOUSLENGTHS";
const std::string FieldIOSIONlib::ATTRNAME_NUMMODES = "NUMMODESPERDIR";
const std::string FieldIOSIONlib::ATTRNAME_POINTSTYPE = "POINTSTYPE";
const std::string FieldIOSIONlib::ATTRNAME_NUMPOINTS = "NUMPOINTSPERDIR";
  
const std::string FieldIOSIONlib::ATTRVALUE_MIXORDER = "MIXORDER";
const std::string FieldIOSIONlib::ATTRVALUE_UNIORDER = "UNIORDER";

sion_int32 FieldIOSIONlib::block_size = -1;
sion_int64 FieldIOSIONlib::chunk_size = -1;
std::string FieldIOSIONlib::read_mode = "";
std::string FieldIOSIONlib::write_mode = "";

  
FieldIOSIONlib::FieldIOSIONlib(LibUtilities::CommSharedPtr pComm,
                         bool sharedFilesystem)
    : FieldIO(pComm, sharedFilesystem)
{
}


void FieldIOSIONlib::v_Init(const LibUtilities::SessionReaderSharedPtr session)
{
    block_size = -1;
    if (session->DefinesSolverInfo("IOBlockSize"))
    {
        std::string int_str = session->GetSolverInfo("IOBlockSize");
        std::istringstream(int_str) >> block_size;
    }

    NekDouble blocks_per_chunk = 1.0;
    if (session->DefinesSolverInfo("IOBlocksPerChunk"))
    {
        std::string int_str = session->GetSolverInfo("IOBlocksPerChunk");
        std::istringstream(int_str) >> blocks_per_chunk;
    }
    chunk_size = (sion_int64) (block_size*blocks_per_chunk);

    read_mode = "br";
    if (session->DefinesSolverInfo("IOReadMode"))
    {
        read_mode = session->GetSolverInfo("IOReadMode");
    }

    write_mode = "bw";
    if (session->DefinesSolverInfo("IOWriteMode"))
    {
        write_mode = session->GetSolverInfo("IOWriteMode");
    }
}


void FieldIOSIONlib::v_InitFromBenchmarker(const LibUtilities::IOSettingsSharedPtr iosettings)
{
    LibUtilities::IOSettings::iterator it;
    
    block_size = -1;
    it = iosettings->find("IOBlockSize");
    if (iosettings->end() != it)
    {
        std::istringstream(it->second) >> block_size;
    }

    NekDouble blocks_per_chunk = 1.0;
    it = iosettings->find("IOBlocksPerChunk");
    if (iosettings->end() != it)
    {
        std::istringstream(it->second) >> blocks_per_chunk;
    }
    chunk_size = (sion_int64) (block_size*blocks_per_chunk);

    read_mode = "br";
    it = iosettings->find("IOReadMode");
    if (iosettings->end() != it)
    {
        read_mode = it->second;
    }

    write_mode = "bw";
    it = iosettings->find("IOWriteMode");
    if (iosettings->end() != it)
    {
        write_mode = it->second;
    }
}

  
/**
 * Open a SIONlib file for writing by constructing and returning
 * a pointer to a SIONlib object, and also outputting the format version.
 *
 * @param outFilename  name of output file.
 * @return Pointer to SIONlib object.
 */
SIONlib::SIONFile *FieldIOSIONlib::OpenFileForWriting(const std::string &outFilename)
{
    CommMpi *commMpiPtr =  (CommMpi*) m_comm.get();
    SIONlib::SIONFile *fp = new SIONlib::SIONFile(outFilename, write_mode,
        1, chunk_size, block_size,
        m_comm->GetRank(), commMpiPtr->GetComm(), commMpiPtr->GetComm());

    ASSERTL0(NULL != fp,
             "Unable to construct SIONFile object for " + outFilename);

    SIONlib::SIONFile &f = *fp;
    f.open();
    
    ASSERTL0(-1 != f.getReturnCode(),
             "Unable to open SIONlib file " + outFilename);

    f << FORMAT_VERSION;
    
    return fp;
}


  
void FieldIOSIONlib::WriteRawDataToBinaryVector(std::vector<NekByte> &v, const NekByte* dptr, const size_t dsize)
{
    std::vector<NekByte> v2(dsize);
    std::copy(dptr, dptr+dsize, v2.begin());

    v.reserve(v.size() + v2.size());
    v.insert(v.end(), v2.begin(), v2.end());
}

void FieldIOSIONlib::WriteStringToBinaryVector(std::vector<NekByte> &v, const std::string &s)
{
    size_t ssize = s.size();
    WriteRawDataToBinaryVector(v, (const NekByte*) &ssize, sizeof(size_t));
    WriteRawDataToBinaryVector(v, (const NekByte*) s.c_str(), ssize);  
}

template<class T>   
void FieldIOSIONlib::WriteDatumToBinaryVector(std::vector<NekByte> &v, const T d)
{
    WriteRawDataToBinaryVector(v, (const NekByte*) &d, sizeof(T));
}

template<class T>   
void FieldIOSIONlib::WriteDataToBinaryVector(std::vector<NekByte> &v, const std::vector<T> &d)
{
    size_t dsize = d.size();
    WriteRawDataToBinaryVector(v, (const NekByte*) &dsize, sizeof(size_t));
    for (auto it = d.begin(); it != d.end(); ++it)
    {
        WriteRawDataToBinaryVector(v, (const NekByte*) &(*it), sizeof(T));
    }
}

  
unsigned long FieldIOSIONlib::v_Write(const std::string &outFilePrefix,
    std::vector<FieldDefinitionsSharedPtr> &fielddefs,
    std::vector<std::vector<NekDouble> > &fielddata,
    const FieldMetaDataMap &fieldmetadatamap,
    const bool backup)
{
    // We make a number of assumptions in this code:
    //   1. All element ids have the same type: unsigned int
    //   2. All elements within a given field have the same number of values
    //   3. All element values have the same type, NekDouble

    ASSERTL1(fielddefs.size() == fielddata.size(),
             "fielddefs and fielddata have incompatible lengths.");

    size_t nFields = fielddefs.size();
    size_t nMaxFields = nFields;

    int homDim = -1;
    int varOrder = 0;

    for (size_t i = 0; i < nFields; ++i)
    {
        if (!fielddefs[i]->m_uniOrder)
        {
            varOrder = 1;
            break;
        }
    }

    m_comm->AllReduce(varOrder, LibUtilities::ReduceMax);
    m_comm->AllReduce(nMaxFields, LibUtilities::ReduceMax);

    SetUpOutput(outFilePrefix, false, backup);

    // Each MPI process iterates through its fields and outputs field
    // definitions and field data to the SIONlib file.
    SIONlib::SIONFile* fp = OpenFileForWriting(outFilePrefix);
    ASSERTL0(NULL != fp, "Cannot open SIONlib file.");
    SIONlib::SIONFile &f = *fp; 

    if (f.isModeCollectiveMerge())
    {
        f << (size_t) m_comm->GetRank();
    }
    f << nFields;

    std::vector<NekByte> attrVector;
    
    for (size_t i = 0; i < nFields; ++i)
    {
        ASSERTL1(fielddata[i].size() > 0,
            "fielddata vector must contain at least one value.");
        ASSERTL1(fielddata[i].size() ==
            fielddefs[i]->m_fields.size()*CheckFieldDefinition(fielddefs[i]),
            "fielddata vector has invalid size.");

        FieldDefinitionsSharedPtr def = fielddefs[i];

        size_t nAttrs = (def->m_numHomogeneousDir ? 5 : 4);
        if (def->m_pointsDef) nAttrs++;
        if (def->m_numPointsDef) nAttrs++;
        
        WriteDatumToBinaryVector(attrVector, nAttrs);
	
        // FIELDS
        ////////////////////////////////////////////////////////////////
        WriteStringToBinaryVector(attrVector, ATTRNAME_FIELDS);
        
        size_t nAttrFields = def->m_fields.size();
        WriteDatumToBinaryVector(attrVector, nAttrFields);

        if (nAttrFields > 0)
        {
            for (size_t j = 0; j < nAttrFields; ++j)
            { 
                WriteStringToBinaryVector(attrVector, def->m_fields[j]);
            }
        }
        ////////////////////////////////////////////////////////////////
                
        // BASIS
        ////////////////////////////////////////////////////////////////
        WriteStringToBinaryVector(attrVector, ATTRNAME_BASIS);

	WriteDataToBinaryVector(attrVector, def->m_basis);
        ////////////////////////////////////////////////////////////////
        
        // SHAPE
        ////////////////////////////////////////////////////////////////
        WriteStringToBinaryVector(attrVector, ATTRNAME_SHAPE);

        std::stringstream shapeStringStream;
        shapeStringStream << ShapeTypeMap[def->m_shapeType];

        if (def->m_numHomogeneousDir > 0)
        {
            if (homDim == -1)
            {
                homDim = def->m_numHomogeneousDir;
            }

            ASSERTL1(homDim == def->m_numHomogeneousDir,
                "HDF5 does not support variable homogeneous directions in "
                "the same file.");

            shapeStringStream << "-HomogenousExp"
                << def->m_numHomogeneousDir << "D";
        }

        if (def->m_homoStrips)
        {
            shapeStringStream << "-Strips";
        }
        
        WriteStringToBinaryVector(attrVector, shapeStringStream.str());
        ////////////////////////////////////////////////////////////////
        
        // Determine HOMOGENEOUS attributes
        if (def->m_numHomogeneousDir)
        {
            // HOMOGENEOUSLENGTHS
            ////////////////////////////////////////////////////////////////
	    WriteStringToBinaryVector(attrVector, ATTRNAME_HOMOLENS);
	    WriteDataToBinaryVector(attrVector, def->m_homogeneousLengths);
            ////////////////////////////////////////////////////////////////

            // homo IDs are streamed from separate arrays
            ////////////////////////////////////////////////////////////////
	    WriteDataToBinaryVector(attrVector, def->m_homogeneousYIDs);
	    WriteDataToBinaryVector(attrVector, def->m_homogeneousZIDs);
	    WriteDataToBinaryVector(attrVector, def->m_homogeneousSIDs);
            ////////////////////////////////////////////////////////////////
        }
        
        // NUMMODESPERDIR
        ////////////////////////////////////////////////////////////////
        WriteStringToBinaryVector(attrVector, ATTRNAME_NUMMODES);

        std::vector<unsigned int>& numModes = def->m_numModes;
        size_t nNumModes = def->m_basis.size();
        ASSERTL1(nNumModes == def->m_numModes.size(),
             "basis and numModes have incompatible lengths.");
        
        if (def->m_uniOrder && !varOrder)
        {
	    WriteStringToBinaryVector(attrVector, ATTRVALUE_UNIORDER);
        }
        else
        {
	    WriteStringToBinaryVector(attrVector, ATTRVALUE_MIXORDER);
        }

        if (nNumModes > 0)
        {
	    WriteDataToBinaryVector(attrVector, def->m_numModes);
        }
        ////////////////////////////////////////////////////////////////

        // POINTSTYPE
        ////////////////////////////////////////////////////////////////
        if (def->m_pointsDef)
        {
	    WriteStringToBinaryVector(attrVector, ATTRNAME_POINTSTYPE);

            stringstream pointsTypeStringStream;
            std::vector<LibUtilities::PointsType>& points = def->m_points;
            size_t nPoints = points.size();
            
            for (size_t j = 0; j < nPoints; ++j)
            {
                if (j > 0)
                {
                    pointsTypeStringStream << ",";
                }
                pointsTypeStringStream << kPointsTypeStr[points[j]];
            }

            WriteStringToBinaryVector(attrVector, pointsTypeStringStream.str());
        }
        ////////////////////////////////////////////////////////////////

        // NUMPOINTSPERDIR
        ////////////////////////////////////////////////////////////////
        if (def->m_numPointsDef)
        {
	    WriteStringToBinaryVector(attrVector, ATTRNAME_NUMPOINTS);

            stringstream numPointsStringStream;
            std::vector<unsigned int>& numPoints = def->m_numPoints;
            size_t nNumPoints = numPoints.size();
            
            for (size_t j = 0; j < nNumPoints; ++j)
            {
                if (j > 0)
                {
                    numPointsStringStream << ",";
                }
                numPointsStringStream << numPoints[j];
            }

            WriteStringToBinaryVector(attrVector, numPointsStringStream.str());
        }
        ////////////////////////////////////////////////////////////////

	// field attributes
	////////////////////////////////////////////////////////////////
        f << attrVector;
	attrVector.clear();
	////////////////////////////////////////////////////////////////

	// elementIDs
        ////////////////////////////////////////////////////////////////
        f << def->m_elementIDs;
        ////////////////////////////////////////////////////////////////

        // data
        ////////////////////////////////////////////////////////////////
        f << fielddata[i];
        ////////////////////////////////////////////////////////////////
    }
    
    for (size_t i = nFields; i < nMaxFields; ++i)
    {
        const std::vector<NekByte> dumAttributes = {0};
	const std::vector<unsigned int> dumElementIds = {0};
	const std::vector<NekDouble> dumData = {0.0};
        f << dumAttributes;
	f << dumElementIds;
	f << dumData;
    }
    
    unsigned long nWritten = f.getBytesWritten();
    f.close();
    
    m_comm->Block();
    return nWritten;
}
  

/**
 * Open a SIONlib file for reading by constructing and returning
 * a pointer to a SIONlib object, and also read in the format version,
 * throwing an error if it is greater than FieldIOSIONlib::FORMAT_VERSION.
 *
 * @param inFilename  name of input file.
 * @return Pointer to SIONlib object.
 */
SIONlib::SIONFile *FieldIOSIONlib::OpenFileForReading(const std::string &inFilename)
{
    CommMpi *commMpiPtr =  (CommMpi*) m_comm.get();
    SIONlib::SIONFile *fp = new SIONlib::SIONFile(inFilename, read_mode,
        1, chunk_size, block_size,
        m_comm->GetRank(), commMpiPtr->GetComm(), commMpiPtr->GetComm());

    ASSERTL0(NULL != fp,
             "Unable to construct SIONFile object for " + inFilename);

    SIONlib::SIONFile &f = *fp;
    f.open();
    
    ASSERTL0(-1 != f.getReturnCode(),
             "Unable to open SIONlib file " + inFilename);

    unsigned int formatVersion = 0;
    f >> formatVersion;
    
    std::stringstream version;
    version << formatVersion;

    ASSERTL0(formatVersion <= FORMAT_VERSION,
             "File " + inFilename + " is at version " + version.str() +
             ", which is higher than supported in this version of Nektar++.");
    
    return fp;
}


unsigned long FieldIOSIONlib::v_Import(const std::string &infilename,
    std::vector<FieldDefinitionsSharedPtr> &fielddefs,
    std::vector<std::vector<NekDouble> > &fielddata,
    FieldMetaDataMap &fieldinfomap,
    const Array<OneD, int> &ElementIDs)  
{
    SIONlib::SIONFile* fp = OpenFileForReading(infilename);
    ASSERTL0(NULL != fp, "Cannot open SIONlib file.");
    SIONlib::SIONFile &f = *fp;
    
    if (f.isModeCollectiveMerge())
    {
        unsigned long bytesWritten = f.getBytesWritten();
        unsigned long bytesRead = 0;
    
        if (bytesWritten > 0)
        {
            // this process is one of the collectors that originally wrote the checkpoint file
            do
	    {
	        size_t rk = 0;
	        f >> rk;
	        bytesRead += sizeof(rk);
	    
	        if (rk != m_comm->GetRank())
	        {
	            // the data stored hereafter originated from process rk
		  
		    // TODO: send the number of fields to rk
	            // TODO: for each field do the following,
	            //         pack the data from rk into three buffers
	            //             (i) field attributes, vector<NekByte>
	            //            (ii) element ids, vector<unsigned int>
	            //           (iii) element data, vector<double>
	            //         send each buffer to process rk

	            // TODO: update bytesRead
	        }
	        else
	        {
		    bytesRead += DirectImport(f, fielddefs,fielddata,fieldinfomap,ElementIDs);
	        }
	    }
	    while (bytesWritten > bytesRead);
	}
	else
	{
	    // this process is not a collector, i.e., its data was written to the area of the
	    // checkpoint file reserved for a collector
	  
	    // TODO: wait to receive the field count from an unknown collector
	    // TODO: for each field do the following,
	    //           (i) receive the field attributes
	    //          (ii) receive the element ids
	    //         (iii) receive the element data
	}
    }
    else
    {
        return DirectImport(f, fielddefs,fielddata,fieldinfomap,ElementIDs);
    }
}


unsigned long FieldIOSIONlib::DirectImport(SIONlib::SIONFile &f,
    std::vector<FieldDefinitionsSharedPtr> &fielddefs,
    std::vector<std::vector<NekDouble> > &fielddata,
    FieldMetaDataMap &fieldinfomap,
    const Array<OneD, int> &ElementIDs)  
{
    size_t nFields = 0;
    f >> nFields;
    
    for (size_t i = 0; i < nFields; ++i)
    {
        FieldDefinitionsSharedPtr def =
            MemoryManager<FieldDefinitions>::AllocateSharedPtr();
        
        std::string attrName, attrValue;
        size_t nAttrs, nAttrFields;

        f >> nAttrs;
        
        for (size_t j = 0; j < nAttrs; ++j)
        {         
            f >> attrName;

            if (attrName == ATTRNAME_FIELDS)
            {
                f >> nAttrFields;
                
                def->m_fields.resize(nAttrFields);

                for (size_t k = 0; k < nAttrFields; ++k)
                {
                    f >> def->m_fields[k];
                }
            }
            else if (attrName == ATTRNAME_BASIS)
            {
                f >> def->m_basis;

                // check the basis is in range
                std::vector<BasisType>::const_iterator bIt  = def->m_basis.begin();
                std::vector<BasisType>::const_iterator bEnd = def->m_basis.end();
                for (; bIt != bEnd; ++bIt)
                {
                    BasisType bt = *bIt;
                    ASSERTL0(bt >= 0 && bt < SIZE_BasisType,
                        "Unable to correctly parse the basis types.");
                }
            }
            else if (attrName == ATTRNAME_SHAPE)
            {
                f >> attrValue;        

                // check to see if homogeneous expansion and if so
                // strip down the shapeString definition
                size_t loc;
                //---> this finds the first location of 'n'!
                if (attrValue.find("Strips") != string::npos)
                {
                    def->m_homoStrips = true;
                }

                if ((loc = attrValue.find_first_of("-")) != string::npos)
                {
                    if (attrValue.find("Exp1D") != string::npos)
                    {
                        def->m_numHomogeneousDir = 1;
                    }
                    else // HomogeneousExp1D
                    {
                        def->m_numHomogeneousDir = 2;
                    }

                    attrValue.erase(loc, attrValue.length());
                }

                // get the geometrical shape
                bool valid = false;
                for (unsigned int k = 0; k < SIZE_ShapeType; k++)
                {
                    if (ShapeTypeMap[k] == attrValue)
                    {
                        def->m_shapeType = (ShapeType) k;
                        valid = true;
                        break;
                    }
                }

                ASSERTL0(valid,
                    std::string("Unable to correctly parse the shape type: ")
                        .append(attrValue)
                        .c_str());
            }
            else if (attrName == ATTRNAME_HOMOLENS)
            {
                def->m_numHomogeneousDir = true;
          
                f >> def->m_homogeneousLengths;
                
                f >> def->m_homogeneousYIDs;
                f >> def->m_homogeneousZIDs;
                f >> def->m_homogeneousSIDs;
            }
            else if (attrName == ATTRNAME_NUMMODES)
            {
                f >> attrValue;

                if (attrValue == ATTRVALUE_UNIORDER)
                {
                    def->m_uniOrder = true;
                    f >> def->m_numModes;
                }
                else if (attrValue == ATTRVALUE_MIXORDER)
                {
                    f >> def->m_numModes;
                }
                else
                {
                    std::string errstr("Unknown " + attrName + " value: ");
                    errstr += attrValue;
                    ASSERTL1(false, errstr.c_str());
                }
            }
            else if (attrName == ATTRNAME_POINTSTYPE)
            {
                def->m_pointsDef = true;
            
                f >> attrValue;
            
                std::vector<std::string> pointsStrings;
                bool valid = ParseUtils::GenerateVector(
                    attrValue.c_str(), pointsStrings);

                ASSERTL0(valid,
                    "Unable to correctly parse the points types.");
                
                for (std::vector<std::string>::size_type k = 0;
                     k < pointsStrings.size();
                     ++k)
                {
                    valid = false;
                    for (unsigned int l = 0; l < SIZE_PointsType; ++l)
                    {
                        if (kPointsTypeStr[l] == pointsStrings[k])
                        {
                            def->m_points.push_back((PointsType)l);
                            valid = true;
                            break;
                        }
                    }

                    ASSERTL0(valid,
                        std::string(
                            "Unable to correctly parse the points type: ")
                            .append(pointsStrings[k])
                            .c_str());
                }
            }
            else if (attrName == ATTRNAME_NUMPOINTS)
            {
                def->m_numPointsDef = true;
            
                f >> attrValue;

                bool valid = ParseUtils::GenerateVector(
                    attrValue.c_str(), def->m_numPoints);
                ASSERTL0(valid,
                    "Unable to correctly parse the number of points.");
            }
            else
            {
                std::string errstr("Unknown field attribute: ");
                errstr += attrName;
                ASSERTL1(false, errstr.c_str());
            }
        }
        
        // element IDs
        ////////////////////////////////////////////////
        f >> def->m_elementIDs;
        ////////////////////////////////////////////////
        
        // element data
        ////////////////////////////////////////////////
        ASSERTL0(fielddata != NullVectorNekDoubleVector,
            "Null fielddata.");
        
        std::vector<NekDouble> data;
        f >> data;
        
        int datasize = CheckFieldDefinition(def);
        ASSERTL0(data.size() == datasize*def->m_fields.size(),
            "Input data is not the same length as header information.");
        ////////////////////////////////////////////////

        fielddefs.push_back(def);
        fielddata.push_back(data);
    }

    unsigned long nRead = f.getBytesRead();
    f.close();

    m_comm->Block();
    
    return nRead;
}

  
DataSourceSharedPtr FieldIOSIONlib::v_ImportFieldMetaData(
    const std::string &filename,
    FieldMetaDataMap &fieldmetadatamap)
{
    DataSourceSharedPtr ans = SIONlibDataSource::create(filename);
    
    return ans;
}


}
}
