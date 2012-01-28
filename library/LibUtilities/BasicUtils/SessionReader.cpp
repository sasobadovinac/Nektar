///////////////////////////////////////////////////////////////////////////////
//
// File SessionReader.cpp
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
// Description: Session reader
//
///////////////////////////////////////////////////////////////////////////////

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include <iostream>
#include <string>
using namespace std;

#include <boost/algorithm/string.hpp>
//#include <boost/algorithm/string/regex.hpp>
#include <tinyxml/tinyxml.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/MeshPartition.h>

namespace Nektar
{
    namespace LibUtilities
    {
        /**
         * @class SessionReader
         *
         * This class provides an interface to Nektar++-specific content in a
         * supplied XML document. It also initialises a Nektar++ session
         * including setting up communication for parallel execution and where
         * necessary partitioning the supplied mesh for running across multiple
         * processes.
         *
         * A session should be initialised at the beginning of a user's
         * application by passing the command-line arguments. This not only
         * allows the SessionReader to extract the name of the XML document to
         * load containing Nektar++ session information, but also supplies the
         * MPI arguments necessary for setting up parallel communication. The
         * SessionReader should be initialised using the #CreateInstance
         * function:
         * @code
         * LibUtilities::SessionReaderSharedPtr vSession
         *          = LibUtilities::SessionReader::CreateInstance(argc, argv);
         * @endcode
         * The instance \c vSession can now be passed to other key Nektar++
         * components during their construction.
         * @note At the end of the user application, it is important to call the
         * #Finalise routine in order to finalise any MPI communication and
         * correctly free resources.
         *
         * The SessionReader class provides streamlined, validated access to
         * session parameters, solver information and functions defined within a
         * Nektar++ XML document. The available routines and their usage is
         * documented below.
         *
         * In the case of solver information properties, the classes to which
         * these parameters are pertinent may register with the SessionReader
         * class the set of valid values for a given property. Such values may
         * also be associated with an enumeration value for more transparent use
         * of the property values in code.
         */

        /**
         * This map of maps stores the list of valid string values for a number
         * of solver information parameters. The top level map connects
         * different parameter names to their list of possible values. The list
         * of possible values is also a map, mapping a valid string to a
         * corresponding enum value.
         *
         * This list is populated through the #RegisterEnumValue static member
         * function which is called statically from various classes to register
         * the valid values for solver info parameters associated with them. The
         * map is therefore fully populated before the SessionReader class is
         * instantiated and a file is read in and parsed.
         */
        EnumMapList SessionReader::m_enums;


        /**
         * List of default values for solver information parameters to be used
         * in the case of them not being provided.
         *
         * This list is populated through the #RegisterDefaultSolverInfo static
         * member variable which is called statically from various classes to
         * register the default value for a given parameter.
         */
        SolverInfoMap SessionReader::m_solverInfoDefaults;


        /**
         * This constructor parses the command-line arguments given to the user
         * application to set up any MPI communication, read supplied XML
         * session files, and partition meshes where necessary.
         *
         * @param   argc        Number of command-line arguments
         * @param   argv        Array of command-line arguments
         */
        SessionReader::SessionReader(int argc, char *argv[])
        {
            ASSERTL0(argc > 1, "No filename argument specified.");

            std::vector<std::string> vFilenames;
            for (unsigned int i = 1; i < argc; ++i)
            {
                std::string vFile = argv[i];

                if (vFile.substr(vFile.find_last_of('.')) == ".xml")
                {
                    vFilenames.push_back(vFile);
                }
            }

            m_filename = vFilenames[0];
            m_sessionName = m_filename.substr(0, m_filename.find_last_of('.'));
            m_xmlDoc = MergeDoc(vFilenames);

            // Create communicator
            CreateComm(argc, argv, m_filename);
        }


        SessionReader::SessionReader(int argc, char *argv[], const std::vector<std::string> &pFilenames, const CommSharedPtr &pComm)
        {
            ASSERTL0(pFilenames.size() > 0, "No filenames specified.");

            m_filename = pFilenames[0];
            m_sessionName = m_filename.substr(0, m_filename.find_last_of('.'));
            m_xmlDoc = MergeDoc(pFilenames);

            // Create communicator
            if (!pComm.get())
            {
                CreateComm(argc, argv, m_filename);
            }
            else
            {
                m_comm = pComm;
            }
        }


        /**
         *
         */
        SessionReader::~SessionReader()
        {

        }


        /**
         * Performs the main initialisation of the object. The XML file provided
         * on the command-line is loaded and any mesh partitioning is done. The
         * resulting process-specific XML file (containing the process's
         * geometry partition) is then reloaded and parsed.
         */
        void SessionReader::InitSession()
        {
            // Partition mesh
            PartitionMesh();

            // Parse the XML data in #m_xmlDoc
            ParseDocument();
        }


        /**
         *
         */
        TiXmlDocument& SessionReader::GetDocument()
        {
            ASSERTL1(m_xmlDoc, "XML Document not defined.");
            return *m_xmlDoc;
        }


        /**
         * The single parameter specifies a path to the requested element in a
         * similar format to the filesystem path. Given the following XML:
         * @code
         * <NEKTAR>
         *   <CONDITIONS>
         *     <PARAMETERS>
         *     ...
         *     </PARAMETERS>
         *   </CONDITIONS>
         * </NEKTAR>
         * @endcode
         * the PARAMETERS element would be retrieved by requesting the path:
         * @code
         * Nektar/Conditions/Parameters
         * @endcode
         * @note Paths are case-insensitive.
         *
         * @param   pPath       Path to requested element.
         * @returns Direct pointer to requested XML Element.
         */
        TiXmlElement* SessionReader::GetElement(const string& pPath)
        {
            std::string vPath = boost::to_upper_copy(pPath);
            std::vector<std::string> strs;
            boost::split(strs, vPath, boost::is_any_of("\\/ "));
            ASSERTL0(strs.size() > 0, "No path given in XML element request.");

            TiXmlElement* vReturn = m_xmlDoc->FirstChildElement(strs[0].c_str());
            ASSERTL0(vReturn, std::string("Cannot find element '")
                              + strs[0] + std::string("'."));
            for (int i = 1; i < strs.size(); ++i)
            {
                vReturn = vReturn->FirstChildElement(strs[i].c_str());
                ASSERTL0(vReturn, std::string("Cannot find element '")
                                  + strs[i] + std::string("'."));
            }
            return vReturn;
        }


        /**
         *
         */
        bool SessionReader::DefinesElement(const std::string &pPath) const
        {
            std::string vPath = boost::to_upper_copy(pPath);
            std::vector<std::string> strs;
            boost::split(strs, vPath, boost::is_any_of("\\/ "));
            ASSERTL0(strs.size() > 0, "No path given in XML element request.");

            TiXmlElement* vReturn = m_xmlDoc->FirstChildElement(strs[0].c_str());
            ASSERTL0(vReturn, std::string("Cannot find element '")
                              + strs[0] + std::string("'."));
            for (int i = 1; i < strs.size(); ++i)
            {
                vReturn = vReturn->FirstChildElement(strs[i].c_str());
                if (!vReturn) return false;
            }
            return true;
        }


        /**
         *
         */
        const std::string& SessionReader::GetFilename() const
        {
            return m_filename;
        }


        /**
         *
         */
        const std::string& SessionReader::GetSessionName() const
        {
            return m_sessionName;
        }


        /**
         *
         */
        CommSharedPtr& SessionReader::GetComm()
        {
            return m_comm;
        }


        /**
         * This routine finalises any parallel communication.
         *
         * @note This routine should be called at the very end of a users
         * application.
         */
        void SessionReader::Finalise()
        {
            m_comm->Finalise();
        }


        /**
         *
         */
        bool SessionReader::DefinesParameter(const std::string& pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            ParameterMap::const_iterator paramMapIter = m_parameters.find(vName);
            return (paramMapIter != m_parameters.end());
        }


        /**
         * If the parameter is not defined, termination occurs. Therefore, the
         * parameters existence should be tested for using #DefinesParameter
         * before calling this function.
         *
         * @param   pName       The name of a floating-point parameter.
         * @returns The value of the floating-point parameter.
         */
        const NekDouble& SessionReader::GetParameter(const std::string& pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            ParameterMap::const_iterator paramMapIter = m_parameters.find(vName);

            ASSERTL0(paramMapIter != m_parameters.end(),
                (std::string("Unable to find requested parameter: ") + pName).c_str());

            return paramMapIter->second;
        }


        /**
         *
         */
        void SessionReader::LoadParameter(const std::string &pName, int &pVar) const
        {
            std::string vName = boost::to_upper_copy(pName);
            ParameterMap::const_iterator paramMapIter = m_parameters.find(vName);
            ASSERTL0(paramMapIter != m_parameters.end(),
                    "Required parameter '" + pName + "' not specified in session.");
            pVar = (int)floor(paramMapIter->second);
        }


        /**
         *
         */
        void SessionReader::LoadParameter(const std::string &pName, int &pVar, const int &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            ParameterMap::const_iterator paramMapIter = m_parameters.find(vName);
            if(paramMapIter != m_parameters.end())
            {
                pVar = (int)floor(paramMapIter->second);
            }
            else
            {
                pVar  = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::LoadParameter(const std::string &pName, NekDouble& pVar) const
        {
            std::string vName = boost::to_upper_copy(pName);
            ParameterMap::const_iterator paramMapIter = m_parameters.find(vName);
            ASSERTL0(paramMapIter != m_parameters.end(),
                    "Required parameter '" + pName + "' not specified in session.");
            pVar = paramMapIter->second;
        }


        /**
         *
         */
        void SessionReader::LoadParameter(const std::string &pName, NekDouble &pVar, const NekDouble &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            ParameterMap::const_iterator paramMapIter = m_parameters.find(vName);
            if(paramMapIter != m_parameters.end())
            {
                pVar = paramMapIter->second;
            }
            else
            {
                pVar  = pDefault;
            }
        }



        /**
         *
         */
        void SessionReader::SetParameter(const std::string &pName, int &pVar) 
        {
            std::string vName = boost::to_upper_copy(pName);
            m_parameters[vName] = pVar;
        }


        /**
         *
         */
        void SessionReader::SetParameter(const std::string &pName, NekDouble& pVar) 
        {
            std::string vName = boost::to_upper_copy(pName);
            m_parameters[vName] = pVar;
        }



        /**
         *
         */
        bool SessionReader::DefinesSolverInfo(const std::string &pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            SolverInfoMap::const_iterator solverInfoMapIter = m_solverInfo.find(vName);
            return (solverInfoMapIter != m_solverInfo.end());
        }


        /**
         *
         */
        const std::string& SessionReader::GetSolverInfo(const std::string &pProperty) const
        {
            std::string vProperty = boost::to_upper_copy(pProperty);
            SolverInfoMap::const_iterator slvIter = m_solverInfo.find(vProperty);

            ASSERTL1(slvIter != m_solverInfo.end(),
                (std::string("Unable to find requested property: ") + pProperty).c_str());

            return slvIter->second;
        }


        /**
         *
         */
        void SessionReader::LoadSolverInfo(const std::string &pName, std::string &pVar, const std::string &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            SolverInfoMap::const_iterator solverInfoMapIter = m_solverInfo.find(vName);
            if(solverInfoMapIter != m_solverInfo.end())
            {
                pVar = solverInfoMapIter->second;
            }
            else
            {
                pVar  = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::MatchSolverInfo(const std::string &pName, const std::string &pTrueVal, bool &pVar, const bool &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            SolverInfoMap::const_iterator solverInfoMapIter = m_solverInfo.find(vName);
            if(solverInfoMapIter != m_solverInfo.end())
            {
                pVar = (NoCaseStringCompare(solverInfoMapIter->second, pTrueVal) == 0);
            }
            else
            {
                pVar  = pDefault;
            }
        }


        /**
         *
         */
        bool SessionReader::MatchSolverInfo(const std::string &pName, const std::string &pTrueVal) const
        {
            if (DefinesSolverInfo(pName))
            {
                std::string vName = boost::to_upper_copy(pName);
                SolverInfoMap::const_iterator solverInfoMapIter = m_solverInfo.find(vName);
                if(solverInfoMapIter != m_solverInfo.end())
                {
                    return true;
                }
            }
            return false;
        }


        /**
         *
         */
        bool SessionReader::DefinesGeometricInfo(const std::string &pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            GeometricInfoMap::const_iterator geometricInfoMapIter = m_geometricInfo.find(vName);
            return (geometricInfoMapIter != m_geometricInfo.end());
        }


        /**
         *
         */
        void SessionReader::LoadGeometricInfo(const std::string &pName,
                                std::string &pVar, const std::string &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            GeometricInfoMap::const_iterator geometricInfoMapIter = m_geometricInfo.find(vName);
            if(geometricInfoMapIter != m_geometricInfo.end())
            {
                pVar = geometricInfoMapIter->second;
            }
            else
            {
                pVar  = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::LoadGeometricInfo(const std::string &pName, bool &pVar,
                                const bool &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            GeometricInfoMap::const_iterator geometricInfoMapIter = m_geometricInfo.find(vName);
            if(geometricInfoMapIter != m_geometricInfo.end())
            {
                if (geometricInfoMapIter->second == "TRUE")
                {
                    pVar = true;
                }
                else
                {
                    pVar = false;
                }
            }
            else
            {
                pVar  = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::LoadGeometricInfo(const std::string &pName, NekDouble &pVar,
                                const NekDouble &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            GeometricInfoMap::const_iterator geometricInfoMapIter = m_geometricInfo.find(vName);
            if(geometricInfoMapIter != m_geometricInfo.end())
            {
                pVar = std::atoi(geometricInfoMapIter->second.c_str());
            }
            else
            {
                pVar  = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::MatchGeometricInfo(const std::string &pName,
                                const std::string &pTrueVal, bool &pVar,
                                const bool &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            GeometricInfoMap::const_iterator geometricInfoMapIter = m_geometricInfo.find(vName);
            if(geometricInfoMapIter != m_geometricInfo.end())
            {
                pVar = (NoCaseStringCompare(geometricInfoMapIter->second, pTrueVal) == 0);
            }
            else
            {
                pVar  = pDefault;
            }
        }


        /**
         *
         */
        const std::string& SessionReader::GetVariable(const unsigned int &idx) const
        {
            ASSERTL0(idx < m_variables.size(), "Variable index out of range.");
            return m_variables[idx];
        }


        /**
         *
         */
        std::vector<std::string> SessionReader::GetVariables() const
        {
            return m_variables;
        }


        /**
         *
         */
        bool SessionReader::DefinesFunction(const std::string &pName) const
        {
            FunctionMap::const_iterator it1;
            std::string vName = boost::to_upper_copy(pName);

            if ((it1 = m_functions.find(vName)) != m_functions.end())
            {
                return true;
            }
            return false;
        }


        /**
         *
         */
        bool SessionReader::DefinesFunction(const std::string &pName, const std::string &pVariable) const
        {
            FunctionMap::const_iterator it1;
            EquationMap::const_iterator it2;
            std::string vName = boost::to_upper_copy(pName);

            if ((it1 = m_functions.find(vName)) != m_functions.end()
                    && (it2 = it1->second.m_expressions.find(pVariable))
                            != it1->second.m_expressions.end())
            {
                return true;
            }
            return false;
        }


        /**
         *
         */
        EquationSharedPtr SessionReader::GetFunction(const std::string &pName, const std::string &pVariable) const
        {
            FunctionMap::const_iterator it1;
            EquationMap::const_iterator it2;
            std::string vName = boost::to_upper_copy(pName);

            ASSERTL0((it1 = m_functions.find(vName)) != m_functions.end(),
                     std::string("No such function '") + pName
                     + std::string("' has been defined in the session file."));
            ASSERTL0((it2 = it1->second.m_expressions.find(pVariable)) != it1->second.m_expressions.end(),
                     std::string("No such variable '") + pVariable
                     + std::string("' defined for function '") + pName
                     + std::string("' in session file."));
            return it2->second;
        }


        /**
         *
         */
        EquationSharedPtr SessionReader::GetFunction(const std::string &pName, const unsigned int &pVar) const
        {
            ASSERTL0(pVar < m_variables.size(), "Variable index out of range.");
            return GetFunction(pName, m_variables[pVar]);
        }


        /**
         *
         */
        enum FunctionType SessionReader::GetFunctionType(const std::string &pName) const
        {
            FunctionMap::const_iterator it1;
            std::string vName = boost::to_upper_copy(pName);

            it1 = m_functions.find(vName);
            ASSERTL0 (it1 != m_functions.end(),
                      std::string("Function '") + pName
                      + std::string("' not found."));
            return it1->second.m_type;
        }


        /**
         *
         */
        std::string SessionReader::GetFunctionFilename(const std::string &pName) const
        {
            FunctionMap::const_iterator it1;
            std::string vName = boost::to_upper_copy(pName);

            it1 = m_functions.find(vName);
            ASSERTL0 (it1 != m_functions.end(),
                      std::string("Function '") + pName
                      + std::string("' not found."));
            return it1->second.m_filename;
        }


        /**
         *
         */
        bool SessionReader::DefinesTag(const std::string &pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            TagMap::const_iterator vTagIterator = m_tags.find(vName);
            return (vTagIterator != m_tags.end());
        }


        /**
         *
         */
        void SessionReader::SetTag(const std::string &pName, const std::string &pValue)
        {
            std::string vName = boost::to_upper_copy(pName);
            m_tags[vName] = pValue;
        }


        /**
         *
         */
        const std::string &SessionReader::GetTag(const std::string& pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            TagMap::const_iterator vTagIterator = m_tags.find(vName);
            ASSERTL0(vTagIterator != m_tags.end(),
                     "Requested tag does not exist.");
            return vTagIterator->second;
        }


        /**
         *
         */
        void SessionReader::SubstituteExpressions(std::string& pExpr)
        {
            ExpressionMap::iterator exprIter;
            for (exprIter = m_expressions.begin(); exprIter != m_expressions.end(); ++exprIter)
            {
                //boost::regex re("\b" + exprIter->first + "\b");
                //boost::replace_all_regex(pExpr, re,
                //        std::string("(") + exprIter->second + std::string(")"));
                boost::replace_all(pExpr, exprIter->first, exprIter->second);
            }
        }


        /**
         *
         */
        TiXmlDocument *SessionReader::MergeDoc(const std::vector<std::string> &pFilenames) const
        {
            ASSERTL0(pFilenames.size() > 0, "No filenames for merging.");

            // Read the first document
            TiXmlDocument *vMainDoc = new TiXmlDocument(pFilenames[0]);
            ASSERTL0(vMainDoc, "Failed to create XML document object.");
            bool loadOkay = vMainDoc->LoadFile();
            ASSERTL0(loadOkay, std::string("Unable to load file: ") +
                    pFilenames[0] + ". Check XML standards compliance. Error on line: "
                    + boost::lexical_cast<std::string>(vMainDoc->Row()) + ": " + std::string(vMainDoc->ErrorDesc()));
            TiXmlHandle vMainHandle(vMainDoc);
            TiXmlElement* vMainNektar = vMainHandle.FirstChildElement("NEKTAR").Element();

            // Read all subsequent XML documents.
            // For each element within the NEKTAR tag, use it to replace the
            // version already present in the loaded XML data.
            for (int i = 1; i < pFilenames.size(); ++i)
            {
                TiXmlDocument vTempDoc (pFilenames[i]);
                loadOkay = vTempDoc.LoadFile();
                ASSERTL0(loadOkay, std::string("Unable to load file: ") +
                    pFilenames[i] + ". Check XML standards compliance. Error on line: "
                    + boost::lexical_cast<std::string>(vTempDoc.Row()));

                TiXmlHandle docHandle(&vTempDoc);
                TiXmlElement* vTempNektar;
                vTempNektar = docHandle.FirstChildElement("NEKTAR").Element();
                ASSERTL0(vTempNektar, "Unable to find NEKTAR tag in file.");
                TiXmlElement* p = vTempNektar->FirstChildElement();

                while (p)
                {
                    TiXmlElement * vMainEntry = vMainNektar->FirstChildElement(p->Value());
                    TiXmlElement * q = new TiXmlElement(*p);
                    if (vMainEntry)
                    {
                        vMainNektar->RemoveChild(vMainEntry);
                    }
                    vMainNektar->LinkEndChild(q);
                    p = p->NextSiblingElement();
                }
            }

            return vMainDoc;
        }


        /**
         *
         */
        void SessionReader::ParseDocument()
        {
            // Check we actually have a document loaded.
            ASSERTL0(m_xmlDoc, "No XML document loaded.");

            // Look for all data in CONDITIONS block.
            TiXmlHandle docHandle(m_xmlDoc);
            TiXmlElement* e;
            e = docHandle.FirstChildElement("NEKTAR").FirstChildElement("CONDITIONS").Element();
            ASSERTL0(e, "Unable to find CONDITIONS tag in file.");

            // Read the various sections of the document
            ReadParameters(e);
            ReadSolverInfo(e);
            ReadExpressions(e);
            ReadVariables (e);
            ReadFunctions (e);

            e = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();

            ReadGeometricInfo(e);
        }


        /**
         *
         */
        void SessionReader::CreateComm(int &argc, char* argv[], const std::string &pFilename)
        {
            if (argc == 0)
            {
                m_comm = GetCommFactory().CreateInstance("Serial", 0, 0);
            }
            else
            {
                TiXmlHandle docHandle(m_xmlDoc);
                TiXmlElement* e;
                e = docHandle.FirstChildElement("NEKTAR").FirstChildElement("CONDITIONS").Element();
                ASSERTL0(e, "Unable to find CONDITIONS tag in file.");
                ReadSolverInfo(e);

                string vCommModule("Serial");
                if (DefinesSolverInfo("Communication"))
                {
                    vCommModule = GetSolverInfo("Communication");
                }
                else if (GetCommFactory().ModuleExists("ParallelMPI"))
                {
                    vCommModule = "ParallelMPI";
                }

                m_comm = GetCommFactory().CreateInstance(vCommModule,argc,argv);

                //If running in parallel change the default global sys soln type
                if (m_comm->GetSize() > 1)
                {
                    m_solverInfoDefaults["GLOBALSYSSOLN"] = "IterativeStaticCond";
                }
            }
        }


        /**
         *
         */
        void SessionReader::PartitionMesh()
        {
            ASSERTL0(m_comm.get(), "Communication not initialised.");
            if (m_comm->GetSize() > 1)
            {
                if (m_comm->GetRank() == 0)
                {
                    SessionReaderSharedPtr vSession = GetSharedThisPtr();
                    MeshPartitionSharedPtr vPartitioner = MemoryManager<MeshPartition>::AllocateSharedPtr(vSession);
                    vPartitioner->PartitionMesh(m_comm->GetSize());
                    vPartitioner->WritePartitions(vSession);
                }

                m_comm->Block();

                m_sessionName += "_P" + boost::lexical_cast<std::string>(m_comm->GetRank());
                m_filename = m_sessionName + ".xml";

                delete m_xmlDoc;
                m_xmlDoc = new TiXmlDocument(m_filename);
                ASSERTL0(m_xmlDoc, "Failed to create XML document object.");

                bool loadOkay = m_xmlDoc->LoadFile();
                ASSERTL0(loadOkay, std::string("Unable to load file: ") +
                        m_filename + ". Check XML standards compliance. Error on line: "
                        + boost::lexical_cast<std::string>(m_xmlDoc->Row()));
            }
        }


        /**
         *
         */
        void SessionReader::ReadParameters(TiXmlElement *conditions)
        {
            m_parameters.clear();

            TiXmlElement *parametersElement = conditions->FirstChildElement("PARAMETERS");

            // See if we have parameters defined.  They are optional so we go on if not.
            if (parametersElement)
            {
                TiXmlElement *parameter = parametersElement->FirstChildElement("P");
                LibUtilities::ExpressionEvaluator expEvaluator;
                ParameterMap caseSensitiveParameters;

                // Multiple nodes will only occur if there is a comment in between
                // definitions.
                while (parameter)
                {
                    TiXmlNode *node = parameter->FirstChild();

                    while (node && node->Type() != TiXmlNode::TEXT)
                    {
                        node = node->NextSibling();
                    }

                    if (node)
                    {
                        // Format is "paramName = value"
                        std::string line = node->ToText()->Value();
                        std::string lhs;
                        std::string rhs;

                        try {
                            /// Pull out lhs and rhs and eliminate any spaces.
                            int beg = line.find_first_not_of(" ");
                            int end = line.find_first_of("=");
                            // Check for no parameter name
                            if (beg == end) throw 1;
                            // Check for no parameter value
                            if (end != line.find_last_of("=")) throw 1;
                            // Check for no equals sign
                            if (end == std::string::npos) throw 1;

                            lhs = line.substr(line.find_first_not_of(" "), end-beg);
                            lhs = lhs.substr(0, lhs.find_last_not_of(" ")+1);

                            rhs = line.substr(line.find_last_of("=")+1);
                            rhs = rhs.substr(rhs.find_first_not_of(" "));
                            rhs = rhs.substr(0, rhs.find_last_not_of(" ")+1);
                        }
                        catch (...)
                        {
                            ASSERTL0(false, "Syntax error. "
                                    "File: '" + m_filename + "', line: "
                                    + boost::lexical_cast<string>(node->Row()));
                        }

                        /// We want the list of parameters to have their RHS evaluated,
                        /// so we use the expression evaluator to do the dirty work.
                        if (!lhs.empty() && !rhs.empty())
                        {
                            NekDouble value=0.0;
                            try
                            {
                                expEvaluator.DefineFunction("", rhs);
                                value =  expEvaluator.Evaluate();
                            }
                            catch (const std::runtime_error &)
                            {
                                ASSERTL0(false, "Error evaluating parameter expression '" + rhs + "'." );
                            }
                            expEvaluator.SetParameter(lhs, value);
                            caseSensitiveParameters[lhs] = value;
                            boost::to_upper(lhs);
                            m_parameters[lhs] = value;
                        }
                    }

                    parameter = parameter->NextSiblingElement();
                }

                try
                {
                    // Set ourselves up for evaluation later.
                    Equation::SetConstParameters(caseSensitiveParameters);
                }
                catch (const std::runtime_error&)
                {
                    // Attempted to set parameters more than once, but we let
                    // this go.
                }
            }
        }


        /**
         *
         */
        void SessionReader::ReadSolverInfo(TiXmlElement *conditions)
        {
            m_solverInfo.clear();
            m_solverInfo = m_solverInfoDefaults;

            TiXmlElement *solverInfoElement = conditions->FirstChildElement("SOLVERINFO");

            if (solverInfoElement)
            {
                TiXmlElement *solverInfo = solverInfoElement->FirstChildElement("I");

                while (solverInfo)
                {
                    // read the property name
                    ASSERTL0(solverInfo->Attribute("PROPERTY"),
                            "Missing PROPERTY attribute in solver info section. "
                            "File: '" + m_filename + "', line: "
                            + boost::lexical_cast<string>(solverInfo->Row()));
                    std::string solverProperty = solverInfo->Attribute("PROPERTY");
                    ASSERTL0(!solverProperty.empty(),
                            "Solver info properties must have a non-empty name. "
                            "File: '" + m_filename + "', line: "
                            + boost::lexical_cast<string>(solverInfo->Row()));

                    // make sure that solver property is capitalised
                    std::string solverPropertyUpper = boost::to_upper_copy(solverProperty);

                    // read the value
                    ASSERTL0(solverInfo->Attribute("VALUE"),
                            "Missing VALUE attribute in solver info section. "
                            "File: '" + m_filename + "', line: "
                            + boost::lexical_cast<string>(solverInfo->Row()));
                    std::string solverValue    = solverInfo->Attribute("VALUE");
                    ASSERTL0(!solverValue.empty(),
                            "Solver info properties must have a non-empty value. "
                            "File: '" + m_filename + "', line: "
                            + boost::lexical_cast<string>(solverInfo->Row()));

                    EnumMapList::const_iterator propIt = m_enums.find(solverPropertyUpper);
                    if (propIt != m_enums.end())
                    {
                        EnumMap::const_iterator valIt = propIt->second.find(solverValue);
                        ASSERTL0(valIt != propIt->second.end(),
                                "Value '" + solverValue + "' is not valid for property '"
                               + solverProperty + "'. File: '" + m_filename
                               + "', line: "
                               + boost::lexical_cast<string>(solverInfo->Row()));
                    }

                    // Set Variable
                    m_solverInfo[solverPropertyUpper] = solverValue;
                    solverInfo = solverInfo->NextSiblingElement("I");
                }
            }

            if (m_comm && m_comm->GetSize() > 1)
            {
                ASSERTL0 (m_solverInfo["GLOBALSYSSOLN"] == "IterativeFull"
                    || m_solverInfo["GLOBALSYSSOLN"] == "IterativeStaticCond",
                    "An iterative solver must be used when run in parallel.");
            }
        }


        /**
         *
         */
        void SessionReader::ReadGeometricInfo(TiXmlElement *geometry)
        {
            m_geometricInfo.clear();

            TiXmlElement *geometricInfoElement = geometry->FirstChildElement("GEOMINFO");

            if (geometricInfoElement)
            {
                TiXmlElement *geometricInfo = geometricInfoElement->FirstChildElement("I");

                while (geometricInfo)
                {
                    ASSERTL0(geometricInfo->Attribute("PROPERTY"),
                            "Missing PROPERTY attribute in geometric info section. "
                            "File: '" + m_filename + "', line: "
                            + boost::lexical_cast<string>(geometricInfo->Row()));
                    std::string geometricProperty = geometricInfo->Attribute("PROPERTY");
                    ASSERTL0(!geometricProperty.empty(),
                            "Geometric info properties must have a non-empty name. "
                            "File: '" + m_filename + "', line: "
                            + boost::lexical_cast<string>(geometricInfo->Row()));

                    // make sure that geometric property is capitalised
                    boost::to_upper(geometricProperty);

                    // check the property has not already been defined
                    GeometricInfoMap::iterator geometricInfoIter = m_geometricInfo.find(geometricProperty);
                    ASSERTL0(geometricInfoIter == m_geometricInfo.end(),
                             (std::string("geometricInfo value: ") + geometricProperty
                              + std::string(" already specified.")).c_str());

                    // read the property value
                    ASSERTL0(geometricInfo->Attribute("VALUE"),
                            "Missing VALUE attribute in geometric info section. "
                            "File: '" + m_filename + ", line: "
                            + boost::lexical_cast<string>(geometricInfo->Row()));
                    std::string geometricValue    = geometricInfo->Attribute("VALUE");
                    ASSERTL0(!geometricValue.empty(),
                            "Geometric info properties must have a non-empty value. "
                            "File: '" + m_filename + ", line: "
                            + boost::lexical_cast<string>(geometricInfo->Row()));

                    // Set Variable
                    m_geometricInfo[geometricProperty] = geometricValue;
                    geometricInfo = geometricInfo->NextSiblingElement("I");
                }
            }
        }


        /**
         *
         */
        void SessionReader::ReadExpressions(TiXmlElement *conditions)
        {
            m_expressions.clear();

            TiXmlElement *expressionsElement = conditions->FirstChildElement("EXPRESSIONS");

            if (expressionsElement)
            {
                TiXmlElement *expr = expressionsElement->FirstChildElement("E");

                while (expr)
                {
                    ASSERTL0(expr->Attribute("NAME"),
                             "Missing NAME attribute in expression definition. "
                             "File: '" + m_filename + "', line: "
                             + boost::lexical_cast<std::string>(expr->Row()));
                    std::string nameString = expr->Attribute("NAME");
                    ASSERTL0(!nameString.empty(),
                             "Expressions must have a non-empty name. "
                             "File: '" + m_filename + "', line: "
                             + boost::lexical_cast<std::string>(expr->Row()));

                    ASSERTL0(expr->Attribute("VALUE"),
                             "Missing VALUE attribute in expression definition. "
                             "File: '" + m_filename + "', line: "
                             + boost::lexical_cast<std::string>(expr->Row()));
                    std::string valString = expr->Attribute("VALUE");
                    ASSERTL0(!valString.empty(),
                             "Expressions must have a non-empty value. "
                             "File: '" + m_filename + "', line: "
                             + boost::lexical_cast<std::string>(expr->Row()));

                    ExpressionMap::iterator exprIter
                                            = m_expressions.find(nameString);
                    ASSERTL0(exprIter == m_expressions.end(),
                             std::string("Expression '") + nameString
                             + std::string("' already specified."));

                    m_expressions[nameString] = valString;
                    expr = expr->NextSiblingElement("E");
                }
            }
        }


        /**
         *
         */
        void SessionReader::ReadVariables(TiXmlElement *conditions)
        {
            m_variables.clear();

            TiXmlElement *variablesElement = conditions->FirstChildElement("VARIABLES");

            int varIndex = 0;   // Current index, should be zero-based.

            // See if we have parameters defined.  They are optional so we go on if not.
            if (variablesElement)
            {
                TiXmlElement *variableElement = variablesElement->FirstChildElement("V");

                // Sequential counter for the composite numbers.
                int nextVariableNumber = -1;

                while (variableElement)
                {
                    /// All elements are of the form: "<V ID="#"> name = value </V>", with
                    /// ? being the element type.

                    nextVariableNumber++;

                    int indx;
                    int err = variableElement->QueryIntAttribute("ID", &indx);
                    ASSERTL0(err == TIXML_SUCCESS,
                            "Variables must have a unique ID number attribute. "
                            "File: '" + m_filename + "', line: "
                            + boost::lexical_cast<string>(variableElement->Row()));
                    ASSERTL0(indx == nextVariableNumber,
                            "ID numbers for variables must begin with zero and be sequential. "
                            "File: '" + m_filename + "', line: "
                            + boost::lexical_cast<string>(variableElement->Row()));

                    TiXmlNode* variableChild = variableElement->FirstChild();
                    // This is primarily to skip comments that may be present.
                    // Comments appear as nodes just like elements.
                    // We are specifically looking for text in the body
                    // of the definition.
                    while(variableChild && variableChild->Type() != TiXmlNode::TEXT)
                    {
                        variableChild = variableChild->NextSibling();
                    }

                    ASSERTL0(variableChild,
                            "Unable to read variable definition body for variable with ID "
                            + boost::lexical_cast<string>(indx) + ". "
                            " File: '" + m_filename + "', line: "
                            + boost::lexical_cast<string>(variableElement->Row()));
                    std::string variableName = variableChild->ToText()->ValueStr();

                    std::istringstream variableStrm(variableName);
                    variableStrm >> variableName;

                    ASSERTL0(std::find(m_variables.begin(), m_variables.end(), variableName) == m_variables.end(),
                            "Variable with ID " + boost::lexical_cast<string>(indx)
                            + " has already been defined. "
                            " File: '" + m_filename + "', line: "
                            + boost::lexical_cast<string>(variableElement->Row()));

                    m_variables.push_back(variableName);

                    variableElement = variableElement->NextSiblingElement("V");
                }

                ASSERTL0(nextVariableNumber > -1, "Number of variables must be greater than zero.");
            }
        }


        /**
         *
         */
        void SessionReader::ReadFunctions(TiXmlElement *conditions)
        {
            m_functions.clear();

            // Scan through conditions section looking for functions.
            TiXmlElement *function = conditions->FirstChildElement("FUNCTION");
            while (function)
            {
                // Every function must have a NAME attribute
                ASSERTL0(function->Attribute("NAME"),
                         "Functions must have a NAME attribute defined. "
                         "File: '" + m_filename + "', line: "
                         + boost::lexical_cast<std::string>(function->Row()));
                std::string functionStr = function->Attribute("NAME");
                ASSERTL0(!functionStr.empty(),
                         "Functions must have a non-empty name. "
                         "File: '" + m_filename + "', line: "
                         + boost::lexical_cast<string>(function->Row()));

                // Store function names in uppercase to remain case-insensitive.
                boost::to_upper(functionStr);

                // Retrieve first entry (variable, or file)
                TiXmlElement *variable  = function->FirstChildElement();

                // Create new function structure with default type of none.
                FunctionDefinition functionDef;
                functionDef.m_type = eFunctionTypeNone;

                // Process all entries in the function block
                while (variable)
                {
                    std::string conditionType = variable->Value();

                    // Expressions are denoted by E
                    if (conditionType == "E")
                    {
                        // Ensure we haven't already found a file to read.
                        ASSERTL0(functionDef.m_type != eFunctionTypeFile,
                               "Cannot mix expressions and files in function.");
                        functionDef.m_type = eFunctionTypeExpression;

                        // Expression must have a VAR and VALUE.
                        ASSERTL0(variable->Attribute("VAR"),
                                 "Attribute VAR expected for function '"
                                 + functionStr + "'.");
                        std::string variableStr = variable->Attribute("VAR");

                        ASSERTL0(variable->Attribute("VALUE"),
                                 "Attribute VALUE expected for function '"
                                 + functionStr + "'.");
                        std::string fcnStr      = variable->Attribute("VALUE");

                        ASSERTL0(!fcnStr.empty(),
                                 (std::string("Expression for var: ")
                                 + variableStr
                                 + std::string(" must be specified.")).c_str());

                        SubstituteExpressions(fcnStr);

                        // Check it has not already been defined
                        EquationMap::iterator fcnsIter
                                = functionDef.m_expressions.find(variableStr);

                        ASSERTL0(fcnsIter == functionDef.m_expressions.end(),
                                "Error setting expression '" + variableStr
                                + "' in function '" + functionStr + "'. "
                                "Expression has already been defined.");

                        // Add variable
                        functionDef.m_expressions[variableStr]
                            = MemoryManager<Equation>::AllocateSharedPtr(fcnStr);
                    }

                    // Files are denoted by F
                    else if (conditionType == "F")
                    {
                        // Ensure we haven't already read expressions
                        ASSERTL0(functionDef.m_type != eFunctionTypeExpression,
                               "Cannot mix expressions and files in function.");
                        functionDef.m_type = eFunctionTypeFile;

                        // A file must specify the FILE attribute
                        ASSERTL0(variable->Attribute("FILE"),
                                 "Attribute FILE expected for function '"
                                 + functionStr + "'.");
                        std::string filenameStr = variable->Attribute("FILE");

                        ASSERTL0(!filenameStr.empty(),
                                 "A filename must be specified for the FILE "
                                 "attribute of function '" + functionStr
                                 + "'.");

                        // set the filename for the function structure
                        functionDef.m_filename = filenameStr;
                    }

                    // Nothing else supported so throw an error
                    else
                    {
                        ASSERTL0(false,
                                "Identifier " + conditionType + " in function "
                                + std::string(function->Attribute("NAME"))
                                + " is not recognised. "
                                "File: '" + m_filename + "', line: "
                                + boost::lexical_cast<string>(variable->Row()));
                    }
                    variable = variable->NextSiblingElement();
                }
                // Add function definition to map
                m_functions[functionStr] = functionDef;
                function = function->NextSiblingElement("FUNCTION");
            }
        }


        /**
         *
         */
        int SessionReader::NoCaseStringCompare(const std::string & s1, const std::string& s2) const
        {
            std::string::const_iterator it1=s1.begin();
            std::string::const_iterator it2=s2.begin();

            //stop when either string's end has been reached
            while ( (it1!=s1.end()) && (it2!=s2.end()) )
            {
                if(::toupper(*it1) != ::toupper(*it2)) //letters differ?
                {
                    // return -1 to indicate smaller than, 1 otherwise
                    return (::toupper(*it1)  < ::toupper(*it2)) ? -1 : 1;
                }

                //proceed to the next character in each string
                ++it1;
                ++it2;
            }

            size_t size1=s1.size();
            size_t size2=s2.size();// cache lengths

            //return -1,0 or 1 according to strings' lengths
            if (size1==size2)
            {
                return 0;
            }

            return (size1 < size2) ? -1 : 1;
        }
    }
}
