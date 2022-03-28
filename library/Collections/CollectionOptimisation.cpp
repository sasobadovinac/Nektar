///////////////////////////////////////////////////////////////////////////////
//
// File: CollectionOptimisation.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Collection top class definition
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <Collections/CollectionOptimisation.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/Timer.h>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <tinyxml.h>

using namespace std;

namespace Nektar
{
namespace Collections
{

// static manager for Operator ImplementationMap
map<OpImpTimingKey,OperatorImpMap> CollectionOptimisation::m_opImpMap;

CollectionOptimisation::CollectionOptimisation(
        LibUtilities::SessionReaderSharedPtr pSession,
        const int shapedim,
        ImplementationType defaultType):
    m_shapeDim(shapedim)
{
    map<ElmtOrder, ImplementationType> defaults, defaultsPhysDeriv,
        defaultsHelmholtz;
    bool verbose  = (pSession.get()) &&
                    (pSession->DefinesCmdLineArgument("verbose")) &&
                    (pSession->GetComm()->GetRank() == 0);

    m_autotune    = false;
    m_maxCollSize = 0;
    m_defaultType = defaultType == eNoImpType ? eIterPerExp : defaultType;

    map<string, LibUtilities::ShapeType> elTypes;
    elTypes["S"] = LibUtilities::eSegment;
    elTypes["T"] = LibUtilities::eTriangle;
    elTypes["Q"] = LibUtilities::eQuadrilateral;
    elTypes["A"] = LibUtilities::eTetrahedron;
    elTypes["P"] = LibUtilities::ePyramid;
    elTypes["R"] = LibUtilities::ePrism;
    elTypes["H"] = LibUtilities::eHexahedron;

    // Set defaults for all element types.
    for (auto &it2 : elTypes)
    {
        defaults          [ElmtOrder(it2.second, -1)] = m_defaultType;
        defaultsPhysDeriv [ElmtOrder(it2.second, -1)] = m_defaultType;
        defaultsHelmholtz [ElmtOrder(it2.second, -1)] = m_defaultType;
    }

    if (defaultType == eNoImpType)
    {
        for (auto &it2 : elTypes)
        {
            // use Nocollection for Phys Deriv 
            defaultsPhysDeriv [ElmtOrder(it2.second, -1)] = eNoCollection;

            // Use IterPerExp 
            defaultsHelmholtz[ElmtOrder(it2.second, -1)] = eMatrixFree;
        }
    }

    map<string, OperatorType> opTypes;
    for (int i = 0; i < SIZE_OperatorType; ++i)
    {
        opTypes[OperatorTypeMap[i]] = (OperatorType)i;
        switch ((OperatorType)i)
        {
        case eHelmholtz:
            m_global[(OperatorType)i] = defaultsHelmholtz;
            break;
        case ePhysDeriv:
            m_global[(OperatorType)i] = defaultsPhysDeriv;
            break;
        default:
            m_global[(OperatorType)i] = defaults;
        }
    }

    map<string, ImplementationType> impTypes;
    for (int i = 0; i < SIZE_ImplementationType; ++i)
    {
        impTypes[ImplementationTypeMap[i]] = (ImplementationType)i;
    }

    // turn off file reader if dummy pointer is given or if default
    // option is passed and by default calling argument.
    if ((defaultType == eNoImpType)&&(pSession.get()))
    {
        TiXmlDocument &doc = pSession->GetDocument();
        TiXmlHandle docHandle(&doc);
        TiXmlElement *master = docHandle.FirstChildElement("NEKTAR").Element();
        ASSERTL0(master, "Unable to find NEKTAR tag in file.");
        bool WriteFullCollections = false;
        
        TiXmlElement *xmlCol = master->FirstChildElement("COLLECTIONS");

        // Check if user has specified some options
        if (xmlCol)
        {
            // Set the maxsize and default implementation type if provided
            const char *maxSize = xmlCol->Attribute("MAXSIZE");
            m_maxCollSize = (maxSize ? atoi(maxSize) : 0);

            const char *defaultImpl = xmlCol->Attribute("DEFAULT");
            m_defaultType = defaultType;

            // If user has specified a default impl type or  autotuning
            // and set this default across all operators.
            if (defaultImpl)
            {
                const std::string collinfo = string(defaultImpl);
                m_autotune = boost::iequals(collinfo, "auto");

                if (!m_autotune)
                {
                    bool collectionFound{false};
                    for(int i = 1; i < Collections::SIZE_ImplementationType; ++i)
                    {
                        if(boost::iequals(collinfo,
                                Collections::ImplementationTypeMap[i]))
                        {
                            m_defaultType = (Collections::ImplementationType) i;
                            collectionFound = true;
                            break;
                        }
                    }

                    ASSERTL0(collectionFound,
                        "Unknown default collection scheme: "+collinfo);

                    defaults.clear();
                    // Override default types
                    for (auto &it2 : elTypes)
                    {
                        defaults[ElmtOrder(it2.second, -1)] = m_defaultType;
                    }

                    for (int i = 0; i < SIZE_OperatorType; ++i)
                    {
                        m_global[(OperatorType)i] = defaults;
                    }
                }
            }            
            const char *write = xmlCol->Attribute("WRITE");
            if(write&&boost::iequals(write,"true"))
            {
                WriteFullCollections = true; 
            }

            // Now process operator-specific implementation selections
            ReadCollOps(xmlCol,m_global,verbose);

            // Print out operator map
            if (verbose)
            {
                if(WriteFullCollections)
                {
                    for (auto &mIt : m_global)
                    {
                        cout << "Operator " << OperatorTypeMap[mIt.first]
                             << ":" << endl;

                        for (auto &eIt : mIt.second)
                        {
                            cout << "- "
                                 << LibUtilities::ShapeTypeMap[eIt.first.first]
                                 << " order " << eIt.first.second << " -> "
                                 << ImplementationTypeMap[eIt.second] << endl;
                        }
                    }

                }

            }
        }
    }
}

void  CollectionOptimisation::ReadCollOps(TiXmlElement *xmlCol,
                                          GlobalOpMap &global, bool verbose)
{
    bool verboseHeader = true;
    map<string, LibUtilities::ShapeType> elTypes;
    elTypes["S"] = LibUtilities::eSegment;
    elTypes["T"] = LibUtilities::eTriangle;
    elTypes["Q"] = LibUtilities::eQuadrilateral;
    elTypes["A"] = LibUtilities::eTetrahedron;
    elTypes["P"] = LibUtilities::ePyramid;
    elTypes["R"] = LibUtilities::ePrism;
    elTypes["H"] = LibUtilities::eHexahedron;

    map<string, OperatorType> opTypes;
    for (int i = 0; i < SIZE_OperatorType; ++i)
    {
        opTypes[OperatorTypeMap[i]] = (OperatorType)i;
    }

    map<string, ImplementationType> impTypes;
    for (int i = 0; i < SIZE_ImplementationType; ++i)
    {
        impTypes[ImplementationTypeMap[i]] = (ImplementationType)i;
    }
    

    TiXmlElement *elmt = xmlCol->FirstChildElement();
    while (elmt)
    {
        string tagname = elmt->ValueStr();
        
        ASSERTL0(boost::iequals(tagname, "OPERATOR"),
                 "Only OPERATOR tags are supported inside the "
                 "COLLECTIONS tag.");
        
        const char *attr = elmt->Attribute("TYPE");
        ASSERTL0(attr, "Missing TYPE in OPERATOR tag.");
        string opType(attr);

        ASSERTL0(opTypes.count(opType) > 0,
                 "Unknown OPERATOR type " + opType + ".");
        
        OperatorType ot = opTypes[opType];
        
        
        TiXmlElement *elmt2 = elmt->FirstChildElement();

        map<int, pair<int,std::string>> verboseWrite; 
        while (elmt2)
        {
            string tagname = elmt2->ValueStr();
            ASSERTL0(boost::iequals(tagname, "ELEMENT"),
                     "Only ELEMENT tags are supported inside the "
                     "OPERATOR tag.");
            
            const char *attr = elmt2->Attribute("TYPE");
            ASSERTL0(attr, "Missing TYPE in ELEMENT tag.");
            
            string elType(attr);
            auto it2 = elTypes.find(elType);
            ASSERTL0(it2 != elTypes.end(),
                     "Unknown element type "+elType+" in ELEMENT "
                     "tag");
            
            const char *attr2 = elmt2->Attribute("IMPTYPE");
            ASSERTL0(attr2, "Missing IMPTYPE in ELEMENT tag.");
            string impType(attr2);
            ASSERTL0(impTypes.count(impType) > 0,
                     "Unknown IMPTYPE type " + impType + ".");
            
            const char *attr3 = elmt2->Attribute("ORDER");
            ASSERTL0(attr3, "Missing ORDER in ELEMENT tag.");
            string order(attr3);
            
            // load details relevant to this shape dimension. 
            if(LibUtilities::ShapeTypeDimMap[it2->second] == m_shapeDim)
            {
                if (order == "*")
                {
                    global[ot][ElmtOrder(it2->second, -1)]
                        = impTypes[impType];
                    
                    if(verbose)
                    {
                        verboseWrite[it2->second] =
                            pair<int,std::string>(-1,impType);
                    }
                }
                else
                {
                    vector<unsigned int> orders;
                    ParseUtils::GenerateSeqVector(order, orders);
                    
                    for (int i = 0; i < orders.size(); ++i)
                    {
                        global[ot][ElmtOrder(it2->second, orders[i])]
                            = impTypes[impType];
                        
                        if(verbose)
                        {
                            verboseWrite[it2->second] =
                                pair<int,std::string>(orders[i],impType);

                        }
                    }
                }
            }
            
            elmt2 = elmt2->NextSiblingElement();
        }
        
        if(verboseWrite.size())
        {
            if(verboseHeader)
            {
                cout << "Collection settings from file:  " << endl;
                verboseHeader = false;
            }
            
            cout << "\t Operator " << OperatorTypeMap[ot]
                 << ":" << endl;
            
            for(auto &it: verboseWrite)
            {
                cout << "\t - "
                     << LibUtilities::ShapeTypeMap[it.first]
                     << " order " << it.second.first << " -> "
                     << it.second.second << endl;
            }
            
        }
        
        elmt = elmt->NextSiblingElement();
    }
}
    
OperatorImpMap  CollectionOptimisation::GetOperatorImpMap(
        StdRegions::StdExpansionSharedPtr pExp)
{
    OperatorImpMap ret;
    ElmtOrder searchKey(pExp->DetShapeType(),
                        pExp->GetBasisNumModes(0));
    ElmtOrder defSearch(pExp->DetShapeType(), -1);

    for (auto &it : m_global)
    {
        ImplementationType impType;

        auto it2 = it.second.find(searchKey);

        if (it2 == it.second.end())
        {
            it2 = it.second.find(defSearch);
            if (it2 == it.second.end())
            {
                // Shouldn't be able to reach here.
                impType = eNoCollection;
            }
            else
            {
                impType = it2->second;
            }
        }
        else
        {
            impType = it2->second;
        }

        ret[it.first] = impType;
    }

    return ret;
}

OperatorImpMap CollectionOptimisation::SetWithTimings(
        vector<StdRegions::StdExpansionSharedPtr> pCollExp,
        OperatorImpMap &impTypes,
        bool verbose )
{
    boost::ignore_unused(impTypes);

    OperatorImpMap ret;

    StdRegions::StdExpansionSharedPtr pExp = pCollExp[0];

    // check to see if already defined for this expansion
    OpImpTimingKey OpKey(pExp,pCollExp.size(),pExp->GetNumBases());
    if(m_opImpMap.count(OpKey) != 0)
    {
        ret = m_opImpMap[OpKey];
        return ret;
    }

    int maxsize = pCollExp.size()*max(pExp->GetNcoeffs(),pExp->GetTotPoints());
    Array<OneD, NekDouble> inarray(maxsize,1.0);
    Array<OneD, NekDouble> outarray1(maxsize);
    Array<OneD, NekDouble> outarray2(maxsize);
    Array<OneD, NekDouble> outarray3(maxsize);

    LibUtilities::Timer t;

    if(verbose)
    {
        cout << "Collection Implemenation for "
             << LibUtilities::ShapeTypeMap[pExp->DetShapeType()] << " ( ";
        for(int i = 0; i < pExp->GetNumBases(); ++i)
        {
            cout << pExp->GetBasis(i)->GetNumModes() <<" ";
        }
        cout << ")" <<  " for ngeoms = " << pCollExp.size() << endl;
    }
    // set  up an array of collections
    CollectionVector coll;
    
    StdRegions::ConstFactorMap factors; // required for helmholtz operator
    factors[StdRegions::eFactorLambda] = 1.5; 
    
    for(int imp = 1; imp < SIZE_ImplementationType; ++imp)
    {
        ImplementationType impType = (ImplementationType)imp;
        OperatorImpMap impTypes;
        for (int i = 0; i < SIZE_OperatorType; ++i)
        {
            OperatorType opType = (OperatorType)i;
            OperatorKey opKey(pCollExp[0]->DetShapeType(), opType, impType,
                              pCollExp[0]->IsNodalNonTensorialExp());

            if (GetOperatorFactory().ModuleExists(opKey))
            {
                impTypes[opType] = impType;
            }
        }

        Collection collLoc(pCollExp,impTypes);
        for (int i = 0; i < SIZE_OperatorType; ++i)
        {
            collLoc.Initialise((OperatorType)i, factors);
        }
        coll.push_back(collLoc);
    }

    // Determine the number of tests to do in one second
    Array<OneD, int> Ntest(SIZE_OperatorType);
    for(int i = 0; i < SIZE_OperatorType; ++i)
    {
        OperatorType OpType = (OperatorType)i;

        t.Start();

        coll[0].ApplyOperator(OpType,    inarray,
                              outarray1, outarray2,
                              outarray3);
        t.Stop();

        NekDouble oneTest = t.TimePerTest(1);

        Ntest[i] = max((int)(0.25/oneTest),1);
    }

    Array<OneD, NekDouble> timing(SIZE_ImplementationType);


    if(verbose)
    {
        cout << "\t " << "   Op.    "  << ":\t"
             << "opt. Impl."  << "\t (IterLocExp,  IterStdExp,  "
            "StdMat,     SumFac,      MatrixFree)" << endl;
    }

    // loop over all operators and determine fastest implementation
    for(int i = 0; i < SIZE_OperatorType; ++i)
    {
        OperatorType OpType = (OperatorType)i;

        // call collection implementation in thorugh ExpList.
        for (int imp = 0; imp < coll.size(); ++imp)
        {
            if (coll[imp].HasOperator(OpType))
            {
                t.Start();
                for(int n = 0; n < Ntest[i]; ++n)
                {
                    coll[imp].ApplyOperator(OpType,    inarray,
                                            outarray1, outarray2,
                                            outarray3);
                }
                t.Stop();
                timing[imp] = t.TimePerTest(Ntest[i]);
            }
            else
            {
                timing[imp] = 1000.0;
            }
        }
        // determine optimal implementation. Note +1 to
        // remove NoImplementationType flag
        int minImp = Vmath::Imin(coll.size(),timing,1)+1;

        if(verbose)
        {
            cout << "\t " << OperatorTypeMap1[i] << ": \t"
                 << ImplementationTypeMap1[minImp] << "\t (";
            for(int j = 0; j < coll.size(); ++j)
            {
                if (timing[j] > 999.0)
                {
                    cout << "     --    ";
                }
                else
                {
                    cout << timing[j] ;
                }
                if(j != coll.size()-1)
                {
                    cout <<", ";
                }
            }
            cout << ")" << endl;
        }

        // could reset global map if reusing  method?
        //m_global[OpType][pExp->DetShapeType()] = (ImplementationType)minImp;
        // set up new map
        ret[OpType] = (ImplementationType)minImp;
    }

    // store map for use by another expansion.
    m_opImpMap[OpKey] = ret;
    return ret;
}

void CollectionOptimisation::UpdateOptFile(std::string sessName,
                                           LibUtilities::CommSharedPtr &comm)
{
    std::string outname = sessName + ".opt"; 
    
    TiXmlDocument doc;
    TiXmlElement *root;
    TiXmlElement *xmlCol = new TiXmlElement("COLLECTIONS");
    GlobalOpMap  global;
    int rank = comm->GetRank();
    int nprocs = comm->GetSize(); 
    if(rank == 0)
    {
        if(!doc.LoadFile(outname)) // set up new file
        {
            TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "utf-8", "");
            doc.LinkEndChild(decl);
            root = new TiXmlElement("NEKTAR");
            doc.LinkEndChild(root);
            root->LinkEndChild(xmlCol);
        }
        else  // load file and read operator information
        {
            root   = doc.FirstChildElement("NEKTAR");
            xmlCol = root->FirstChildElement("COLLECTIONS");
            
            bool verbose = false;
            ReadCollOps(xmlCol,global,verbose);
        }
    }
        
    // update global with m_opImpMap info
    map<LibUtilities::ShapeType, int> ShapeMaxSize; 
    for(auto &opimp : m_opImpMap)
    {
        bool updateShape = true; 
        LibUtilities::ShapeType shape = opimp.first.GetShapeType();
        
        // check to see if already added this shapes details but with
        // a larger collection and if so do not update.
        if(ShapeMaxSize.count(shape))
        {
            int ngeoms = opimp.first.GetNGeoms();
            if(ngeoms > ShapeMaxSize[shape])
            {
                ShapeMaxSize[shape] = ngeoms;
            }
            else
            {
                updateShape = false;
            }
        }

        if(updateShape)
        {
            for(auto &op : opimp.second)
            {
                global[op.first][ElmtOrder(shape,-1)] = op.second;
            }
        }
    }
    
    // share 
    if(nprocs)
    {
        // loop over operators
        for(auto &op : global)
        {
            // check to see which shapes are defined in this proc
            Array<OneD, int>  ElmtImp(LibUtilities::SIZE_ShapeType,-1);
            Array<OneD, bool> ElmtDef(LibUtilities::SIZE_ShapeType,false);
            for(auto &el : op.second)
            {
                ElmtImp[el.first.first] = el.second;
                ElmtDef[el.first.first] = true; 
            }

            comm->AllReduce(ElmtImp,LibUtilities::ReduceMax);

            // loop over elements and update if not already defined
            if(rank == 0)
            {
                for(int i = 1; i < LibUtilities::SIZE_ShapeType; ++i)
                {
                    if((ElmtImp[i] != -1)&&(ElmtDef[i] == false))
                    {
                        global[op.first]
                            [ElmtOrder((LibUtilities::ShapeType) i,-1)]
                            = (ImplementationType) ElmtImp[i]; 
                    }
                }
            }
        }
    }
    
    // Update Collection section with global data on root 
    if(rank == 0)
    {
        xmlCol->Clear();
        
        map<LibUtilities::ShapeType,string>
                ShapeLetMap = {{LibUtilities::eSegment,      "S"},
                               {LibUtilities::eTriangle,     "T"},
                               {LibUtilities::eQuadrilateral,"Q"},
                               {LibUtilities::eTetrahedron,  "A"},
                               {LibUtilities::ePyramid,      "P"},
                               {LibUtilities::ePrism,        "R"},
                               {LibUtilities::eHexahedron,   "H"}};
        
        for(auto &op : global)
        {
            TiXmlElement* ColOp = new TiXmlElement("OPERATOR");
            xmlCol->LinkEndChild(ColOp);
            ColOp->SetAttribute("TYPE", OperatorTypeMap[op.first]);
            
            for(auto &el : op.second)
            {
                TiXmlElement* ElmtOp = new TiXmlElement("ELEMENT");
                ColOp->LinkEndChild(ElmtOp);
                
                ElmtOp->SetAttribute("TYPE",ShapeLetMap[el.first.first]);
                ElmtOp->SetAttribute("ORDER","*");
                ElmtOp->SetAttribute("IMPTYPE",ImplementationTypeMap[el.second]);
            }
        }
        
        doc.SaveFile(outname);
    }
}

}
}
