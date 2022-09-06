////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessCombine.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
//  Description: Combine two mesh files into one output file.
//
////////////////////////////////////////////////////////////////////////////////

#include "ProcessCombine.h"
#include <LibUtilities/BasicUtils/Timer.h>
#include <NekMesh/MeshElements/Element.h>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

using namespace std;
using namespace Nektar::NekMesh;

namespace Nektar
{
namespace NekMesh
{

ModuleKey ProcessCombine::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "combine"), ProcessCombine::create,
        "Combine two mesh files into one output file.");

ProcessCombine::ProcessCombine(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["file"] = ConfigOption(
        false, "", "Second mesh file to be combined with the input mesh file.");
}

ProcessCombine::~ProcessCombine()
{
}

void ProcessCombine::Process()
{
    ASSERTL0(m_config["file"].beenSet,
             "The 'file' module parameter must be set.")

    std::string fname = m_config["file"].as<std::string>();

    // Check to see if filename exists.
    if (!boost::filesystem::exists(fname))
    {
        m_log(FATAL) << "Unable to read file: '" << fname << "'" << std::endl;
    }

    // Process the second mesh input
    ModuleKey module;
    module.first = eInputModule;

    vector<string> tmp1;
    boost::split(tmp1, fname, boost::is_any_of(":"));
    if (tmp1.size() == 1)
    {
        int dot    = tmp1[0].find_last_of('.') + 1;
        string ext = tmp1[0].substr(dot, tmp1[0].length() - dot);

        if (ext == "gz")
        {
            string tmp = tmp1[0].substr(0, tmp1[0].length() - 3);
            dot        = tmp.find_last_of('.') + 1;
            ext        = tmp.substr(dot, tmp.length() - dot);
        }

        module.second = ext;
        tmp1.push_back("infile=" + tmp1[0]);
    }
    else
    {
        module.second = tmp1[1];
        tmp1.push_back("infile=" + tmp1[0]);
    }

    // Copy communicator
    MeshSharedPtr mesh = std::make_shared<Mesh>();
    mesh->m_comm       = m_mesh->m_comm;

    ModuleSharedPtr mod = GetModuleFactory().CreateInstance(module, mesh);
    mod->SetLogger(m_log);
    mod->GetLogger().SetPrefix("ProcessCombine");

    // Set options for this module.
    for (int j = 1; j < tmp1.size(); ++j)
    {
        vector<string> tmp2;
        boost::split(tmp2, tmp1[j], boost::is_any_of("="));

        if (tmp2.size() == 1)
        {
            mod->RegisterConfig(tmp2[0]);
        }
        else if (tmp2.size() == 2)
        {
            mod->RegisterConfig(tmp2[0], tmp2[1]);
        }
        else
        {
            NEKERROR(ErrorUtil::efatal, "ERROR: Invalid module configuration: "
                                        "format is either :arg or :arg=val");
        }
    }

    // Ensure configuration options have been set.
    mod->SetDefaults();

    // Run second input module
    mod->Process();
    MeshSharedPtr mesh2 = mod->GetMesh();

    // Check dimensions of both meshes are the same
    ASSERTL0(m_mesh->m_expDim == mesh2->m_expDim,
             "The expansion dimensions of the meshes being combined must be "
             "the same.")

    // Renumber vertices and copy in to m_vertexSet
    int vid = m_mesh->m_vertexSet.size();
    for (auto &elmt : mesh2->m_element[m_mesh->m_expDim])
    {
        for (int j = 0; j < elmt->GetVertexCount(); ++j)
        {
            pair<NodeSet::iterator, bool> testIns =
                m_mesh->m_vertexSet.insert(elmt->GetVertex(j));

            if (testIns.second)
            {
                (*testIns.first)->m_id = vid++;
            }

            elmt->SetVertex(j, *testIns.first);
        }
    }

    // Renumber edges and copy in to m_edgeSet
    int eid = m_mesh->m_edgeSet.size();
    for (auto &elmt : mesh2->m_element[m_mesh->m_expDim])
    {
        for (int j = 0; j < elmt->GetEdgeCount(); ++j)
        {
            EdgeSharedPtr ed = elmt->GetEdge(j);
            pair<EdgeSet::iterator, bool> testIns =
                m_mesh->m_edgeSet.insert(ed);

            if (testIns.second)
            {
                EdgeSharedPtr ed2 = *testIns.first;
                ed2->m_id         = eid++;
                // ed2->m_elLink.push_back(
                //        pair<ElementSharedPtr, int>(elmt, j));
            }
            else
            {
                EdgeSharedPtr e2 = *(testIns.first);
                elmt->SetEdge(j, e2);
                if (e2->m_edgeNodes.size() == 0 && ed->m_edgeNodes.size() > 0)
                {
                    e2->m_curveType = ed->m_curveType;
                    e2->m_edgeNodes = ed->m_edgeNodes;

                    // Reverse nodes if appropriate.
                    if (e2->m_n1->m_id != ed->m_n1->m_id)
                    {
                        reverse(e2->m_edgeNodes.begin(), e2->m_edgeNodes.end());
                    }
                }

                // if (ed->m_parentCAD)
                //{
                //    e2->m_parentCAD = ed->m_parentCAD;
                //}

                // Update edge to element map.
                // e2->m_elLink.push_back(
                //        pair<ElementSharedPtr, int>(elmt, j));
            }

            elmt->SetEdge(j, *testIns.first);
        }
    }

    // @TODO: Create links for 1D elements?

    // Renumber faces and copy in to m_faceSet
    int fid = m_mesh->m_faceSet.size();
    for (auto &elmt : mesh2->m_element[m_mesh->m_expDim])
    {
        for (int j = 0; j < elmt->GetFaceCount(); ++j)
        {
            pair<FaceSet::iterator, bool> testIns =
                m_mesh->m_faceSet.insert(elmt->GetFace(j));

            if (testIns.second)
            {
                (*(testIns.first))->m_id = fid++;
                // Update face to element map.
                //(*(testIns.first))->m_elLink.push_back(
                //        pair<ElementSharedPtr,int>(elmt,j));
            }
            else
            {
                // Update face to element map.
                //(*(testIns.first))->m_elLink.push_back(
                //        pair<ElementSharedPtr,int>(elmt,j));
            }

            elmt->SetFace(j, *testIns.first);
        }
    }

    // @TODO: Create links for 2D elements?

    // Renumber elements and copy in to m_element
    auto &elmts = mesh2->m_element;
    for (int d = 0; d < 4; ++d)
    {
        int numEl  = m_mesh->m_element[d].size();
        auto elVec = elmts[d];
        for (auto &el : elVec)
        {
            el->SetId(numEl++);
            m_mesh->m_element[d].emplace_back(el);
        }
    }

    // Renumber composites and copy in to m_composite
    std::map<int, int> compRenumber;
    for (int d = 0; d <= m_mesh->m_expDim; ++d)
    {
        vector<ElementSharedPtr> &elmt = m_mesh->m_element[d];
        for (int i = elmt.size() - mesh2->m_element[d].size(); i < elmt.size();
             ++i) // for number of elements added by mesh 2
        {
            CompositeMap::iterator it;
            unsigned int tagid = elmt[i]->GetTagList()[0];

            auto findKey = compRenumber.find(elmt[i]->GetTagList()[0]);
            if (findKey != compRenumber.end())
            {
                tagid = findKey->second;
            }
            else
            {
                // Checks if composite tag is already defined in mesh 1
                while (m_mesh->m_composite.find(tagid) !=
                       m_mesh->m_composite.end())
                {
                    tagid++;
                    // Checks if trying to replace the tag with a tag already
                    // defined in mesh 2
                    if (mesh2->m_composite.find(tagid) !=
                        mesh2->m_composite.end())
                    {
                        tagid++;
                    }
                }

                compRenumber[elmt[i]->GetTagList()[0]] = tagid;
            }

            it = m_mesh->m_composite.find(tagid);
            if (it == m_mesh->m_composite.end())
            {
                CompositeSharedPtr tmp =
                    std::shared_ptr<Composite>(new Composite());
                pair<CompositeMap::iterator, bool> testIns;
                tmp->m_id  = tagid;
                tmp->m_tag = elmt[i]->GetTag();
                if (mesh2->m_faceLabels.count(tmp->m_id) != 0)
                {
                    tmp->m_label = mesh2->m_faceLabels[tmp->m_id];
                }

                testIns = m_mesh->m_composite.insert(
                    pair<unsigned int, CompositeSharedPtr>(tagid, tmp));
                it = testIns.first;
            }

            elmt[i]->GetTagList()[0] = tagid;
            it->second->m_items.push_back(elmt[i]);
        }
    }

    if (!compRenumber.empty())
    {
        m_log << "Duplicate composite IDs from mesh 1 detected in mesh 2."
              << endl;
        m_log << "These will be remapped in the output file to:" << endl;

        for (auto &cIt : compRenumber)
        {
            if (cIt.first != cIt.second)
            {
                m_log << "- C[" << cIt.first << "] => "
                      << "C[" << cIt.second << "]" << endl;
            }
        }
    }
}
} // namespace NekMesh
} // namespace Nektar
