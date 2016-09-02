////////////////////////////////////////////////////////////////////////////////
//
//  File: InputCAD.cpp
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
//  Description: create mesh from cad using mesh utils
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

#include <boost/filesystem.hpp>

#include <NekMeshUtils/MeshElements/Element.h>

#include <NekMeshUtils/CADSystem/CADSystem.h>
#include <NekMeshUtils/Octree/Octree.h>
#include <NekMeshUtils/SurfaceMeshing/SurfaceMesh.h>
#include <NekMeshUtils/BLMeshing/BLMesh.h>
#include <NekMeshUtils/TetMeshing/TetMesh.h>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Communication/CommSerial.h>

#include "InputCAD.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey InputCAD::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eInputModule, "mcf"),
    InputCAD::create,
    "Reads CAD geometry and will generate the mesh file.");

/**
 * @brief Set up InputCAD object.
 */
InputCAD::InputCAD(MeshSharedPtr m) : InputModule(m)
{
}

InputCAD::~InputCAD()
{
}

void InputCAD::Process()
{
    vector<string> filename;
    filename.push_back(m_config["infile"].as<string>());
    string fn = filename[0].substr(0, filename[0].find("."));

    LibUtilities::SessionReaderSharedPtr pSession =
        LibUtilities::SessionReader::CreateInstance(0, NULL, filename);

    // these parameters must be defined for any mesh generation to work
    pSession->LoadParameter("MinDelta", m_minDelta);
    pSession->LoadParameter("MaxDelta", m_maxDelta);
    pSession->LoadParameter("EPS", m_eps);
    pSession->LoadParameter("Order", m_order);
    m_CADName = pSession->GetSolverInfo("CADFile");

    if (pSession->DefinesSolverInfo("MeshType"))
    {
        if (pSession->GetSolverInfo("MeshType") == "BL")
        {
            m_makeBL = true;
            pSession->LoadParameter("BLThick", m_blthick);
        }
        else
        {
            m_makeBL = false;
        }
    }
    else
    {
        m_makeBL = false;
    }

    if (pSession->DefinesSolverInfo("WriteOctree"))
    {
        m_writeoctree = pSession->GetSolverInfo("WriteOctree") == "TRUE";
    }

    CADSystemSharedPtr m_cad =
        MemoryManager<CADSystem>::AllocateSharedPtr(m_CADName);

    if (m_mesh->m_verbose)
    {
        cout << "Building mesh for: " << m_CADName << endl;
    }

    m_mesh->m_hasCAD = true;
    m_mesh->m_CADId = m_CADName;

    ASSERTL0(m_cad->LoadCAD(), "Failed to load CAD");


    vector<int> bs = m_cad->GetBoundarySurfs();

    vector<unsigned int> symsurfs;
    vector<unsigned int> blsurfs, blsurfst;
    if (m_makeBL)
    {
        string sym, bl;
        bl = pSession->GetSolverInfo("BLSurfs");
        ParseUtils::GenerateSeqVector(bl.c_str(), blsurfst);
        sort(blsurfst.begin(), blsurfst.end());
        ASSERTL0(blsurfst.size() > 0,
                 "No surfaces selected to make boundary layer on");
        for(int i = 0; i < blsurfst.size(); i++)
        {
            bool add = true;
            for(int j = 0; j < bs.size(); j++)
            {
                if(bs[j] == blsurfst[i])
                {
                    add = false;
                    break;
                }
            }
            if(add)
            {
                blsurfs.push_back(blsurfst[i]);
            }
        }
    }

    if (m_mesh->m_verbose)
    {
        cout << "With parameters:" << endl;
        cout << "\tmin delta: " << m_minDelta << endl
             << "\tmax delta: " << m_maxDelta << endl
             << "\tesp: " << m_eps << endl
             << "\torder: " << m_order << endl;
        m_cad->Report();
    }

    if (m_makeBL && m_mesh->m_verbose)
    {

        cout << "\tWill make boundary layers on surfs: ";
        for (int i = 0; i < blsurfs.size(); i++)
        {
            cout << blsurfs[i] << " ";
        }
        cout << endl << "\tWith the symmetry planes: ";
        for (int i = 0; i < symsurfs.size(); i++)
        {
            cout << symsurfs[i] << " ";
        }
        cout << endl << "\tWith thickness " << m_blthick << endl;
    }

    // create octree
    OctreeSharedPtr m_octree = MemoryManager<Octree>::AllocateSharedPtr(
        m_cad, m_mesh->m_verbose, m_minDelta, m_maxDelta, m_eps);

    if(pSession->DefinesSolverInfo("SourcePoints"))
    {
        ASSERTL0(boost::filesystem::exists(pSession->GetSolverInfo("SourcePoints").c_str()),
                 "sourcepoints file does not exist");
        vector<vector<NekDouble> > points;
        ifstream file;
        file.open(pSession->GetSolverInfo("SourcePoints").c_str());
        string line;

        while (getline(file, line))
        {
            vector<NekDouble> point(3);
            stringstream s(line);
            s >> point[0] >> point[1] >> point[2];
            points.push_back(point);
        }
        NekDouble sp;
        pSession->LoadParameter("SPSize", sp);
        m_octree->SetSourcePoints(points, sp);

    }

    if (pSession->DefinesSolverInfo("UserDefinedSpacing"))
    {
        string udsName = pSession->GetSolverInfo("UserDefinedSpacing");
        ASSERTL0(boost::filesystem::exists(udsName.c_str()),
                 "UserDefinedSpacing file does not exist");

        m_octree->SetUDSFile(udsName);
    }

    m_octree->Build();

    if (m_writeoctree)
    {
        MeshSharedPtr oct = boost::shared_ptr<Mesh>(new Mesh());
        oct->m_expDim     = 3;
        oct->m_spaceDim   = 3;
        oct->m_nummode    = 2;

        m_octree->GetOctreeMesh(oct);

        ModuleSharedPtr mod = GetModuleFactory().CreateInstance(
            ModuleKey(eOutputModule, "xml"), oct);
        mod->RegisterConfig("outfile", fn + "_oct.xml");
        mod->ProcessVertices();
        mod->ProcessEdges();
        mod->ProcessFaces();
        mod->ProcessElements();
        mod->ProcessComposites();
        mod->Process();
    }

    m_mesh->m_expDim   = 3;
    m_mesh->m_spaceDim = 3;
    m_mesh->m_nummode = m_order + 1;
    if (m_makeBL)
    {
        m_mesh->m_numcomp = 2; // prisms and tets
    }
    else
    {
        m_mesh->m_numcomp = 1; // just tets
    }

    //create surface mesh
    m_mesh->m_expDim--; //just to make it easier to surface mesh for now
    SurfaceMeshSharedPtr m_surfacemesh = MemoryManager<SurfaceMesh>::
                AllocateSharedPtr(m_mesh, m_cad, m_octree);

    m_surfacemesh->Mesh();

    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    m_surfacemesh->Report();

    EdgeSet::iterator eit;
    int count = 0;
    for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
    {
        if((*eit)->m_elLink.size() != 2)
        {
            count++;
        }
    }

    //curve cad taging is working and tested to this point

    if (count > 0)
    {
        cerr << "Error: mesh contains unconnected edges and is not valid"
             << endl;
        abort();
    }

    map<int, FaceSharedPtr> surftopriface;
    // map of surface element id to opposite prism
    // face for psudo surface in tetmesh

    TetMeshSharedPtr m_tet;
    if (m_makeBL)
    {
        BLMeshSharedPtr m_blmesh = MemoryManager<BLMesh>::AllocateSharedPtr(
                                        m_cad, m_mesh, blsurfs, m_blthick);

        m_blmesh->Mesh();

        //m_mesh->m_element[2].clear();
        //m_mesh->m_expDim = 3;
        /*ClearElementLinks();
        ProcessVertices();
        ProcessEdges();
        ProcessFaces();
        ProcessElements();
        ProcessComposites();
        return;*/

        m_surfacemesh->Remesh(m_blmesh);

        /*m_mesh->m_nummode = 2;
        m_mesh->m_expDim = 3;
        m_mesh->m_element[2].clear();
        ClearElementLinks();
        ProcessVertices();
        ProcessEdges();
        ProcessFaces();
        ProcessElements();
        ProcessComposites();
        return;*/

        //create tet mesh
        m_mesh->m_expDim = 3;

        m_tet = MemoryManager<TetMesh>::AllocateSharedPtr(
            m_mesh, m_octree, m_blmesh);
    }
    else
    {
        m_mesh->m_expDim = 3;
        m_tet = MemoryManager<TetMesh>::AllocateSharedPtr(m_mesh, m_octree);

    }

    m_tet->Mesh();

    //m_mesh->m_element[2].clear();

    ClearElementLinks();
    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    FaceSet::iterator fit;
    count = 0;
    for(fit = m_mesh->m_faceSet.begin(); fit != m_mesh->m_faceSet.end(); fit++)
    {
        if((*fit)->m_elLink.size() != 2)
        {
            count++;
        }
    }

    if (count - m_mesh->m_element[2].size() > 0)
    {
        cerr << "Error: mesh contains unconnected faces and is not valid "
             << count - m_mesh->m_element[2].size()
             << endl;
        abort();
    }

    //return;

    m_surfacemesh->HOSurf();

    if (m_mesh->m_verbose)
    {
        cout << endl;
        cout << m_mesh->m_element[3].size() << endl;
    }

}
}
}
