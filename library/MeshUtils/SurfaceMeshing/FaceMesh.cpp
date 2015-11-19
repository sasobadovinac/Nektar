////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMesh.cpp
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
//  Description: surfacemesh object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <limits>
#include <MeshUtils/SurfaceMeshing/FaceMesh.h>
#include <MeshUtils/ExtLibInterface/TriangleInterface.h>

#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;
namespace Nektar
{
namespace MeshUtils
{

void FaceMesh::Mesh()
{
    Stretching();

    OrientateCurves();

    int numPoints = 0;
    for(int i = 0; i < orderedLoops.size(); i++)
    {
        numPoints+=orderedLoops[i].size();
    }

    stringstream ss;
    ss << "3 points required for triangulation, " << numPoints << " in loop" << endl;
    ss << "curves: ";
    for(int i = 0; i < m_edgeloops.size(); i++)
    {
        for(int j = 0; j < m_edgeloops[i].edges.size(); j++)
        {
            ss << m_edgeloops[i].edges[j]->GetID() << " ";
        }
    }

    ASSERTL0(numPoints > 2, ss.str());

    //create interface to triangle thirdparty library
    TriangleInterfaceSharedPtr pplanemesh =
        MemoryManager<TriangleInterface>::AllocateSharedPtr();

    vector<Array<OneD, NekDouble> > centers;
    for(int i = 0; i < m_edgeloops.size(); i++)
    {
        centers.push_back(m_edgeloops[i].center);
    }

    pplanemesh->Assign(orderedLoops, centers, m_id, asr/pasr);

    pplanemesh->Mesh();

    pplanemesh->Extract(m_connec);

    bool repeat = true;
    int meshcounter = 1;

    while (repeat)
    {
        repeat = Validate();
        if(!repeat)
        {
            break;
        }
        m_connec.clear();
        pplanemesh->AssignStiener(m_stienerpoints);
        pplanemesh->Mesh();
        pplanemesh->Extract(m_connec);
        meshcounter++;
    }

    BuildLocalMesh();

    OptimiseLocalMesh();

    //clear local element links
    EdgeSet::iterator eit;
    for(eit = m_localEdges.begin(); eit != m_localEdges.end(); eit++)
    {
        (*eit)->m_elLink.clear();
    }

    //make new elements and add to list from list of nodes and connectivity from triangle
    for(int i = 0; i < m_localElements.size(); i++)
    {
        vector<EdgeSharedPtr> e  = m_localElements[i]->GetEdgeList();
        for(int j = 0; j < e.size(); j++)
        {
            e[j]->m_elLink.clear();
        }
        m_mesh->m_element[m_mesh->m_expDim].push_back(m_localElements[i]);
    }

    cout << "\r                                                                                             ";
    cout << scientific << "\r\t\tFace " << m_id << endl
         << "\t\t\tNodes: " << m_localNodes.size() << endl
         << "\t\t\tEdges: " << m_localEdges.size() << endl
         << "\t\t\tTriangles: " << m_localElements.size() << endl
         << endl;
}

void FaceMesh::OptimiseLocalMesh()
{
    DiagonalSwap();

    Smoothing();

    DiagonalSwap();

    Smoothing();
}

void FaceMesh::Smoothing()
{
    EdgeSet::iterator eit;
    NodeSet::iterator nit;

    map<int, vector<EdgeSharedPtr> > connectingedges;

    map<int, vector<ElementSharedPtr> > connectingelements;

    for(eit = m_localEdges.begin(); eit != m_localEdges.end(); eit++)
    {
        connectingedges[(*eit)->m_n1->m_id].push_back(*eit);
        connectingedges[(*eit)->m_n2->m_id].push_back(*eit);
    }

    for(int i = 0; i < m_localElements.size(); i++)
    {
        vector<NodeSharedPtr> v = m_localElements[i]->GetVertexList();
        for(int j = 0; j < 3; j++)
        {
            connectingelements[v[j]->m_id].push_back(m_localElements[i]);
        }
    }

    //perform 8 runs of elastic relaxation based on the octree
    for(int q = 0; q < 4; q++)
    {
        if(m_mesh->m_verbose)
            cout << "\t\t Elastic relaxation run: " << q+1 << endl;

        for(nit = m_localNodes.begin(); nit != m_localNodes.end(); nit++)
        {
            vector<int> c = (*nit)->GetListCADCurve();
            if(c.size()>0) //node is on curve so skip
                continue;

            vector<NodeSharedPtr> connodes; //this can be real nodes or dummy nodes depending on the system

            vector<EdgeSharedPtr> edges = connectingedges[(*nit)->m_id];
            vector<ElementSharedPtr> els = connectingelements[(*nit)->m_id];

            vector<NodeSharedPtr> nodesystem;
            vector<NekDouble> lamp;

            for(int i = 0; i < edges.size(); i++)
            {
                vector<NekDouble> lambda;

                NodeSharedPtr J;
                if(*nit == edges[i]->m_n1)
                    J = edges[i]->m_n2;
                else if(*nit == edges[i]->m_n2)
                    J = edges[i]->m_n1;
                else
                    ASSERTL0(false,"could not find node");

                Array<OneD, NekDouble> ui = (*nit)->GetCADSurf(m_id);
                Array<OneD, NekDouble> uj = J->GetCADSurf(m_id);

                for(int j = 0; j < els.size(); j++)
                {
                    vector<NodeSharedPtr> v = els[j]->GetVertexList();
                    if(v[0] == J || v[1] == J || v[2] == J)
                        continue; //elememt is adjacent to J therefore no intersection on IJ

                    //need to find other edge
                    EdgeSharedPtr LK;
                    vector<EdgeSharedPtr> es = els[j]->GetEdgeList();
                    if(!(es[0]->m_n1 == *nit || es[0]->m_n2 == *nit))
                    {
                        LK = es[0];
                    }
                    else
                    {
                        if(!(es[1]->m_n1 == *nit || es[1]->m_n2 == *nit))
                        {
                            LK = es[1];
                        }
                        else
                        {
                            if(!(es[2]->m_n1 == *nit || es[2]->m_n2 == *nit))
                            {
                                LK = es[2];
                            }
                            else
                            {
                                ASSERTL0(false,"failed to find edge");
                            }
                        }
                    }

                    Array<OneD, NekDouble> uk = LK->m_n1->GetCADSurf(m_id);
                    Array<OneD, NekDouble> ul = LK->m_n2->GetCADSurf(m_id);



                    Array<OneD, NekDouble> n(2);
                    n[0] = -1.0*(uk[1] - ul[1]);
                    n[1] = uk[0] - ul[0];
                    n[0] = n[0] / sqrt(n[0]*n[0] + n[1]*n[1]);
                    n[1] = n[1] / sqrt(n[0]*n[0] + n[1]*n[1]);

                    NekDouble lam = -1.0*(ui[0]*uk[0]+ui[1]*uk[1]) /
                                    ((uj[0]-ui[0])*n[0] +
                                     (uj[1]-ui[1])*n[1]);
                    if(!(lam < 0) && !(lam > 1))
                        lambda.push_back(lam);
                }

                if(lambda.size() > 0)
                {
                    sort(lambda.begin(),lambda.end());
                    //make a new dummy node based on the system
                    Array<OneD, NekDouble> ud(2);
                    ud[0] = ui[0] + lambda[0]*(uj[0] - ui[0]);
                    ud[1] = ui[1] + lambda[0]*(uj[1] - ui[1]);
                    Array<OneD, NekDouble> locd = m_cadsurf->P(ud);
                    NodeSharedPtr dn = boost::shared_ptr<Node>(new Node(0,locd[0],locd[1],locd[2]));
                    dn->SetCADSurf(m_id,ud);
                    nodesystem.push_back(dn);

                    lamp.push_back(lambda[0]);
                }
                else
                {
                    nodesystem.push_back(J);
                    lamp.push_back(1.0);
                }
            }

            //we want to move the node to the centroid of the system,
            //then perfrom one newton interation to move the node according to deltaspec
            Array<OneD, NekDouble> ui(2);
            ui[0]=0.0; ui[1]=0.0;

            DNekMat f(2,1,0.0);
            DNekMat J(2,2,0.0);
            for(int i = 0; i < nodesystem.size(); i++)
            {
                Array<OneD, NekDouble> uj = nodesystem[i]->GetCADSurf(m_id);
                ui[0]+=uj[0]/nodesystem.size();
                ui[1]+=uj[1]/nodesystem.size();
            }

            Array<OneD, NekDouble> Xi = m_cadsurf->P(ui);
            Array<OneD, NekDouble> ri = m_cadsurf->D1(ui);

            for(int i = 0; i < nodesystem.size();i++)
            {
                Array<OneD, NekDouble> uj = nodesystem[i]->GetCADSurf(m_id);
                Array<OneD, NekDouble> Xj = nodesystem[i]->GetLoc();

                NekDouble ljstar = m_octree->Query(Xj);
                NekDouble lj = sqrt((Xi[0]-Xj[0])*(Xi[0]-Xj[0]) +
                                    (Xi[1]-Xj[1])*(Xi[1]-Xj[1]) +
                                    (Xi[2]-Xj[2])*(Xi[2]-Xj[2]));

                NekDouble Lj = sqrt((ui[0]-uj[0])*(ui[0]-uj[0]) +
                                    (ui[1]-uj[1])*(ui[1]-uj[1]));

                f(0,0) += (lj - ljstar)*(ui[0] - uj[0])/Lj;
                f(1,0) += (lj - ljstar)*(ui[1] - uj[1])/Lj;

                J(0,0) += (ri[3]*(Xi[0]-Xj[0]) + ri[4]*(Xi[1]-Xj[1]) + ri[5]*(Xi[2]-Xj[2]))/lj/Lj*(ui[0]-uj[0]) + (lj-ljstar)/Lj/Lj*(Lj - (ui[0]-uj[0])*(ui[0]-uj[0])/Lj);

                J(1,1) += (ri[6]*(Xi[0]-Xj[0]) + ri[7]*(Xi[1]-Xj[1]) + ri[8]*(Xi[2]-Xj[2]))/lj/Lj*(ui[1]-uj[1]) + (lj-ljstar)/Lj/Lj*(Lj - (ui[1]-uj[1])*(ui[1]-uj[1])/Lj);

                J(1,0) += (ri[3]*(Xi[0]-Xj[0]) + ri[4]*(Xi[1]-Xj[1]) + ri[5]*(Xi[2]-Xj[2]))/lj/Lj*(ui[1]-uj[1]) - (lj-ljstar)/Lj/Lj*(ui[0]-uj[0])*(ui[1]-uj[1])/Lj;

                J(0,1) += (ri[6]*(Xi[0]-Xj[0]) + ri[7]*(Xi[1]-Xj[1]) + ri[8]*(Xi[2]-Xj[2]))/lj/Lj*(ui[0]-uj[0]) - (lj-ljstar)/Lj/Lj*(ui[0]-uj[0])*(ui[1]-uj[1])/Lj;
            }

            NekDouble fmagbefore = sqrt(f(0,0)*f(0,0)+f(1,0)*f(1,0));

            J.Invert();

            Array<OneD, NekDouble> bounds = m_cadsurf->GetBounds();

            DNekMat U = J*f;
            Array<OneD, NekDouble> uvn(2);

            uvn[0] = ui[0]; //- U(0,0);
            uvn[1] = ui[1]; //- U(1,0);

            if(!(uvn[0] < bounds[0] ||
                       uvn[0] > bounds[1] ||
                       uvn[1] < bounds[2] ||
                       uvn[1] > bounds[3]))
            {

                /*f(0,0) = 0; f(1,0) = 0;
                Xi = m_cadsurf->P(uvn);
                for(int i = 0; i < nodesystem.size();i++)
                {
                    Array<OneD, NekDouble> uj = nodesystem[i]->GetCADSurf(m_id);
                    Array<OneD, NekDouble> Xj = nodesystem[i]->GetLoc();

                    NekDouble ljstar = m_octree->Query(Xj);
                    NekDouble lj = sqrt((Xi[0]-Xj[0])*(Xi[0]-Xj[0]) +
                                        (Xi[1]-Xj[1])*(Xi[1]-Xj[1]) +
                                        (Xi[2]-Xj[2])*(Xi[2]-Xj[2]));

                    NekDouble Lj = sqrt((uvn[0]-uj[0])*(uvn[0]-uj[0]) +
                                        (uvn[1]-uj[1])*(uvn[1]-uj[1]));

                    f(0,0) += (lj - ljstar)*(uvn[0] - uj[0])/Lj;
                    f(1,0) += (lj - ljstar)*(uvn[1] - uj[1])/Lj;
                }

                NekDouble fmagafter = sqrt(f(0,0)*f(0,0)+f(1,0)*f(1,0));

                if(fmagafter > fmagbefore)
                    cout << "not optimsing" << endl;*/

                Array<OneD, NekDouble> l2 = m_cadsurf->P(uvn);
                (*nit)->Move(l2,m_id,uvn);
            }
        }
    }
}

void FaceMesh::DiagonalSwap()
{
    map<int, int> idealConnec;
    map<int, int> actualConnec;
    map<int, vector<EdgeSharedPtr> > nodetoedge;
    //figure out ideal node count and actual node count
    EdgeSet::iterator eit;
    for(eit = m_localEdges.begin(); eit != m_localEdges.end(); eit++)
    {
        nodetoedge[(*eit)->m_n1->m_id].push_back(*eit);
        nodetoedge[(*eit)->m_n2->m_id].push_back(*eit);
    }
    NodeSet::iterator nit;
    for(nit = m_localNodes.begin(); nit != m_localNodes.end(); nit++)
    {
        if((*nit)->GetListCADCurve().size() == 0)
        {
            //node is interior
            idealConnec[(*nit)->m_id] = 6;
        }
        else
        {
            //need to identify the two other nodes on the boundary to find
            //interior angle
            vector<NodeSharedPtr> ns;
            vector<EdgeSharedPtr> e = nodetoedge[(*nit)->m_id];
            for(int i = 0; i < e.size(); i++)
            {
                if(e[i]->CADCurveID == -1)
                    continue; //the linking nodes are not going to exist on interior edges

                if(e[i]->m_n1 == (*nit))
                    ns.push_back(e[i]->m_n2);
                else
                    ns.push_back(e[i]->m_n1);
            }
            ASSERTL0(ns.size() == 2,"failed to find 2 nodes in the angle system");

            idealConnec[(*nit)->m_id] = ceil((*nit)->Angle(ns[0],ns[1])/3.142*3) + 1;
        }
    }
    for(nit = m_localNodes.begin(); nit != m_localNodes.end(); nit++)
    {
        actualConnec[(*nit)->m_id] = nodetoedge[(*nit)->m_id].size();
    }

    //edgeswapping fun times
    //perfrom edge swap based on node defect and then angle
    for(int q = 0; q < 4; q++)
    {
        if(m_mesh->m_verbose)
        {
            cout << "\t\t Edge swap ";
            if(q<2)
            {
                cout << "defect run: " << q+1;
            }
            else
            {
                cout << "angle run: " << q+1-2;
            }
        }

        int edgesStart = m_localEdges.size();
        EdgeSet edges = m_localEdges;
        m_localEdges.clear();

        int swappedEdges = 0;

        EdgeSet::iterator it;

        for(it = edges.begin(); it != edges.end(); it++)
        {
            EdgeSharedPtr e = *it;

            if(e->m_elLink.size() != 2)
            {
                m_localEdges.insert(e);
                continue;
            }

            ElementSharedPtr tri1 = e->m_elLink[0].first;
            ElementSharedPtr tri2 = e->m_elLink[1].first;

            NodeSharedPtr n1 = e->m_n1;
            NodeSharedPtr n2 = e->m_n2;

            vector<NodeSharedPtr> nt = tri1->GetVertexList();

            //identify node a,b,c,d of the swapping
            NodeSharedPtr A,B,C,D;
            if(nt[0] != n1 && nt[0] != n2)
            {
                C = nt[0];
                B = nt[1];
                A = nt[2];
            }
            else if(nt[1] != n1 && nt[1] != n2)
            {
                C = nt[1];
                B = nt[2];
                A = nt[0];
            }
            else if(nt[2] != n1 && nt[2] != n2)
            {
                C = nt[2];
                B = nt[0];
                A = nt[1];
            }
            else
            {
                ASSERTL0(false,"failed to identify verticies in tri1");
            }

            nt = tri2->GetVertexList();

            if(nt[0] != n1 && nt[0] != n2)
            {
                D = nt[0];
            }
            else if(nt[1] != n1 && nt[1] != n2)
            {
                D = nt[1];
            }
            else if(nt[2] != n1 && nt[2] != n2)
            {
                D = nt[2];
            }
            else
            {
                ASSERTL0(false,"failed to identify verticies in tri2");
            }

            //determine signed area of alternate config
            Array<OneD, NekDouble> ai,bi,ci,di;
            ai = A->GetCADSurf(m_id);
            bi = B->GetCADSurf(m_id);
            ci = C->GetCADSurf(m_id);
            di = D->GetCADSurf(m_id);

            NekDouble CDA, CBD;

            CDA = 0.5*(-di[0]*ci[1] + ai[0]*ci[1] + ci[0]*di[1] - ai[0]*di[1] -
                            ci[0]*ai[1] + di[0]*ai[1]);

            CBD = 0.5*(-bi[0]*ci[1] + di[0]*ci[1] + ci[0]*bi[1] - di[0]*bi[1] -
                            ci[0]*di[1] + bi[0]*di[1]);

            //if signed area of the swapping triangles is less than zero
            //that configuration is invalid and swap cannot be performed
            if(!(CDA > 0.001 && CBD > 0.001))
            {
                m_localEdges.insert(e);
                continue;
            }

            bool swap = false; //assume do not swap

            if(q<2)
            {
                int nodedefectbefore = 0;
                nodedefectbefore += abs(actualConnec[A->m_id] - idealConnec[A->m_id]);
                nodedefectbefore += abs(actualConnec[B->m_id] - idealConnec[B->m_id]);
                nodedefectbefore += abs(actualConnec[C->m_id] - idealConnec[C->m_id]);
                nodedefectbefore += abs(actualConnec[D->m_id] - idealConnec[D->m_id]);

                int nodedefectafter = 0;
                nodedefectafter += abs(actualConnec[A->m_id] -1 - idealConnec[A->m_id]);
                nodedefectafter += abs(actualConnec[B->m_id] -1 - idealConnec[B->m_id]);
                nodedefectafter += abs(actualConnec[C->m_id] +1 - idealConnec[C->m_id]);
                nodedefectafter += abs(actualConnec[D->m_id] +1 - idealConnec[D->m_id]);

                if(nodedefectafter < nodedefectbefore)
                {
                    swap = true;
                }

            }
            else
            {
                NekDouble minanglebefore = C->Angle(A,B);
                minanglebefore = min(minanglebefore, A->Angle(B,C));
                minanglebefore = min(minanglebefore, B->Angle(A,C));
                minanglebefore = min(minanglebefore, B->Angle(A,D));
                minanglebefore = min(minanglebefore, A->Angle(B,D));
                minanglebefore = min(minanglebefore, D->Angle(A,B));

                NekDouble minangleafter = C->Angle(B,D);
                minangleafter = min(minangleafter, D->Angle(B,C));
                minangleafter = min(minangleafter, B->Angle(C,D));
                minangleafter = min(minangleafter, C->Angle(A,D));
                minangleafter = min(minangleafter, A->Angle(C,D));
                minangleafter = min(minangleafter, D->Angle(A,C));

                if(minangleafter > minanglebefore)
                {
                    swap = true;
                }
            }

            if(swap)
            {
                actualConnec[A->m_id]--;
                actualConnec[B->m_id]--;
                actualConnec[C->m_id]++;
                actualConnec[D->m_id]++;

                //make the 4 other edges
                EdgeSharedPtr CA, AD, DB, BC, CAt, ADt, DBt, BCt;
                CAt = boost::shared_ptr<Edge>(new Edge(C,A));
                ADt = boost::shared_ptr<Edge>(new Edge(A,D));
                DBt = boost::shared_ptr<Edge>(new Edge(D,B));
                BCt = boost::shared_ptr<Edge>(new Edge(B,C));

                vector<EdgeSharedPtr> es = tri1->GetEdgeList();
                for(int i = 0; i < 3; i++)
                {
                    if(es[i] == CAt)
                    {
                        CA = es[i];
                    }
                    if(es[i] == BCt)
                    {
                        BC = es[i];
                    }
                }
                es = tri2->GetEdgeList();
                for(int i = 0; i < 3; i++)
                {
                    if(es[i] == DBt)
                    {
                        DB = es[i];
                    }
                    if(es[i] == ADt)
                    {
                        AD = es[i];
                    }
                }

                //now sort out links for the 4 edges surrounding the patch
                vector<pair<ElementSharedPtr, int> > links;

                links = CA->m_elLink;
                CA->m_elLink.clear();
                for(int i = 0; i < links.size(); i++)
                {
                    if(links[i].first->GetId() == tri1->GetId())
                        continue;
                    CA->m_elLink.push_back(links[i]);
                }

                links = BC->m_elLink;
                BC->m_elLink.clear();
                for(int i = 0; i < links.size(); i++)
                {
                    if(links[i].first->GetId() == tri1->GetId())
                        continue;
                    BC->m_elLink.push_back(links[i]);
                }

                links = AD->m_elLink;
                AD->m_elLink.clear();
                for(int i = 0; i < links.size(); i++)
                {
                    if(links[i].first->GetId() == tri2->GetId())
                        continue;
                    AD->m_elLink.push_back(links[i]);
                }

                links = DB->m_elLink;
                DB->m_elLink.clear();
                for(int i = 0; i < links.size(); i++)
                {
                    if(links[i].first->GetId() == tri2->GetId())
                        continue;
                    DB->m_elLink.push_back(links[i]);
                }

                EdgeSharedPtr newe = boost::shared_ptr<Edge>(new Edge(C,D));

                vector<NodeSharedPtr> t1,t2;
                t1.push_back(B); t1.push_back(D); t1.push_back(C);
                t2.push_back(A); t2.push_back(C); t2.push_back(D);

                ElmtConfig conf(LibUtilities::eTriangle,1,false,false);
                vector<int> tags;
                tags.push_back(m_id);

                int id1 = tri1->GetId();
                int id2 = tri2->GetId();

                ElementSharedPtr ntri1 = GetElementFactory().
                            CreateInstance(LibUtilities::eTriangle,
                                           conf,t1,tags);
                ElementSharedPtr ntri2 = GetElementFactory().
                            CreateInstance(LibUtilities::eTriangle,
                                           conf,t2,tags);

                ntri1->SetId(id1);
                ntri2->SetId(id2);

                vector<EdgeSharedPtr> t1es = ntri1->GetEdgeList();
                for(int i = 0; i < 3; i++)
                {
                    if(t1es[i] == DB)
                    {
                        ntri1->SetEdge(i,DB);
                        DB->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri1,i));
                    }
                    else if(t1es[i] == BC)
                    {
                        ntri1->SetEdge(i,BC);
                        BC->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri1,i));
                    }
                    else if(t1es[i] == newe)
                    {
                        ntri1->SetEdge(i,newe);
                        newe->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri1,i));
                    }
                    else
                    {
                        ASSERTL0(false,"weird edge in new tri 1");
                    }
                }
                vector<EdgeSharedPtr> t2es = ntri2->GetEdgeList();
                for(int i = 0; i < 3; i++)
                {
                    if(t2es[i] == CA)
                    {
                        ntri2->SetEdge(i,CA);
                        CA->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri2,i));
                    }
                    else if(t2es[i] == AD)
                    {
                        ntri2->SetEdge(i,AD);
                        AD->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri2,i));
                    }
                    else if(t2es[i] == newe)
                    {
                        ntri2->SetEdge(i,newe);
                        newe->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri2,i));
                    }
                    else
                    {
                        ASSERTL0(false,"weird edge in new tri 2");
                    }
                }

                newe->CADSurfID.push_back(m_id);
                m_localEdges.insert(newe);

                m_localElements[id1] = ntri1;
                m_localElements[id2] = ntri2;

                swappedEdges++;

            }
            else
            {
                m_localEdges.insert(e);
            }
        }

        ASSERTL0(m_localEdges.size() == edgesStart, "mismatch edge count");

        if(m_mesh->m_verbose)
            cout << ".\tEdges swapped: " << swappedEdges << endl;
    }
}

void FaceMesh::BuildLocalMesh()
{
    /*************************
    // build a local set of nodes edges and elemenets for optimstaion prior to putting them into m_mesh
    */

    for(int i = 0; i < m_connec.size(); i++)
    {
        ElmtConfig conf(LibUtilities::eTriangle, 1, false, false);

        vector<int> tags;
        tags.push_back(m_id);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
                                LibUtilities::eTriangle, conf, m_connec[i], tags);

        vector<NodeSharedPtr> nods = E->GetVertexList();
        for(int j = 0; j < nods.size(); j++)
        {
            //nodes are already unique some will insert some wont
            m_localNodes.insert(nods[j]);
        }
        vector<EdgeSharedPtr> edgs = E->GetEdgeList();
        for(int j = 0; j < edgs.size(); j++)
        {
            //look for edge in m_mesh edgeset from curves
            EdgeSet::iterator s = m_mesh->m_edgeSet.find(edgs[j]);
            if(!(s == m_mesh->m_edgeSet.end()))
            {
                edgs[j] = *s;
                E->SetEdge(j, edgs[j]);
            }

            pair<EdgeSet::iterator, bool> test = m_localEdges.insert(edgs[j]);

            if(test.second)
            {
                (*test.first)->m_elLink.push_back(pair<ElementSharedPtr,int>(E,j));
                (*test.first)->CADSurfID.push_back(m_id);
            }
            else
            {
                E->SetEdge(j, *test.first);
                (*test.first)->m_elLink.push_back(pair<ElementSharedPtr,int>(E,j));
            }
        }
        E->SetId(i);
        m_localElements.push_back(E);
    }
}

void FaceMesh::Stretching()
{
    //define a sampling and calculate the aspect ratio of the paramter plane
    asr = 0.0;
    Array<OneD, NekDouble> bnds = m_cadsurf->GetBounds();
    pasr = (bnds[1] - bnds[0])/
           (bnds[3] - bnds[2]);

    Array<TwoD, Array<OneD,NekDouble> > stretch(40,40);

    NekDouble du = (bnds[1]-bnds[0])/(40-1);
    NekDouble dv = (bnds[3]-bnds[2])/(40-1);

    for(int i = 0; i < 40; i++)
    {
        for(int j = 0; j < 40; j++)
        {
            Array<OneD, NekDouble> uv(2);
            uv[0] = bnds[0] + i*du;
            uv[1] = bnds[2] + j*dv;
            if(i==40-1) uv[0]=bnds[1];
            if(j==40-1) uv[1]=bnds[3];
            stretch[i][j]=m_cadsurf->P(uv);
        }
    }

    int ct = 0;

    for(int i = 0; i < 40-1; i++)
    {
        for(int j = 0; j < 40-1; j++)
        {
            NekDouble ru = sqrt((stretch[i][j][0]-stretch[i+1][j][0])*
                                (stretch[i][j][0]-stretch[i+1][j][0])+
                                (stretch[i][j][1]-stretch[i+1][j][1])*
                                (stretch[i][j][1]-stretch[i+1][j][1])+
                                (stretch[i][j][2]-stretch[i+1][j][2])*
                                (stretch[i][j][2]-stretch[i+1][j][2]));
            NekDouble rv = sqrt((stretch[i][j][0]-stretch[i][j+1][0])*
                                (stretch[i][j][0]-stretch[i][j+1][0])+
                                (stretch[i][j][1]-stretch[i][j+1][1])*
                                (stretch[i][j][1]-stretch[i][j+1][1])+
                                (stretch[i][j][2]-stretch[i][j+1][2])*
                                (stretch[i][j][2]-stretch[i][j+1][2]));

            if(rv < 1E-8)
                continue;

            asr += ru/rv;
            ct++;
        }
    }

    asr/=ct;
}

bool FaceMesh::Validate()
{
    //check all edges in the current mesh for length against the octree
    //if the octree is not conformed to add a new point inside the triangle
    //if no new points are added meshing can stop
    int pointBefore = m_stienerpoints.size();
    for(int i = 0; i < m_connec.size(); i++)
    {
        Array<OneD, NekDouble> triDelta(3);

        Array<OneD, NekDouble> r(3);

        r[0]=m_connec[i][0]->Distance(m_connec[i][1]);
        r[1]=m_connec[i][1]->Distance(m_connec[i][2]);
        r[2]=m_connec[i][2]->Distance(m_connec[i][0]);

        triDelta[0] = m_octree->Query(m_connec[i][0]->GetLoc());
        triDelta[1] = m_octree->Query(m_connec[i][1]->GetLoc());
        triDelta[2] = m_octree->Query(m_connec[i][2]->GetLoc());

        int numValid = 0;

        if(r[0] < triDelta[0])
            numValid++;

        if(r[1] < triDelta[1])
            numValid++;

        if(r[2] < triDelta[2])
            numValid++;

        if(numValid != 3)
        {
            Array<OneD,NekDouble> ainfo,binfo,cinfo;
            ainfo = m_connec[i][0]->GetCADSurf(m_id);
            binfo = m_connec[i][1]->GetCADSurf(m_id);
            cinfo = m_connec[i][2]->GetCADSurf(m_id);

            Array<OneD, NekDouble> uvc(2);
            uvc[0] = (ainfo[0]+binfo[0]+cinfo[0])/3.0;
            uvc[1] = (ainfo[1]+binfo[1]+cinfo[1])/3.0;
            AddNewPoint(uvc);
        }
    }

    if(m_stienerpoints.size() == pointBefore)
    {
        return false;
    }
    else
    {
        return true;
    }
}

void FaceMesh::AddNewPoint(Array<OneD, NekDouble> uv)
{
    //adds a new point but checks that there are no other points nearby first
    Array<OneD, NekDouble> np = m_cadsurf->P(uv);
    NekDouble npDelta = m_octree->Query(np);

    NodeSharedPtr n = boost::shared_ptr<Node>(new Node(m_mesh->m_numNodes++,
                                                       np[0],np[1],np[2]));

    bool add = true;

    for(int i = 0; i < orderedLoops.size(); i++)
    {
        for(int j = 0; j < orderedLoops[i].size(); j++)
        {
            NekDouble r = orderedLoops[i][j]->Distance(n);

            if(r<npDelta/1.414)
            {
                add = false;
                break;
            }
        }
    }

    if(add)
    {
        for(int i = 0; i < m_stienerpoints.size(); i++)
        {
            NekDouble r = m_stienerpoints[i]->Distance(n);

            if(r<npDelta/1.414)
            {
                add = false;
                break;
            }
        }
    }

    if(add)
    {
        n->SetCADSurf(m_id,uv);
        m_stienerpoints.push_back(n);
    }
}

void FaceMesh::OrientateCurves()
{
    //create list of bounding loop nodes
    for(int i = 0; i < m_edgeloops.size(); i++)
    {
        vector<NodeSharedPtr> cE;
        for(int j = 0; j < m_edgeloops[i].edges.size(); j++)
        {
            int cid = m_edgeloops[i].edges[j]->GetID();
            vector<NodeSharedPtr> edgePoints = m_curvemeshes[cid]->GetMeshPoints();

            int numPoints = m_curvemeshes[cid]->GetNumPoints();

            if(m_edgeloops[i].edgeo[j] == 0)
            {
                for(int k = 0; k < numPoints-1; k++)
                {
                    cE.push_back(edgePoints[k]);
                }
            }
            else
            {
                for(int k = numPoints-1; k >0; k--)
                {
                    cE.push_back(edgePoints[k]);
                }
            }
        }
        orderedLoops.push_back(cE);
    }

    //loops made need to orientate on which is biggest and define holes

    for(int i = 0; i < orderedLoops.size(); i++)
    {
        int half = int(orderedLoops[i].size()/2) - 1;

        NodeSharedPtr n1,n2,nh;

        n1 = orderedLoops[i][0];
        n2 = orderedLoops[i][1];
        nh = orderedLoops[i][half];

        Array<OneD,NekDouble> n1info,n2info,nhinfo;
        n1info = n1->GetCADSurf(m_id);
        n2info = n2->GetCADSurf(m_id);
        nhinfo = nh->GetCADSurf(m_id);

        NekDouble ua = (100.0*n1info[0]+
                        100.0*n2info[0]+
                        1.0* nhinfo[0])/201.0 ;
        NekDouble va = (100.0*n1info[1]+
                        100.0*n2info[1]+
                        1.0* nhinfo[1])/201.0 ;

        Array<OneD, NekDouble> tmp(2);
        tmp[0]=ua;
        tmp[1]=va;
        m_edgeloops[i].center=tmp;
    }

    for(int i = 0; i < orderedLoops.size(); i++)
    {
        NekDouble area=0.0;
        for(int j = 0; j < orderedLoops[i].size()-1; j++)
        {
            Array<OneD,NekDouble> n1info,n2info;
            n1info = orderedLoops[i][j]->GetCADSurf(m_id);
            n2info = orderedLoops[i][j+1]->GetCADSurf(m_id);

            area += -n2info[1]*(n2info[0]-n1info[0])
                    +n1info[0]*(n2info[1]-n1info[1]);
        }
        area*=0.5;
        m_edgeloops[i].area = area;
    }

    int ct=0;

    do
    {
        ct=0;
        for(int i = 0; i < m_edgeloops.size()-1; i++)
        {
            if(fabs(m_edgeloops[i].area)<fabs(m_edgeloops[i+1].area))
            {
                //swap
                vector<NodeSharedPtr> orderedlooptmp = orderedLoops[i];
                EdgeLoop edgeLoopstmp = m_edgeloops[i];

                orderedLoops[i]=orderedLoops[i+1];
                m_edgeloops[i]=m_edgeloops[i+1];

                orderedLoops[i+1]=orderedlooptmp;
                m_edgeloops[i+1]=edgeLoopstmp;

                ct+=1;
            }
        }

    }while(ct>0);

    if(m_edgeloops[0].area<0) //reverse the first uvLoop
    {
        vector<NodeSharedPtr> tmp = orderedLoops[0];
        reverse(tmp.begin(), tmp.end());
        orderedLoops[0]=tmp;
    }

    for(int i = 1; i < orderedLoops.size(); i++)
    {
        if(m_edgeloops[i].area>0) //reverse the loop
        {
            vector<NodeSharedPtr> tmp = orderedLoops[i];
            reverse(tmp.begin(), tmp.end());
            orderedLoops[i]=tmp;
        }
    }
}

}
}
