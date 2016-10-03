////////////////////////////////////////////////////////////////////////////////
//
//  File: ElUtil.cpp
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
//  Description: Calculate jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include "ElUtil.h"
#include "ProcessVarOpti.h"

#include <LibUtilities/Foundations/ManagerAccess.h>

#include <Kokkos_Core.hpp>

using namespace std;

namespace Nektar
{
namespace Utilities
{

boost::mutex mtx2;

ElUtil::ElUtil(ElementSharedPtr e, DerivUtilSharedPtr d,
               ResidualSharedPtr r, int n)
{
    m_el = e;
    derivUtil = d;
    res = r;
    m_mode = n;
    m_dim = m_el->GetDim();
    vector<NodeSharedPtr> ns;
    m_el->GetCurvedNodes(ns);
    nodes.resize(ns.size());
    for (int i = 0; i < ns.size(); ++i)
    {
        nodes[i].resize(m_dim);
        nodes[i][0] = &ns[i]->m_x;

        if (m_dim >= 2)
        {
            nodes[i][1] = &ns[i]->m_y;
        }

        if (m_dim >= 3)
        {
            nodes[i][2] = &ns[i]->m_z;
        }
    }
    maps = MappingIdealToRef();
}

vector<Array<OneD, NekDouble> > ElUtil::MappingIdealToRef()
{
    //need to make ideal element out of old element
    /*ElmtConfig ec = m_el->GetConf();
    ec.m_order  = 1;
    ec.m_faceNodes = false;
    ec.m_volumeNodes = false;
    ec.m_reorient = false;

    ElementSharedPtr E = GetElementFactory().CreateInstance(
                            ec.m_e, ec, m_el->GetVertexList(),
                            m_el->GetTagList());

    SpatialDomains::GeometrySharedPtr    geom = E->GetGeom(m_dim);
    geom->FillGeom();
    StdRegions::StdExpansionSharedPtr    chi  = geom->GetXmap();*/

    vector<Array<OneD, NekDouble> > ret;

    if(m_el->GetConf().m_e == LibUtilities::eQuadrilateral)
    {
        ASSERTL0(false,"Not coded");
        /*vector<Array<OneD, NekDouble> > xy;
        for(int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(2);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xy.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> u = b[0]->GetZ();
        Array<OneD, NekDouble> v = b[1]->GetZ();

        for(int j = 0; j < b[1]->GetNumPoints(); j++)
        {
            for(int i = 0; i < b[0]->GetNumPoints(); i++)
            {
                NekDouble a1 = 0.5*(1.0-u[i]), a2 = 0.5*(1.0+u[i]);
                NekDouble b1 = 0.5*(1.0-v[j]), b2 = 0.5*(1.0+v[j]);
                DNekMat dxdz(2,2,1.0,eFULL);

                dxdz(0,0) = 0.5*(-b1*xy[0][0] + b1*xy[1][0] + b2*xy[2][0] - b2*xy[3][0]);
                dxdz(1,0) = 0.5*(-b1*xy[0][1] + b1*xy[1][1] + b2*xy[2][1] - b2*xy[3][1]);

                dxdz(0,1) = 0.5*(-a1*xy[0][0] - a2*xy[1][0] + a2*xy[2][0] + a1*xy[3][0]);
                dxdz(1,1) = 0.5*(-a1*xy[0][1] - a2*xy[1][1] + a2*xy[2][1] + a1*xy[3][1]);

                NekDouble det = 1.0/(dxdz(0,0)*dxdz(1,1) - dxdz(1,0)*dxdz(0,1));

                dxdz.Invert();
                Array<OneD, NekDouble> r(9,0.0);
                r[0] = dxdz(0,0);
                r[1] = dxdz(1,0);
                r[3] = dxdz(0,1);
                r[4] = dxdz(1,1);
                ret.push_back(r);
            }
        }*/
    }
    else if(m_el->GetConf().m_e == LibUtilities::eTriangle)
    {
        DNekMat J(2,2,0.0);
        J(0,0) = (*nodes[1][0] - *nodes[0][0]);
        J(1,0) = (*nodes[1][1] - *nodes[0][1]);
        J(0,1) = (*nodes[2][0] - *nodes[0][0]);
        J(1,1) = (*nodes[2][1] - *nodes[0][1]);

        J.Invert();

        DNekMat R(2,2,0.0);
        R(0,0) = 2.0;
        R(1,1) = 2.0;

        J = J * R;

        for(int i = 0 ; i < derivUtil->ptsHigh; i++)
        {
            Array<OneD, NekDouble> r(10,0.0); //store det in 10th entry

            r[9] = 1.0 / (J(0,0) * J(1,1) - J(0,1) * J(1,0));
            r[0] = J(0,0);
            r[1] = J(1,0);
            r[2] = 0.0;
            r[3] = J(0,1);
            r[4] = J(1,1);
            r[5] = 0.0;
            r[6] = 0.0;
            r[7] = 0.0;
            r[8] = 0.0;
            ret.push_back(r);
        }
    }
    else if(m_el->GetConf().m_e == LibUtilities::eTetrahedron)
    {
        DNekMat J(3,3,0.0);
        J(0,0) = (*nodes[1][0] - *nodes[0][0]);
        J(1,0) = (*nodes[1][1] - *nodes[0][1]);
        J(2,0) = (*nodes[1][2] - *nodes[0][2]);
        J(0,1) = (*nodes[2][0] - *nodes[0][0]);
        J(1,1) = (*nodes[2][1] - *nodes[0][1]);
        J(2,1) = (*nodes[2][2] - *nodes[0][2]);
        J(0,2) = (*nodes[3][0] - *nodes[0][0]);
        J(1,2) = (*nodes[3][1] - *nodes[0][1]);
        J(2,2) = (*nodes[3][2] - *nodes[0][2]);

        J.Invert();

        DNekMat R(3,3,0.0);
        R(0,0) = 2.0;
        R(1,1) = 2.0;
        R(2,2) = 2.0;

        J = J * R;

        for(int i = 0 ; i < derivUtil->ptsHigh; i++)
        {
            Array<OneD, NekDouble> r(10,0.0); //store det in 10th entry

            r[9] = 1.0/(J(0,0)*(J(1,1)*J(2,2)-J(2,1)*J(1,2))
                       -J(0,1)*(J(1,0)*J(2,2)-J(2,0)*J(1,2))
                       +J(0,2)*(J(1,0)*J(2,1)-J(2,0)*J(1,1)));

            r[0] = J(0,0);
            r[1] = J(1,0);
            r[2] = J(2,0);
            r[3] = J(0,1);
            r[4] = J(1,1);
            r[5] = J(2,1);
            r[6] = J(0,2);
            r[7] = J(1,2);
            r[8] = J(2,2);
            ret.push_back(r);
        }
    }
    else if(m_el->GetConf().m_e == LibUtilities::ePrism)
    {
        ASSERTL0(false, "not coded");
        /*vector<Array<OneD, NekDouble> > xyz;
        for(int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(3);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xyz.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> eta1 = b[0]->GetZ();
        Array<OneD, NekDouble> eta2 = b[1]->GetZ();
        Array<OneD, NekDouble> eta3 = b[2]->GetZ();

        for(int k = 0; k < b[2]->GetNumPoints(); k++)
        {

            for(int j = 0; j < b[1]->GetNumPoints(); j++)
            {
                for(int i = 0; i < b[0]->GetNumPoints(); i++)
                {
                    NekDouble xi1 = 0.5*(1+eta1[i])*(1-eta3[k])-1.0;
                    NekDouble a1 = 0.5*(1-xi1),     a2 = 0.5*(1+xi1);
                    NekDouble b1 = 0.5*(1-eta2[j]), b2 = 0.5*(1+eta2[j]);
                    NekDouble c1 = 0.5*(1-eta3[k]), c2 = 0.5*(1+eta3[k]);

                    DNekMat dxdz(3,3,1.0,eFULL);

                    dxdz(0,0) = 0.5*(-b1*xyz[0][0] + b1*xyz[1][0] + b2*xyz[2][0] - b2*xyz[3][0]);
                    dxdz(1,0) = 0.5*(-b1*xyz[0][1] + b1*xyz[1][1] + b2*xyz[2][1] - b2*xyz[3][1]);
                    dxdz(2,0) = 0.5*(-b1*xyz[0][2] + b1*xyz[1][2] + b2*xyz[2][2] - b2*xyz[3][2]);

                    dxdz(0,1) = 0.5*((a2-c1)*xyz[0][0] - a2*xyz[1][0] + a2*xyz[2][0] + (c1-a2)*xyz[3][0] - c2*xyz[4][0] + c2*xyz[5][0]);
                    dxdz(1,1) = 0.5*((a2-c1)*xyz[0][1] - a2*xyz[1][1] + a2*xyz[2][1] + (c1-a2)*xyz[3][1] - c2*xyz[4][1] + c2*xyz[5][1]);
                    dxdz(2,1) = 0.5*((a2-c1)*xyz[0][2] - a2*xyz[1][2] + a2*xyz[2][2] + (c1-a2)*xyz[3][2] - c2*xyz[4][2] + c2*xyz[5][2]);

                    dxdz(0,2) = 0.5*(-b1*xyz[0][0] - b2*xyz[3][0] + b1*xyz[4][0] + b2*xyz[5][0]);
                    dxdz(1,2) = 0.5*(-b1*xyz[0][1] - b2*xyz[3][1] + b1*xyz[4][1] + b2*xyz[5][1]);
                    dxdz(2,2) = 0.5*(-b1*xyz[0][2] - b2*xyz[3][2] + b1*xyz[4][2] + b2*xyz[5][2]);

                    dxdz.Invert();
                    Array<OneD, NekDouble> r(9,0.0);
                    r[0] = dxdz(0,0);
                    r[1] = dxdz(1,0);
                    r[3] = dxdz(0,1);
                    r[4] = dxdz(1,1);
                    r[2] = dxdz(2,0);
                    r[5] = dxdz(2,1);
                    r[6] = dxdz(0,2);
                    r[7] = dxdz(1,2);
                    r[8] = dxdz(2,2);
                    ret.push_back(r);
                }
            }
        }*/
    }
    else
    {
        ASSERTL0(false,"not coded");
    }

    return ret;
}



void ElUtil::Evaluate()
{
    NekDouble mx = -1.0 * DBL_MAX;
    NekDouble mn =  DBL_MAX;

    if(m_dim == 2)
    {
        int nodes_size = nodes.size();

        NekMatrix<NekDouble> derivUtil_VdmDL_0 = derivUtil->VdmDL[0];
        NekMatrix<NekDouble> derivUtil_VdmDL_1 = derivUtil->VdmDL[1];

        NekVector<NekDouble> X(nodes_size),Y(nodes_size),Z(nodes_size);

        for(int j = 0; j < nodes_size; j++)
        {
            X(j) = *nodes[j][0];
            Y(j) = *nodes[j][1];
        }

        NekVector<NekDouble> x1i(nodes_size),y1i(nodes_size),
                             x2i(nodes_size),y2i(nodes_size);

        x1i = derivUtil_VdmDL_0*X;
        y1i = derivUtil_VdmDL_0*Y;
        x2i = derivUtil_VdmDL_1*X;
        y2i = derivUtil_VdmDL_1*Y;

        for(int j = 0; j < nodes_size; j++)
        {
            NekDouble jacDet = x1i(j) * y2i(j) - x2i(j)*y1i(j);
            //mx = max(mx,jacDet);
            mx = (mx < jacDet ? jacDet : mx);
            //mn = min(mn,jacDet);
            mn = (mn > jacDet ? jacDet : mn);
        }
    }
    else if(m_dim == 3)
    {
        typedef Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace> range_policy;
        typedef Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace> range_policy_host;
        typedef Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> team_policy;
        typedef Kokkos::TeamPolicy<Kokkos::DefaultHostExecutionSpace> team_policy_host;
        typedef Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type  member_type;
        typedef Kokkos::TeamPolicy<Kokkos::DefaultHostExecutionSpace>::member_type  member_type_host;
        typedef Kokkos::MemoryTraits<Kokkos::RandomAccess> random_memory;
        //typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> random_memory;

        //int res_startInv = res->startInv;
        //NekDouble res_worstJac = res->worstJac;
        int nodes_size = nodes.size();
           

        // initialising and copying node coordinates X Y Z onto the GPU
        NekVector<NekDouble> Xn(nodes_size),Yn(nodes_size),Zn(nodes_size);
        //std::vector<std::vector<NekDouble *> > nodes
        for(int j = 0; j < nodes_size; j++)
        {
            Xn(j) = *nodes[j][0];
            Yn(j) = *nodes[j][1];
            Zn(j) = *nodes[j][2];
        }
        Kokkos::View<double*> X("X",nodes_size);
        typename Kokkos::View< double*>::HostMirror h_X = Kokkos::create_mirror_view(X);
        Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> t_X( Xn.GetRawPtr() ,nodes_size);
        h_X = t_X;
        Kokkos::View<double*> Y("Y",nodes_size);
        typename Kokkos::View< double*>::HostMirror h_Y = Kokkos::create_mirror_view(Y);
        Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> t_Y( Yn.GetRawPtr() ,nodes_size);
        h_Y = t_Y;
        Kokkos::View<double*> Z("Z",nodes_size);
        typename Kokkos::View< double*>::HostMirror h_Z = Kokkos::create_mirror_view(Z);
        Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> t_Z( Zn.GetRawPtr() ,nodes_size);
        h_Z = t_Z;
        
        Kokkos::deep_copy(X,h_X);
        Kokkos::deep_copy(Y,h_Y);
        Kokkos::deep_copy(Z,h_Z);



        // initialising and copying Vandermonde coefficients onto the GPU
        NekMatrix<NekDouble> derivUtil_VdmDL_0 = derivUtil->VdmDL[0];
        NekMatrix<NekDouble> derivUtil_VdmDL_1 = derivUtil->VdmDL[1];
        NekMatrix<NekDouble> derivUtil_VdmDL_2 = derivUtil->VdmDL[2];

        //Kokkos::View<double**> VdmDL_0("VdmDL_0",nodes_size,nodes_size);
        Kokkos::View<double**,random_memory> VdmDL_0("VdmDL_0",nodes_size,nodes_size);
        typename Kokkos::View< const double**>::HostMirror h_VdmDL_0 = Kokkos::create_mirror_view(VdmDL_0);
        Kokkos::View<double**,random_memory> VdmDL_1("VdmDL_1",nodes_size,nodes_size);
        typename Kokkos::View< const double**>::HostMirror h_VdmDL_1 = Kokkos::create_mirror_view(VdmDL_1);        
        Kokkos::View<double**,random_memory> VdmDL_2("VdmDL_2",nodes_size,nodes_size);
        typename Kokkos::View< const double**>::HostMirror h_VdmDL_2 = Kokkos::create_mirror_view(VdmDL_2);
        //Kokkos::View<double**,Kokkos::LayoutLeft, Kokkos::DefaultHostExecutionSpace> t_VdmDL_0( derivUtil_VdmDL_0.GetRawPtr() ,nodes_size, nodes_size);
        //h_VdmDL_0 = t_VdmDL_0;     // easy version, but will result in wrong layout
        //Kokkos::View<double**,Kokkos::LayoutLeft, Kokkos::DefaultHostExecutionSpace> t_VdmDL_1( derivUtil_VdmDL_1.GetRawPtr() ,nodes_size, nodes_size);
        //h_VdmDL_1 = t_VdmDL_1;     // easy version, but will result in wrong layout
        //Kokkos::View<double**,Kokkos::LayoutLeft, Kokkos::DefaultHostExecutionSpace> t_VdmDL_2( derivUtil_VdmDL_2.GetRawPtr() ,nodes_size, nodes_size);
        //h_VdmDL_2 = t_VdmDL_2;     // easy version, but will result in wrong layout
        
        int N2 = nodes_size;
        int M2 = nodes_size;
        Kokkos::parallel_for(team_policy_host(N2,M2), KOKKOS_LAMBDA (const member_type_host& teamMember)
        {         
            const int i = teamMember.league_rank();
            Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, M2 ), [&] (const int j)
            {                
                h_VdmDL_0(i,j) = derivUtil_VdmDL_0(i,j);
                h_VdmDL_1(i,j) = derivUtil_VdmDL_1(i,j);
                h_VdmDL_2(i,j) = derivUtil_VdmDL_2(i,j);                
            }); 
        });      
       
        Kokkos::deep_copy(VdmDL_0,h_VdmDL_0);
        Kokkos::deep_copy(VdmDL_1,h_VdmDL_1);
        Kokkos::deep_copy(VdmDL_2,h_VdmDL_2);


        // do the matrix vector multiplication on the GPU
        Kokkos::View<double*> x1i("x1i",nodes_size);
        Kokkos::View<double*> y1i("y1i",nodes_size);
        Kokkos::View<double*> z1i("z1i",nodes_size);
        Kokkos::View<double*> x2i("x2i",nodes_size);
        Kokkos::View<double*> y2i("y2i",nodes_size);
        Kokkos::View<double*> z2i("z2i",nodes_size);
        Kokkos::View<double*> x3i("x3i",nodes_size);
        Kokkos::View<double*> y3i("y3i",nodes_size);
        Kokkos::View<double*> z3i("z3i",nodes_size);
        
        int N = nodes_size;
        int M = nodes_size;
        /*Kokkos::parallel_for( team_policy( N , Kokkos::AUTO ), KOKKOS_LAMBDA ( const member_type& teamMember)
        {
            const int i = teamMember.league_rank();
            y1i(i) = 0.0;           
            Kokkos::parallel_reduce( Kokkos::TeamThreadRange( teamMember, M ), [&] (const int j, double &update )
            {
                update += VdmDL_0( i , j ) * Y( j );
            }, y1i(i));      
        }); // not working yet */
        Kokkos::parallel_for( range_policy( 0 , N ), KOKKOS_LAMBDA ( const int i)
        {
            x1i[i] = 0;
            y1i[i] = 0;
            z1i[i] = 0;
            x2i[i] = 0;
            y2i[i] = 0;
            z2i[i] = 0;
            x3i[i] = 0;
            y3i[i] = 0;
            z3i[i] = 0;           
            for (int j = 0; j < M; ++j)
            {
                x1i(i) += VdmDL_0( i , j ) * X( j );
                y1i(i) += VdmDL_0( i , j ) * Y( j );
                z1i(i) += VdmDL_0( i , j ) * Z( j );
                x2i(i) += VdmDL_1( i , j ) * X( j );
                y2i(i) += VdmDL_1( i , j ) * Y( j );
                z2i(i) += VdmDL_1( i , j ) * Z( j );
                x3i(i) += VdmDL_2( i , j ) * X( j );
                y3i(i) += VdmDL_2( i , j ) * Y( j );
                z3i(i) += VdmDL_2( i , j ) * Z( j );
            }      
        });     

        
        // pure GPU version =====================================================
        // calculate the Jacobian       
        Kokkos::View<double[3][3]> dxdz("dxdz");
        Kokkos::View<double*> jacDet("jacDet", nodes_size);
        Kokkos::parallel_for(range_policy(0,nodes_size), KOKKOS_LAMBDA (const int i)
        {
            dxdz(0,0) = x1i(i);
            dxdz(0,1) = y1i(i);
            dxdz(0,2) = z1i(i);
            dxdz(1,0) = x2i(i);
            dxdz(1,1) = y2i(i);
            dxdz(1,2) = z2i(i);
            dxdz(2,0) = x3i(i);
            dxdz(2,1) = y3i(i);
            dxdz(2,2) = z3i(i);

            jacDet(i) = dxdz(0,0)*(dxdz(1,1)*dxdz(2,2)-dxdz(2,1)*dxdz(1,2))
                              -dxdz(0,1)*(dxdz(1,0)*dxdz(2,2)-dxdz(2,0)*dxdz(1,2))
                              +dxdz(0,2)*(dxdz(1,0)*dxdz(2,1)-dxdz(2,0)*dxdz(1,1));             
        });
        
        // CPU-Version with GPU reduce =================================================
        /*
        Kokkos::View<double*[9],Kokkos::DefaultHostExecutionSpace> deriv("deriv",nodes_size);
        Kokkos::View<double*[9],Kokkos::DefaultHostExecutionSpace> h_deriv("deriv",nodes_size);        
        
        Kokkos::parallel_for(range_policy_host(0,nodes_size), KOKKOS_LAMBDA (const int i)
        {         
            h_deriv(i,0) = h2_x1i[i];
            h_deriv(i,1) = h2_y1i[i];
            h_deriv(i,2) = h2_z1i[i];
            h_deriv(i,3) = h2_x2i[i];
            h_deriv(i,4) = h2_y2i[i];
            h_deriv(i,5) = h2_z2i[i];
            h_deriv(i,6) = h2_x3i[i];
            h_deriv(i,7) = h2_y3i[i];
            h_deriv(i,8) = h2_z3i[i];
        });

        Kokkos::deep_copy(deriv,h_deriv);
        
        Kokkos::View<double*[3][3],Kokkos::DefaultHostExecutionSpace> dxdz("dxdz", nodes_size);
        Kokkos::View<double*,Kokkos::DefaultExecutionSpace> jacDet("jacDet", nodes_size);
        typename Kokkos::View< double*>::HostMirror h_jacDet = Kokkos::create_mirror_view(jacDet);

        Kokkos::parallel_for(range_policy_host(0,nodes_size), KOKKOS_LAMBDA (const int j)
        {
            dxdz(j,0,0) = deriv(j,0);
            dxdz(j,0,1) = deriv(j,1);
            dxdz(j,0,2) = deriv(j,2);
            dxdz(j,1,0) = deriv(j,3);
            dxdz(j,1,1) = deriv(j,4);
            dxdz(j,1,2) = deriv(j,5);
            dxdz(j,2,0) = deriv(j,6);
            dxdz(j,2,1) = deriv(j,7);
            dxdz(j,2,2) = deriv(j,8);

            h_jacDet(j) = dxdz(j,0,0)*(dxdz(j,1,1)*dxdz(j,2,2)-dxdz(j,2,1)*dxdz(j,1,2))
                              -dxdz(j,0,1)*(dxdz(j,1,0)*dxdz(j,2,2)-dxdz(j,2,0)*dxdz(j,1,2))
                              +dxdz(j,0,2)*(dxdz(j,1,0)*dxdz(j,2,1)-dxdz(j,2,0)*dxdz(j,1,1));
        });
        Kokkos::deep_copy(jacDet,h_jacDet); // ========================================*/
        
        // resume
        MaxFunctor <double> mxfunctor(jacDet);
        Kokkos::parallel_reduce(range_policy(0, nodes_size) , mxfunctor, mx);
        MinFunctor <double> mnfunctor(jacDet);
        Kokkos::parallel_reduce(range_policy(0, nodes_size) , mnfunctor, mn);

        //mx = double(mx);
        //mn = double(mn);
        /*for(int j = 0; j < nodes_size; j++)
        {
            //  mx = max(mx,jacDet);
            mx = (mx < jacDet(j) ? jacDet(j) : mx);
            //  mn = min(mn,jacDet);
            mn = (mn > jacDet(j) ? jacDet(j) : mn);
        }*/
    }

    mtx2.lock();
    if(mn < 0)
    {
        res->startInv++;
    }
    //  res->worstJac = min(res->worstJac,mn/mx);
    res->worstJac = (res->worstJac > mn/mx ? mn/mx : res->worstJac);

    //res->worstJac = res_worstJac;
    //res->startInv = res_startInv;
    mtx2.unlock();

    //mtx2.lock();
    maps = MappingIdealToRef();
    //mtx2.unlock();

    minJac = mn;
    //printf("minJac: %f\n", minJac);
    scaledJac = mn/mx;
    //printf("scaledJac: %f\n", scaledJac);

}

ElUtilJob* ElUtil::GetJob()
{
    return new ElUtilJob(this);
}

}
}
