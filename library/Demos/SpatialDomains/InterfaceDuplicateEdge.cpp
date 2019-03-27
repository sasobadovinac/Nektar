#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/Interface.h>
#include <MultiRegions/DisContField2D.h>

using namespace std;
using namespace Nektar;

double SearchForPoint(NekDouble xs[2], SpatialDomains::SegGeomSharedPtr &seg)
{



    Array<OneD, NekDouble> xi(1, 0.0);
    const NekDouble c1 = 1e-4, c2 = 0.9;

    auto xmap = seg->GetXmap();
    int nq = xmap->GetTotPoints();

    Array<OneD, NekDouble> x(nq), y(nq);
    xmap->BwdTrans(seg->GetCoeffs(0), x);
    xmap->BwdTrans(seg->GetCoeffs(1), y);

    Array<OneD, NekDouble> xder(nq), yder(nq);
    xmap->PhysDeriv(x, xder);
    xmap->PhysDeriv(y, yder);

    Array<OneD, NekDouble> xder2(nq), yder2(nq);
    xmap->PhysDeriv(xder, xder2);
    xmap->PhysDeriv(yder, yder2);

    bool opt_succeed = false;

    for (int i = 0; i < 1000; ++i)
    {

        // Compute f(x_k) and its derivatives
        NekDouble xc = xmap->PhysEvaluate(xi, x);
        NekDouble yc = xmap->PhysEvaluate(xi, y);

        NekDouble xc_der = xmap->PhysEvaluate(xi, xder);
        NekDouble yc_der = xmap->PhysEvaluate(xi, yder);

        NekDouble xc_der2 = xmap->PhysEvaluate(xi, xder2);
        NekDouble yc_der2 = xmap->PhysEvaluate(xi, yder2);

        NekDouble fx = (xc - xs[0])*(xc - xs[0]) + (yc - xs[1])*(yc - xs[1]);
        NekDouble fxp = 2.0 * (xc_der * (xc - xs[0]) + yc_der * (yc - xs[1]));
        NekDouble fxp2 = 2.0 * (xc_der2 * (xc - xs[0]) + xc_der * xc_der +
yc_der2 * (yc - xs[1]) + yc_der * yc_der);

        std::cout <<"Gradient descent iteration = " << i << "\t xi = " << xi[0]
<< "\t fx = " << fx << "\t grad = " << fxp << "\t hess = " << fxp2 << std::endl;

        // Check for convergence
        if (fx < 1e-16)
        {
            opt_succeed = true;
            break;
        }


        NekDouble gamma = 1.0;
        bool conv = false;

        // Search direction: quasi-Newton
        NekDouble pk = - fxp / fxp2;

        int l =0;
        // Backtracking line search
        while (gamma > 1e-10)
        {
            cout << "\tLine search iteration: " << l <<"\t gamma = " << gamma <<
endl; l++; Array<OneD, NekDouble> xi_pk(1); xi_pk[0] = xi[0] + pk * gamma;

            if (xi_pk[0] < -1.0 || xi_pk[0] > 1.0)
            {
                gamma /= 2.0;
                continue;
            }

            NekDouble xc_pk = xmap->PhysEvaluate(xi_pk, x);
            NekDouble yc_pk = xmap->PhysEvaluate(xi_pk, y);

            NekDouble xc_der_pk = xmap->PhysEvaluate(xi_pk, xder);
            NekDouble yc_der_pk = xmap->PhysEvaluate(xi_pk, yder);

            NekDouble fx_pk = (xc_pk - xs[0])*(xc_pk - xs[0]) + (yc_pk -
xs[1])*(yc_pk - xs[1]); NekDouble fxp_pk = 2.0 * (xc_der_pk * (xc_pk - xs[0]) +
yc_der_pk * (yc_pk - xs[1]));

            // Check Wolfe conditions
            if (fx_pk <= fx + c1 * gamma * pk * fxp && -pk * fxp_pk <= - c2 * pk
* fxp)
            {
                conv = true;
                break;
            }

            gamma /= 2.0;
        }

        if (!conv)
        {
            opt_succeed = false;
            break;
        }

        xi[0] += gamma * pk;
    }

    if (opt_succeed)
    {
        return xi[0];
    }
    else
    {
        return std::numeric_limits<NekDouble>::max();
    }
}

double rand_float( double low, double high ) {
    return ((double)rand() * (high - low)) / (double)RAND_MAX + low;
}


int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr    session;
    LibUtilities::FieldIOSharedPtr          fld;
    SpatialDomains::MeshGraphSharedPtr      graph;
    SpatialDomains::Interfaces              interfaceCollection;
    MultiRegions::DisContField2DSharedPtr   field;
//  LibUtilities::EquationSharedPtr         icond, ex_sol;
//  StdRegions::ConstFactorMap              factors;

    session = LibUtilities::SessionReader::CreateInstance(argc, argv);
    graph = SpatialDomains::MeshGraph::Read(session);

    interfaceCollection = SpatialDomains::Interfaces(session, graph);
    auto interfaces = interfaceCollection.GetInterfaces();

    //Split interfaces
    for (auto &interface : interfaces)
    {
        auto interfaceProperties = interface.second;
        interfaceProperties->SeparateGraph(graph);
    }

    fld = LibUtilities::FieldIO::CreateDefault(session);

    /*
    string       sessionName = session->GetSessionName();
    string       outFile     = sessionName + ".fld";
    unsigned int nSteps      = session->GetParameter("NumSteps");
    NekDouble    delta_t     = session->GetParameter("TimeStep");
    NekDouble    epsilon     = session->GetParameter("epsilon" );
    */

    //field = MemoryManager<MultiRegions::DisContField2D> ::AllocateSharedPtr(session, graph, session->GetVariable(0));
    field = MemoryManager<MultiRegions::DisContField2D> ::AllocateSharedPtr();
    // Write geometry
    std::string filename = "out.xml";
    graph->WriteGeometry(filename, true);

    // TESTS for SearchForPoint and the MoveDomain function
#if 0
    auto seg = graph->GetSegGeom(4);
    seg->FillGeom();

    std::map<NekDouble, NekDouble> foundPoints;

    for (int i = 0; i < 1000; i++)
    {
        NekDouble random = rand_float(1.0, 0.5);
        NekDouble xs[2] = {0.5, random};

        auto foundPoint = SearchForPoint(xs, seg);

        foundPoints[random] = foundPoint;
    }

    std::ofstream myfile;
    myfile.open ("example.csv");
    myfile << "Random number, Found point"<< endl;
    for (auto i : foundPoints)
    {
        myfile << i.first <<"," << i.second << endl;

    }
    myfile.close();

    // MOVE DOMAIN
    std::set<int> seenVerts, seenEdges;

    for (auto &comp : interfaceProperties->GetMovingDomain())
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            auto newGeom = graph->GetGeometry2D(geom->GetGlobalID());
            for (int i = 0; i < newGeom->GetNumVerts(); ++i)
            {
                PointGeomSharedPtr vert = newGeom->GetVertex(i);

                if (seenVerts.find(vert->GetGlobalID()) != seenVerts.end())
                {
                    continue;
                }

                (*vert)(1) += 10.0;
                seenVerts.insert(vert->GetGlobalID());

            }

            for (int i = 0; i < newGeom->GetNumEdges(); ++i)
            {
                SegGeomSharedPtr edge =
std::static_pointer_cast<SegGeom>(newGeom->GetEdge(i));

                // move curve points
                if (seenEdges.find(edge->GetGlobalID()) != seenEdges.end())
                {
                    continue;
                }

                CurveSharedPtr curve = edge->GetCurve();
                if (!curve)
                {
                    continue;
                }

                for (auto &pt : curve->m_points)
                {
                    (*pt)(1) += 10.0;
                }

                seenEdges.insert(edge->GetGlobalID());
            }
        }
    }

}
#endif
}
