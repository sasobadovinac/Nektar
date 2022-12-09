///////////////////////////////////////////////////////////////////////////////
//
// File: PostProcessing.cpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Kernel/kernel.h>
#include <MultiRegions/ExpList.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(argc, argv);

    int j, e;
    if (argc != 2)
    {
        fprintf(stderr, "Usage: Post-processor meshfile \n");
        exit(1);
    }
    /************************* Mesh File ***********************/
    // Read the mesh
    SpatialDomains::MeshGraphSharedPtr graph1D =
        SpatialDomains::MeshGraph::Read(vSession);

    /***********************************************************/

    // Construct an object from the class ExpList.
    // This is the class which represents a multi-elemental expansion
    // This object can be constructed based on the input mesh
    MultiRegions::ExpListSharedPtr u =
        MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(vSession,
                                                                graph1D);

    /***********************************************************/
    // Get the number of elements
    int nelements = u->GetExpSize();
    NekDouble h   = 1.0 / u->GetExpSize(); // mesh spacing

    // Get the number of modes
    Array<OneD, LibUtilities::BasisSharedPtr> base = (u->GetExp(0))->GetBase();
    int order                                      = base[0]->GetNumModes();

    /****** Read in the data file and update the coefficients ***/
    // Contains the local polynomial modes for each element
    Array<OneD, NekDouble> uhat_pre = u->UpdateCoeffs();

    ifstream inFile;
    inFile.open("uhatPre.txt");
    double x;
    int index = 0;
    while (inFile >> x)
    {
        uhat_pre[index] = x;
        index++;
    }
    inFile.close();

    /************** Define the kernel object *******************/

    LibUtilities::KernelSharedPtr post_kernel =
        MemoryManager<LibUtilities::Kernel>::AllocateSharedPtr(order);
    post_kernel->UpdateKernelBreaks(h);

    /******************* Post-Processing ***********************/
    // Construct the appropriate post-processed multi-element expansion object
    int post_order = 2 * order;

    // reset expansiosn to be of new order
    graph1D->SetExpansionInfosToPolyOrder(post_order);

    MultiRegions::ExpListSharedPtr u_post =
        MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(vSession,
                                                                graph1D);

    // The post-processing for the first element
    // Step1: Define the element ID on which the post-processing is done
    e = 0;

    // Step2: Calculate the points used for evaluating the post-processed
    // solution
    StdRegions::StdExpansionSharedPtr elmExp = u_post->GetExp(e);
    int eval_npoints                         = elmExp->GetTotPoints();
    Array<OneD, NekDouble> eval_points(eval_npoints);
    elmExp->GetCoords(eval_points);

    // Step3: Call the appropriate post-processing function
    Array<OneD, NekDouble> ustar_elm(eval_npoints);
    u->PostProcess(post_kernel, eval_points, ustar_elm, h, e);

    // The Post-Processing for the entire domain
    // Step1: Calculate the points used for evaluating the post-processed
    // solution
    int eval_npoints2 = u_post->GetNpoints();
    Array<OneD, NekDouble> eval_points2(eval_npoints2);
    u_post->GetCoords(eval_points2);

    // Step2: Call the appropriate post-processing function
    Array<OneD, NekDouble> ustar(eval_npoints2);
    u->PostProcess(post_kernel, eval_points2, ustar, h);

    // Calculate the post-processed coefficients for the entire domain
    Array<OneD, NekDouble> uhat_post = u_post->UpdateCoeffs();
    u_post->FwdTrans(ustar, uhat_post);

    ofstream postCoeff;
    postCoeff.open("postCoeff.txt");
    cout << "Printing the uhat_post: "
         << endl; // These coefficients are different than my previous result,
                  // given the ustars are the same
    for (e = 0; e < nelements; e++)
    {
        for (j = 0; j < post_order; j++)
        {
            cout << uhat_post[e * post_order + j] << "  ";
            postCoeff << uhat_post[e * post_order + j] << endl;
        }
        cout << endl;
    }
}
