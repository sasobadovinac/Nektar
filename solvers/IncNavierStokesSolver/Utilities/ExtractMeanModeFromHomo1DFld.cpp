///////////////////////////////////////////////////////////////////////////////
//
// File: ExtractMeanModeFromHomo1DFld.cpp
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

#include <SpatialDomains/MeshGraph.h> // for FieldDefinitions, etc
#include <StdRegions/StdTriExp.h>
#include <cstdio>
#include <cstdlib>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr,
                "Usage: ExtractmeanModeFromHomo1DFld fieldfile outfield\n");
        exit(1);
    }

    int i       = 0;
    int k       = 0;
    int n       = 0;
    int nz      = 0;
    int ncoeffs = 0;

    //----------------------------------------------
    // Import fieldfile.
    string fieldfile(argv[argc - 2]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble>> fielddata;
    LibUtilities::Import(fieldfile, fielddef, fielddata);
    //----------------------------------------------

    vector<vector<NekDouble>> combineddata;
    vector<LibUtilities::FieldDefinitionsSharedPtr> newfielddef;

    //----------------------------------------------
    // put mean data consecutively
    for (i = 0; i < fielddata.size(); ++i)
    {
        ASSERTL0(fielddef[i]->m_numHomogeneousDir == 1,
                 "Expected fieldfile to have one homogeneous direction");

        if (fielddef[i]->m_homogeneousZIDs[0] != 0)
        {
            continue;
        }
        else
        {
            nz = fielddef[i]->m_homogeneousZIDs.size();

            fielddef[i]->m_numHomogeneousDir = 0;
            fielddef[i]->m_basis.resize(2);
            newfielddef.push_back(fielddef[i]);

            // Determine the number of coefficients per element
            switch (fielddef[i]->m_shapeType)
            {
                case LibUtilities::eTriangle:
                    ncoeffs = LibUtilities::StdTriData::getNumberOfCoefficients(
                        fielddef[i]->m_numModes[0], fielddef[i]->m_numModes[1]);
                    break;
                case LibUtilities::eQuadrilateral:
                    ncoeffs =
                        fielddef[i]->m_numModes[0] * fielddef[i]->m_numModes[1];
                    break;
                default:
                    ASSERTL0(false, "Shape not recognised");
                    break;
            }

            vector<NekDouble> newdata;
            auto vec_iter = fielddata[i].begin();

            for (k = 0; k < fielddef[i]->m_fields.size(); ++k)
            {
                // copy data from each field into consecutive order
                for (n = 0; n < fielddef[i]->m_elementIDs.size(); ++n)
                {
                    // put zero mode into newdata
                    newdata.insert(newdata.end(), vec_iter, vec_iter + ncoeffs);
                    vec_iter += nz * ncoeffs;
                }
            }
            combineddata.push_back(newdata);
        }
    }
    //----------------------------------------------

    //-----------------------------------------------
    // Write out datafile.
    LibUtilities::Write(argv[argc - 1], newfielddef, combineddata);
    //-----------------------------------------------

    return 0;
}
