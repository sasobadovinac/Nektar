///////////////////////////////////////////////////////////////////////////////
//
// File LSIACPostProcessor.h
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
// Description: LSIACPostProcessor definition
//
///////////////////////////////////////////////////////////////////////////////
#ifndef LSIACPOSTPROC_H
#define LSIACPOSTPROC_H

// Note TOLERENCE should always be smaller than TOLERENCE_MESH_COMP
#define TOLERENCE 1e-9
#define TOLERENCE_MESH_COMP 1e-8

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ContField1D.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace Nektar;

namespace Nektar
{
namespace LSIAC
{
/// This class help import all nektar base classes.
class LSIACPostProcessor
{
private:
protected:
public:
    void printNekArray(const Array<OneD, NekDouble> &ar) const;
    void printNekArray(const vector<NekDouble> &ar) const;
    void printNekArray(const vector<int> &ar) const;
    void writeNekArray(vector<int> &ar, string filename) const;
    void writeNekArray(vector<NekDouble> &ar, string filename) const;
    void writeNekArray(Array<OneD, NekDouble> &ar, string filename) const;
    void writeNekArrayBin(Array<OneD, NekDouble> &ar, string filename) const;
    void readNekArray(Array<OneD, NekDouble> &ar, string filename) const;
    void readNekArray(vector<NekDouble> &ar, string filename) const;
    void readNekArray(vector<int> &ar, string filename) const;
    void printGraphArray(const Array<OneD, NekDouble> &test, NekDouble down,
                         NekDouble up, NekDouble increment = 1.0) const;
};
} // namespace LSIAC
} // namespace Nektar
#endif
