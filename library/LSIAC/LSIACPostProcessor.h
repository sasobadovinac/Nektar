#ifndef LSIACPOSTPROC_H
#define LSIACPOSTPROC_H

//#define NEKTAR_USING_LAPACK
//#define NEKTAR_USING_BLAS
// Note TOLERENCE should always be smaller than TOLERENCE_MESH_COMP
#define TOLERENCE 1e-8
#define TOLERENCE_MESH_COMP 1e-9
//#define TOLERENCE 1e-13
//#define TOLERENCE_MESH_COMP 1e-12
// TO test SpeedLSIAC V2
//#define TOLERENCE 1e-15
//#define TOLERENCE_MESH_COMP 1e-13

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
    //	void readNekArray(vector<NekDouble> &ar, string filename) const;
    void readNekArray(Array<OneD, NekDouble> &ar, string filename) const;
    void readNekArray(vector<NekDouble> &ar, string filename) const;
    void readNekArray(vector<int> &ar, string filename) const;
    void printGraphArray(const Array<OneD, NekDouble> &test, NekDouble down,
                         NekDouble up, NekDouble increment = 1.0) const;
};
} // namespace LSIAC
} // namespace Nektar
#endif
