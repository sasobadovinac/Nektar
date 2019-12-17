#pragma once

#include "LSIACPostProcessor.h"
#include <iostream>
#include <string>
#include <vector>
using namespace std;

/// This is base class to handle different type of meshes.

/** Each subclass can handle one particular format of the mesh.
        This class allows users to add their own mesh format in the future
   purpose.
*/
namespace Nektar
{
namespace LSIAC
{

class HandleMesh : public LSIACPostProcessor
{
protected:
public:
};

} // namespace LSIAC
} // namespace Nektar
