#pragma once

#include "FilterType.h"
#include "LSIACPostProcessor.h"
#include <iostream>
#include <vector>
using namespace std;

/// This subclasses of the class can be used to postprocess any point on a given
/// mesh.
/** This class and its subclasses enable users to post process different points
   the mesh. These classes contain methods which are used by the end user of the
   tool.
*/
namespace Nektar
{
namespace LSIAC
{
class Smoothie : public LSIACPostProcessor
{
private:
protected:
public:
};
} // namespace LSIAC
} // namespace Nektar
