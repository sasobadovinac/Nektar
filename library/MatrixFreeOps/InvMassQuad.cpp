#include "InvMass.h"

namespace Nektar
{
namespace MatrixFree
{

std::string __register_InvMass_Quad = GetOperatorFactory().
    RegisterCreatorFunction(
    std::string("InvMass_Quad_Regular"), &InvMassQuad<>::Create);

std::string __register_InvMass_Quad_Deformed = GetOperatorFactory().
    RegisterCreatorFunction(
    std::string("InvMass_Quad_Deformed"), &InvMassQuad<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
