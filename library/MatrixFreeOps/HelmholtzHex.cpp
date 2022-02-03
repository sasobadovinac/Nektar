#include "Helmholtz.h"

namespace Nektar
{
namespace MatrixFree
{
std::string __register_Helmholtz_Hex = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Hex_Regular"), &HelmholtzHex<>::Create);

std::string __register_Helmholtz_Hex_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Hex_Deformed"), &HelmholtzHex<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
