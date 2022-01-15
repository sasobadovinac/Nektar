#include "Helmholtz.h"

namespace Nektar
{
namespace MatrixFree
{
std::string __register_Helmholtz_Pyr = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Pyr_Regular"), &HelmholtzPyr<>::Create);

std::string __register_Helmholtz_Pyr_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Pyr_Deformed"), &HelmholtzPyr<true>::Create);
} // namespace MatrixFree
} // namespace Nektar
