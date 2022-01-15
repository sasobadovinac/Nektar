#include "Helmholtz.h"

namespace Nektar
{
namespace MatrixFree
{
std::string __register_Helmholtz_Prism = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Prism_Regular"), &HelmholtzPrism<>::Create);

std::string __register_Helmholtz_Prism_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Prism_Deformed"), &HelmholtzPrism<true>::Create);

} // namespace MatrixFree
} // namespace Nektar
