#include "Helmholtz.h"

namespace Nektar
{
namespace MatrixFree
{
std::string __register_Helmholtz_Tri = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tri_Regular"), &HelmholtzTri<>::Create);

std::string __register_Helmholtz_Tri_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tri_Deformed"), &HelmholtzTri<true>::Create);
    
} // namespace MatrixFree
} // namespace Nektar
