#include "Helmholtz.h"

namespace Nektar
{
namespace MatrixFree
{
    std::string __register_Helmholtz_Tet = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tet_Regular"), &HelmholtzTet<>::Create);

std::string __register_Helmholtz_Tet_Deformed = GetOperatorFactory().RegisterCreatorFunction(
    std::string("Helmholtz_Tet_Deformed"), &HelmholtzTet<true>::Create);
    
} // namespace MatrixFree
} // namespace Nektar
