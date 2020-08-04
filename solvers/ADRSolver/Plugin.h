#ifndef adrsolver_plugin
#define adrsolver_plugin

#include <vector>
#include <SolverUtils/Plugin.h>

#include <iostream>

namespace Nektar
{
namespace SolverUtils
{

class ADRSolverPlugin : public SolverPlugin
{
public:
    ADRSolverPlugin() = default;

    std::string Name() override final
    {
        return "ADRSolver";
    }
};

extern "C" ADRSolverPlugin plugin;

}
}

#endif
