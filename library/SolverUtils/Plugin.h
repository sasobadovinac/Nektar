#ifndef nektar_plugin_h
#define nektar_plugin_h

#include <functional>
#include <memory>
#include <iostream>

typedef std::function<void()> RegisterFunction;

namespace Nektar
{
namespace SolverUtils
{

class SolverPlugin
{
public:
    SolverPlugin() = default;
    virtual ~SolverPlugin() = default;

    virtual std::string Name() = 0;

    virtual void Initialise()
    {
        // Call initialisation functions
        for (auto &func : m_funcs)
        {
            func();
        }
    }

    bool RegisterCallback(RegisterFunction fn)
    {
        m_funcs.push_back(fn);
        return true;
    }

protected:
    std::vector<RegisterFunction> m_funcs;
};

std::shared_ptr<SolverPlugin> LoadSolverPlugin(std::string name);

}
}

#endif
