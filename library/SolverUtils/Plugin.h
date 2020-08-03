#ifndef nektar_plugin_h
#define nektar_plugin_h

#include <boost/shared_ptr.hpp>

#include <functional>
#include <memory>
#include <iostream>

typedef std::function<void()> RegisterFunction;

namespace Nektar
{
namespace SolverUtils
{

class SolverPluginAPI
{
public:
    SolverPluginAPI() = default;
    virtual ~SolverPluginAPI() = default;

    virtual std::string Name() = 0;

    virtual void Initialise()
    {
        std::cout << "called initialise" << std::endl;

        // Call initialisation functions
        for (auto &func : m_funcs)
        {
            func();
        }
    }

    bool RegisterCallback(RegisterFunction fn)
    {
        std::cout << "called register" << std::endl;
        m_funcs.push_back(fn);
        return true;
    }

protected:
    std::vector<RegisterFunction> m_funcs;
};

boost::shared_ptr<SolverPluginAPI> LoadSolverPlugin(std::string name);

}
}

#endif
