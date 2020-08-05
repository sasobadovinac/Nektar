///////////////////////////////////////////////////////////////////////////////
//
// File: Plugin.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Python wrapper for SolverPlugin base class.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <SolverUtils/Plugin.h>

using namespace Nektar;
using namespace Nektar::SolverUtils;

/**
 * @brief SolverPlugin wrapper to handle virtual function calls.
 */
struct SolverPluginWrap : public SolverPlugin, public py::wrapper<SolverPlugin>
{
    /**
     * @brief Constructor, which is identical to NekMesh::SolverUtils::SolverPlugin.
     */
    SolverPluginWrap()
         : SolverPlugin(), py::wrapper<SolverPlugin>()
    {
    }

    /**
     * @brief Concrete implementation of the SolverPlugin::Name function.
     */
    std::string Name()
    {
        this->get_override("Name")();
    }
};

struct SolverPluginWrapConverter
{
    SolverPluginWrapConverter()
    {
        // An important bit of code which will allow shared_ptr<SolverPlugin>
        // to be interpreted as something that boost::python recognises,
        // otherwise classes constructed from the factory will not
        // work from Python.
        py::objects::class_value_wrapper<
            std::shared_ptr<SolverPlugin>,
            py::objects::make_ptr_instance<
                SolverPlugin,
                py::objects::pointer_holder<std::shared_ptr<SolverPlugin>,
                                            SolverPlugin>>
            >();
    }
};


void export_Plugin()
{
    py::implicitly_convertible<std::shared_ptr<SolverPluginWrap>,
                               std::shared_ptr<SolverPlugin>>();

    py::class_<SolverPluginWrap, boost::noncopyable>("Plugin")
        .def("Name", py::pure_virtual(&SolverPlugin::Name))
        ;

    py::def("LoadSolverPlugin", &LoadSolverPlugin);

    SolverPluginWrapConverter();
}
