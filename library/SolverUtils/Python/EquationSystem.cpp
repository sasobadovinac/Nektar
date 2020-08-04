///////////////////////////////////////////////////////////////////////////////
//
// File: EquationSystem.cpp
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
// Description: Python wrapper for Module.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <SolverUtils/EquationSystem.h>

using namespace Nektar;
using namespace Nektar::SolverUtils;
using namespace Nektar::LibUtilities;
using namespace Nektar::SpatialDomains;
/**
 * @brief EquationSystem wrapper to handle virtual function calls.
 */
struct EquationSystemWrap : public EquationSystem, public py::wrapper<EquationSystem>
{
    /**
     * @brief Constructor, which is identical
     * to NekMesh::SolverUtils::EquationSystem.
     *
     * @param pSession  given session
     * @param pGraph    given graph
     */
    EquationSystemWrap(
         const SessionReaderSharedPtr& pSession,
         const MeshGraphSharedPtr& pGraph)
         : EquationSystem(pSession, pGraph), py::wrapper<EquationSystem>()
    {
    }

    /**
     * @brief Concrete implementation of the EquationSystem::DoSolve function.
     */
    void DoSolve()
    {
        this->get_override("v_DoSolve")();
    }

    /**
     * @brief Concrete implementation of the EquationSystem::v_DoSolve function.
     */
    void v_DoSolve()
    {
        this->get_override("v_DoSolve")();
    }

};

struct EquationSystemWrapConverter
{
    EquationSystemWrapConverter()
    {
        // An important bit of code which will allow shared_ptr<EquationSystem>
        // to be interpreted as something that boost::python recognises,
        // otherwise classes constructed from the factory will not
        // work from Python.
        py::objects::class_value_wrapper<
            std::shared_ptr<EquationSystem>,
            py::objects::make_ptr_instance<
                EquationSystem,
                py::objects::pointer_holder<std::shared_ptr<EquationSystem>,
                                            EquationSystem>>
            >();
    }
};

/**
 * @brief Lightweight wrapper for EquationSystem creation.
 *
 * @param  args[0]  The name of the EquationSystem
 * @param  args[1]  SessionReader
 * @param  args[2]  MeshGraph
 */
EquationSystemSharedPtr CreateEquationSystem(py::tuple args, py::dict kwargs)
{
    std::string            name    = py::extract<std::string>(args[0]);
    SessionReaderSharedPtr session = py::extract<SessionReaderSharedPtr>(args[1]);
    MeshGraphSharedPtr     graph   = py::extract<MeshGraphSharedPtr>(args[2]);
    EquationSystemSharedPtr eqsystem =
           GetEquationSystemFactory().CreateInstance(name, session, graph);

    return eqsystem;
}

// /**
//  * @brief Wrapper for subclasses of the EquationSystem class, which
//  * can then be inhereted from Python classes.
//  */
// struct PythonEquationSystemClass
// {
//     PythonEquationSystemClass()
//     {
//         py::class_<EquationSystemWrap,
//                    std::shared_ptr<EquationSystemWrap>,
//                    py::bases<EquationSystem>,
//                    boost::noncopyable>(
//                        "EquationSystemWrap",
//                        py::init<SessionReaderSharedPtr, MeshGraphSharedPtr>()
//             )
//             .def("DoSolve", py::pure_virtual(&EquationSystem::DoSolve))
//             .def("Create", py::raw_function(CreateEquationSystem))
//             .staticmethod("Create")
//             ;
//
//         EquationSystemWrapConverter();
//     }
// };


void export_EquationSystem()
{
    py::implicitly_convertible<std::shared_ptr<EquationSystemWrap>,
                           std::shared_ptr<EquationSystem>>();
    // Export EquationSystem
    py::class_<EquationSystemWrap,
               std::shared_ptr<EquationSystemWrap>,
               boost::noncopyable>(
            "EquationSystem",
            py::init<SessionReaderSharedPtr, MeshGraphSharedPtr>())
        .def("DoSolve", py::pure_virtual(&EquationSystem::DoSolve))
        .def("Create", py::raw_function(CreateEquationSystem))
        .staticmethod("Create")
        ;
    EquationSystemWrapConverter();
}
