///////////////////////////////////////////////////////////////////////////////
//
// File: VarCoeffEntry.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
// Description: Python wrapper for VarCoeffEntry.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>
#include <StdRegions/StdRegions.hpp>

#include <type_traits>

using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;

#if PY_MAJOR_VERSION == 2
void CapsuleDestructor(void *ptr)
{
    VarCoeffEntry *tmp = (VarCoeffEntry *)ptr;
    delete tmp;
}
#else
void CapsuleDestructor(PyObject *ptr)
{
    VarCoeffEntry *tmp = (VarCoeffEntry *)PyCapsule_GetPointer(ptr, 0);
    delete tmp;
}
#endif

struct VarCoeffEntryToPython
{
    static PyObject *convert(VarCoeffEntry const &entry)
    {
        // Create a Python capsule to hold a pointer that contains a lightweight
        // copy of arr. That way we guarantee Python will still have access to
        // the memory allocated inside arr even if arr is deallocated in C++.
#if PY_MAJOR_VERSION == 2
        py::object capsule(py::handle<>(PyCObject_FromVoidPtr(
            new VarCoeffEntry(entry), CapsuleDestructor)));
#else
        py::object capsule(py::handle<>(
            PyCapsule_New(new VarCoeffEntry(entry), 0,
                          (PyCapsule_Destructor)&CapsuleDestructor)));
#endif
        PyObject *tmp =
            py::incref(np::from_data(entry.GetValue().data(),
                                     np::dtype::get_builtin<NekDouble>(),
                                     py::make_tuple(entry.GetValue().size()),
                                     py::make_tuple(sizeof(NekDouble)), capsule)
                           .ptr());

        return tmp;
    }
};

struct PythonToVarCoeffEntry
{
    PythonToVarCoeffEntry()
    {
        py::converter::registry::push_back(&convertible, &construct,
                                           py::type_id<VarCoeffEntry>());
    }
    static void *convertible(PyObject *objPtr)
    {
        try
        {
            py::object obj((py::handle<>(py::borrowed(objPtr))));
            np::ndarray array = py::extract<np::ndarray>(obj);

            // Check data types match
            np::dtype dtype = np::dtype::get_builtin<
                typename boost::remove_const<NekDouble>::type>();
            if (dtype != array.get_dtype())
            {
                return 0;
            }

            // Check shape is 1D
            if (array.get_nd() != 1)
            {
                return 0;
            }
        }
        catch (boost::python::error_already_set &)
        {
            py::handle_exception();
            PyErr_Clear();
            return 0;
        }

        return objPtr;
    }

    static void decrement(void *objPtr)
    {
        if (!Py_IsInitialized())
        {
            // In deinitialisation phase, reference counters are not terribly
            // robust; decremementing counters here can lead to segfaults during
            // program exit in some cases.
            return;
        }

        // Otherwise decrement reference counter.
        py::decref((PyObject *)objPtr);
    }

    static void construct(PyObject *objPtr,
                          py::converter::rvalue_from_python_stage1_data *data)
    {
        // This has to be a _borrowed_ reference, otherwise at the end of this
        // scope it seems memory gets deallocated
        py::object obj((py::handle<>(py::borrowed(objPtr))));
        Array<OneD, NekDouble> array = py::extract<Array<OneD, NekDouble>>(obj);

        void *storage =
            ((py::converter::rvalue_from_python_storage<VarCoeffEntry> *)data)
                ->storage.bytes;
        data->convertible = storage;

        new (storage) VarCoeffEntry(array);

        py::incref(objPtr);
    }
};

void export_VarCoeffEntry()
{
    py::to_python_converter<VarCoeffEntry, VarCoeffEntryToPython>();
    PythonToVarCoeffEntry();
}
