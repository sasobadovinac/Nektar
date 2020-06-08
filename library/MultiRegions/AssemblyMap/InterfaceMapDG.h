///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyCommDG.cpp
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
// Description: Communication for interfaces
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_INTERFACEMAPDG_H
#define NEKTAR_INTERFACEMAPDG_H

#include <MultiRegions/ExpList.h>
#include <SpatialDomains/Interface.h>

namespace Nektar
{
namespace MultiRegions
{

class InterfaceExchange
{

public:
    /// Default destructor
    MULTI_REGIONS_EXPORT virtual ~InterfaceExchange() = default;

    /// Default constructor
    MULTI_REGIONS_EXPORT InterfaceExchange(
        const SpatialDomains::InterfaceBaseShPtr &interface,
        const int &rank)
        : m_interface(interface),
        m_rank(rank)
    {
    }
private:
    SpatialDomains::InterfaceBaseShPtr m_interface;
    int m_rank;
};

typedef std::shared_ptr<InterfaceExchange>  InterfaceExchangeSharedPtr;

class InterfaceMapDG
{
public:
    /// Default destructor
    MULTI_REGIONS_EXPORT ~InterfaceMapDG() = default;

    // Constructor for interface communication
    MULTI_REGIONS_EXPORT InterfaceMapDG(
        const SpatialDomains::InterfacesSharedPtr &interfaces,
        const ExpList &locExp, const ExpListSharedPtr &trace,
        const std::map<int, int> geomIdToTraceId);

private:
    std::vector<InterfaceExchangeSharedPtr> m_exchange;
};

typedef std::shared_ptr<InterfaceMapDG>  InterfaceMapDGSharedPtr;


} // namespace MultiRegions
} // namespace Nektar

#endif