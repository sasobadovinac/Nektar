///////////////////////////////////////////////////////////////////////////////
//
// File: FSIBody.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Abstract base class for fluid-structure interaction bodies.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/FSI/FSIBody.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>

using namespace std;

namespace Nektar
{

FSIBodyFactory& GetFSIBodyFactory()
{
    static FSIBodyFactory instance;
    return instance;
}

FSIBody::FSIBody(
        const LibUtilities::SessionReaderSharedPtr         &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr>  &pFields)
    : m_session(pSession)
{
}

/**
 * 
 */
void FSIBody::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const std::map<std::string, std::string>          &pParams)
{
    // Read boundary parameter
    auto it = pParams.find("Boundary");
    ASSERTL0(it != pParams.end(),     "Missing parameter 'Boundary'.");
    ASSERTL0(it->second.length() > 0, "Empty parameter 'Boundary'.");
    std::string boundaryString = it->second;

    // Parse the boundary regions into a list.
    std::string::size_type firstInd =
                            boundaryString.find_first_of('[') + 1;
    std::string::size_type lastInd =
                            boundaryString.find_last_of(']') - 1;

    ASSERTL0(firstInd <= lastInd,
            (std::string("Error reading boundary region definition:") +
             boundaryString).c_str());

    std::string indString =
            boundaryString.substr(firstInd, lastInd - firstInd + 1);
    std::vector<unsigned int> boundaryRegionsIdList;
    bool parseGood = ParseUtils::GenerateSeqVector(indString,
                                                   boundaryRegionsIdList);
    ASSERTL0(parseGood && !boundaryRegionsIdList.empty(),
             (std::string("Unable to read boundary regions index "
              "range for FSIBody: ") + indString).c_str());

    // determine what boundary regions need to be considered
    unsigned int numBoundaryRegions =
                        pFields[0]->GetBndConditions().num_elements();
    m_boundaryRegionIsInList.insert(m_boundaryRegionIsInList.end(),
                                    numBoundaryRegions, 0);

    SpatialDomains::BoundaryConditions bcs(m_session,
                                            pFields[0]->GetGraph());
    const SpatialDomains::BoundaryRegionCollection &bregions =
                                            bcs.GetBoundaryRegions();

    int cnt = 0;
    for (auto &it : bregions)
    {
        if ( std::find(boundaryRegionsIdList.begin(),
                       boundaryRegionsIdList.end(), it.first) !=
                boundaryRegionsIdList.end() )
        {
            m_boundaryRegionIsInList[cnt] = 1;
        }
        cnt++;
    }
}

}
