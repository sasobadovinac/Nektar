////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessSurfCurvature.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Computes wss field.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include "ProcessSurfCurvature.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LocalRegions/Expansion1D.h>
#include <MultiRegions/ExpList.h>

using namespace std;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessSurfCurv::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "surfcurvature"), ProcessSurfCurv::create,
    "Computes surface curvature.");

ProcessSurfCurv::ProcessSurfCurv(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessSurfCurv::~ProcessSurfCurv()
{
}

void ProcessSurfCurv::Process(po::variables_map &vm)
{
    m_f->SetUpExp(vm);

    ASSERTL0(m_f->m_exp[0]->GetShapeDimension() == 1,"Currently this method is "
             "only set up for one-dimensional surfaces");
    ASSERTL0(m_f->m_graph->GetSpaceDimension() == 2,"Currently this method is "
             "only set up for two-dimensional meshes");

    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    int npoints = m_f->m_exp[0]->GetTotPoints();
    int nel = m_f->m_exp[0]->GetNumElmts();

    
    Array<OneD, NekDouble> x(npoints),y(npoints);
    Array<OneD, NekDouble> xd(npoints), yd(npoints);
    Array<OneD, NekDouble> xdd(npoints), ydd(npoints);
    Array<OneD, NekDouble> curvature(npoints), tmp;
    
    m_f->m_exp[0]->GetCoords(x,y);

    int cnt = 0; 
    for(int e = 0; e < nel; ++e)
    {
        LocalRegions::Expansion1DSharedPtr elmt = std::dynamic_pointer_cast
            <LocalRegions::Expansion1D>(m_f->m_exp[0]->GetExp(e));
        
        elmt->PhysTensorDeriv(x+cnt, tmp = xd+cnt);
        elmt->PhysTensorDeriv(xd+cnt,tmp = xdd+cnt);

        elmt->PhysTensorDeriv(y+cnt, tmp = yd+cnt);
        elmt->PhysTensorDeriv(yd+cnt,tmp = ydd+cnt);

        int nelmtpts = elmt->GetTotPoints();
        
        for(int i = 0; i < nelmtpts; ++i)
        {
            curvature[cnt + i] = (ydd[cnt+i]*xd[cnt+i] - yd[cnt+i]*xdd[cnt+i])
                *pow(xd[cnt+i]*xd[cnt+i] + yd[cnt+1]*yd[cnt+1],-1.5);
            curvature[cnt +i] = fabs(curvature[cnt+1]);
        }
        cnt += nelmtpts; 
    }

    m_f->m_exp[0]->FwdTransLocalElmt(curvature,m_f->m_exp[0]->UpdateCoeffs());
    
    if(m_f->m_variables.size())
    {
        m_f->m_variables[0] = "Curvature";
    }
    else
    {
        m_f->m_variables.push_back("Curvature");
    }

}


} // namespace FieldUtils
} // namespace Nektar
