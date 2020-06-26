////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessLSIAC.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
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
//  Description: postprocesses the field quantities using the LSIAC filter.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include "ProcessLSIAC.h"
#include "ProcessMapping.h"
#include <GlobalMapping/Mapping.h>
#include <algorithm>
#include <sstream>
#include <tinyxml.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

// Add necessary LSIAC modules.
#include <LSIAC/CentralBSplines.h>
#include <LSIAC/GeneralBSplines.h>
#include <LSIAC/HandleNekMesh1D.h>
#include <LSIAC/HandleNekMesh2D.h>
#include <LSIAC/HandleNekMesh3D.h>
#include <LSIAC/SmoothieSIAC1D.cpp>
#include <LSIAC/SmoothieSIAC2D.cpp>
#include <LSIAC/SmoothieSIAC3D.cpp>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessLSIAC::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "LSIAC"), ProcessLSIAC::create,
    "Computes L-SIAC on the field quatities.");

ProcessLSIAC::ProcessLSIAC(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["configxml"] = ConfigOption(
        false, "-1", "Configuration flile for the L-SIAC paramters");
}

ProcessLSIAC::~ProcessLSIAC()
{
}

vector<string> ProcessLSIAC::getVecFromStr(string str)
{
    std::istringstream ss(str);
    vector<string> strV;
    std::string token;
    while (std::getline(ss, token, ','))
    {
        token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
        strV.push_back(token);
    }
    return strV;
}

Array<OneD, NekDouble> ProcessLSIAC::getDirFromStr(string str)
{
    Array<OneD, NekDouble> dir(3, 0.0);
    vector<string> strV = this->getVecFromStr(str);
    for (int i = 0; i < strV.size(); i++)
    {
        string temp = strV[i];
        dir[i]      = (atof(temp.c_str()));
    }
    return dir;
}

string ProcessLSIAC::getStrFromNode(TiXmlNode *pChild, string str)
{
    TiXmlNode *pNode;
    TiXmlElement *pElem;
    string rtString = "";
    pNode           = pChild->FirstChild(str);
    if (pNode != NULL)
    {
        pElem = pNode->ToElement();
        if (pElem->GetText() != NULL)
        {
            rtString = pElem->GetText();
        }
    }
    return rtString;
}

bool ProcessLSIAC::getIntVFromVals(vector<string> InVar,
                                   vector<int> &inVarIndices)
{
    for (int i = 0; i < InVar.size(); i++)
    {
        for (int j = 0; j < m_f->m_variables.size(); j++)
        {
            if (InVar[i].compare(m_f->m_variables[j]) == 0)
            {
                inVarIndices[i] = j;
            }
        }
    }

    for (int i = 0; i < inVarIndices.size(); i++)
    {
        if (inVarIndices[i] == -1)
            return false;
    }
    return true;
}

void ProcessLSIAC::ApplyLSIAC(LSIACParams set1)
{
    int i;
    int expdim  = m_f->m_graph->GetMeshDimension();
    int nfields = m_f->m_variables.size();

    // int addfields = (m_spacedim == 2) ? 1 : 3;
    int addfields = set1.outVar.size();

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    int npoints                             = m_f->m_exp[0]->GetNpoints();
    GlobalMapping::MappingSharedPtr mapping = ProcessMapping::GetMapping(m_f);

    LSIAC::HandleNekMesh *HNM = NULL;
    LSIAC::SmoothieSIAC *sm   = NULL;

    switch (expdim)
    {
        case 1:
            HNM = new LSIAC::HandleNekMesh1D(m_f);
            if (set1.derivative == 0)
            {
                sm = new LSIAC::SmoothieSIAC1D(
                    static_cast<LSIAC::SIACUtilities::FilterType>(
                        set1.filterType),
                    HNM, set1.order, 0.1);
            }
            else
            {
                sm = new LSIAC::SmoothieSIAC1D(
                    static_cast<LSIAC::SIACUtilities::FilterType>(
                        set1.filterType),
                    HNM, set1.order, 0.1, set1.derivative);
            }
            break;
        case 2:
            HNM = new LSIAC::HandleNekMesh2D(m_f);
            if (set1.derivative == 0)
            {
                sm = new LSIAC::SmoothieSIAC2D(
                    static_cast<LSIAC::SIACUtilities::FilterType>(
                        set1.filterType),
                    HNM, set1.order, 0.1);
            }
            else
            {
                sm = new LSIAC::SmoothieSIAC2D(
                    static_cast<LSIAC::SIACUtilities::FilterType>(
                        set1.filterType),
                    HNM, set1.order, 0.1, set1.derivative);
            }
            break;
        case 3:
            HNM = new LSIAC::HandleNekMesh3D(m_f);
            if (set1.derivative == 0)
            {
                sm = new LSIAC::SmoothieSIAC3D(
                    static_cast<LSIAC::SIACUtilities::FilterType>(
                        set1.filterType),
                    HNM, set1.order, 0.1);
            }
            else
            {
                sm = new LSIAC::SmoothieSIAC3D(
                    static_cast<LSIAC::SIACUtilities::FilterType>(
                        set1.filterType),
                    HNM, set1.order, 0.1, set1.derivative);
            }
            break;
    }
    // count the number of fields.
    // Apply LSIAC on all the fields. - for now only 1 field.
    //
    // Loop through all the quadrature points.
    // Apply L-SIAC at all the quadrature points.
    // Add them into the Expansion and apply write to an fld file.

    // get quadrature points.
    // Initializing
    vector<Array<OneD, NekDouble>> LSIAC_Fields;
    HNM->CalculateDynamicScaling();
	if (expdim==3)
	{
    	HNM->LoadExpListIntoRTree();
	}

    for (int i_f = 0; i_f < nfields; i_f++)
    {
        HNM->m_Arrays.push_back(m_f->m_exp[i_f]->GetPhys());
    }
    for (int i = 0; i < set1.outVar.size(); i++)
    {
        m_f->m_variables.push_back(set1.outVar[i]);
        LSIAC_Fields.push_back(Array<OneD, NekDouble>(npoints));
    }

    Array<OneD, NekDouble> xc0(npoints, 0.0), xc1(npoints, 0.0),
        xc2(npoints, 0.0);
    m_f->m_exp[0]->GetCoords(xc0, xc1, xc2);
    Array<OneD, NekDouble> direction(3, 0.0), glCoords(3, 0.0);
    direction = set1.direction;
    NekDouble dynScaling, valY, valZ;

    vector<int> outVarVal(set1.inVar.size(), -1);
    if (!getIntVFromVals(set1.inVar, outVarVal))
    {
        cout << "One or more of the input variables are not defined" << endl;
    }
    // for (int i_f=0;i_f<nfields; i_f++)
    for (int i_f = 0; i_f < set1.outVar.size(); i_f++)
    {
        for (int q = 0; q < npoints; q++)
        {
            glCoords[0] = xc0[q];
            glCoords[1] = xc1[q];
            glCoords[2] = xc2[q];
            dynScaling  = HNM->GetDynamicScaling(glCoords);
            dynScaling =
                (dynScaling > set1.maxChLen) ? set1.maxChLen : dynScaling;
            sm->EvaluateAt(xc0[q], xc1[q], xc2[q], LSIAC_Fields[i_f][q], valY,
                           valZ, direction, dynScaling, outVarVal[i_f]);
        }
    }

    vector<MultiRegions::ExpListSharedPtr> Exp(addfields);
    for (i = 0; i < addfields; ++i)
    {
        int n  = i;
        Exp[n] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
        Vmath::Vcopy(npoints, LSIAC_Fields[i], 1, Exp[n]->UpdatePhys(), 1);
        Exp[n]->FwdTrans_IterPerExp(LSIAC_Fields[i], Exp[n]->UpdateCoeffs());
    }

    int totalFldAlreadyPresent = m_f->m_exp.size();
    for (i = 0; i < addfields; ++i)
    {
        m_f->m_exp.insert(m_f->m_exp.begin() + totalFldAlreadyPresent + i,
                          Exp[i]);
    }
}

void ProcessLSIAC::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);
    string configFilename = m_config["configxml"].as<string>();
    TiXmlDocument doc(configFilename);
    bool loadOkay = doc.LoadFile();
    if (loadOkay)
    {
        TiXmlHandle hDoc(&doc);
        TiXmlElement *pElem;
        TiXmlNode *pChild;
        TiXmlHandle hRoot(0);

        hRoot = hDoc.FirstChild("LSIAC");
        for (pChild = hRoot.FirstChild("SET").Node(); pChild != 0;
             pChild = pChild->NextSibling("SET"))
        {
            LSIACParams set1;
            pElem = pChild->FirstChild("INVARIABLES")->ToElement();
            vector<string> varIn = this->getVecFromStr(pElem->GetText());
            set1.inVar           = varIn;

            pElem = pChild->FirstChild("OUTVARIABLES")->ToElement();
            vector<string> varOut = this->getVecFromStr(pElem->GetText());
            set1.outVar           = varOut;

            pElem = pChild->FirstChild("DIRECTION")->ToElement();
            Array<OneD, NekDouble> direction =
                this->getDirFromStr(pElem->GetText());
            set1.direction = direction;

            string str      = this->getStrFromNode(pChild, "DERIVATIVE");
            int derivative  = (!str.empty()) ? atoi(str.c_str()) : 0;
            set1.derivative = derivative;

            str             = this->getStrFromNode(pChild, "ORDER");
            int filterOrder = (!str.empty()) ? atoi(str.c_str()) : 2;
            set1.order      = filterOrder;

            // Advanced option not exposed to user yet.
            str                = this->getStrFromNode(pChild, "MAXCHLENGTH");
            double maxChLength = (!str.empty()) ? atof(str.c_str()) : 1e10;
            set1.maxChLen      = maxChLength;

            // Advanced option not exposed to user yet.
            str = this->getStrFromNode(pChild, "TYPE");
            if (str.empty() && derivative == 0)
            {
                set1.filterType = 0; // Simple filter
            }
            else if (str.empty() && derivative != 0)
            {
                set1.filterType = 1;
            }
            else
            {
                set1.filterType = atoi(str.c_str());
            }

            ApplyLSIAC(set1);
        }
    }
    else
    {
        // Load of xml filed.. Give help instructions
        this->printHelpMessage();
    }
}

void ProcessLSIAC::printHelpMessage()
{
    cout << "Error: Unable to load config file. Need to specify configxml file"
         << endl;
    cout << "Usage FieldConvert -LSIAC:configxml=config.xml Input1.xml "
            "Inpput2.fld output.fld"
         << endl;
    cout << "Example of an config.xml file is follows" << endl;
    cout << "<?xml version=\"1.0\" ?>" << endl;
    cout << "<LSIAC>" << endl;
    cout << "<Set>" << endl;
    cout << "<Invariables> u,v,w <\\Invariables>" << endl;
    cout << "<Outvariables> L_u,L_v,L_w <\\Outvariables>" << endl;
    cout << "<Order> 2<\\Order>" << endl;
    cout << "<Derivative> 2<\\Derivative>" << endl;
    cout << "<Direction>0.0,1.0,0.0<\\Direction>" << endl;
    cout << "<MaxChLength>0.5<\\MaxChLength>" << endl;
    cout << "<\\Set>" << endl;
    cout << "<\\LSIAC>" << endl;
}

} // namespace FieldUtils

} // namespace Nektar
