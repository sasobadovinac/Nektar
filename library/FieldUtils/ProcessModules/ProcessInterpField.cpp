////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessInterpField.cpp
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
//  Description: Interpolate one field to another.
//
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>
#include <boost/geometry.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <FieldUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessInterpField.h"

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessInterpField::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "interpfield"), ProcessInterpField::create,
        "Interpolates one field to another, requires fromxml, "
        "fromfld to be defined");

ProcessInterpField::ProcessInterpField(FieldSharedPtr f) : ProcessModule(f)
{

    m_config["fromxml"] = ConfigOption(
        false, "NotSet", "Xml file from which to interpolate field");
    m_config["fromfld"] = ConfigOption(
        false, "NotSet", "Fld file from which to interpolate field");

    m_config["clamptolowervalue"] =
        ConfigOption(false, "-10000000", "Lower bound for interpolation value");
    m_config["clamptouppervalue"] =
        ConfigOption(false, "10000000", "Upper bound for interpolation value");
    m_config["defaultvalue"] =
        ConfigOption(false, "0", "Default value if point is outside domain");
    m_config["realmodetoimag"] =
        ConfigOption(false, "NotSet", "Take fields as sin mode");
}

ProcessInterpField::~ProcessInterpField()
{
}

void ProcessInterpField::Process(po::variables_map &vm)
{
    m_f->SetUpExp(vm);

    FieldSharedPtr fromField = std::shared_ptr<Field>(new Field());

    std::vector<std::string> files;

    // set up session file for from field
    char *argv[] = {const_cast<char *>("FieldConvert"), nullptr};
    ParseUtils::GenerateVector(m_config["fromxml"].as<string>(), files);
    fromField->m_session = LibUtilities::SessionReader::CreateInstance(
        1, argv, files,
        LibUtilities::GetCommFactory().CreateInstance("Serial", 0, 0));

    // Set up range based on min and max of local parallel partition
    LibUtilities::DomainRangeShPtr rng =
        MemoryManager<LibUtilities::DomainRange>::AllocateSharedPtr();

    int numHomoDir = m_f->m_numHomogeneousDir;
    int coordim    = m_f->m_exp[0]->GetCoordim(0) + numHomoDir;
    int npts       = m_f->m_exp[0]->GetTotPoints();
    Array<OneD, Array<OneD, NekDouble>> coords(3);

    for (int i = 0; i < coordim; ++i)
    {
        coords[i] = Array<OneD, NekDouble>(npts);
    }

    for (int i = coordim; i < 3; ++i)
    {
        coords[i] = NullNekDouble1DArray;
    }

    m_f->m_exp[0]->GetCoords(coords[0], coords[1], coords[2]);

    rng->m_checkShape = false;
    switch (coordim)
    {
        case 3:
            rng->m_doZrange = true;
            rng->m_zmin     = Vmath::Vmin(npts, coords[2], 1);
            rng->m_zmax     = Vmath::Vmax(npts, coords[2], 1);
            /* Falls through. */
        case 2:
            rng->m_doYrange = true;
            rng->m_ymin     = Vmath::Vmin(npts, coords[1], 1);
            rng->m_ymax     = Vmath::Vmax(npts, coords[1], 1);
            /* Falls through. */
        case 1:
            rng->m_doXrange = true;
            rng->m_xmin     = Vmath::Vmin(npts, coords[0], 1);
            rng->m_xmax     = Vmath::Vmax(npts, coords[0], 1);
            break;
        default:
            NEKERROR(ErrorUtil::efatal, "coordim should be <= 3");
    }

    // setup rng parameters.
    fromField->m_graph =
        SpatialDomains::MeshGraph::Read(fromField->m_session, rng);

    // Read in local from field partitions
    const SpatialDomains::ExpansionInfoMap &expansions =
        fromField->m_graph->GetExpansionInfo();

    // check for case where no elements are specified on this
    // parallel partition
    if (!expansions.size())
    {
        return;
    }

    Array<OneD, int> ElementGIDs(expansions.size());

    int i = 0;
    for (auto &expIt : expansions)
    {
        ElementGIDs[i++] = expIt.second->m_geomShPtr->GetGlobalID();
    }

    string fromfld = m_config["fromfld"].as<string>();
    m_f->FieldIOForFile(fromfld)->Import(
        fromfld, fromField->m_fielddef, fromField->m_data,
        LibUtilities::NullFieldMetaDataMap, ElementGIDs);

    int fromNumHomoDir = fromField->m_fielddef[0]->m_numHomogeneousDir;
    for (i = 0; i < fromField->m_fielddef.size(); ++i)
    {
        int d1 = fromField->m_fielddef[i]->m_basis.size();
        d1 -= 1;
        if (d1 >= 0 && (fromField->m_fielddef[i]->m_basis[d1] ==
                            LibUtilities::eFourierHalfModeRe ||
                        fromField->m_fielddef[i]->m_basis[d1] ==
                            LibUtilities::eFourierHalfModeIm))
        {
            fromField->m_fielddef[i]->m_homogeneousZIDs[0] += 2;
            fromField->m_fielddef[i]->m_numModes[d1] = 4;
            fromField->m_fielddef[i]->m_basis[d1]    = LibUtilities::eFourier;
        }
    }

    //----------------------------------------------
    // Set up Expansion information to use mode order from field
    fromField->m_graph->SetExpansionInfo(fromField->m_fielddef);

    int nfields = fromField->m_fielddef[0]->m_fields.size();

    fromField->m_exp.resize(nfields);
    fromField->m_exp[0] = fromField->SetUpFirstExpList(fromNumHomoDir, true);

    m_f->m_exp.resize(nfields);

    // declare auxiliary fields.
    for (i = 1; i < nfields; ++i)
    {
        m_f->m_exp[i]       = m_f->AppendExpList(numHomoDir);
        fromField->m_exp[i] = fromField->AppendExpList(fromNumHomoDir);
    }

#if EXPLISTDATA
#else
    // delcare memory
    m_f->m_fieldCoeffs= std::make_shared<NekField<NekDouble,eCoeff>>(m_f->m_exp);
    m_f->m_fieldPhys  = std::make_shared<NekField<NekDouble,ePhys >>(m_f->m_exp);
    fromField->m_fieldCoeffs= std::make_shared<NekField<NekDouble,eCoeff>>(fromField->m_exp);
    fromField->m_fieldPhys  = std::make_shared<NekField<NekDouble,ePhys >>(fromField->m_exp);
#endif
    
    // load field into expansion in fromfield.
    set<int> sinmode;
    if (m_config["realmodetoimag"].as<string>().compare("NotSet"))
    {
        vector<int> value;
        ASSERTL0(ParseUtils::GenerateVector(
                     m_config["realmodetoimag"].as<string>(), value),
                 "Failed to interpret realmodetoimag string");
        for (int j : value)
        {
            sinmode.insert(j);
        }
    }
    for (int j = 0; j < nfields; ++j)
    {
        for (i = 0; i < fromField->m_fielddef.size(); i++)
        {
            fromField->m_exp[j]->ExtractDataToCoeffs(
                fromField->m_fielddef[i], fromField->m_data[i],
                fromField->m_fielddef[0]->m_fields[j],
#if EXPLISTDATA
                fromField->m_exp[j]->UpdateCoeffs());
#else
                fromField->m_fieldCoeffs->UpdateArray1D(j));
#endif
        }
        if (fromNumHomoDir == 1)
        {
            fromField->m_exp[j]->SetWaveSpace(true);
            if (sinmode.count(j))
            {
                int Ncoeff = fromField->m_exp[j]->GetPlane(2)->GetNcoeffs();
#if EXPLISTDATA
                Vmath::Smul(
                    Ncoeff, -1., fromField->m_exp[j]->GetPlane(2)->GetCoeffs(),
                    1, fromField->m_exp[j]->GetPlane(3)->UpdateCoeffs(), 1);
                Vmath::Zero(Ncoeff,
                            fromField->m_exp[j]->GetPlane(2)->UpdateCoeffs(),
                            1);
#else
                Array<OneD, NekDouble> tmp;
                Vmath::Smul(Ncoeff, -1., fromField->m_fieldCoeffs->GetArray1D(j)
                            + 2*Ncoeff, 1,
                            tmp = fromField->m_fieldCoeffs->UpdateArray1D(j)
                            + 3*Ncoeff, 1);
                Vmath::Zero(Ncoeff,tmp = fromField->m_fieldCoeffs->
                            UpdateArray1D(j) + 2*Ncoeff, 1);
#endif
                
            }
        }
#if EXPLISTDATA
        fromField->m_exp[j]->BwdTrans(fromField->m_exp[j]->GetCoeffs(),
                                      fromField->m_exp[j]->UpdatePhys());
#else
        Array<OneD, NekDouble> tmp;
        fromField->m_exp[j]->BwdTrans(fromField->m_fieldCoeffs->GetArray1D(j),
                              tmp = fromField->m_fieldPhys->UpdateArray1D(j));
#endif
    }

    int nq1 = m_f->m_exp[0]->GetTotPoints();

    NekDouble clamp_low = m_config["clamptolowervalue"].as<NekDouble>();
    NekDouble clamp_up  = m_config["clamptouppervalue"].as<NekDouble>();
    NekDouble def_value = m_config["defaultvalue"].as<NekDouble>();

    for (int i = 0; i < nfields; i++)
    {
#if EXPLISTDATA
        for (int j = 0; j < nq1; ++j)
        {
            m_f->m_exp[i]->UpdatePhys()[j] = def_value;
        }
#else
        Vmath::Fill(nq1,def_value,m_f->m_fieldPhys->UpdateArray1D(i),1);
#endif
    }

    Interpolator interp;
    if (m_f->m_verbose && m_f->m_comm->TreatAsRankZero())
    {
        interp.SetProgressCallback(&ProcessInterpField::PrintProgressbar, this);
    }

#if EXPLISTDATA
    interp.Interpolate(fromField->m_exp, m_f->m_exp);
#else
    interp.Interpolate(fromField->m_exp, fromField->m_fieldPhys,
                       m_f->m_exp, m_f->m_fieldPhys);
#endif

    if (m_f->m_verbose && m_f->m_comm->TreatAsRankZero())
    {
        cout << endl;
    }

    for (int i = 0; i < nfields; ++i)
    {
#if EXPLISTDATA
        for (int j = 0; j < nq1; ++j)
        {            
            if (m_f->m_exp[i]->GetPhys()[j] > clamp_up)
            {
                m_f->m_exp[i]->UpdatePhys()[j] = clamp_up;
            }
            else if (m_f->m_exp[i]->GetPhys()[j] < clamp_low)
            {
                m_f->m_exp[i]->UpdatePhys()[j] = clamp_low;
            }
        }
        m_f->m_exp[i]->FwdTransLocalElmt(m_f->m_exp[i]->GetPhys(),
                                         m_f->m_exp[i]->UpdateCoeffs());
#else
        Array<OneD, NekDouble> phys = m_f->m_fieldPhys->UpdateArray1D(i);
        for (int j = 0; j < nq1; ++j)
        {            
            if (phys[j] > clamp_up)
            {
                phys[j] = clamp_up;
            }
            else if (phys[j] < clamp_low)
            {
                phys[j] = clamp_low;
            }
        }
        m_f->m_exp[i]->FwdTransLocalElmt(phys,
                                    m_f->m_fieldCoeffs->UpdateArray1D(i));
#endif
    }
    // save field names
    m_f->m_variables = fromField->m_fielddef[0]->m_fields;
}

void ProcessInterpField::PrintProgressbar(const int position,
                                          const int goal) const
{
    LibUtilities::PrintProgressbar(position, goal, "Interpolating");
}
} // namespace FieldUtils
} // namespace Nektar
