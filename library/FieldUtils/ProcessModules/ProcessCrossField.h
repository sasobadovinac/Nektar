////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessCrossField.h
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
//  Description: Post-processes cross field simulation results.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_PROCESSCROSSFIELD
#define FIELDUTILS_PROCESSCROSSFIELD

#include "../Module.h"

namespace Nektar
{
namespace FieldUtils
{
/**
 * @brief This processing module post-processes cross field simulation results.
 */
class ProcessCrossField : public ProcessModule
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<ProcessCrossField>::AllocateSharedPtr(f);
    }
    static ModuleKey className;

    ProcessCrossField(FieldSharedPtr f);
    virtual ~ProcessCrossField();

    virtual void Process(po::variables_map &vm);

    virtual std::string GetModuleName()
    {
        return "ProcessCrossField";
    }

    virtual std::string GetModuleDescription()
    {
        return "Post-processing cross field simulation results";
    }

    virtual ModulePriority GetModulePriority()
    {
        return eModifyExp;
    }

private:
    // Find all singularities and their number of branches
    vector<pair<Array<OneD, NekDouble>, int>> FindAllSingularities();
    // Find all elements on isocontour lines
    vector<set<int>> FindIsocontourElements();
    // Find singularity if any
    Array<OneD, NekDouble> FindSingularityInElmt(int id);
    // Calculate the number of branches around the singularity
    int CalculateNumberOfBranches(int id, Array<OneD, NekDouble> eta);

    // Space dimension
    int m_dim;
    // Output CSV file for streamlines
    ofstream m_csvfile;
    // Step size for streamline marching
    NekDouble m_step;
};

class Streamline
{
public:
    // Default constructor
    Streamline(FieldSharedPtr f, Array<OneD, NekDouble> singularity,
               NekDouble step, ofstream &csvfile)
        : m_dim(2), m_f(f), m_step(step), m_csvfile(csvfile)
    {
        AddPoint(singularity);
    }

    // Add point to history
    void AddPoint(Array<OneD, NekDouble> point, NekDouble angle = 0.0)
    {
        m_points.push_back(point);
        m_angles.push_back(angle);

        m_csvfile << point[0] << "," << point[1] << endl;
    }

    // Find first point of nearest branch based on angle
    void Initialise(NekDouble &angle);

    // Make one step
    bool Advance();

private:
    // Dimensions
    int m_dim;
    // Field object
    FieldSharedPtr m_f;
    // Output CSV file for streamlines
    ofstream &m_csvfile;
    // History of points
    vector<Array<OneD, NekDouble>> m_points;
    // History of directions
    vector<NekDouble> m_angles;
    // Step size
    NekDouble m_step;
    // Number of pi/2 adjustment rotations of last point
    int m_rot;
    // Negative values for u and v of last point
    pair<bool, bool> m_neg;

    // Coefficients for the Adams-Bashforth method
    static NekDouble AdamsBashforth_coeffs[4][4];
};
}
}

#endif
