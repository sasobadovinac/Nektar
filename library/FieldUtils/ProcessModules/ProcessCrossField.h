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
 * @brief This class handles streamline tracing.
 */
class Streamline
{
public:
    // Default constructor
    Streamline(FieldSharedPtr f, Array<OneD, NekDouble> singularity,
               NekDouble step, bool splitter = false)
        : m_dim(2), m_f(f), m_step(step), m_locked(false), m_splitter(splitter)
    {
        AddPoint(singularity);
    }

    // Merged streamline constructor
    Streamline(vector<Array<OneD, NekDouble>> points)
        : m_dim(2), m_step(0.0), m_locked(true), m_splitter(false)
    {
        for (auto &it : points)
        {
            AddPoint(it);
        }
    }

    // Add point to history
    void AddPoint(Array<OneD, NekDouble> point, NekDouble angle = 0.0)
    {
        m_points.push_back(point);
        m_angles.push_back(angle);
    }

    const Array<OneD, NekDouble> &GetFirstPoint()
    {
        return m_points.front();
    }

    const Array<OneD, NekDouble> &GetLastPoint()
    {
        return m_points.back();
    }

    const vector<Array<OneD, NekDouble>> &GetAllPoints()
    {
        return m_points;
    }

    const NekDouble &GetLastAngle()
    {
        return m_angles.back();
    }

    bool IsLocked()
    {
        return m_locked;
    }

    bool IsSplitter()
    {
        return m_splitter;
    }

    // Find first point of nearest branch based on angle
    void Initialise(NekDouble &angle);

    // Make one step
    bool Advance();

    // Check if it should be merged with another streamline
    int CheckMerge(Streamline sl, NekDouble tol);

    // Merge the 2 streamlines and return the new one
    Streamline MergeWith(Streamline sl);

    // Convert splitter into normal streamline and return the other 2 branches
    vector<Streamline> ConvertSplitter(int n);

    // Write out all points to a CSV file for use in NekMesh
    void WritePoints(ofstream &csvfile);

private:
    // Dimensions
    int m_dim;
    // Field object
    FieldSharedPtr m_f;
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
    // Reached out of domain
    bool m_locked;
    // Streamline used for splitting a collapsed quad
    bool m_splitter;

    // Coefficients for the Adams-Bashforth method
    static NekDouble AdamsBashforth_coeffs[4][4];
};

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
    // Initialise process
    void Initialise();
    // Find all singularities and their number of quadrants
    vector<pair<Array<OneD, NekDouble>, int>> FindAllSingularities();
    // Analyse vertices and find their number of quadrants
    vector<pair<Array<OneD, NekDouble>, int>> AnalyseVertices();
    // Create each streamline starting from singularities and initialise them
    vector<Streamline> CreateStreamlines(
        vector<pair<Array<OneD, NekDouble>, int>> &vertices,
        vector<pair<Array<OneD, NekDouble>, int>> &singularities);
    // Advance streamlines and merge them if necessary
    void AdvanceStreamlines(vector<Streamline> &sls);
    // Finalise process (incl. output)
    void Finalise(vector<Streamline> &sls);

    // Find all elements on isocontour lines
    vector<set<int>> FindIsocontourElements();
    // Find singularity if any
    Array<OneD, NekDouble> FindSingularityInElmt(int id);
    // Calculate the number of branches around the singularity
    int CalculateNumberOfQuadrants(int id, Array<OneD, NekDouble> eta);

    // Space dimension
    int m_dim;
    // Output CSV file for streamlines
    ofstream m_csvfile;
    // Step size for streamline marching
    NekDouble m_step;
    // Tolerance for merging streamlines
    NekDouble m_mergeTol;
};
}
}

#endif
