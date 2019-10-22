////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessPhiFromFile.cpp
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
//  Description: Reads an STL file.
//
////////////////////////////////////////////////////////////////////////////////

#include "ProcessPhiFromFile.h"
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>

using namespace Nektar;
using namespace std;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessPhiFromFile::m_className = {
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "phifile"),
        ProcessPhiFromFile::create,
        "Computes the Phi function from a file, used in IB methods.")
};

/**
 * @brief Set up ProcessPhiFromFile object.
 *
 */
ProcessPhiFromFile::ProcessPhiFromFile(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["scale"] = ConfigOption(false, "NotSet",
                                     "Scale coefficient for Phi.");
    m_config["file"] = ConfigOption(false, "NotSet",
                                    "File with the IB definition.");
}

/**
 *
 */
ProcessPhiFromFile::~ProcessPhiFromFile()
{
}

/**
 *
 */
void ProcessPhiFromFile::Process(po::variables_map &vm)
{
    // Check if required params are defined
    ASSERTL0(m_f->m_graph, "A session file file must be provided before the "
                           "STL file.");

    ASSERTL0(m_config["scale"].as<string>().compare("NotSet") != 0,
             "Need to specify a scale coefficient, scale=value");

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    // Setup the random number generator
    std::random_device dev;
    m_rng   = mt19937(dev());
    m_uDist = uniform_real_distribution<NekDouble>(-1.0, 1.0);
    
    // Read Phi function from the session file...
    if (m_config["file"].as<string>().compare("NotSet") == 0)
    {
        GetPhifromSession();
    }
    // ...or Read STL file and append Phi values to the existing expansions
    else
    {
        STLfile phiFile = ReadSTL(m_config["file"].as<string>());
        GetPhifromSTL(phiFile);
    }
}

/**
 * @brief Read one 3D vector from a STL file, starting from the next line
 * of the input 'ifstream'. Numbers in ifstream are defined as 'float'
 * 
 * @param in
 * @return Array<OneD, NekDouble>
 */
Array<OneD, NekDouble> ProcessPhiFromFile::ReadVector(ifstream &in)
{
    Array<OneD, NekDouble> out(3);
    char buf[4];

    in.read(buf, 4);
    out[0] = *((float*) buf);
    in.read(buf, 4);
    out[1] = *((float*) buf);
    in.read(buf, 4);
    out[2] = *((float*) buf);

    return out;
}

/**
 * @brief Read an STL binary file and returns a struct of type 'STLfile'
 * containing the parsed data
 * 
 * @param filename
 * @return ProcessPhiFromFile::STLfile
 */
ProcessPhiFromFile::STLfile ProcessPhiFromFile::ReadSTL(string filename)
{
    STLfile out;

    // Open file
    ifstream fileStl(filename.c_str(), ios::binary);
    ASSERTL0(fileStl, "An error occurred while trying to open the STL file.")

    // Buffers
    char headerBuf[80];
    char numTriBuf[4];
    char dumpBuf[2];

    // Read header and num of triangles
    fileStl.read(headerBuf, 80);
    fileStl.read(numTriBuf, 4);
    unsigned int numTri = *((unsigned int*) numTriBuf);

    out.header = headerBuf;
    out.numTri = numTri;

    // Read triangle data
    out.triangles = Array<OneD, triangle>(numTri);
    for (uint i = 0; i < numTri; ++i)
    {
        // Read normal vector
        triangle tmpTri;
        tmpTri.normal = ReadVector(fileStl);

        // Read three vertices
        tmpTri.v0 = ReadVector(fileStl);
        tmpTri.v1 = ReadVector(fileStl);
        tmpTri.v2 = ReadVector(fileStl);
        out.triangles[i] = tmpTri;

        // Dump triangle type
        fileStl.read(dumpBuf, 2);
    }

    // Close the file
    fileStl.close();

    return out;
}

/**
 * @brief Smoothing function for the SPM method given a distance value
 * and a scaling coefficient
 * 
 * @param dist
 * @param coeff
 * @return double
 */
double ProcessPhiFromFile::PhiFunction(double dist, double coeff)
{
    return -0.5*(std::tanh(dist/coeff)-1.0);
}

/**
 * @brief 
 * 
 */
void ProcessPhiFromFile::GetPhifromSession()
{
    // Check that 'ShapeFunction' appears in the session file
    ASSERTL0(m_f->m_session->DefinesFunction("ShapeFunction"),
             "If file=file.stl is not supplied as an argument, a "
             "'ShapeFunction' block must be provided in the session "
             "file.")

    // Phi function in session file
    LibUtilities::EquationSharedPtr phiFunction =
        m_f->m_session->GetFunction("ShapeFunction", "Phi");

    // Get info about the domain
    int nPts  = m_f->m_exp[0]->GetNpoints();
    int nVars = m_f->m_variables.size();
    int nStrips;
    m_f->m_session->LoadParameter("Strip_Z", nStrips, 1);

    // Add new variable
    m_f->m_variables.push_back("phi");

    for (int s = 0; s < nStrips; ++s)
    {
        // Get current coords of the point
        Array<OneD, Array<OneD, NekDouble> > coords(3);
        for (int i = 0; i < 3; ++i)
        {
            coords[i] = Array<OneD, NekDouble>(nPts, 0.0);
        }
        m_f->m_exp[s*nVars]->GetCoords(coords[0], coords[1], coords[2]);

        // Append Phi expansion to 'm_f'
        MultiRegions::ExpListSharedPtr Exp;
        Exp = m_f->AppendExpList(m_f->m_numHomogeneousDir);
        phiFunction->Evaluate(coords[0], coords[1], coords[2],
                              Exp->UpdatePhys());
        Exp->FwdTrans_IterPerExp(Exp->GetPhys(), Exp->UpdateCoeffs());

        auto it = m_f->m_exp.begin() + s * (nVars + 1) + nVars;
        m_f->m_exp.insert(it, Exp);
    }
}

/**
 * @brief Assigns to 'm_phi' the corresponding values of Phi
 * 
 * @param file
 */
void ProcessPhiFromFile::GetPhifromSTL(const ProcessPhiFromFile::STLfile &file)
{
    // Get info about the domain
    int nPts  = m_f->m_exp[0]->GetNpoints();
    int nVars = m_f->m_variables.size();

    // Add new variable
    m_f->m_variables.push_back("phi");

    // Number of homogeneous strips
    int nStrips;
    m_f->m_session->LoadParameter("Strip_Z", nStrips, 1);
    
    for (int s = 0; s < nStrips; ++s)
    {
        // Phi array allocation
        Array<OneD, NekDouble> phi(nPts);

        // Get current coords of the point
        Array<OneD, Array<OneD, NekDouble> > coords(3);
        for (int i = 0; i < 3; ++i)
        {
            coords[i] = Array<OneD, NekDouble>(nPts, 0.0);
        }
        m_f->m_exp[s*nVars]->GetCoords(coords[0], coords[1], coords[2]);

        // Parallelisation is highly recommended here
        for (int i = 0; i < nPts; ++i)
        {
            // Get coordinates of each point
            Array<OneD, NekDouble> tmpCoords(3);
            tmpCoords[2] = coords[2][i];
            tmpCoords[1] = coords[1][i];
            tmpCoords[0] = coords[0][i];

            // Find the shortest distance to the body(ies)
            double dist;
            bool inside = IsInterior(file, tmpCoords);
            FindShortestDist(file, tmpCoords, dist);
            if (inside)
            {
                dist = -dist;
            }

            // Get corresponding value of Phi
            phi[i] = PhiFunction(dist, m_config["scale"].as<double>());
        }

        // Append Phi expansion to 'm_f'
        MultiRegions::ExpListSharedPtr Exp;
        Exp = m_f->AppendExpList(m_f->m_numHomogeneousDir);
        Vmath::Vcopy(nPts, phi, 1, Exp->UpdatePhys(), 1);
        Exp->FwdTrans_IterPerExp(phi, Exp->UpdateCoeffs());

        auto it = m_f->m_exp.begin() + s * (nVars + 1) + nVars;
        m_f->m_exp.insert(it, Exp);
    }
}

/**
 * @brief Checks if a ray traced from 'Origin' with direction 'Dvec' hits
 * the triangle defined by 'tri'. Returns the distance to the plane
 * defined by 'tri' in any case. A negative distance means that the hit
 * happend in the direction oposite that of the ray
 * 
 * @param tri
 * @param Origin
 * @param Dvec
 * @param distance
 * @param u
 * @param v
 * @return true
 * @return false
 */
bool ProcessPhiFromFile::CheckHit(const ProcessPhiFromFile::triangle &tri,
                                  const Array<OneD, NekDouble> &Origin,
                                  const Array<OneD, NekDouble> &Dvec,
                                  double &distance, double &u, double &v)
                                
{
    // Edge vectors
    Array<OneD, NekDouble> E1(3);
    Array<OneD, NekDouble> E2(3);
    for (int i = 0; i < 3; ++i)
    {
        E1[i] = tri.v1[i]-tri.v0[i];
        E2[i] = tri.v2[i]-tri.v0[i];
    }

    // If det == 0, ray parallel to triangle
    Array<OneD, NekDouble> Pvec = Cross(Dvec, E2);
    double det = Vmath::Dot(3, Pvec, E1);
    double inv_det = 1.0 / det;
    if (IsZero(det))
    {
        distance = std::numeric_limits<double>::infinity();
        u        = std::numeric_limits<double>::infinity();
        v        = std::numeric_limits<double>::infinity();
        return false;
    }

    // Vector T and parameter u = (0.0, 1.0)
    Array<OneD, NekDouble> Tvec(3);
    for (int i = 0; i < 3; ++i)
    {
        Tvec[i] = Origin[i]-tri.v0[i];
    }
    u = Vmath::Dot(3, Pvec, Tvec) * inv_det;

    // Vector Q and parameter v = (0.0, 1.0)
    Array<OneD, NekDouble> Qvec = Cross(Tvec, E1);
    v = Vmath::Dot(3, Qvec, Dvec) * inv_det;

    // There is a hit if (u,v) coordinates are bounded
    distance = Vmath::Dot(3, Qvec, E2) * inv_det;
    if ((u < 0.0 || u > 1.0) || (v < 0.0 || u+v > 1.0))
    {
        return false;
    }
    else
    {
        return true;
    }
}

/**
 * @brief Returns true if a point is inside the 3D object defined in the
 * STL file. It is based in the idea that a ray traced from inside a closed
 * surface will go through an odd number of surfaces, no matter how complex
 * the geometry is
 * 
 * @param file
 * @param x
 * @return true
 * @return false
 */
bool ProcessPhiFromFile::IsInterior(const STLfile &file,
                                    const Array<OneD, NekDouble> &x)
{
    // Choose a random direction
    Array<OneD, NekDouble> dir(3);
    dir[0]   = m_uDist(m_rng);
    dir[1]   = m_uDist(m_rng);
    dir[2]   = sqrt(1.0 - dir[0]*dir[0] - dir[1]*dir[1]);
    int hits = 0;

    // Stores the distances of the hits with each surface
    // It has to be a dynamic container, it will always be small
    vector<NekDouble> distVec;

    // Check hits with all the triangles
    for (triangle tri : file.triangles)
    {
        double dist;
        double u, v;
        bool hit = CheckHit(tri, x, dir, dist, u, v);

        if (hit && dist > 0.0 &&
            std::find_if(distVec.begin(), distVec.end(),
                [&](double x){ return IsZero(x-dist); }) == distVec.end())
        {
            distVec.push_back(dist);
            hits++;
        }
    }

    // Odd number of hits -> the point lies INSIDE
    if (hits % 2)
    {
        return true;
    }
    // Otherwise, it falls OUTSIDE
    else
    {
        return false;
    }
}

/**
 * @brief Calculates the shortest distance from a point \f[x\f] to the closed
 * body contained in the STL file
 * 
 * @param file
 * @param x
 * @param dist
 */
void ProcessPhiFromFile::FindShortestDist(
                                const ProcessPhiFromFile::STLfile &file,
                                const Array<OneD, NekDouble> &x,
                                double &dist)
{
    // Set 'dist' to an unreal value
    dist = numeric_limits<double>::infinity();

    for (triangle tri : file.triangles)
    {
        double tmpDist;
        double u, v;
        bool hit = CheckHit(tri, x, tri.normal, tmpDist, u, v);
        tmpDist  = abs(tmpDist);

        if (!hit)
        {
            // The minimum has to be in one of the edges
            if (v < 0)   // Edge V0-V1
            {
                tmpDist = Distance2edge(x, tri.v0, tri.v1);
            }
            else if (u < 0)   // Edge V0-V2
            {
                tmpDist = Distance2edge(x, tri.v0, tri.v2);
            }
            else   // Edge V1-V2
            {
                tmpDist = Distance2edge(x, tri.v1, tri.v2);
            }
        }

        // Update 'dist'
        if (tmpDist < dist)
        {
            dist = tmpDist;
        }
    }
}

/**
 * @brief Returns true if the argument is CLOSE to zero. Tuned for
 * the STL parsing module, do not use for other purposes
 * 
 * @param x
 * @return true
 * @return false
 */
bool ProcessPhiFromFile::IsZero(double x)
{
    double EPS = 0.000001;
    return (x > -EPS && x < EPS);
}

Array<OneD, NekDouble> ProcessPhiFromFile::Cross(
                                const Array<OneD, NekDouble> &v0,
                                const Array<OneD, NekDouble> &v1)
{
    Array<OneD, NekDouble> out(3);

    out[0] = v0[1]*v1[2] - v0[2]*v1[1];
    out[1] = v0[2]*v1[0] - v0[0]*v1[2];
    out[2] = v0[0]*v1[1] - v0[1]*v1[0];

    return out;
}

/**
 * @brief Calculates the distance between two n-dimensional points
 * 
 * @param v0
 * @param v1
 * @return double
 */
double ProcessPhiFromFile::Distance2point(
                                const Array<OneD, NekDouble> &v0,
                                const Array<OneD, NekDouble> &v1)
{
    size_t n   = v0.num_elements();
    double out = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        out += (v1[i]-v0[i]) * (v1[i]-v0[i]);
    }
    
    return sqrt(out);
}

/**
 * @brief Determines the shortest distance from a point 'x' to the segment
 * defined by the points 'e1' and 'e2'. Note that this distance may be
 * equal to that to one of the end points
 * 
 * @param x
 * @param e1
 * @param e2
 * @return double
 */
double ProcessPhiFromFile::Distance2edge(
                                const Array<OneD, NekDouble> &x,
                                const Array<OneD, NekDouble> &e1,
                                const Array<OneD, NekDouble> &e2)
{
    size_t n = x.num_elements();
    Array<OneD, NekDouble> e1x(n);
    Array<OneD, NekDouble> e1e2(n);
    for (size_t i = 0; i < n; ++i)
    {
        e1x[i]  = x[i]-e1[i];
        e1e2[i] = e2[i]-e1[i];
    }
    double norm = sqrt(Vmath::Dot(n, e1e2, e1e2));
    for (size_t i = 0; i < n; ++i)
    {
        e1e2[i] /= norm;
    }

    double proj = Vmath::Dot(n, e1x, e1e2);
    if (proj < 0.0)
    {
        return Distance2point(x, e1);
    }
    else if (proj > norm)
    {
        return Distance2point(x, e2);
    }

    Array<OneD, NekDouble> distVec(n);
    for (size_t i = 0; i < n; ++i)
    {
        distVec[i] = e1x[i]-proj*e1e2[i];
    }

    return sqrt(Vmath::Dot(n, distVec, distVec));
}
}
}
