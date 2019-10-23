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

        // Calculate barycentre
        for (int j = 0; j < 3; ++j)
        {
            tmpTri.centre[j] = tmpTri.v0[j]/3.0 +
                               tmpTri.v1[j]/3.0 +
                               tmpTri.v2[j]/3.0;
        }

        // Calculate surface
        Array<OneD, NekDouble> side1(3);
        Array<OneD, NekDouble> side2(3);
        Vmath::Vsub(3, tmpTri.v1, 1, tmpTri.v0, 1, side1, 1);
        Vmath::Vsub(3, tmpTri.v2, 1, tmpTri.v0, 1, side2, 1);
        tmpTri.surf = 0.5*sqrt(
                      Vmath::Dot(3, side1, side1)*Vmath::Dot(3, side2, side2) -
                      Vmath::Dot(3, side1, side2)*Vmath::Dot(3, side1, side2));
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
 * @return NekDouble
 */
NekDouble ProcessPhiFromFile::PhiFunction(double dist, double coeff)
{
    return -0.5*(std::tanh(dist/coeff)-1.0);
}

/**
 * @brief Assigns to 'phi' the values indicated by 'ShapeFunction'
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
        Exp->FwdTrans(Exp->GetPhys(), Exp->UpdateCoeffs());

        auto it = m_f->m_exp.begin() + s * (nVars + 1) + nVars;
        m_f->m_exp.insert(it, Exp);
    }
}

/**
 * @brief Assigns to 'phi' the corresponding values of Phi
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
        // Append Phi expansion to 'm_f'
        MultiRegions::ExpListSharedPtr phi;
        phi = m_f->AppendExpList(m_f->m_numHomogeneousDir);

        // Get current coords of the point
        Array<OneD, Array<OneD, NekDouble> > coords(3);
        for (int i = 0; i < 3; ++i)
        {
            coords[i] = Array<OneD, NekDouble>(nPts, 0.0);
        }
        phi->GetCoords(coords[0], coords[1], coords[2]);

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
            phi->UpdatePhys()[i] = PhiFunction(dist, m_config["scale"].as<double>());
        }

        // Update vector of expansions
        phi->FwdTrans(phi->GetPhys(), phi->UpdateCoeffs());
        auto it = m_f->m_exp.begin() + s * (nVars + 1) + nVars;
        m_f->m_exp.insert(it, phi);
    }
}

/**
 * @brief Checks if a ray traced from 'Origin' with direction 'Dvec' hits
 * the triangle defined by 'tri'. Returns the distance to the plane
 * defined by 'tri' in any case. A negative distance means that the hit
 * happened in the direction oposite that of the ray. Approach to calculate
 * the intersection point found in:
 * 
 * Fast, minimum storage ray/triangle intersection,
 * Tomas Moeller, Ben Trumbore
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
    Vmath::Vsub(3, tri.v1, 1, tri.v0, 1, E1, 1);
    Vmath::Vsub(3, tri.v2, 1, tri.v0, 1, E2, 1);

    // If det == 0, ray parallel to triangle
    Array<OneD, NekDouble> Pvec = Cross(Dvec, E2);
    double det = Vmath::Dot(3, Pvec, E1);
    double inv_det = 1.0 / det;
    if (IsZero(det))
    {
        distance = numeric_limits<double>::infinity();
        u        = numeric_limits<double>::infinity();
        v        = numeric_limits<double>::infinity();
        return false;
    }

    // Vector T and parameter u = (0.0, 1.0)
    Array<OneD, NekDouble> Tvec(3);
    Vmath::Vsub(3, Origin, 1, tri.v0, 1, Tvec, 1);
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
 * STL file. Based on the idea that the solid angle covered by the
 * 3D shape for inner points should be 4*Pi, whereas it is 0 for points
 * lying outside. Formula for the solid angle found in:
 * 
 * The solid angle of a plane triangle
 * A. Van Oosterom and J. Strackee
 * 
 * @param file
 * @param x
 * @return true
 * @return false
 */
bool ProcessPhiFromFile::IsInterior(const STLfile &file,
                                    const Array<OneD, NekDouble> &x)
{
    // Add up the solid angle covered by each triangle
    NekDouble solidAngle = 0.0;
    for (triangle tri : file.triangles)
    {
        // Relative position of triangle vertices
        Array<OneD, NekDouble> v0(3);
        Vmath::Vsub(3, tri.v0, 1, x, 1, v0, 1);
        NekDouble v0mag = sqrt(Vmath::Dot(3, v0, v0));

        Array<OneD, NekDouble> v1(3);
        Vmath::Vsub(3, tri.v1, 1, x, 1, v1, 1);
        NekDouble v1mag = sqrt(Vmath::Dot(3, v1, v1));

        Array<OneD, NekDouble> v2(3);
        Vmath::Vsub(3, tri.v2, 1, x, 1, v2, 1);
        NekDouble v2mag = sqrt(Vmath::Dot(3, v2, v2));

        // Some calculations for the solid angle formula
        Array<OneD, Array<OneD, NekDouble> > tmpMat(3);
        tmpMat[0] = v0;
        tmpMat[1] = v1;
        tmpMat[2] = v2;

        NekDouble num = Det3(tmpMat);
        NekDouble den = v0mag*v1mag*v2mag + Vmath::Dot(3, v0, v1)*v2mag +
                    Vmath::Dot(3, v0, v2)*v1mag + Vmath::Dot(3, v1, v2)*v0mag;
        
        // Solid angle
        solidAngle += 2.0*atan2(num, den);
    }

    // Low values of 'solidAngle' correspond to an EXTERIOR point
    return (solidAngle > 0.1);
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

        if (!hit)
        {
            /* No need to check this if the shape is closed */
            
            // // The minimum has to be in one of the edges
            // if (v < 0)   // Edge V0-V1
            // {
            //     tmpDist = Distance2edge(x, tri.v0, tri.v1);
            // }
            // else if (u < 0)   // Edge V0-V2
            // {
            //     tmpDist = Distance2edge(x, tri.v0, tri.v2);
            // }
            // else   // Edge V1-V2
            // {
            //     tmpDist = Distance2edge(x, tri.v1, tri.v2);
            // }
        }
        else
        {
            // Update 'dist'
            tmpDist = abs(tmpDist);
            if (tmpDist < dist)
            {
                dist = tmpDist;
            }
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
    double EPS = numeric_limits<double>::epsilon();
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
 * @brief Calculates the determinant of a 3x3 matrix
 * 
 * @param mat
 * @return NekDouble
 */
NekDouble ProcessPhiFromFile::Det3(const Array<OneD, Array<OneD, NekDouble> > &mat)
{
    return mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) -
           mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
           mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
}

/**
 * @brief Calculates the distance between two n-dimensional points
 * 
 * @param v0
 * @param v1
 * @return NekDouble
 */
NekDouble ProcessPhiFromFile::Distance2point(const Array<OneD, NekDouble> &v0,
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
 * @return NekDouble
 */
NekDouble ProcessPhiFromFile::Distance2edge(const Array<OneD, NekDouble> &x,
                                         const Array<OneD, NekDouble> &e1,
                                         const Array<OneD, NekDouble> &e2)
{
    size_t n = x.num_elements();
    Array<OneD, NekDouble> e1x(n);
    Array<OneD, NekDouble> e1e2(n);
    Vmath::Vsub(n, x, 1, e1, 1, e1x, 1);
    Vmath::Vsub(n, e2, 1, e1, 1, e1e2, 1);
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
    Vmath::Svtsvtp(n, 1.0, e1x, 1, -proj, e1e2, 1, distVec, 1);

    return sqrt(Vmath::Dot(n, distVec, distVec));
}
}
}
