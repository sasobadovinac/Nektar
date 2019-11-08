#include "Octree.h"
#include <stdexcept>
#include <math.h>

namespace Nektar
{
namespace FieldUtils
{

/**
 * @brief Construct a new octree object
 *
 */
octree::octree() : m_maxPts(0), m_nMshPts(0), m_nNodes(0), m_nLeaves(0),
                   m_maxDepth(0), m_root(nullptr)
{
}

/**
 * @brief Construct a new octree object
 *
 * @param pts
 * @param maxPts
 * @param bounds
 */
octree::octree(const Array<OneD, Array<OneD, NekDouble> > &pts, int maxPts,
               const Array<OneD, NekDouble> &bounds)
{
    // Set some values
    m_maxPts  = maxPts;
    m_nMshPts = pts.num_elements();

    // Create first (root) node
    std::vector<int> indices(m_nMshPts);
    for (int i = 0; i < m_nMshPts; ++i)
    {
        indices[i] = i;
    }
    m_root = std::make_shared<octant>(0, 1, 0, bounds);
    m_root->SetIndices(indices);
    m_nodes.push_back(m_root);

    // Create the tree
    m_root->Subdivide(m_maxPts, pts, m_nodes);

    // Get some data after the tree is computed
    m_maxDepth = 0;
    AdvanceToStats(m_root->GetID());
    m_nNodes = m_nodes.size();

    // Set the pointers to the neighbouring nodes
    SetNeighbours(m_root->GetID());
}

/**
 * @brief Construct a new octree object
 *
 * @param pts
 * @param maxPts
 */
octree::octree(const Array<OneD, Array<OneD, NekDouble> > &pts, int maxPts)
{
    // Find coordinates of the bounding box
    Array<OneD, NekDouble> bounds(6);
    bounds[0] = pts[0][0];
    bounds[1] = pts[0][0];
    bounds[2] = pts[0][1];
    bounds[3] = pts[0][1];
    bounds[4] = pts[0][2];
    bounds[5] = pts[0][2];
    for (int i = 1; i < m_nMshPts; ++i)
    {
        bounds[0] = (bounds[0] < pts[i][0]) ? bounds[0] : pts[i][0];
        bounds[1] = (bounds[1] > pts[i][0]) ? bounds[1] : pts[i][0];
        bounds[2] = (bounds[2] < pts[i][1]) ? bounds[2] : pts[i][1];
        bounds[3] = (bounds[3] > pts[i][1]) ? bounds[3] : pts[i][1];
        bounds[4] = (bounds[4] < pts[i][2]) ? bounds[4] : pts[i][2];
        bounds[5] = (bounds[5] > pts[i][2]) ? bounds[5] : pts[i][2];
    }

    // Add a small margin
    bounds[0] -= fabs(bounds[0])*0.01;
    bounds[1] += fabs(bounds[1])*0.01;
    bounds[2] -= fabs(bounds[2])*0.01;
    bounds[3] += fabs(bounds[3])*0.01;
    bounds[4] -= fabs(bounds[4])*0.01;
    bounds[5] += fabs(bounds[5])*0.01;

    // Call the octree constructor
    octree(pts, maxPts, bounds);
}

/**
 * @brief Given the coordinates 'coords' of a point, returns the leaf octant
 * that contains it. If 'depth' is specified, returs the node that contains the
 * point such that its depth is lower or equal to the one indicated. If the
 * point lies outside the tree, it returns -1
 *
 * @param coords
 * @param depth
 * @return int
 */
int octree::QueryNode(const Array<OneD, NekDouble> &coords, int depth)
{
    int nodeID = -1;

    if (coords.num_elements())
    {
        // First, check if inside the octree
        Array<OneD, NekDouble> bounds = m_root->GetBounds();
        if ((coords[0] >= bounds[0]) && (coords[0] <= bounds[1]) &&
            (coords[1] >= bounds[2]) && (coords[1] <= bounds[3]) &&
            (coords[2] >= bounds[4]) && (coords[2] <= bounds[5]))
        {
            // Initialise 'node'
            octantSharedPtr node = m_root;

            // Keep advancing towards the end of the branch
            while (!node->IsLeaf() && node->GetDepth() < depth)
            {
                int loc = node->GetLocInNode(coords);
                node = node->GetChildren()[loc-1];
            }

            nodeID = node->GetID();
        }
    }
    
    return nodeID;
}

/**
 * @brief Finds the ID of the closest point in 'pts' to the one specified by
 * 'coords'. It also returns the distance between both points in 'distance'
 *
 * @param pts
 * @param coords
 * @param distance
 * @param pointInd
 * @return int
 */
int octree::QueryClosest(const Array<OneD, Array<OneD, NekDouble> > &pts,
                         const Array<OneD, NekDouble> &coords,
                         double &distance, int pointInd)
{
    int index = -1;
    distance  = std::numeric_limits<double>::max();

    if (coords.num_elements())
    {
        // Find the corresponding node to 'coords'
        int nodeInd;
        if (pointInd > 0)
        {
            nodeInd = pointInd;
        }   
        else
        {
            nodeInd = QueryNode(coords);
        }

        // List the indices of all the candidate points
        std::vector<int> indices(m_nodes[nodeInd]->GetIndices());
        for (octantWeakPtr neigh : m_nodes[nodeInd]->GetNeighbours())
        {
            for (int i : neigh.lock()->GetIndices())
            {
                indices.push_back(i);
            }
        }

        // Check the distances with all the nodes
        for (int i : indices)
        {
            double sub = pts[i][0]-coords[0];
            double tmpDistance = sub * sub;
            for (int j = 1; j < 3; ++j)
            {
                sub = pts[i][j]-coords[j];
                tmpDistance += sub * sub;
            }
            tmpDistance = std::sqrt(tmpDistance);
            
            if (distance > tmpDistance)
            {
                distance = tmpDistance;
                index = i;
            }
        }
    }

    return index;
}

/**
 * @brief Returns the indices of the points of the mesh contained in the tree
 *
 * @param nodeID
 * @return std::vector<int>
 */
std::vector<int> octree::QueryPoints(int nodeID)
{
    return m_nodes[nodeID]->GetIndices();
}

/**
 * @brief Returns the IDs of the octants that surround the queried node. First,
 * it finds the neighbouring nodes with the same or lower depth and, if they
 * are not leaf nodes, return all the leaf octants contained in them. This
 * means that, for octants of depth 2, this function returns all the leaf nodes
 * in the tree except those lying inside the queried octant
 *
 * @param nodeID
 * @return std::vector<int>
 */
std::vector<int> octree::QueryNeighbours(int nodeID)
{
    std::vector<int> indices;
    for (const octantWeakPtr node : m_nodes[nodeID]->GetNeighbours())
    {
        indices.push_back(node.lock()->GetID());
    }
    return indices;
}

/**
 * @brief Returns some characteristic values of the tree.
 *
 * @param maxPts
 * @param nPts
 * @param nNodes
 * @param nLeaves
 * @param depth
 */
void octree::GetStats(int &maxPts, int &nPts, int &nNodes,
                      int &nLeaves, int &depth)
{
    maxPts  = m_maxPts;
    nPts    = m_nMshPts;
    nNodes  = m_nNodes;
    nLeaves = m_nLeaves;
    depth   = m_maxDepth;
}

/**
 * @brief Goes through all the nodes of the octree counting the number of
 * octants and the maximum depth reached
 * 
 * @param nodeID 
 */
void octree::AdvanceToStats(int nodeID)
{
    // Update stats if we reached the end of the branch
    if (m_nodes[nodeID]->IsLeaf())
    {
        m_nLeaves++;
        m_maxDepth = (m_maxDepth > m_nodes[nodeID]->GetDepth()) ? m_maxDepth :
                                                m_nodes[nodeID]->GetDepth();
    }
    // In any other case, dig into the tree
    else
    {
        Array<OneD, octantSharedPtr> children = m_nodes[nodeID]->GetChildren();
        for (octantSharedPtr child : children)
        {
            AdvanceToStats(child->GetID());
        }
    }
}

/**
 * @brief Once the nodes of the octree are created, sets their neighbours as
 * explained in 'octree::QueryNeighbours'
 *
 * @param nodeID
 */
void octree::SetNeighbours(int nodeID)
{
    // Array with the different steps
    static int steps[26][3] = {{-1,-1,-1}, {1,0,0}, {1,0,0}, {0,1,0}, {0,1,0},
                               {-1,0,0}, {-1,0,0}, {0,-1,0}, {1,0,0},
                               {-1,-1,1}, {1,0,0}, {1,0,0}, {0,1,0}, {0,1,0},
                               {-1,0,0}, {-1,0,0}, {0,-1,0}, {0,-1,1}, {1,0,0},
                               {1,0,0}, {0,1,0}, {0,1,0}, {-1,0,0}, {-1,0,0},
                               {0,-1,0}, {1,0,0}};

    // Advance to the leaves of the octree
    if (!m_nodes[nodeID]->IsLeaf())
    {
        for (octantSharedPtr child : m_nodes[nodeID]->GetChildren())
        {
            SetNeighbours(child->GetID());
        }
    }
    else
    {
        // delta * steps
        Array<OneD, NekDouble> probeCoords(m_nodes[nodeID]->GetCentre());
        for (int step = 0; step < 26; ++step)
        {
            for (int i = 0; i < 3; ++i)
            {
                probeCoords[i] += m_nodes[nodeID]->GetDelta()*steps[step][i];
            }

            // For each neighbour, find the leaves and add them
            int neighInd = QueryNode(probeCoords, m_nodes[nodeID]->GetDepth());
            if (neighInd > -1)
            {
                std::vector<octantSharedPtr> leaves;
                m_nodes[neighInd]->GetLeaves(leaves);
                m_nodes[nodeID]->AddNeighbours(leaves);
            }
        }
    }
}

/**
 * @brief Construct a new octree::octant object
 *
 */
octree::octant::octant() : m_nPts(-1), m_loc(-1), m_depth(-1), m_id(-1),
                           m_delta(-1), m_centre(3), m_bounds(6), m_isLeaf(true)
{
}

/**
 * @brief Construct a new octree::octant object
 *
 * @param loc
 * @param depth
 * @param id
 * @param bounds
 */
octree::octant::octant(int loc, int depth, int id,
                       const Array<OneD, NekDouble> &bounds) :
                            m_nPts(0), m_loc(loc), m_depth(depth),
                            m_id(id), m_isLeaf(true)
{
    // Check the size of 'bounds'
    if (bounds.num_elements() != 6)
    {
        throw std::out_of_range("Size of bounds must be 6.");
    }

    // If all deltas are not equal, use the largest ones
    double deltaX = bounds[1] - bounds[0];
    double deltaY = bounds[3] - bounds[2];
    double deltaZ = bounds[5] - bounds[4];
    if (deltaX != deltaY || deltaY != deltaZ)
    {
        m_delta = (deltaX > deltaY) ? deltaX : deltaY;
        m_delta = (m_delta > deltaZ) ? m_delta : deltaZ;
    }
    else
    {
        m_delta = deltaX;
    }

    // Fill in the rest of the data
    m_centre = Array<OneD, NekDouble>(3);
    m_bounds = Array<OneD, NekDouble>(6);
    for (int i = 0; i < 3; ++i)
    {
        m_centre[i]     = (bounds[2*i+1] + bounds[2*i])/2.0;
        m_bounds[2*i]   = m_centre[i] - m_delta/2.0;
        m_bounds[2*i+1] = m_centre[i] + m_delta/2.0;
    }
}

/**
 * @brief Construct a new octree::octant object
 *
 * @param loc
 * @param parent
 */
octree::octant::octant(int loc, octant &parent) : m_nPts(0), m_loc(loc),
                                                  m_id(-1), m_isLeaf(true)
{
    // Set depth
    m_depth = parent.GetDepth() + 1;

    // Set delta
    m_delta = parent.GetDelta()/2.0;

    // Set centre
    double centreDX;
    double centreDY;
    double centreDZ;
    switch (loc)
    {
        case 1:  // x-, y-, z-
            centreDX = -m_delta/2.0;
            centreDY = -m_delta/2.0;
            centreDZ = -m_delta/2.0;
            break;
        case 2:  // x+, y-, z-
            centreDX =  m_delta/2.0;
            centreDY = -m_delta/2.0;
            centreDZ = -m_delta/2.0;
            break;
        case 3:  // x+, y+, z-
            centreDX =  m_delta/2.0;
            centreDY =  m_delta/2.0;
            centreDZ = -m_delta/2.0;
            break;
        case 4:  // x-, y+, z-
            centreDX = -m_delta/2.0;
            centreDY =  m_delta/2.0;
            centreDZ = -m_delta/2.0;
            break;
        case 5:  // x-, y-, z+
            centreDX = -m_delta/2.0;
            centreDY = -m_delta/2.0;
            centreDZ =  m_delta/2.0;
            break;
        case 6:  // x+, y-, z+
            centreDX =  m_delta/2.0;
            centreDY = -m_delta/2.0;
            centreDZ =  m_delta/2.0;
            break;
        case 7:  // x+, y+, z+
            centreDX =  m_delta/2.0;
            centreDY =  m_delta/2.0;
            centreDZ =  m_delta/2.0;
            break;
        case 8:  // x-, y+, z+
            centreDX = -m_delta/2.0;
            centreDY =  m_delta/2.0;
            centreDZ =  m_delta/2.0;
            break;
        default:
            throw std::out_of_range("Loc must be in the range (1,8).");
    }
    Array<OneD, NekDouble> pCentre = parent.GetCentre();
    m_centre = Array<OneD, NekDouble>(3);
    m_centre[0] = pCentre[0] + centreDX;
    m_centre[1] = pCentre[1] + centreDY;
    m_centre[2] = pCentre[2] + centreDZ;

    // Set bounds
    m_bounds = Array<OneD, NekDouble>(6);
    for (int i = 0; i < 3; ++i)
    {
        m_bounds[2*i]   = m_centre[i] - m_delta/2.0;
        m_bounds[2*i+1] = m_centre[i] + m_delta/2.0;
    }
}

/**
 * @brief Updates 'leaves' so that it contains all the leaf nodes belonging to
 * the octant
 *
 * @param leaves
 */
void octree::octant::GetLeaves(std::vector<octantSharedPtr>& leaves)
{
    if (m_isLeaf)
    {
        leaves.push_back(shared_from_this());
    }
    else
    {
        for (octantSharedPtr child : m_children)
        {
            child->GetLeaves(leaves);
        }
    }
}

/**
 * @brief Sets the values of 'm_pointInd' to those in 'indices'
 *
 * @param indices
 */
void octree::octant::SetIndices(const std::vector<int> &indices)
{
    for (int i : indices)
    {
        m_pointInd.push_back(i);
    }
    m_nPts = indices.size();
}

/**
 * @brief Adds to 'm_neighbours' the octants that are not already in the list
 *
 * @param neighbours
 */
void octree::octant::AddNeighbours(const std::vector<octantSharedPtr> &neighbours)
{
    for (const octantSharedPtr neighbour : neighbours)
    {
        bool equal = false;
        for (const octantWeakPtr neigh: m_neighbours)
        {
            if (neigh.lock()->GetID() == neighbour->GetID())
            {
                equal = true;
                break;
            }
        }
        if (!equal)
        {
            m_neighbours.push_back(neighbour);
        }
    }
}

/**
 * @brief Adds to 'm_pointInd' the IDs of the points in 'pts' that fall inside
 * the octant
 *
 * @param pts
 * @param indices
 */
void octree::octant::AddPoints(const Array<OneD, Array<OneD, NekDouble> > &pts,
                               const std::vector<int> &indices)
{
    for (int i : indices)
    {
        // Check if the point is inside the node
        Array<OneD, NekDouble> pt = pts[i];
        if ((pt[0] < m_bounds[0]) || (pt[0] > m_bounds[1]))
        {
            continue;
        }
        if ((pt[1] < m_bounds[2]) || (pt[1] > m_bounds[3]))
        {
            continue;
        }
        if ((pt[2] < m_bounds[4]) || (pt[2] > m_bounds[5]))
        {
            continue;
        }

        // If so, add it to the list
        m_nPts++;
        m_pointInd.push_back(i);

        // Flag it as a leaf node
        m_isLeaf = true;
    }
}

/**
 * @brief Recursively divides the octant into 8 children and fills the leaf
 * nodes with their corresponding points. Does NOT add neighbours
 *
 * @param maxPts
 * @param pts
 * @param nodes
 */
void octree::octant::Subdivide(int maxPts,
                               const Array<OneD, Array<OneD, NekDouble> > &pts,
                               std::vector<octantSharedPtr> &nodes)
{
    // For a non-leaf node
    if (m_nPts > maxPts)
    {
        // Create and fill children RECURSIVELY
        m_children = Array<OneD, octantSharedPtr>(8);
        for (int i = 0; i < 8; ++i)
        {
            octantSharedPtr newChild =
                std::make_shared<octant>(i+1, *shared_from_this());
            newChild->AddPoints(pts, m_pointInd);
            newChild->SetID(nodes.size());  // ID's start from 0

            // Add it to the list
            m_children[i] = newChild;
            nodes.push_back(newChild);

            // Keep dividing
            newChild->Subdivide(maxPts, pts, nodes);  // Recursion
        }

        // Not a leaf node anymore
        m_pointInd.clear();
        m_nPts   = 0;
        m_isLeaf = false;
    }    
}

/**
 * @brief Returns the position inside an octant in the range (1-8). The name
 * convention is as follows: let \f[\Delta\f] be half the length of the
 * father octant side, and \f[x_c,y_c,z_c\f] be the coordinates of the centre
 * of the father node. Then, position 1 corresponds to \f[x=x_c-\Delta\f],
 * \f[y=y_c-\Delta\f] and \f[y=y_c-\Delta\f]. The next positions are obtained
 * by rotating counter-clockwise around the Z axis and then, making the same
 * rotation for \f[z=z_c+\Delta\f]
 *
 * @param coords
 * @return int
 */
int octree::octant::GetLocInNode(const Array<OneD, NekDouble> &coords)
{
    // Different positions as bits in 'posByte'
    // MSB <==> LSB
    unsigned char posByte;

    if (coords[0] <= m_centre[0])  // x-
    {
        posByte = 153;   //0b10011001;
    }
    else                           // x+
    {
        posByte = 102;   //0b01100110;
    }
    if (coords[1] <= m_centre[1])  // y-
    {
        posByte &= 51;   //0b00110011;
    }
    else                           // y+
    {
        posByte &= 204;  //0b11001100;
    }
    if (coords[2] <= m_centre[2])  // z-
    {
        posByte &= 15;   //0b00001111;
    }
    else                           // z+
    {
        posByte &= 240;  //0b11110000;
    }

    // Transform into a position in the range (1,8)
    int position = 1;
    while (posByte > 1)
    {
        posByte = posByte >> 1;
        position++;
    }

    return position;
}
}
}