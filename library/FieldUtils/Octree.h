#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <vector>
#include <memory>
#include <limits>

namespace Nektar
{
namespace FieldUtils
{

class octree
{
    class octant;
    using octantSharedPtr = std::shared_ptr<octant>;
    using octantWeakPtr   = std::weak_ptr<octant>;

public:
    // Empty constructor
    octree();
    // Constructor with bounds
    octree(const Array<OneD, Array<OneD, NekDouble> > &pts, int maxPts,
           const Array<OneD, NekDouble> &bounds);
    // Constructor without bounds
    octree(const Array<OneD, Array<OneD, NekDouble> > &pts, int maxPts);

    // Returns the ID of the leaf node that contains 'coords'
    int QueryNode(const Array<OneD, NekDouble> &coords,
                  int depth=std::numeric_limits<int>::max());
    // Returns the ID of the point in 'pts' closest to 'coords'
    int QueryClosest(const Array<OneD, Array<OneD, NekDouble> > &pts,
                     const Array<OneD, NekDouble> &coords, double &distance,
                     int pointInd=-1);

    // Returns the ID of the points inside the node 'nodeID'
    std::vector<int> QueryPoints(int nodeID);
    // Returns the IDs of the leaf nodes neighbouring 'nodeID'
    std::vector<int> QueryNeighbours(int nodeID);

    int QueryNpoints(int nodeID) { return m_nodes[nodeID]->GetNpoints(); }
    int QueryLocation(int nodeID) { return m_nodes[nodeID]->GetLoc(); }
    int QueryDepth(int nodeID) { return m_nodes[nodeID]->GetDepth(); }
    int QueryDelta(int nodeID) { return m_nodes[nodeID]->GetDelta(); }

    // Get some statistics about the octree
    void GetStats(int &maxPts, int &nPts, int &nNodes,
                  int &nLeaves, int &depth);

private:
    class octant : public std::enable_shared_from_this<octant>
    {
    public:
        // Empty constructor
        octant();
        // Constructor
        octant(int loc, int depth, int id, const Array<OneD, NekDouble> &bounds);
        // Copy constructor
        octant(int loc, octant &parent);

        int GetNpoints() const { return m_nPts; }
        int GetLoc() const { return m_loc; }
        int GetDepth() const { return m_depth; }
        double GetDelta() const { return m_delta; }
        int GetID() const { return m_id; }
        bool IsLeaf() const { return m_isLeaf; }
        Array<OneD, NekDouble>& GetCentre() { return m_centre; }
        Array<OneD, NekDouble>& GetBounds() { return m_bounds; }
        std::vector<int>& GetIndices() { return m_pointInd; }
        Array<OneD, octantSharedPtr>& GetChildren() { return m_children; }
        std::vector<octantWeakPtr>& GetNeighbours() { return m_neighbours; }

        void SetID(int id) { m_id = id; }

        void GetLeaves(std::vector<octantSharedPtr>& leaves);
        void SetIndices(const std::vector<int> &indices);
        void AddNeighbours(const std::vector<octantSharedPtr> &neighbours);
        void AddPoints(const Array<OneD, Array<OneD, NekDouble> > &pts,
                       const std::vector<int> &indices);
        void Subdivide(int maxPts, const Array<OneD, Array<OneD, NekDouble> > &pts,
                       std::vector<octantSharedPtr> &nodes);
        int GetLocInNode(const Array<OneD, NekDouble> &coords);

    private:
        /// Number of points in the octant
        int m_nPts;
        /// Position in the father octant (1-8)
        int m_loc;
        /// Depth in octree (root is 1)
        int m_depth;
        /// ID of the octant (index in 'm_nodes' vector)
        int m_id;
        /// Length of the octant side
        double m_delta;
        /// Coordinates of the centre of the octant
        Array<OneD, NekDouble> m_centre;
        /// Min/max coordinates of the octant (x, y, z)
        Array<OneD, NekDouble> m_bounds;
        /// Indices of the points comprised by the octant
        std::vector<int> m_pointInd;
        /// True for a leaf octant, false otherwise
        bool m_isLeaf;

        /// Vector of pointers to the children octants
        Array<OneD, octantSharedPtr> m_children;
        /// Vector of pointers to the neighbouring octants
        std::vector<octantWeakPtr> m_neighbours;
    };

    /// Max points allowed in each octant
    int m_maxPts;
    /// Number of points in the mesh
    int m_nMshPts;
    /// Number of octants in the tree
    int m_nNodes;
    /// Number of leaf nodes in the tree
    int m_nLeaves;
    /// Maximum depth of the tree
    int m_maxDepth;

    /// First node of the tree, linked to the rest
    octantSharedPtr m_root;
    /// Vector of pointers to every octant in the tree
    std::vector<octantSharedPtr> m_nodes;

    void AdvanceToStats(int nodeID);
    void SetNeighbours(int nodeID);
};
}
}