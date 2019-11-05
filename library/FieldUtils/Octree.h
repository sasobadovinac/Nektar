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
    octree();
    octree(const Array<OneD, Array<OneD, NekDouble> > &pts, int maxPts,
           const Array<OneD, NekDouble> &bounds);
    octree(const Array<OneD, Array<OneD, NekDouble> > &pts, int maxPts);

    int QueryNode(const Array<OneD, NekDouble> &coords,
                  int depth=std::numeric_limits<int>::max());
    int QueryClosest(const Array<OneD, Array<OneD, NekDouble> > &pts,
                     const Array<OneD, NekDouble> &coords, double &distance,
                     int pointInd=-1);

    std::vector<int> QueryPoints(int nodeID);
    std::vector<int> QueryNeighbours(int nodeID);

    int QueryNpoints(int nodeID) { return m_nodes[nodeID]->GetNpoints(); }
    int QueryLocation(int nodeID) { return m_nodes[nodeID]->GetLoc(); }
    int QueryDepth(int nodeID) { return m_nodes[nodeID]->GetDepth(); }
    int QueryDelta(int nodeID) { return m_nodes[nodeID]->GetDelta(); }

    void GetStats(int &maxPts, int &nPts, int &nNodes,
                  int &nLeaves, int &depth);

private:
    class octant : public std::enable_shared_from_this<octant>
    {
    public:
        octant();
        octant(int loc, int depth, int id, const Array<OneD, NekDouble> &bounds);
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
        int m_nPts;
        int m_loc;
        int m_depth;
        int m_id;
        double m_delta;
        Array<OneD, NekDouble> m_centre;
        Array<OneD, NekDouble> m_bounds;  // x, y, z
        std::vector<int> m_pointInd;
        bool m_isLeaf;

        Array<OneD, octantSharedPtr> m_children;
        std::vector<octantWeakPtr> m_neighbours;
    };

    int m_maxPts;
    int m_nMshPts;
    int m_nNodes;
    int m_nLeaves;
    int m_maxDepth;

    octantSharedPtr m_root;  // First node of the tree, linked to the rest
    std::vector<octantSharedPtr> m_nodes;

    void AdvanceToStats(int nodeID);
    void SetNeighbours(int nodeID);
};
}
}