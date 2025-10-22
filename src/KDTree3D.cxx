#include "../include/KDTree3D.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <stdexcept>


// Utility function to calculate squared distance
double KDTree3D::squaredDistance(const Point& a, double x, double y, double z) {
    return (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y) + (a.z - z) * (a.z - z);
}

// Constructor to build the KD-tree from a list of points
KDTree3D::KDTree3D(const std::vector<Point>& points) {
    std::vector<Point> pts(points.begin(), points.end());
    root = buildKDTree(pts);
}

// Recursive function to build the KD-tree
KDNode* KDTree3D::buildKDTree(std::vector<Point>& points, int depth) {
    if (points.empty()) return nullptr;

    int axis = depth % 3;

    std::sort(points.begin(), points.end(), [axis](const Point& a, const Point& b) {
        return (axis == 0) ? a.x < b.x : (axis == 1) ? a.y < b.y : a.z < b.z;
    });

    size_t median = points.size() / 2;
    KDNode* node = new KDNode(points[median]);

    std::vector<Point> leftPoints(points.begin(), points.begin() + median);
    std::vector<Point> rightPoints(points.begin() + median + 1, points.end());

    node->left = buildKDTree(leftPoints, depth + 1);
    node->right = buildKDTree(rightPoints, depth + 1);

    return node;
}

// Recursive function to find the nearest neighbor in the KD-tree
void KDTree3D::findNearest(KDNode* node, double x, double y, double z, int depth,
                           KDNode*& best, double& bestDistSq) const {
    if (!node) return;

    double distSq = squaredDistance(node->point, x, y, z);
    if (distSq < bestDistSq) {
        bestDistSq = distSq;
        best = node;
    }

    int axis = depth % 3;
    double diff = (axis == 0) ? x - node->point.x : (axis == 1) ? y - node->point.y : z - node->point.z;
    KDNode* near = (diff < 0) ? node->left : node->right;
    KDNode* far  = (diff < 0) ? node->right : node->left;

    findNearest(near, x, y, z, depth + 1, best, bestDistSq);

    if (diff * diff < bestDistSq) {
        findNearest(far, x, y, z, depth + 1, best, bestDistSq);
    }
}

// Public function to find the closest point to the query (x, y, z)
Point KDTree3D::find_closest_point(double x, double y, double z) const {
    KDNode* best = nullptr;
    double bestDistSq = std::numeric_limits<double>::max();
    findNearest(root, x, y, z, 0, best, bestDistSq);

    if (!best) {
        throw std::runtime_error("KD-tree search failed.");
    }

    return best->point;  // Return the closest point directly
}



void KDTree3D::findKNN(KDNode* node, double x, double y, double z, int depth, int k,
    std::priority_queue<std::pair<double,const KDNode*>>& heap) const
{
if (!node) return;

const double distSq = squaredDistance(node->point, x, y, z);
if ((int)heap.size() < k) {
heap.emplace(distSq, node);
} else if (distSq < heap.top().first) {
heap.pop();
heap.emplace(distSq, node);
}

const int axis = depth % 3;
const double diff = (axis == 0) ? x - node->point.x
  : (axis == 1) ? y - node->point.y
                : z - node->point.z;

KDNode* near = (diff < 0) ? node->left : node->right;
KDNode* far  = (diff < 0) ? node->right : node->left;

// Explore nearer side first
findKNN(near, x, y, z, depth + 1, k, heap);

// Prune with splitting plane distance
const double planeDistSq = diff * diff;
if ((int)heap.size() < k || planeDistSq < heap.top().first) {
findKNN(far, x, y, z, depth + 1, k, heap);
}
}

std::vector<Point> KDTree3D::k_closest_points(double x, double y, double z, int k) const
{
std::vector<Point> out;
if (!root || k <= 0) return out;

std::priority_queue<std::pair<double,const KDNode*>> heap; // max-heap
findKNN(root, x, y, z, 0, k, heap);

out.reserve(heap.size());
// Extract from heap (largest first) then reverse to get nearest-first
while (!heap.empty()) {
out.push_back(heap.top().second->point);
heap.pop();
}
std::reverse(out.begin(), out.end());
return out;
}

