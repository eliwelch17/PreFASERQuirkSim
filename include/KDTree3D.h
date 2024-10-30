#ifndef KDTREE3D_H
#define KDTREE3D_H

#include <vector>

// Structure to represent a 3D point with force components
struct Point {
    double x, y, z;
    double Fx, Fy, Fz;
    Point(double x_, double y_, double z_, double Fx_, double Fy_, double Fz_)
        : x(x_), y(y_), z(z_), Fx(Fx_), Fy(Fy_), Fz(Fz_) {}
};

// Node structure for KD-tree
struct KDNode {
    Point point;
    KDNode* left;
    KDNode* right;

    KDNode(const Point& pt) : point(pt), left(nullptr), right(nullptr) {}
};

// KD-tree class for building and querying 3D points
class KDTree3D {
public:
    // Constructor to build a KD-tree from a list of points
    KDTree3D(const std::vector<Point>& points);

    // Function to find the closest point to the query point (x, y, z)
    Point find_closest_point(double x, double y, double z) const;

private:
    KDNode* root;

    // Recursive function to build the KD-tree
    KDNode* buildKDTree(std::vector<Point>& points, int depth = 0);

    // Recursive nearest-neighbor search
    void findNearest(KDNode* node, double x, double y, double z, int depth,
                     KDNode*& best, double& bestDistSq) const;

    // Utility to calculate squared distance between two points
    static double squaredDistance(const Point& a, double x, double y, double z);
};

#endif  // KDTREE3D_H
