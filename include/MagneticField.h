#pragma once

#include "KDTree3D.h"

#include <memory>
#include <string>
#include <vector>

struct Magnet {
    std::vector<Point> points;
    std::unique_ptr<KDTree3D> kd_tree;
    Magnet(const std::vector<Point>& pts);
    double x_min, x_max, y_min, y_max, z_min, z_max;
};

extern std::vector<Magnet> all_magnets;

bool load_csv(const std::string& filepath, std::vector<Point>& points);
std::vector<double> BctAdv(const std::vector<Magnet>& magnets, double x, double y, double z);
void initializeFieldMaps();
std::vector<double> Bct(double x, double y, double z);
std::vector<long double> BctLong(long double x, long double y, long double z);
