#include "../include/MagneticField.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <unistd.h>

using namespace std;

Magnet::Magnet(const std::vector<Point> &pts) : points(pts), kd_tree(nullptr) {
    if (!points.empty()) {
        x_min = x_max = points[0].x;
        y_min = y_max = points[0].y;
        z_min = z_max = points[0].z;
        for (const auto &pt : points) {
            if (pt.x < x_min) x_min = pt.x; if (pt.x > x_max) x_max = pt.x;
            if (pt.y < y_min) y_min = pt.y; if (pt.y > y_max) y_max = pt.y;
            if (pt.z < z_min) z_min = pt.z; if (pt.z > z_max) z_max = pt.z;
        }
        kd_tree = std::make_unique<KDTree3D>(points); // build here
    } else {
        x_min = x_max = y_min = y_max = z_min = z_max = 0.0;
    }
}

std::vector<Magnet> all_magnets;

bool load_csv(const std::string &filepath, std::vector<Point> &points)
{
    std::ifstream file(filepath);
    if (!file.is_open())
    {
        std::cerr << "Failed to open " << filepath << "\n";
        return false;
    }

    std::string line;
    if (!std::getline(file, line))
    { // Skip header line
        std::cerr << "Empty file: " << filepath << "\n";
        return false;
    }

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string item;
        std::vector<double> row;

        while (std::getline(ss, item, ','))
        {
            row.push_back(std::stod(item));
        }

        if (row.size() < 6)
            continue; // Ensure at least 6 values
        // Create a Point and add to points vector
        points.emplace_back(
            row[0] * 1e6, // x (converted to micrometers)
            row[1] * 1e6, // y (converted to micrometers)
            row[2] * 1e6, // z (converted to micrometers)
            row[3],       // Fx
            row[4],       // Fy
            row[5]        // Fz
        );
    }

    file.close();
    return !points.empty();
}




#include <limits>

static inline bool in_bounds_of(const Magnet &m, double x, double y, double z) {
    return !(x < m.x_min || x > m.x_max ||
             y < m.y_min || y > m.y_max ||
             z < m.z_min || z > m.z_max);
}

/*
std::vector<long double> BctAdv(const std::vector<Magnet> &magnets,
                                double x, double y, double z)
{
    static thread_local size_t last_idx = std::numeric_limits<size_t>::max();

    size_t idx = magnets.size();
    if (last_idx < magnets.size() && in_bounds_of(magnets[last_idx], x, y, z)) {
        idx = last_idx;
    } else {
        for (size_t i = 0; i < magnets.size(); ++i) {
            if (in_bounds_of(magnets[i], x, y, z)) { idx = i; break; }
        }
        if (idx == magnets.size()) return {0.0L, 0.0L, 0.0L};
        last_idx = idx;
    }

    const Magnet &m = magnets[idx];
    if (!m.kd_tree) return {0.0L, 0.0L, 0.0L};

    const Point p = m.kd_tree->find_closest_point(x, y, z);
    return { (long double)p.Fx, (long double)p.Fy, (long double)p.Fz };
}
*/

std::vector<double> BctAdv(const std::vector<Magnet> &magnets,
    double x, double y, double z)
{
static thread_local size_t last_idx = std::numeric_limits<size_t>::max();

auto in_bounds_of = [](const Magnet &m, double X, double Y, double Z)->bool {
return !(X < m.x_min || X > m.x_max ||
Y < m.y_min || Y > m.y_max ||
Z < m.z_min || Z > m.z_max);
};

size_t idx = magnets.size();
if (last_idx < magnets.size() && in_bounds_of(magnets[last_idx], x, y, z)) {
idx = last_idx;
} else {
for (size_t i = 0; i < magnets.size(); ++i) {
if (in_bounds_of(magnets[i], x, y, z)) { idx = i; break; }
}
if (idx == magnets.size()) return {0.0L, 0.0L, 0.0L};
last_idx = idx;
}

const Magnet &m = magnets[idx];
if (!m.kd_tree) return {0.0L, 0.0L, 0.0L};

constexpr int K = 8;            
constexpr double EPS = 1e-9L;

const auto neigh = m.kd_tree->k_closest_points(x, y, z, K);
if (neigh.empty()) return {0.0L, 0.0L, 0.0L};

double wsum = 0.0L, bx = 0.0L, by = 0.0L, bz = 0.0L;
for (const auto &pt : neigh) {
    const double dx = (double)x - (double)pt.x;
    const double dy = (double)y - (double)pt.y;
    const double dz = (double)z - (double)pt.z;
    const double d2 = dx*dx + dy*dy + dz*dz;
    const double w  = 1.0L / std::sqrt(d2 + EPS);   // IDW p=1
    bx += w * (double)pt.Fx;
    by += w * (double)pt.Fy;
    bz += w * (double)pt.Fz;
    wsum += w;
}
if (wsum == 0.0L) return {0.0L, 0.0L, 0.0L};
return { bx/wsum, by/wsum, bz/wsum };
}

void initializeFieldMaps()
{
    std::cout << "Initializing magnetic fields...\n";

    std::string export_dir = "./ExportedMagneticFields";
    std::vector<std::string> magnet_types = {"MainDipole", "D1", "D2", "InnerQuad", "RevInnerQuad"};
    const int max_retries = 5;  
    const int retry_delay_seconds = 2;

    // Define expected magnet counts for each type
    std::map<std::string, int> expected_counts = {
        {"MainDipole", 6},
        {"D1", 6},
        {"D2", 1},
        {"InnerQuad", 2},
        {"RevInnerQuad", 2}
    };

    // First pass: validate ALL files exist
    std::cout << "Validating magnetic field files...\n";
    for (const auto &magnet_type : magnet_types)
    {
        int expected = expected_counts[magnet_type];
        for (int idx = 1; idx <= expected; ++idx)
        {
            std::stringstream ss;
            ss << export_dir << "/" << magnet_type << "/" << magnet_type << "_magnet_" << idx << ".csv";
            std::string filepath = ss.str();
            
            if (!std::filesystem::exists(filepath))
            {
                std::cerr << "FATAL ERROR: Required file missing: " << filepath << std::endl;
                std::cerr << "Cannot proceed without all magnetic field files." << std::endl;
                exit(1);
            }
        }
    }
    std::cout << "All magnetic field files validated.\n";

    // Second pass: load ALL files sequentially with retries
    for (const auto &magnet_type : magnet_types)
    {
        int expected = expected_counts[magnet_type];
        
        for (int magnet_idx = 1; magnet_idx <= expected; ++magnet_idx)
        {
            std::stringstream ss;
            ss << export_dir << "/" << magnet_type << "/" << magnet_type << "_magnet_" << magnet_idx << ".csv";
            std::string filepath = ss.str();

            std::vector<Point> points;
            bool loaded = false;

            // Retry loading up to max_retries times
            for (int attempt = 1; attempt <= max_retries; ++attempt)
            {
                points.clear(); // Clear points before each attempt
                
                if (load_csv(filepath, points))
                {
                    Magnet magnet(points);
                    all_magnets.push_back(std::move(magnet));
                    std::cout << "Successfully loaded and initialized KD-tree for " << filepath << "\n";
                    loaded = true;
                    break;
                }
                else
                {
                    std::cerr << "Failed to load " << filepath << " (Attempt " << attempt << "/" << max_retries << ")\n";
                    if (attempt < max_retries)
                    {
                        std::cerr << "Retrying in " << retry_delay_seconds << " seconds...\n";
                        sleep(retry_delay_seconds);
                    }
                }
            }

            if (!loaded)
            {
                std::cerr << "FATAL ERROR: Failed to load " << filepath << " after " << max_retries << " attempts.\n";
                std::cerr << "Magnetic field configuration is incomplete. Exiting.\n";
                exit(1);
            }
        }
    }
    
    std::cout << "All magnetic fields loaded successfully. Total magnets: " << all_magnets.size() << "\n";
}
std::vector<double> Bct(double x, double y, double z)
{
    // B field function of location region
    if (((sqrt(x * x + y * y) < 0.06e6) && (abs(z - 72.287e6) < 12.365e6)) ||
        ((((sqrt((x - 0.093e6) * (x - 0.093e6) + y * y) < 0.04e6) || (sqrt((x + 0.093e6) * (x + 0.093e6) + y * y) < 0.04e6)) && (abs(z - 158.2e6) < 4.725e6))))
    {
        if ((sqrt(x * x + y * y) < 0.06e6) && (abs(z - 72.287e6) < 12.365e6))
            return {0, 3.5, 0};
        if (((sqrt((x - 0.093e6) * (x - 0.093e6) + y * y) < 0.04e6) || (sqrt((x + 0.093e6) * (x + 0.093e6) + y * y) < 0.04e6)) && (abs(z - 158.2e6) < 4.725e6))
            return {0, -3.5, 0};
    }
    return {0, 0, 0};
}

std::vector<long double> BctLong(long double x, long double y, long double z)
{
    // B long function

    if (((std::sqrt(x * x + y * y) < 0.06e6L) && (std::abs(z - 72.287e6L) < 12.365e6L)) ||
        ((((std::sqrt((x - 0.093e6L) * (x - 0.093e6L) + y * y) < 0.04e6L) ||
           (std::sqrt((x + 0.093e6L) * (x + 0.093e6L) + y * y) < 0.04e6L)) &&
          (std::abs(z - 158.2e6L) < 4.725e6L))))
    {
        if ((std::sqrt(x * x + y * y) < 0.06e6L) && (std::abs(z - 72.287e6L) < 12.365e6L))
        {
            return {0.0L, 3.5L, 0.0L};
        }
        if (((std::sqrt((x - 0.093e6L) * (x - 0.093e6L) + y * y) < 0.04e6L) ||
             (std::sqrt((x + 0.093e6L) * (x + 0.093e6L) + y * y) < 0.04e6L)) &&
            (std::abs(z - 158.2e6L) < 4.725e6L))
        {
            return {0.0L, -3.5L, 0.0L};
        }
    }
    return {0.0L, 0.0L, 0.0L};
}
