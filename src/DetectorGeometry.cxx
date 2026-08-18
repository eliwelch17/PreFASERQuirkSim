#include "../include/DetectorGeometry.h"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

// TAN piece 1: only this piece was included in the original Mathematica implementaton
static bool InTanPiece1(double x, double y, double z)
{
    if (abs(z - 140.5e6) >= 0.5e6)
        return false;
    if (abs(x) >= (0.094 / 2) * 1e6)
        return false;
    if (abs(y - (0.605 / 2 - 0.067) * 1e6) >= (0.605 / 2) * 1e6)
        return false;

    return true;
}

// https://lss.fnal.gov/archive/test-fn/0000/fermilab-fn-0732.pdf
static bool InTanPiece2(double x, double y, double z)
{
    if (abs(z - 141.75e6) >= 1.75e6)
        return false;
    if (abs(x) >= 0.130e6)
        return false;
    if (abs(y) >= 0.105e6)
        return false;

    const double hole_r = 0.025e6;
    const double hole_x = 0.08e6;
    if (sqrt((x - hole_x) * (x - hole_x) + y * y) < hole_r)
        return false;
    if (sqrt((x + hole_x) * (x + hole_x) + y * y) < hole_r)
        return false;
    if (InTanPiece1(x, y, z))
        return false;

    return true;
}

static bool InTanCopper(double x, double y, double z)
{
    return InTanPiece1(x, y, z) || InTanPiece2(x, y, z);
}

int Loct(double x, double y, double z)
{
    // determine location region of quirks
    if ((((sqrt(x * x + y * y) > 0.017e6) && (abs(z - 19.9e6) < 0.9e6)) ||
         InTanCopper(x, y, z) ||
         (abs(z - 385.0e6) < 5.0e6) || (abs(z - 432.3e6) < 42.3e6)))
    {
        // TAS
        if ((sqrt(x * x + y * y) > 0.017e6) && (abs(z - 19.9e6) < 0.9e6))
            return 1;
        // TAN
        if (InTanCopper(x, y, z))
            return 2;
        // CONCRETE
        if (abs(z - 385.0e6) < 5.0e6)
            return 3;
        // ROCK
        if (abs(z - 432.3e6) < 42.3e6)
            return 4;
    }
    return 0;
}

static inline void add_z_breakpoint(std::vector<double> &bp, double z, double z0, double z1)
{
    if (z > z0 && z < z1)
        bp.push_back(z);
}

double fraction_in_loct_com(int loct_code, double z0_um, double z1_um,
                                   double bx, double by, double bz)
{
    if (z1_um <= z0_um)
        return 0.0;
    if (loct_code != 1 && loct_code != 2)
        return 1.0;

    const double kx = bx / bz;
    const double ky = by / bz;
    const double k_perp = std::hypot(kx, ky);

    std::vector<double> bp;
    bp.reserve(24);
    bp.push_back(z0_um);

    if (loct_code == 1)
    {
        constexpr double R = 0.017e6;
        if (k_perp > 0.0)
            add_z_breakpoint(bp, R / k_perp, z0_um, z1_um);
    }
    else if (loct_code == 2)
    {
        const double piece1_half_x = (0.094 / 2) * 1e6;
        const double piece1_half_y = (0.605 / 2) * 1e6;
        const double piece1_y_off = (0.605 / 2 - 0.067) * 1e6;
        const double piece2_half_x = 0.130e6;
        const double piece2_half_y = 0.105e6;
        const double hole_r = 0.025e6;
        const double hole_x = 0.08e6;

        auto add_linear = [&](double k, double val) {
            if (std::abs(k) > 0.0)
                add_z_breakpoint(bp, val / k, z0_um, z1_um);
        };
        add_z_breakpoint(bp, 141.0e6, z0_um, z1_um);
        add_linear(kx, piece1_half_x);
        add_linear(kx, -piece1_half_x);
        add_linear(ky, piece1_y_off + piece1_half_y);
        add_linear(ky, piece1_y_off - piece1_half_y);
        add_linear(kx, piece2_half_x);
        add_linear(kx, -piece2_half_x);
        add_linear(ky, piece2_half_y);
        add_linear(ky, -piece2_half_y);

        auto add_circle = [&](double xc) {
            const double A = kx * kx + ky * ky;
            if (A <= 0.0)
                return;
            const double B = -2.0 * kx * xc;
            const double C = xc * xc - hole_r * hole_r;
            const double D = B * B - 4.0 * A * C;
            if (D < 0.0)
                return;
            const double sd = std::sqrt(D);
            add_z_breakpoint(bp, (-B - sd) / (2.0 * A), z0_um, z1_um);
            add_z_breakpoint(bp, (-B + sd) / (2.0 * A), z0_um, z1_um);
        };
        add_circle(hole_x);
        add_circle(-hole_x);
    }

    bp.push_back(z1_um);
    std::sort(bp.begin(), bp.end());

    double in_len = 0.0;
    for (size_t i = 0; i + 1 < bp.size(); ++i)
    {
        const double za = bp[i];
        const double zb = bp[i + 1];
        if (zb <= za)
            continue;
        const double zm = 0.5 * (za + zb);
        if (Loct(kx * zm, ky * zm, zm) == loct_code)
            in_len += (zb - za);
    }
    return in_len / (z1_um - z0_um);
}


int Layer(double x, double y, double z)
{
    // determine scintilaltor layers
    if (abs(x) < 0.15e6 && abs(y) < 0.15e6 && abs(z - 480.01e6) < 0.01e6)
        return 1;
    if (abs(x) < 0.2e6 && abs(y) < 0.2e6 && abs(z - 481.555e6) < 0.005e6)
        return 2;
    if (abs(x) < 0.15e6 && abs(y) < 0.15e6 && abs(z - 484.18e6) < 0.01e6)
        return 3;
    return 0;
}
