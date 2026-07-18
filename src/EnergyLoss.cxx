#include "../include/EnergyLoss.h"

#include <cmath>

using namespace std;

const double ZCu = 29;
const double ZoACu = ZCu / 63.546;
const double rhoCu = 8.960;
const double I0Cu = 322.0e-6;
const double aCu = 0.14339;
const double kCu = 2.9044;
const double x0Cu = -0.0254;
const double x1Cu = 3.2792;
const double CbarCu = 4.4190;
const double d0Cu = 0.08;

const double ZCc = 8.56;
const double ZoACc = 0.50274;
const double rhoCc = 2.300;
const double I0Cc = 135.2e-6;
const double aCc = 0.07515;
const double kCc = 3.5467;
const double x0Cc = 0.1301;
const double x1Cc = 3.0466;
const double CbarCc = 3.9464;
const double d0Cc = 0.00;

const double ZRock = 11;
const double ZoARock = 0.5;
const double rhoRock = 2.650;
const double I0Rock = 136.4e-6;
const double aRock = 0.08301;
const double kRock = 3.4120;
const double x0Rock = 0.0492;
const double x1Rock = 3.0549;
const double CbarRock = 3.7738;
const double d0Rock = 0.00;

//-------------------- dE/dx funcitons ---------------------

/*
    The variables m (in GeV), z (in e), and \[Beta] (in light speed C) represent the mass, the electric charge number, and the velocity of the quirk particle, respectively;
    \[Delta]x (in cm) represents the distance traveled during this time step; The function EoxCu/Cc/RockGaus[m,z,\[Beta],\[Delta]x] (in 10^-16 GeV^2) provides the de/dx 
    value, which is normally distributed with EoxCu/Cc/RockAll[m,z,\[Beta]] (in 10^-16 GeV^2) as the mean and \[Sigma]Cu/Cc/Rock[m,z,\[Beta],\[Delta]x] (in 10^-16 GeV^2) as the standard deviation.
*/

static inline long double horner_desc(const long double *c, int n, long double x) {
    //for faster high power calcs
    long double y = c[0];
    for (int i = 1; i < n; ++i) y = y * x + c[i];
    return y;
}

static inline double sample_standard_normal(std::mt19937 &gen) {
    static thread_local std::normal_distribution<double> N01(0.0, 1.0);
    return N01(gen);
}

double Gamma(double beta)
{
    return 1.0 / sqrt(1.0 - beta * beta);
}

double Emax(double m, double beta)
{
    return (2 * 0.511 * beta * beta * pow(Gamma(beta), 2)) /
           (1 + 2 * Gamma(beta) * (0.511e-3 / m) + pow(0.511e-3 / m, 2));
}

double Xi(double z, double ZoA, double rho, double beta)
{
    return 0.5 * 0.307075 * z * z * ZoA * rho / (beta * beta);
}

double Delta(double x, double x0, double x1, double Cbar, double a, double k, double d0)
{
    if (x >= x1)
    {
        return 2 * log(10) * x - Cbar;
    }
    else if (x >= x0)
    {
        return 2 * log(10) * x - Cbar + a * pow(x1 - x, k);
    }
    else
    {
        return d0 * pow(10, 2 * (x - x0));
    }
}

double Eox(double m, double z, double ZoA, double rho, double I0, double x0, double x1, double Cbar, double a, double k, double d0, double beta)
{

    return 2 * Xi(z, ZoA, rho, beta) * (0.5 * log((2 * 0.511 * beta * beta * pow(Gamma(beta), 2) * Emax(m, beta)) / pow(I0, 2)) - beta * beta - Delta(log10(beta * Gamma(beta)), x0, x1, Cbar, a, k, d0) / 2);
}

double EoxCuAll(double mq, int param, double beta)
{
    // de/dx copper
    // EoxCuAll with recasting variables to long double
    long double mq_ld = static_cast<long double>(mq);
    long double beta_ld = static_cast<long double>(beta);
    long double scalar;
    long double ZCu_eff = 1.0L;

    if (beta_ld <= 0.00226L)
    {
        scalar = 37597.30061169589L * beta_ld;
    }
    else if (beta_ld >= 0.06345454545454546L)
    {
        scalar = static_cast<long double>(Eox(static_cast<double>(mq_ld), static_cast<double>(ZCu_eff), ZoACu, rhoCu, I0Cu, x0Cu, x1Cu, CbarCu, aCu, kCu, d0Cu, static_cast<double>(beta_ld))) * 0.197L;
    }
    else
    {
        static const long double cu_desc[15] = {
            5.982648430710585e19L,
        -5.714216718459019e19L,
            2.496247979943447e19L,
        -6.596425471576249e18L,
            1.174756123051605e18L,
        -1.486143034696763e17L,
            1.370510669479835e16L,
        -9.294312877895761e14L,
            4.608165819030902e13L,
        -1.6320734910018296e12L,
            3.916027090833348e10L,
        -5.582148886272709e8L,
            2.6438521703577857e6L,
            33909.5576744394L,
            0.34289008715835223L
        };
        scalar = horner_desc(cu_desc, 15, beta_ld);
    }

    return static_cast<double>(scalar);
}

double EoxCcAll(double mq, int param, double beta)
{
    // de/dx concrete
    long double mq_ld = static_cast<long double>(mq);
    long double beta_ld = static_cast<long double>(beta);
    long double scalar;
    long double ZCc_eff = 1.0L;

    if (beta_ld <= 0.00226L)
    {
        scalar = 30376.62714732753L * beta_ld;
    }
    else if (beta_ld >= 0.06181818181818182L)
    {
        scalar = static_cast<long double>(Eox(static_cast<double>(mq_ld), static_cast<double>(ZCc_eff), ZoACc, rhoCc, I0Cc, x0Cc, x1Cc, CbarCc, aCc, kCc, d0Cc, static_cast<double>(beta_ld))) * 0.197L;
    }
    else
    {
        static const long double cc_desc[15] = {
            6.737738168445225e19L,
        -6.442273849826215e19L,
            2.819535817843903e19L,
        -7.470735431733813e18L,
            1.3350858443896056e18L,
        -1.695911162356092e17L,
            1.5707786950520584e16L,
        -1.0693146651849746e15L,
            5.310954532975338e13L,
        -1.8750617184047375e12L,
            4.441507306806682e10L,
        -6.152256215129855e8L,
            2.8807310255188923e6L,
            26377.245398907748L,
            0.3717416673986129L
        };
        scalar = horner_desc(cc_desc, 15, beta_ld);
    }

    return static_cast<double>(scalar);
}

double EoxRockAll(double mq, int param, double beta)
{
    // de/dx rock
    long double mq_ld = static_cast<long double>(mq);
    long double beta_ld = static_cast<long double>(beta);
    long double scalar;
    long double ZRock_eff = 1.0L;

    if (beta_ld <= 0.00226L)
    {
        scalar = 28340.291807946152L * beta_ld;
    }
    else if (beta_ld >= 0.05745454545454545L)
    {
        scalar = static_cast<long double>(Eox(static_cast<double>(mq_ld), static_cast<double>(ZRock_eff), ZoARock, rhoRock, I0Rock, x0Rock, x1Rock, CbarRock, aRock, kRock, d0Rock, static_cast<double>(beta_ld))) * 0.197L;
    }
    else
    {
        static const long double rock_desc[15] = {
            5.686153048402872e19L,
        -5.34284980622759e19L,
            2.30021301170228e19L,
        -6.003163825933505e18L,
            1.0585274314396489e18L,
        -1.3296907090779789e17L,
            1.2214854741842604e16L,
        -8.278484331611662e14L,
            4.113520947417071e13L,
        -1.462059229558537e12L,
            3.513441645396524e10L,
        -4.976688103007817e8L,
            2.350400267498349e6L,
            25065.834574193614L,
            0.3044358906484703L
        };
        scalar = horner_desc(rock_desc, 15, beta_ld);
    }

    return static_cast<double>(scalar);
}

// really the de/dx Gaus functions should use truncated distributions

/*double truncated_normal(double mean, double std_dev) {
    boost::math::normal_distribution<> normal_dist(mean, std_dev);
    boost::random::uniform_real_distribution<> uniform_dist(0.0, 1.0);

    double lower_cdf = boost::math::cdf(normal_dist, 0);
    double upper_cdf = boost::math::cdf(normal_dist, std::numeric_limits<double>::infinity());

    double u = uniform_dist(gen) * (upper_cdf - lower_cdf) + lower_cdf;

    // inverse CDF corresponding to the generated uniform random number
    return boost::math::quantile(normal_dist, u);
}*/


double EoxGaus(double m, double z, double ZoA, double rho, double I0, double x0, double x1, double Cbar, double a, double k, double d0, double beta, double delta_x, std::function<double(double, int, double)> EoxAllFunc, std::mt19937 &gen)
{
    double z_eff = 1;
    double mean = EoxAllFunc(m, z_eff, beta);
    double std_dev = 0.197 * std::sqrt(Xi(z_eff, ZoA, rho, beta) * delta_x * Emax(m, beta) * (1 - beta * beta / 2)) / delta_x;
    double z_draw = sample_standard_normal(gen);
    return mean + std_dev * z_draw;
}

// de/dx gausfor different materials
double EoxCuGaus(double mq, int param, const double beta, double dx, std::mt19937 &gen)
{
    return EoxGaus(mq, ZCu, ZoACu, rhoCu, I0Cu, x0Cu, x1Cu, CbarCu, aCu, kCu, d0Cu, beta, dx, EoxCuAll, gen);
}

double EoxCcGaus(double mq, int param, const double beta, double dx, std::mt19937 &gen)
{
    return EoxGaus(mq, ZCc, ZoACc, rhoCc, I0Cc, x0Cc, x1Cc, CbarCc, aCc, kCc, d0Cc, beta, dx, EoxCcAll, gen);
}

double EoxRockGaus(double mq, int param, const double beta, double dx, std::mt19937 &gen)
{
    return EoxGaus(mq, ZRock, ZoARock, rhoRock, I0Rock, x0Rock, x1Rock, CbarRock, aRock, kRock, d0Rock, beta, dx, EoxRockAll, gen);
}

std::vector<double> EoxCu(double mq, int param, const std::vector<double> &v)
{
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i)
    {
        result[i] = Eox(mq, ZCu, ZoACu, rhoCu, I0Cu, x0Cu, x1Cu, CbarCu, aCu, kCu, d0Cu, v[i]);
    }
    return result;
}

std::vector<double> EoxCc(double mq, int param, const std::vector<double> &v)
{
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i)
    {
        result[i] = Eox(mq, ZCc, ZoACc, rhoCc, I0Cc, x0Cc, x1Cc, CbarCc, aCc, kCc, d0Cc, v[i]);
    }
    return result;
}

std::vector<double> EoxRock(double mq, int param, const std::vector<double> &v)
{
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i)
    {
        result[i] = Eox(mq, ZRock, ZoARock, rhoRock, I0Rock, x0Rock, x1Rock, CbarRock, aRock, kRock, d0Rock, v[i]);
    }
    return result;
}
