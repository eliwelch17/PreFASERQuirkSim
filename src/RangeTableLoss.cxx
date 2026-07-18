#include "../include/RangeTableLoss.h"

#include "../include/DetectorGeometry.h"
#include "../include/VectorUtils.h"

#include "../rangeTables/range_50.h"
#include "../rangeTables/range_65.h"
#include "../rangeTables/range_75.h"
#include "../rangeTables/range_100.h"
#include "../rangeTables/range_125.h"
#include "../rangeTables/range_150.h"
#include "../rangeTables/range_175.h"
#include "../rangeTables/range_200.h"
#include "../rangeTables/range_225.h"
#include "../rangeTables/range_250.h"
#include "../rangeTables/range_275.h"
#include "../rangeTables/range_300.h"
#include "../rangeTables/range_325.h"
#include "../rangeTables/range_350.h"
#include "../rangeTables/range_375.h"
#include "../rangeTables/range_400.h"
#include "../rangeTables/range_425.h"
#include "../rangeTables/range_450.h"
#include "../rangeTables/range_475.h"
#include "../rangeTables/range_500.h"
#include "../rangeTables/range_525.h"
#include "../rangeTables/range_550.h"
#include "../rangeTables/range_575.h"
#include "../rangeTables/range_600.h"
#include "../rangeTables/range_625.h"
#include "../rangeTables/range_650.h"
#include "../rangeTables/range_675.h"
#include "../rangeTables/range_700.h"

#include <algorithm>
#include <cmath>
#include <iostream>

// -----------------------------------------------------------------------------------------------------------------------------
// Range-table based ionization-loss correction for the skipped 0 -> front segment (used when we start the simulation near `back`)
// -----------------------------------------------------------------------------------------------------------------------------

namespace RangeTables {

struct TableView {
  int n;
  const double* beta;
  const double* r_cu_um;
  const double* r_cc_um;
  const double* r_rock_um;
};

#define RT_RET(M) \
  return TableView{ range_m##M::RANGE_N, range_m##M::beta_grid, range_m##M::range_Cu_um, range_m##M::range_Cc_um, range_m##M::range_Rock_um };

static inline TableView get_table_for_mass(int mq_gev_int) {
  switch (mq_gev_int) {
    case 50:  RT_RET(50);
    case 65:  RT_RET(65);
    case 75:  RT_RET(75);
    case 100: RT_RET(100);
    case 125: RT_RET(125);
    case 150: RT_RET(150);
    case 175: RT_RET(175);
    case 200: RT_RET(200);
    case 225: RT_RET(225);
    case 250: RT_RET(250);
    case 275: RT_RET(275);
    case 300: RT_RET(300);
    case 325: RT_RET(325);
    case 350: RT_RET(350);
    case 375: RT_RET(375);
    case 400: RT_RET(400);
    case 425: RT_RET(425);
    case 450: RT_RET(450);
    case 475: RT_RET(475);
    case 500: RT_RET(500);
    case 525: RT_RET(525);
    case 550: RT_RET(550);
    case 575: RT_RET(575);
    case 600: RT_RET(600);
    case 625: RT_RET(625);
    case 650: RT_RET(650);
    case 675: RT_RET(675);
    case 700: RT_RET(700);
    default:
      // Unsupported mass: return empty view.
      return TableView{0, nullptr, nullptr, nullptr, nullptr};
  }
}

#undef RT_RET

static inline double clamp01(double x) {
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

static inline double beta_from_pE(const std::vector<double>& p, double E) {
  const double p2 = DotProduct(p, p);
  if (E <= 0.0) return 0.0;
  return std::sqrt(p2) / E;
}

static inline double energy_from_beta(double mq, double beta) {
  beta = clamp01(beta);
  if (beta >= 1.0) beta = 1.0 - 1e-15;
  const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
  return mq * gamma;
}

static inline double p_mag_from_beta(double mq, double beta) {
  beta = clamp01(beta);
  if (beta >= 1.0) beta = 1.0 - 1e-15;
  const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
  return mq * beta * gamma;
}

static inline double interp_monotone(const double* x, const double* y, int n, double xq) {
  // Linear interpolation on strictly increasing x.
  if (n <= 0) return 0.0;
  if (xq <= x[0]) return y[0];
  if (xq >= x[n - 1]) return y[n - 1];
  const double* it = std::lower_bound(x, x + n, xq);
  int hi = int(it - x);
  int lo = hi - 1;
  const double x0 = x[lo], x1 = x[hi];
  const double y0 = y[lo], y1 = y[hi];
  const double t = (xq - x0) / (x1 - x0);
  return y0 + t * (y1 - y0);
}

static inline double invert_R_to_beta(const TableView& tab, const double* R_um, double Rq_um) {
  // Given monotone increasing R(beta), return beta such that R(beta)=Rq_um.
  if (tab.n <= 0) return 0.0;
  if (Rq_um <= R_um[0]) return tab.beta[0];
  if (Rq_um >= R_um[tab.n - 1]) return tab.beta[tab.n - 1];
  const double* it = std::lower_bound(R_um, R_um + tab.n, Rq_um);
  int hi = int(it - R_um);
  int lo = hi - 1;
  const double R0 = R_um[lo], R1 = R_um[hi];
  const double b0 = tab.beta[lo], b1 = tab.beta[hi];
  const double t = (Rq_um - R0) / (R1 - R0);
  return b0 + t * (b1 - b0);
}

enum class Material { Cu, Cc, Rock };

static inline const double* range_array_for(Material m, const TableView& tab) {
  switch (m) {
    case Material::Cu:   return tab.r_cu_um;
    case Material::Cc:   return tab.r_cc_um;
    case Material::Rock: return tab.r_rock_um;
  }
  return nullptr;
}

static inline bool apply_material_range_step(
    const TableView& tab,
    Material mat,
    double L_um,
    double mq,
    double beta_min,
    std::vector<double>& p,
    double& E)
{
  if (L_um <= 0.0) return true;
  const double beta = beta_from_pE(p, E);
  if (beta < beta_min) return false;

  const double* R_um = range_array_for(mat, tab);
  if (!R_um) return false;

  const double R_old = interp_monotone(tab.beta, R_um, tab.n, beta);
  const double R_new = R_old - L_um;
  if (R_new <= 0.0) return false;

  const double beta_new = invert_R_to_beta(tab, R_um, R_new);
  if (beta_new < beta_min) return false;

  const double p_new_mag = p_mag_from_beta(mq, beta_new);
  const double p_old_mag = std::sqrt(DotProduct(p, p));
  if (p_old_mag <= 0.0) return false;

  p = MultiplyVector(p, p_new_mag / p_old_mag);
  E = energy_from_beta(mq, beta_new);
  return true;
}

// Apply a common ionization loss scaling to the pair momenta

static inline bool apply_material_range_step_common(const TableView& tab,Material mat,double L_um,double mq,double beta_min,std::vector<double>& p1,double& E1,std::vector<double>& p2,double& E2)
{
  if (L_um <= 0.0) return true;

  const double beta1 = beta_from_pE(p1, E1);
  const double beta2 = beta_from_pE(p2, E2);
  const double beta_rep = 0.5 * (beta1 + beta2);
  if (beta_rep < beta_min) return false;

  const double* R_um = range_array_for(mat, tab);
  if (!R_um) return false;

  const double R_old = interp_monotone(tab.beta, R_um, tab.n, beta_rep);
  const double R_new = R_old - L_um;
  if (R_new <= 0.0) return false;

  const double beta_new = invert_R_to_beta(tab, R_um, R_new);
  if (beta_new < beta_min) return false;

  const double p_rep_old = p_mag_from_beta(mq, beta_rep);
  const double p_rep_new = p_mag_from_beta(mq, beta_new);
  if (!(p_rep_old > 0.0)) return false;
  const double f = p_rep_new / p_rep_old;

  p1 = MultiplyVector(p1, f);
  p2 = MultiplyVector(p2, f);
  E1 = std::sqrt(mq * mq + DotProduct(p1, p1));
  E2 = std::sqrt(mq * mq + DotProduct(p2, p2));
  return true;
}

static inline double osc_path_factor(double k) {
  // Average over one period of sqrt(1 + k^2 cos^2 theta).
  // Small-k series: 1 + k^2/4 - 3k^4/64 + O(k^6)
  k = std::abs(k);
  const double PI = std::acos(-1.0);
  if (k < 0.2) {
    const double k2 = k * k;
    return 1.0 + 0.25 * k2 - (3.0 / 64.0) * k2 * k2;
  }
  // Simpson on [0, 2 pi]
  const int N = 400; // even
  const double a = 0.0;
  const double b = 2.0 * PI;
  const double h = (b - a) / N;
  auto f = [&](double t) {
    const double c = std::cos(t);
    return std::sqrt(1.0 + (k * k) * (c * c));
  };
  double s = f(a) + f(b);
  for (int i = 1; i < N; ++i) {
    s += (i % 2 ? 4.0 : 2.0) * f(a + i * h);
  }
  const double integral = (h / 3.0) * s;
  return integral / (2.0 * PI);
}

bool apply_range_table_loss(int mq_int,double mq,double Lambda_eV,double z_end_um,double distance_period_m,  const std::vector<double>& Beta_pair,double beta_min,std::vector<double>& p1,double& E1,std::vector<double>& p2,double& E2,double& t_arrival_ns_out)
{
    //distance_period is in meters
  
  if (z_end_um <= 0.0) return true;
  if (beta_min <= 0.0) beta_min = 0.1;
  t_arrival_ns_out = 0.0;

  const TableView tab = get_table_for_mass(mq_int);
  if (tab.n <= 0) {
    std::cerr << "No range table for mq=" << mq_int << " GeV\n";
    return false;
  }

  // Require forward-going pair direction for the COM ray
  const double Beta_mag = std::sqrt(DotProduct(Beta_pair, Beta_pair));
  if (!(Beta_mag > 0.0)) return false;
  const double Beta_z = Beta_pair[2];
  if (Beta_z <= 0.0) return false;

  // Oscillation amplitude from COM energy :
  // rmax_um = (KE_total / Lambda_GeV^2) * hbarc_GeV_um ; A_each = rmax/2
  const double Lambda_GeV = Lambda_eV * 1e-9;
  if (!(Lambda_GeV > 0.0)) return false;

  const std::vector<double> p1_COM = BoostToCOM(p1, E1, Beta_pair);
  const std::vector<double> p2_COM = BoostToCOM(p2, E2, Beta_pair);
  const double E1_COM = std::sqrt(mq * mq + DotProduct(p1_COM, p1_COM));
  const double E2_COM = std::sqrt(mq * mq + DotProduct(p2_COM, p2_COM));
  const double KE_total = (E1_COM - mq) + (E2_COM - mq);
  if (!(KE_total >= 0.0)) return false;

  constexpr double hbarc_GeV_um = 1.9732697e-10;
  const double rmax_um = (KE_total / (Lambda_GeV * Lambda_GeV)) * hbarc_GeV_um;
  const double A_each_um = 0.5 * rmax_um;

  // Full sin  wavelength along the COM path:
  // distance_period_m is closest approach distance, so wavelength_sin = 2 * distance_period
  const double lambda_osc_um = 2.0 * distance_period_m * 1e6;
  if (!(lambda_osc_um > 0.0)) return false;

  const double PI = std::acos(-1.0);
  const double k = (2.0 * PI * A_each_um) / lambda_osc_um;
  const double osc_factor = osc_path_factor(k);

  // Convert z thickness to slant length along COM ray: L = dz / cosθz, cosθz = Beta_z/|Beta|
  const double cos_theta_z = Beta_z / Beta_mag; // Beta_z > 0 enforced above
  if (!(cos_theta_z > 0.0)) return false;

  struct Slab { double z0, z1; int loct; Material mat; };
  const Slab slabs[] = {
    // Match Loct() definitions
    {19.0e6, 20.8e6, 1, Material::Cu},        // TAS copper window
    {140.0e6, 141.0e6, 2, Material::Cu},      // TAN copper (two 25 mm holes at y=±80 mm)
    {380.0e6, 390.0e6, 3, Material::Cc},      // concrete
    {390.0e6, 474.6e6, 4, Material::Rock},    // rock
  };

  // Coarse arrival time estimate to z=front (ns):
  // We want "arrival at a given z plane", so use dt = dz / (Beta_z * c).
  // Approximation: treat Beta_z constant in vacuum, and across concrete/rock use a trapezoid in 1/Beta_z
  // using the Beta_z before/after applying the corresponding ionization-loss step.

  const double c_um_per_ns = 300000.0;
  auto safe_Bz = [&](double bz) -> double {
    if (bz < 1e-6) return 1e-6;
    return bz;
  };
  auto beta_pair_z = [&](const std::vector<double>& pa, double Ea,
                         const std::vector<double>& pb, double Eb) -> double {
    const double Et = Ea + Eb;
    if (!(Et > 0.0)) return 0.0;
    return (pa[2] + pb[2]) / Et;
  };

  // Region boundaries (micrometers) consistent with Loct()
  const double z_conc_start_um = 380.0e6;
  const double z_conc_end_um   = 390.0e6;
  const double z_rock_start_um = 390.0e6;
  const double z_rock_end_um   = 474.6e6;

  // Beta_z checkpoints (computed from current p/E during the loss pass)
  const double Bz_init = Beta_pair[2];
  double Bz_after_conc = Bz_init;
  double Bz_after_rock = Bz_init;
  bool saw_conc = false;
  bool saw_rock = false;
  double Bz_in_conc = Bz_init, Bz_out_conc = Bz_init;
  double Bz_in_rock = Bz_init, Bz_out_rock = Bz_init;

  for (const auto& s : slabs) {
    if (s.z0 >= z_end_um) break;
    const double z_start = std::max(0.0, s.z0);
    const double z_end = std::min(z_end_um, s.z1);
    if (z_end <= z_start) continue;

    double frac = 1.0;
    // For TAS/TAN
    if (s.loct == 1 || s.loct == 2) {
      frac = fraction_in_loct_com(s.loct, z_start, z_end, Beta_pair[0], Beta_pair[1], Beta_z);
    }
    if (frac <= 0.0) continue;

    const double dz_eff_um = (z_end - z_start) * frac;
    const double L_slant_um = (dz_eff_um / cos_theta_z);
    const double L_path_um = L_slant_um * osc_factor;

    if (s.mat == Material::Cc && !saw_conc) {
      saw_conc = true;
      Bz_in_conc = beta_pair_z(p1, E1, p2, E2);
    }
    if (s.mat == Material::Rock && !saw_rock) {
      saw_rock = true;
      Bz_in_rock = beta_pair_z(p1, E1, p2, E2);
    }

    if (!apply_material_range_step_common(tab, s.mat, L_path_um, mq, beta_min, p1, E1, p2, E2)) return false;

    if (s.mat == Material::Cc) {
      Bz_out_conc = beta_pair_z(p1, E1, p2, E2);
      Bz_after_conc = Bz_out_conc;
    }
    if (s.mat == Material::Rock) {
      Bz_out_rock = beta_pair_z(p1, E1, p2, E2);
      Bz_after_rock = Bz_out_rock;
    }
  }

  // Build dz segments to front (um)
  const double zf = z_end_um;
  // segment A: 0 -> min(front, 380m)
  const double dz_A = std::max(0.0, std::min(zf, z_conc_start_um) - 0.0);
  // segment B: 380m -> min(front, 390m)  (concrete window)
  const double dz_B = std::max(0.0, std::min(zf, z_conc_end_um) - z_conc_start_um);
  // segment C: 390m -> min(front, 474.6m)  (rock window)
  const double dz_C = std::max(0.0, std::min(zf, z_rock_end_um) - z_rock_start_um);
  // segment D: beyond 474.6m (vacuum)
  const double dz_D = std::max(0.0, zf - z_rock_end_um);

  // A: vacuum-ish, constant initial Beta_z
  const double dt_A = (dz_A / c_um_per_ns) * (1.0 / safe_Bz(Bz_init));

  // B: concrete, trapezoid in 1/Beta_z 
  double dt_B = 0.0;
  if (dz_B > 0.0) {
    const double Bzin = safe_Bz(saw_conc ? Bz_in_conc : Bz_init);
    const double Bzout = safe_Bz(saw_conc ? Bz_out_conc : Bz_after_conc);
    dt_B = (dz_B / c_um_per_ns) * 0.5 * ((1.0 / Bzin) + (1.0 / Bzout));
  }

  // C: rock, trapezoid in 1/Beta_z 
  double dt_C = 0.0;
  if (dz_C > 0.0) {
    const double Bzin = safe_Bz(saw_rock ? Bz_in_rock : Bz_after_conc);
    const double Bzout = safe_Bz(saw_rock ? Bz_out_rock : Bz_after_rock);
    dt_C = (dz_C / c_um_per_ns) * 0.5 * ((1.0 / Bzin) + (1.0 / Bzout));
  }

  // D: vacuum beyond rock end (unlikely this ever used..), constant Beta_z at end of rock
  const double dt_D = (dz_D / c_um_per_ns) * (1.0 / safe_Bz(Bz_after_rock));

  t_arrival_ns_out = dt_A + dt_B + dt_C + dt_D;

  return true;
}

// Apply range-table loss only over a z-span [z_start_um, z_end_um] (with z_end > z_start).
// This lets us "phase map" in stages without double-counting loss.
bool apply_range_table_loss_zspan(int mq_int,double mq,double Lambda_eV,double z_start_um,double z_end_um,double distance_period_m,   const std::vector<double>& Beta_pair,double beta_min,std::vector<double>& p1,double& E1,std::vector<double>& p2,double& E2)
{
    // closest-approach distance (meters)
  if (z_end_um <= z_start_um) return true;
  if (beta_min <= 0.0) beta_min = 0.1;

  const TableView tab = get_table_for_mass(mq_int);
  if (tab.n <= 0) return false;

  const double Beta_mag = std::sqrt(DotProduct(Beta_pair, Beta_pair));
  const double Beta_z = Beta_pair[2];
  if (!(Beta_mag > 0.0) || !(Beta_z > 0.0)) return false;

  // amplitude/osc factor (same as in apply_range_table_loss)
  const double Lambda_GeV = Lambda_eV * 1e-9;
  if (!(Lambda_GeV > 0.0)) return false;

  const std::vector<double> p1_COM = BoostToCOM(p1, E1, Beta_pair);
  const std::vector<double> p2_COM = BoostToCOM(p2, E2, Beta_pair);
  const double E1_COM = std::sqrt(mq * mq + DotProduct(p1_COM, p1_COM));
  const double E2_COM = std::sqrt(mq * mq + DotProduct(p2_COM, p2_COM));
  const double KE_total = (E1_COM - mq) + (E2_COM - mq);
  constexpr double hbarc_GeV_um = 1.9732697e-10;
  const double rmax_um = (KE_total / (Lambda_GeV * Lambda_GeV)) * hbarc_GeV_um;
  const double A_each_um = 0.5 * rmax_um;

  const double lambda_osc_um = 2.0 * distance_period_m * 1e6; // wavelength_sin = 2*distance_period
  if (!(lambda_osc_um > 0.0)) return false;
  const double PI = std::acos(-1.0);
  const double k = (2.0 * PI * A_each_um) / lambda_osc_um;
  const double osc_factor = osc_path_factor(k);

  const double cos_theta_z = Beta_z / Beta_mag;
  if (!(cos_theta_z > 0.0)) return false;

  struct Slab { double z0, z1; int loct; Material mat; };
  const Slab slabs[] = {
    {19.0e6, 20.8e6, 1, Material::Cu},
    {140.0e6, 141.0e6, 2, Material::Cu},      // TAN copper (two 25 mm holes at y=±80 mm)
    {380.0e6, 390.0e6, 3, Material::Cc},
    {390.0e6, 474.6e6, 4, Material::Rock},
  };

  for (const auto& s : slabs) {
    const double z0 = std::max(z_start_um, s.z0);
    const double z1 = std::min(z_end_um, s.z1);
    if (z1 <= z0) continue;

    double frac = 1.0;
    if (s.loct == 1 || s.loct == 2)
      frac = fraction_in_loct_com(s.loct, z0, z1, Beta_pair[0], Beta_pair[1], Beta_z);
    if (frac <= 0.0) continue;

    const double dz_eff_um = (z1 - z0) * frac;
    const double L_slant_um = (dz_eff_um / cos_theta_z);
    const double L_path_um = L_slant_um * osc_factor;

    if (!apply_material_range_step(tab, s.mat, L_path_um, mq, beta_min, p1, E1)) return false;
    if (!apply_material_range_step(tab, s.mat, L_path_um, mq, beta_min, p2, E2)) return false;
  }

  return true;
}

} // namespace RangeTables
