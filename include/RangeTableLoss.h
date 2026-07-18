#pragma once

#include <vector>

namespace RangeTables {

bool apply_range_table_loss(int mq_int, double mq, double Lambda_eV, double z_end_um,
                            double distance_period_m, const std::vector<double>& Beta_pair,
                            double beta_min, std::vector<double>& p1, double& E1,
                            std::vector<double>& p2, double& E2, double& t_arrival_ns_out);
bool apply_range_table_loss_zspan(int mq_int, double mq, double Lambda_eV,
                                  double z_start_um, double z_end_um,
                                  double distance_period_m,
                                  const std::vector<double>& Beta_pair, double beta_min,
                                  std::vector<double>& p1, double& E1,
                                  std::vector<double>& p2, double& E2);

} // namespace RangeTables
