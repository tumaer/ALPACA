//===------------------- communication_statistics.h -----------------------===//
//
//                                 ALPACA
//
// Part of ALPACA, under the GNU General Public License as published by
// the Free Software Foundation version 3.
// SPDX-License-Identifier: GPL-3.0-only
//
// If using this code in an academic setting, please cite the following:
// @article{hoppe2022parallel,
//  title={A parallel modular computing environment for three-dimensional
//  multiresolution simulations of compressible flows},
//  author={Hoppe, Nils and Adami, Stefan and Adams, Nikolaus A},
//  journal={Computer Methods in Applied Mechanics and Engineering},
//  volume={391},
//  pages={114486},
//  year={2022},
//  publisher={Elsevier}
// }
//
//===----------------------------------------------------------------------===//
#ifndef COMMUNICATION_STATISTICS_H
#define COMMUNICATION_STATISTICS_H

#include <string>

/**
 * @brief The CommunicationStatistics struct gathers information about the MPI
 * Communication for proper logging output.
 */
struct CommunicationStatistics {
public:
  static long no_jump_halos_recv_;
  static long no_jump_halos_send_;
  static long jump_halos_recv_;
  static long jump_halos_send_;
  static long balance_send_;
  static long balance_recv_;
  static long average_level_send_;
  static long average_level_recv_;
};

std::string SummedCommunicationStatisticsString();

#endif /* COMMUNICATION_STATISTICS_H */
