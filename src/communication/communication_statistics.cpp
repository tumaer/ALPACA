//===------------------ communication_statistics.cpp ----------------------===//
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
#include "communication_statistics.h"
#include "mpi_utilities.h"
#include <mpi.h>

long CommunicationStatistics::no_jump_halos_recv_ = 0;
long CommunicationStatistics::no_jump_halos_send_ = 0;
long CommunicationStatistics::jump_halos_recv_ = 0;
long CommunicationStatistics::jump_halos_send_ = 0;
long CommunicationStatistics::balance_send_ = 0;
long CommunicationStatistics::balance_recv_ = 0;
long CommunicationStatistics::average_level_send_ = 0;
long CommunicationStatistics::average_level_recv_ = 0;

/**
 * @brief Sums the MPI statistic over all MPI ranks and gives a string
 * respresentation of the result.
 * @return .
 */
std::string SummedCommunicationStatisticsString() {

  // Logging Stats
  std::string statistics;
  statistics.append(" Ranks: " + std::to_string(MpiUtilities::NumberOfRanks()) +
                    " | ");

  long global_statistic;
  MPI_Allreduce(&CommunicationStatistics::balance_send_, &global_statistic, 1,
                MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  statistics.append(" Balance Send: " + std::to_string(global_statistic) +
                    " | ");

  MPI_Allreduce(&CommunicationStatistics::balance_recv_, &global_statistic, 1,
                MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  statistics.append(" Balance Recv: " + std::to_string(global_statistic) +
                    " | ");

  MPI_Allreduce(&CommunicationStatistics::no_jump_halos_send_,
                &global_statistic, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  statistics.append(" No Jump Halos Send: " + std::to_string(global_statistic) +
                    " | ");

  MPI_Allreduce(&CommunicationStatistics::no_jump_halos_recv_,
                &global_statistic, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  statistics.append(" No Jump Halos Recv: " + std::to_string(global_statistic) +
                    " | ");

  MPI_Allreduce(&CommunicationStatistics::jump_halos_send_, &global_statistic,
                1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  statistics.append(" Jump Halos Send: " + std::to_string(global_statistic) +
                    " | ");

  MPI_Allreduce(&CommunicationStatistics::jump_halos_recv_, &global_statistic,
                1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  statistics.append(" Jump Halos Recv: " + std::to_string(global_statistic) +
                    " | ");

  MPI_Allreduce(&CommunicationStatistics::average_level_send_,
                &global_statistic, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  statistics.append(" Proj.lvl-send: " + std::to_string(global_statistic) +
                    " | ");

  MPI_Allreduce(&CommunicationStatistics::average_level_recv_,
                &global_statistic, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  statistics.append(" Proj.lvl-recv: " + std::to_string(global_statistic) +
                    " | ");
  return statistics;
}
