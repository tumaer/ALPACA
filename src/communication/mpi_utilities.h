//===------------------------ mpi_utilities.h -----------------------------===//
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
#ifndef MPI_UTILITIES_H
#define MPI_UTILITIES_H

#include <mpi.h>
#include <numeric>
#include <vector>

namespace MpiUtilities {

/**
 * @brief Reduces a bool across MPI ranks.
 * @param input local bool.
 * @param operation The redcution operation. Default: or.
 * @return Global results.
 */
inline bool GloballyReducedBool(bool const input,
                                MPI_Op const operation = MPI_LOR) {
  bool result = input;
  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_CXX_BOOL, operation,
                MPI_COMM_WORLD);
  return result;
}

/**
 * @brief Gives the rank id from MPI directly as int. Avoids handle creation,
 * e.g. for const members in initializer list.
 * @return Rank id.
 */
inline int MyRankId() {
  int rank_id = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
  return rank_id;
}

/**
 * @brief Indicates whether the invoking rank is the master rank, i.e. rank 0.
 * @return True if invoking is master, false otherwise.
 */
inline bool MasterRank() {
  int rank_id = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
  return rank_id == 0;
}

/**
 * @brief Gives the number of ranks in the MPI communicator "MPI_COMM_WORLD".
 * Avoids handle creation, e.g. for const members in initializer list.
 * @return Communicator Size which is the number of ranks.
 */
inline int NumberOfRanks() {
  int communicator_size = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &communicator_size);
  return communicator_size;
}

/**
 * @brief Gives the MPI_TAG_UB. Avoids handle creation, e.g. for const members
 * in initializer list.
 * @return MPI_TAG_UB
 */
inline int MpiTagUb() {
  int *tag_ub;
  int flag;
  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &flag);
  return *tag_ub;
}

/**
 * @brief Wrapper function to collect data from a local (= on one MPI rank)
 * vector into a large global (= data of all ranks) one via a gatherv operation.
 * @param local_data The data present at this rank.
 * @param type The MPI datatype to be used in the gather calls.
 * @param number_of_ranks The number of ranks in the communicator.
 * @param global_data Vector holding the collected data from all ranks (indirect
 * return parameter).
 * @tparam Data type.
 * @note Does not perform sanity checks. If the template type and the MPI
 * datatype do not match the results will be corrupted. Uses MPI_COMM_WORLD as
 * communicator. Overrides the provided global_data array.
 */
template <class T>
void LocalToGlobalData(std::vector<T> const &local_data,
                       MPI_Datatype const type, int const number_of_ranks,
                       std::vector<T> &global_data) {
  int length = local_data.size(); // Must be int due to MPI standard.
  std::vector<int> all_lengths(number_of_ranks);
  MPI_Allgather(&length, 1, MPI_INT, all_lengths.data(), 1, MPI_INT,
                MPI_COMM_WORLD);

  std::vector<int> offsets(number_of_ranks);
  int insert_key = 0;
  for (int i = 0; i < number_of_ranks; ++i) {
    offsets[i] = insert_key;
    insert_key += all_lengths[i];
  }

  global_data.resize(
      std::accumulate(all_lengths.begin(), all_lengths.end(), 0));
  MPI_Allgatherv(local_data.data(), length, type, global_data.data(),
                 all_lengths.data(), offsets.data(), type, MPI_COMM_WORLD);
}
} // namespace MpiUtilities

#endif // MPI_UTILITIES_H
