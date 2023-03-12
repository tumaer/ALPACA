//===----------------------- bit_operations.h -----------------------------===//
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
#ifndef BIT_OPERATIONS_H
#define BIT_OPERATIONS_H

#include <bitset>
#include <cstddef>

namespace BitOperations {

/**
 * @brief Performs a right circular shift on the given bitset.
 * @param argument The bitset to rotate.
 * @param rotations The amount of rotations to perfrom.
 * @return Rotated bitset.
 * @note If the number of rotations is larger than the twice arguments's size it
 * return an bitset with all bits false. This is different to C++20 rotr. The
 * difference is unintentional. Do not rely on it.
 */
template <std::size_t N>
std::bitset<N> RightCircularShift(std::bitset<N> const argument,
                                  std::size_t const rotations) {
  return std::bitset<N>((argument >> rotations) |
                        (argument << (N - rotations)));
}

/**
 * @brief Extracts a subset of bits from a larger bitset. The sub-bitset holds
 * entries from start until start + N of the input bitset.
 * @param long_bitset The bitset from which a subset is to be extracted.
 * @param start The postion in long_bitset from which the subset starts.
 * @tparam N length of the subbitset.
 * @tparam M length of the longer bitset.
 */
template <std::size_t N, std::size_t M>
std::bitset<N> SubBitset(std::bitset<M> const long_bitset, std::size_t start) {
  static_assert(N <= M,
                "Sub-bitset may not be larger the original input bitset.");
  std::bitset<N> subbitset;
  for (std::size_t i = 0; i < N; ++i) {
    subbitset[i] = long_bitset[start + i];
  }
  return subbitset;
}

} // namespace BitOperations

#endif // BIT_OPERATIONS_H
