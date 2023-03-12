//===------------------ datatype_for_mpi_definition.h ---------------------===//
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
#ifndef DATATYPE_FOR_MPI_DEFINITION_H
#define DATATYPE_FOR_MPI_DEFINITION_H

#include <type_traits>

/**
 * @brief Identifier to obtain the correct MPI_Datatype during MPI calls via the
 * CommunicationTypes proxy.
 */
enum class DatatypeForMpi : unsigned short { Double = 0, Byte = 1 };

/**
 * @brief Converts a DatatypeForMpi identifier to a (C++11 standard compliant,
 * i. e. positive) array index. "DTI = Datatype To Index".
 * @param d The Mpi Datatype identifier.
 * @return Index to be used in Arrays.
 */
constexpr std::underlying_type<DatatypeForMpi>::type
DTI(DatatypeForMpi const d) {
  return static_cast<typename std::underlying_type<DatatypeForMpi>::type>(d);
}

#endif // DATATYPE_FOR_MPI_DEFINITION_H
