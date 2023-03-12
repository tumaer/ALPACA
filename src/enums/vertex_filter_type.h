//===---------------------- vertex_filter_type.h --------------------------===//
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
#ifndef VERTEX_FILTER_TYPE_H
#define VERTEX_FILTER_TYPE_H

#include <string>

/**
 * @brief Identifier for the mode of output vertex filter. Type 'mpi' indicates
 * that overlapping ids are sorted out via MPI-communication and 'FinestLevel'
 * indicates that the ids correspond to the IDs a mesh would have if all nodes
 * would lie on the finest level. 'Off' represents the same approach as for the
 * 'Mpi' type but without the mpi-communication for the vertex filter.
 */
enum class VertexFilterType { Off, Mpi, FinestLevel };

/**
 * @brief Converts the vertex filter type into an appropriate string.
 * @param filter VertexFilterType identifier.
 * @return  String for the filter type.
 */
inline std::string VertexFilterTypeToString(VertexFilterType const filter) {
  switch (filter) {
  case VertexFilterType::Off:
    return "Off";
  case VertexFilterType::Mpi:
    return "Mpi";
  case VertexFilterType::FinestLevel:
    return "Finest Level";
  default:
    return "ERROR: This vertex filter is not (yet) defined!";
  }
}

#endif // VERTEX_FILTER_TYPE_H
