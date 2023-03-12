//===--------------- instantiation_multiresolution.h ----------------------===//
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
#ifndef INITIALIZATION_MULTIRESOLUTION_H
#define INITIALIZATION_MULTIRESOLUTION_H

#include "input_output/input_reader.h"
#include "multiresolution/multiresolution.h"
#include "topology/topology_manager.h"

/**
 * @brief Defines all instantiation functions required for the multiresolution.
 */
namespace Instantiation {

// Instantiation function of the multiresolution
Multiresolution
InstantiateMultiresolution(InputReader const &input_reader,
                           TopologyManager const &topology_manager);
} // namespace Instantiation

#endif // INITIALIZATION_MULTIRESOLUTION_H
