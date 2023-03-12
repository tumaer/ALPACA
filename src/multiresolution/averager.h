//===--------------------------- averager.h -------------------------------===//
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
#ifndef AVERAGER_H
#define AVERAGER_H

#include "communication/communication_manager.h"
#include "topology/topology_manager.h"
#include "topology/tree.h"

/**
 * @brief Handles the propagation of information from fine to coarse nodes.
 */
class Averager {

private:
  TopologyManager const &topology_;
  CommunicationManager &communicator_;
  Tree &tree_;

public:
  Averager() = delete;
  explicit Averager(TopologyManager const &topology,
                    CommunicationManager &communicator, Tree &tree);
  ~Averager() = default;
  Averager(Averager const &) = delete;
  Averager &operator=(Averager const &) = delete;
  Averager(Averager &&) = delete;
  Averager &operator=(Averager &&) = delete;

  void AverageMaterial(
      std::vector<unsigned int> const &child_levels_descending) const;

  void AverageParameters(
      std::vector<unsigned int> const child_levels_descending) const;

  void AverageInterfaceTags(std::vector<unsigned int> const
                                &levels_with_updated_parents_descending) const;
};

#endif // AVERAGER_H
