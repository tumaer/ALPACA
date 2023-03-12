//===------------------ vortex_stretching_output.h ------------------------===//
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
#ifndef VORTEX_STRETCHING_OUTPUT_H
#define VORTEX_STRETCHING_OUTPUT_H

#include "input_output/output_writer/output_quantity.h"
#include "topology/node.h"

/**
 * @brief The VortexStretchingOutput class handles the output of the vortex
 * stretching term of the vorticity transport equation into a output file
 * (currently Xdmf + HDF5). VortexStretchingOutput must not change any data.
 */
class VortexStretchingOutput : public OutputQuantity {

private:
  // Compute functions required from base class
  void
  DoComputeCellData(Node const &node, std::vector<double> &cell_data,
                    unsigned long long int &cell_data_counter) const override;
  void DoComputeDebugCellData(Node const &node, std::vector<double> &cell_data,
                              unsigned long long int &cell_data_counter,
                              MaterialName const material) const override;

public:
  VortexStretchingOutput() = delete;
  explicit VortexStretchingOutput(UnitHandler const &unit_handler,
                                  MaterialManager const &material_manager,
                                  std::string const &quantity_name,
                                  std::array<bool, 3> const output_flags);
  virtual ~VortexStretchingOutput() = default;
  VortexStretchingOutput(VortexStretchingOutput const &) = delete;
  VortexStretchingOutput &operator=(VortexStretchingOutput const &) = delete;
  VortexStretchingOutput(VortexStretchingOutput &&) = delete;
  VortexStretchingOutput &operator=(VortexStretchingOutput &&) = delete;
};

#endif // VORTEX_STRETCHING_OUTPUT_H
