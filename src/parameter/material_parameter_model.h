//===--------------------- material_parameter_model.h ---------------------===//
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
#ifndef MATERIAL_PARAMETER_MODEL_H
#define MATERIAL_PARAMETER_MODEL_H

#include "block_definitions/block.h"
#include "materials/material_definitions.h"

/**
 * @brief The MaterialPropertyModel class defines an interface for different
 * discretizations of a model to compute material and material pairing
 * parameters depending on other materials or physical parameter. It works a
 * proxy to define specific models for different material properties.
 */
class MaterialParameterModel {

protected:
  // required functions needed in derived classes
  virtual void DoUpdateParameter(Block &block,
                                 double const cell_size) const = 0;
  virtual void DoUpdateParameter(
      Block &block, double const cell_size,
      std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()],
      std::int8_t const material_sign) const = 0;

  // protected default constructor (can only be called from derived classes)
  explicit MaterialParameterModel() = default;

public:
  virtual ~MaterialParameterModel() = default;
  MaterialParameterModel(MaterialParameterModel const &) = delete;
  MaterialParameterModel &operator=(MaterialParameterModel const &) = delete;
  MaterialParameterModel(MaterialParameterModel &&) = delete;
  MaterialParameterModel &operator=(MaterialParameterModel &&) = delete;

  /**
   * @brief Computes the desired parameter based on the given parameter for a
   * complete block.
   * @param block Block on which parameters are calculated (indirect return).
   * @param cell_size Size of the cell of given block.
   */
  void UpdateParameter(Block &block, double const cell_size) const {
    DoUpdateParameter(block, cell_size);
  }

  /**
   * @brief Computes the desired parameter based on the given parameter of an
   * interface-containing block.
   * @param mat_block Pair of MaterialName and corresponding block on which
   * parameters are computed (indirect return).
   * @param cell_size Size of the cell of given block.
   * @param interface_tags Interface tags to specify location of the interface
   * on the given block.
   * @param material_sign Sign of the material for identification on interface
   * tags.
   */
  void UpdateParameter(
      Block &block, double const cell_size,
      std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()],
      std::int8_t const material_sign) const {
    DoUpdateParameter(block, cell_size, interface_tags, material_sign);
  }
};

#endif // MATERIAL_PARAMETER_MODEL_H
