//===------------------- interface_parameter_model.h ----------------------===//
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
#ifndef INTERFACE_PARAMETER_MODEL_H
#define INTERFACE_PARAMETER_MODEL_H

#include "topology/node.h"

/**
 * @brief The InterfacePropertyModel class defines an interface for different
 * discretizations of a model to compute material and material pairing
 * parameters depending on other materials or physical parameter. It works a
 * proxy to define specific models for different material properties.
 */
class InterfaceParameterModel {

protected:
  // required functions needed in derived classes
  virtual void DoUpdateParameter(Node &node) const = 0;

  // protected default constructor (can only be called from derived classes)
  explicit InterfaceParameterModel() = default;

public:
  virtual ~InterfaceParameterModel() = default;
  InterfaceParameterModel(InterfaceParameterModel const &) = delete;
  InterfaceParameterModel &operator=(InterfaceParameterModel const &) = delete;
  InterfaceParameterModel(InterfaceParameterModel &&) = delete;
  InterfaceParameterModel &operator=(InterfaceParameterModel &&) = delete;

  /**
   * @brief Computes the desired parameter based on the given parameter of an
   * interface-containing block.
   * @param node Node on which the parameter is calculated.
   */
  void UpdateParameter(Node &node) const { DoUpdateParameter(node); }
};

#endif // INTERFACE_PARAMETER_MODEL_H
