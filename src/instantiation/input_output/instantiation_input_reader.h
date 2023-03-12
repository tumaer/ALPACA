//===----------------- instantiation_input_reader.h -----------------------===//
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
#ifndef INSTANTIATION_INPUT_READER_H
#define INSTANTIATION_INPUT_READER_H

#include "input_output/input_reader.h"

/**
 * @brief Defines all instantiation functions required for the input reader.
 */
namespace Instantiation {
// Instantiation function for the input reader
InputReader InstantiateInputReader(std::string const &input_filename);
} // namespace Instantiation

#endif // INSTANTIATION_INPUT_READER_H
