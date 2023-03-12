//===------------------- instantiation_log_writer.h -----------------------===//
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
#ifndef INSTANTIATION_LOG_WRITER_H
#define INSTANTIATION_LOG_WRITER_H

#include "input_output/log_writer/log_writer.h"

/**
 * @brief Defines all instantiation functions required for the Logger.
 */
namespace Instantiation {

LogWriter &InstantiateLogWriter(bool const on_master_rank);

}

#endif // INSTANTIATION_LOG_WRITER_H
