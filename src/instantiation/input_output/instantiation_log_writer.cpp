//===------------------ instantiation_log_writer.cpp ----------------------===//
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
#include "instantiation/input_output/instantiation_log_writer.h"
#include "input_output/log_writer/log_writer.h"
#include <memory>
#include <sstream>
#include <stdexcept>

namespace {
/**
 * @brief Instantiates a 'talking' logger on the master rank and 'silent' ones
 * on all other ranks.
 * @param on_master_rank Indicator wether the function is executed on the master
 * rank.
 * @note Throws ('dies') if a logger could not be instantiated.
 */
LogWriter &CreateLogWriterOrDieTrying(bool const on_master_rank) {
  try {
    std::unique_ptr<std::stringstream> terminal_stream(nullptr);
    std::unique_ptr<std::stringstream> file_stream(nullptr);
    if (on_master_rank) {
      terminal_stream = std::make_unique<std::stringstream>(std::ios_base::out);
      file_stream = std::make_unique<std::stringstream>(std::ios_base::out);
    }
    return LogWriter::Instance(std::move(terminal_stream),
                               std::move(file_stream));
  } catch (std::exception const &exception) {
    std::string message = "Could not instantiate a Logger. Reason: " +
                          std::string(exception.what());
    throw std::logic_error(message);
  }
}
} // namespace

namespace Instantiation {
/**
 * @brief Instantiates the logger. Only the logger. Only the logger on the
 * master rank is instantiated such that it writes to terminal and/or logfile.
 * @param on_master_rank Indicator wether the function is executed on the master
 * rank.
 */
LogWriter &InstantiateLogWriter(bool const on_master_rank) {
  LogWriter &logger = CreateLogWriterOrDieTrying(on_master_rank);
  logger.WelcomeMessage();
  logger.LogMessage("Logger initialised");
  logger.Flush();
  return logger;
}
} // namespace Instantiation
