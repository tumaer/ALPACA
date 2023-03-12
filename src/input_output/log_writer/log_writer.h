//===-------------------------- log_writer.h ------------------------------===//
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
#ifndef LOG_WRITER_H
#define LOG_WRITER_H

#include "input_output/log_writer/log_writer_implementation.h"
#include <filesystem>

/**
 * @brief A light-weight logger to write output to the terminal and to a log
 * file. Messages are internally buffered and only written to external terminal
 * or file on command.
 * @note Singleton.
 */
class LogWriter {

  LogWriterImplementation implementation_;
  explicit LogWriter(std::unique_ptr<std::stringstream> &&terminal_output,
                     std::unique_ptr<std::stringstream> &&file_output);

public:
  // Singelton "Constructor":
  static LogWriter &
  Instance(std::unique_ptr<std::stringstream> &&terminal_output = nullptr,
           std::unique_ptr<std::stringstream> &&file_output = nullptr);

  // Singeltons may never call these methods.
  LogWriter() = delete;
  ~LogWriter() = default;
  LogWriter(LogWriter const &) = delete;
  LogWriter &operator=(LogWriter const &) = delete;
  LogWriter(LogWriter &&) = delete;
  LogWriter &operator=(LogWriter &&) = delete;

  void WelcomeMessage();
  void Flush();
  void FlushToFile();
  void FlushToTerminal();
  void SetLogfile(std::filesystem::path const &logfile);
  void RunningAlpaca(double const percentage, bool const fast_forward = false);
  void LogMessage(std::string const &message);
  void BufferMessage(std::string const &message);
  void LogBufferedMessages();
  void LogBreakLine();
};

#endif // LOG_WRITER_H
