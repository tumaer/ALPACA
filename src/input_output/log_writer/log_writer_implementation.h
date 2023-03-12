//===------------------- log_writer_implementation.h ----------------------===//
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
#include <filesystem>
#include <memory>
#include <sstream>

/**
 * @brief This class provides the functions to actually write logs to terminal
 * or file output. Logged messages are buffered internally until they are
 * actively flushed.
 * @note This class must not be used in production code. The LogWriter class is
 * provided for this purpose.
 */
class LogWriterImplementation {

  std::unique_ptr<std::stringstream> terminal_output_;
  std::unique_ptr<std::stringstream> file_output_;

  std::unique_ptr<std::filesystem::path> logfile_;

  std::string delayed_log_;

  void InsertMessageInAllStreams(std::string const &message);

public:
  explicit LogWriterImplementation(
      std::unique_ptr<std::stringstream> &&terminal_output = nullptr,
      std::unique_ptr<std::stringstream> &&file_output = nullptr);
  ~LogWriterImplementation() = default;
  LogWriterImplementation(LogWriterImplementation const &) = delete;
  LogWriterImplementation &operator=(LogWriterImplementation const &) = delete;
  LogWriterImplementation(LogWriterImplementation &&) = delete;
  LogWriterImplementation &operator=(LogWriterImplementation &&) = delete;

  void WelcomeMessage();
  void Flush();
  void FlushToFile();
  void FlushToTerminal();
  void SetLogfile(std::filesystem::path const &logfile);
  void RunningAlpaca(double const percentage, bool const fast_forward);
  void LogMessage(std::string const &message);
  void BufferMessage(std::string const &message);
  void LogBufferedMessages();
  void LogBreakLine();

  std::unique_ptr<std::stringstream> SwapOutTerminalOutputStream(
      std::unique_ptr<std::stringstream> &&new_in_old_out);
  std::unique_ptr<std::stringstream>
  SwapOutFileOutputStream(std::unique_ptr<std::stringstream> &&new_in_old_out);
};
