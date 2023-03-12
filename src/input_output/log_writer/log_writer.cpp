//===------------------------- log_writer.cpp -----------------------------===//
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
#include "input_output/log_writer/log_writer.h"

/**
 * @brief Default constructor using the provided streams to buffer the messages
 * and reduced amount of external terminal of logfile writes.
 * @param terminal_output stream to use for terminal logging.
 * @param file_output stream to use for logfile logging.
 * @note Takes ownership over the provided streams.
 */
LogWriter::LogWriter(std::unique_ptr<std::stringstream> &&terminal_output,
                     std::unique_ptr<std::stringstream> &&file_output)
    : implementation_(std::move(terminal_output), std::move(file_output)) {
  // Empty besides initializer list.
}

/**
 * @brief This function is used to get the Logger. If no Logger exists yet it is
 * created, otherwise the existing logger is passed back. "Singleton
 * Constructor"
 * @param terminal_output stream to use for terminal logging.
 * @param file_output stream to use for logfile logging.
 * @return The logger instance.
 */
LogWriter &
LogWriter::Instance(std::unique_ptr<std::stringstream> &&terminal_output,
                    std::unique_ptr<std::stringstream> &&file_output) {
  static LogWriter instance(std::move(terminal_output), std::move(file_output));
  return instance;
}

/**
 * @brief Logs the welcome (Greeting) message.
 */
void LogWriter::WelcomeMessage() { implementation_.WelcomeMessage(); }

/**
 * @brief Triggers writing of logged messages to terminal.
 */
void LogWriter::FlushToTerminal() { implementation_.FlushToTerminal(); }

/**
 * @brief Triggers writing of logged messages to file.
 */
void LogWriter::FlushToFile() { implementation_.FlushToFile(); }

/**
 * @brief Triggers writing of logged messages to terminal and/or file.
 */
void LogWriter::Flush() { implementation_.Flush(); }

/**
 * @brief Redirects logs to the specified file and logs the conducted change.
 * @param logfile New file (path) used as logfile.
 */
void LogWriter::SetLogfile(std::filesystem::path const &logfile) {
  implementation_.SetLogfile(logfile);
}

/**
 * @brief Logs the ascii Alpaca. The position of the Alpaca implies a progress
 * bar.
 * @param percentage Indicator how much progress is to be indicated.
 * @param fast_forward Flag whether the progress 'jumps' to the given position.
 */
void LogWriter::RunningAlpaca(double const percentage,
                              bool const fast_forward) {
  implementation_.RunningAlpaca(percentage, fast_forward);
}

/**
 * @brief Logs (and formats) the provided message
 * @param message The (unformatted) message.
 */
void LogWriter::LogMessage(std::string const &message) {
  implementation_.LogMessage(message);
}

/**
 * @brief Buffers the provided message without formatting it, so it can be
 * logged (and formatted) later.
 * @param message The (unformatted) message.
 */
void LogWriter::BufferMessage(std::string const &message) {
  implementation_.BufferMessage(message);
}

/**
 * @brief Logs the so far buffered messages. Thereby the buffered messages are
 * concatenated and treated as one message.
 */
void LogWriter::LogBufferedMessages() { implementation_.LogBufferedMessages(); }

/**
 * @brief Logs a full line of starts to indicated a sections.
 */
void LogWriter::LogBreakLine() { implementation_.LogBreakLine(); }
