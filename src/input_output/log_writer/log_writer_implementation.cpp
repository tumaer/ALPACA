//===----------------- log_writer_implementation.cpp ----------------------===//
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
#include "input_output/log_writer/log_writer_implementation.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

/**
 * @brief Collection strings that re-occur during the logging.
 */
namespace MessageBlocks {
// Some time in the future compiles will support constexpr  string, until then
// we stick with boring normal const
std::string const break_line = "|**********************************************"
                               "**************************************|\n";
namespace Welcome {
std::string const ascii_alpaca = break_line +
                                 "|*  \\\\                                     "
                                 "                                         *|\n"
                                 "|*  l '>                                     "
                                 "                                       *|\n"
                                 "|*  | |                                      "
                                 "                                       *|\n"
                                 "|*  | |                                      "
                                 "                                       *|\n"
                                 "|*  | alpaca~                                "
                                 "                                       *|\n"
                                 "|*  ||    ||                                 "
                                 "                                       *|\n"
                                 "|*  ''    ''                                 "
                                 "                                       *|\n"
                                 "|*                                           "
                                 "                                       *|\n";
std::string const greetings_message_terminal =
    "|*  THE AGE OF ALPACA HAS COME                                            "
    "          *|\n";
std::string const greetings_message_file =
    "|* TERMINAL OUTPUT HAS GONE - THE ALPACA LOGBOOK HAS COME                 "
    "          *|\n";
std::string const end_of_ascii_alpaca_frame =
    "|*                                                                        "
    "          *|\n" +
    break_line;
} // namespace Welcome
} // namespace MessageBlocks

/**
 * @brief Collection of constants used in the logging.
 */
namespace Constants {
constexpr std::size_t line_width = 80;
}

namespace {
std::vector<std::string> FormatMessage(std::string const &message) {
  unsigned int message_size = message.size();
  unsigned int number_of_lines =
      std::ceil(double(message_size) / double(Constants::line_width + 1));

  std::vector<std::string> final_lines;

  for (unsigned int i = 0; i < number_of_lines; ++i) {
    if (message.size() <= i * Constants::line_width + Constants::line_width) {
      final_lines.emplace_back(message.begin() + (i * Constants::line_width),
                               message.end());
    } else {
      final_lines.emplace_back(message.begin() + (i * Constants::line_width),
                               message.begin() + (i * Constants::line_width) +
                                   Constants::line_width);
    }
  }

  for (auto &line : final_lines) {
    if (line.size() < Constants::line_width) { //-1
      int difference = Constants::line_width - line.size();
      line.insert(line.end(), difference, ' '); //+1
    }
    line.insert(line.end(), ' ');
    line.insert(line.end(), '*'); // insert empty space
    line.insert(line.end(), '|');
    line.insert(0, "|* ");
  }

  return final_lines;
}

std::vector<std::string>
SplitMessageAccordingToLinebreaks(std::string const &message) {
  std::vector<std::string> line_broken_strings;
  if (message.find('\n') != std::string::npos) {
    // Instantiate vector containing split strings if needed
    std::string line;
    std::stringstream stream(message);
    while (std::getline(stream, line, '\n')) {
      line_broken_strings.push_back(line);
    }
  } else {
    line_broken_strings.push_back(message);
  }
  return line_broken_strings;
}

} // namespace

/**
 * @brief Constructor using the input streams for writing of output to terminal
 * and file.
 * @param terminal_output Stream to be used to buffer formatted mesages destined
 * to end in terminal output.
 * @param file_output Stream to be used to buffer formatted mesages destined to
 * end in logfile output.
 * @note Takes over ownership of the provided streams.
 */
LogWriterImplementation::LogWriterImplementation(
    std::unique_ptr<std::stringstream> &&terminal_output,
    std::unique_ptr<std::stringstream> &&file_output)
    : terminal_output_(std::move(terminal_output)),
      file_output_(std::move(file_output)) {
  // Empty besides initializer list
}

/**
 * @brief Swaps the terminal output stream with the provided one. Takes
 * ownership of the provided stream and returns the old one.
 * @param new_output_stream Stream to be used to buffer formatted mesages
 * destined to end in terminal output.
 * @return Returns the previously used (and potentially filled) terminal stream.
 */
std::unique_ptr<std::stringstream>
LogWriterImplementation::SwapOutTerminalOutputStream(
    std::unique_ptr<std::stringstream> &&new_output_stream) {
  std::unique_ptr<std::stringstream> result(terminal_output_.release());
  terminal_output_ = std::move(new_output_stream);
  return result;
}

/**
 * @brief Swaps the logfile output stream with the provided one. Takes ownership
 * of the provided stream and returns the old one.
 * @param new_output_stream Stream to be used to buffer formatted mesages
 * destined to end in logfile output.
 * @return Returns the previously used (and potentially filled) logfile stream.
 */
std::unique_ptr<std::stringstream>
LogWriterImplementation::SwapOutFileOutputStream(
    std::unique_ptr<std::stringstream> &&new_output_stream) {
  std::unique_ptr<std::stringstream> result(file_output_.release());
  file_output_ = std::move(new_output_stream);
  return result;
}

/**
 * @brief Flushes the (so far) collected messages to the terminal. Thus external
 * writes are limited to calls to this function.
 * @note This method cannot be unit-tested. By providing nullptr for the stream
 * members no logging text is propagated to the outside, hence the logger is
 * silenced.
 */
void LogWriterImplementation::FlushToTerminal() {
  if (terminal_output_) {
    std::cout << terminal_output_->str();
    std::cout.flush();
    terminal_output_->clear();
    terminal_output_->str(std::string());
  }
}

/**
 * @brief Flushes the (so far) collected messages to the logfile. Thus external
 * writes are limited to calls to this function.
 * @note This method cannot be unit-tested. By providing nullptr for the stream
 * members no logging text is propagated to the outside, hence the logger is
 * silenced.
 */
void LogWriterImplementation::FlushToFile() {
  if (file_output_ && logfile_) {
    std::ofstream file_output_stream(*logfile_, std::ios::app);
    file_output_stream << file_output_->str();
    file_output_stream.flush();
    file_output_stream.close();
    file_output_->clear();
    file_output_->str(std::string());
  }
}

/**
 * @brief Flushes the (so far) collected messages to the terminal and the
 * logfile respectively. Thus external writes are limited to calls to this
 * function.
 * @note This method cannot be unit-tested. By providing nullptr for the stream
 * members no logging text is propagated to the outside, hence the logger is
 * silenced.
 */
void LogWriterImplementation::Flush() {
  FlushToFile();
  FlushToTerminal();
}

/**
 * @brief Logs the (ascii-art) welcome message.
 */
void LogWriterImplementation::WelcomeMessage() {
  if (terminal_output_) {
    *terminal_output_
        << MessageBlocks::Welcome::ascii_alpaca +
               MessageBlocks::Welcome::greetings_message_terminal +
               MessageBlocks::Welcome::end_of_ascii_alpaca_frame;
  }
  if (file_output_) {
    *file_output_ << MessageBlocks::Welcome::ascii_alpaca +
                         MessageBlocks::Welcome::greetings_message_file +
                         MessageBlocks::Welcome::end_of_ascii_alpaca_frame;
  }
}

/**
 * @brief Redirects logs to the specified file and logs the conducted change.
 * @param logfile New file (path) used as logfile.
 */
void LogWriterImplementation::SetLogfile(std::filesystem::path const &logfile) {
  auto const old_name = logfile_ ? logfile_->string() : "Not Set";
  logfile_ = std::make_unique<std::filesystem::path>(logfile);
  LogMessage("Logger: Changed log file from: " + old_name + " to " +
             logfile_->string());
}

/**
 * @brief Inserts a message into the log streams without any checks.
 * @param message The message to be logged.
 */
void LogWriterImplementation::InsertMessageInAllStreams(
    std::string const &message) {
  if (terminal_output_) {
    *terminal_output_ << message;
  }
  if (file_output_) {
    *file_output_ << message;
  }
}

/**
 * @brief Logs the ascii Alpaca in different positions accroding to the given
 * percentage.
 * @param percentage Indicator how far to the right the alpaca is to be
 * positioned. Value as decimal [0,1.0]. Smaller/higher values are rounded to 0
 * or 1.0, respectively.
 * @param fast_forward Flag whether the Alpaca 'jumps' to the given position. If
 * set the jump is indicated by movement lines.
 */
void LogWriterImplementation::RunningAlpaca(double const percentage,
                                            bool const fast_forward) {

  double const cut_percentage = std::min(1.0, std::max(0.0, percentage));
  std::size_t const plot_percentage = static_cast<std::size_t>(
      std::floor(cut_percentage * (Constants::line_width - 10) + 9));

  InsertMessageInAllStreams(MessageBlocks::break_line);
  LogMessage(" ");
  // The minus is the length of the string
  LogMessage(std::string(plot_percentage - 3, ' ') + "\\\\");
  LogMessage(std::string(plot_percentage - 4, fast_forward ? '>' : ' ') +
             " l '>"); // except here, because the snout surpasses to the right
  LogMessage(std::string(plot_percentage - 3, ' ') + "| |");
  LogMessage(std::string(plot_percentage - 4, fast_forward ? '>' : ' ') +
             " | |");
  LogMessage(std::string(plot_percentage - 9, ' ') + "~alpaca |");
  LogMessage(std::string(plot_percentage - 9, fast_forward ? '>' : ' ') +
             " ||    ||");
  LogMessage(std::string(plot_percentage - 9, ' ') + " ''    ''");
  LogMessage(" ");
  InsertMessageInAllStreams(MessageBlocks::break_line);
}

/**
 * @brief Logs the given message. I.e. applies the logger format to the given
 * message and buffers them.
 * @param message The unformatted message.
 */
void LogWriterImplementation::LogMessage(std::string const &message) {

  auto const split_message = SplitMessageAccordingToLinebreaks(message);

  // Loop through all single_messages and call single send message function
  for (auto const &single_message : split_message) {
    std::vector<std::string> formatted(
        FormatMessage(std::move(single_message)));
    for (auto const &line : formatted) {
      if (terminal_output_) {
        *terminal_output_ << std::scientific << std::setprecision(5) << line
                          << std::endl;
      }
      if (file_output_) {
        *file_output_ << std::scientific << std::setprecision(5) << line
                      << std::endl;
      }
    }
  }
}

/**
 * @brief Buffers the given message. Allows, i.e. to concatenate messages
 * without formatting them up-front.
 * @param message The message to be buffered.
 */
void LogWriterImplementation::BufferMessage(std::string const &message) {
  delayed_log_ += message;
}

/**
 * @brief Logs the previously buffered messages. I.e. treats all buffered
 * messages as one (long) message and applies the logging format to this long
 * message.
 */
void LogWriterImplementation::LogBufferedMessages() {
  LogMessage(delayed_log_);
  delayed_log_.clear();
}

/**
 * @brief Logs a line breaker, i.e. a full line of stars.
 */
void LogWriterImplementation::LogBreakLine() {
  InsertMessageInAllStreams(MessageBlocks::break_line);
}
