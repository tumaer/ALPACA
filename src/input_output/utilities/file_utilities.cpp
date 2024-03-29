//===----------------------- file_utilities.cpp ---------------------------===//
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
#include "input_output/utilities/file_utilities.h"

#include <fstream>
#include <sstream>
#include <sys/stat.h>

namespace FileUtilities {

/**
 * @brief Removes the file extension from the given filename.
 * @param filename.
 * @return The name of the file without extension.
 */
std::string RemoveFileExtension(std::string const &filename) {
  return filename.substr(0, filename.find_last_of("."));
}

/**
 * @brief Removes the path from the given filename.
 * @param filename.
 * @return The name of the file without path.
 */
std::string RemoveFilePath(std::string const &filename) {
  return filename.substr(filename.find_last_of("/") + 1);
}

/**
 * @brief Changes the extension of a given file. Can handle filenames with full
 * (absolute or relative) path.
 * @param filename_with_path .
 * @param extension The future extension of the file.
 * @return New filename (with path if input had path).
 */
std::string ChangeFileExtension(std::string const &filename_with_path,
                                std::string const extension) {
  auto const positon_of_last_dot = filename_with_path.find_last_of('.');
  return filename_with_path.substr(0, positon_of_last_dot) + extension;
}

/**
 * @brief Gives the extension of a given file. Can handle filenames with full
 * (absolute or relative) path.
 * @param filename Name of the file.
 * @return Extension of the file (without dot).
 */
std::string GetFileExtension(std::string const &filename) {
  auto const positon_of_last_dot = filename.find_last_of('.');
  return filename.substr(positon_of_last_dot + 1);
}

/**
 * @brief Checks if the given path already exists.
 * @param path The name of the path, whose existence is to be checked.
 * @return "true" if the path exists, otherwise false.
 */
bool CheckIfPathExists(std::string const &path) {
  struct stat info;
  return (stat(path.c_str(), &info) == 0);
}

/**
 * @brief Creates a folder with the name "path".
 * @param path The name of the folder to be created.
 * @return "true" if the folder "path" was successfully created, otherwise
 * false.
 */
bool CreateFolder(std::string const &path) {
  mode_t mode = 0755;
  int ret = mkdir(path.c_str(), mode);
  return (ret == 0);
}

/**
 * @brief Finds and adds the next free integer-counter to the path.
 * $THIS FUNCTION IS NOT THREAD SAFE, i.e. if one thread creates a folder in
 * between calls, the results are different$.
 * @param path The candidate for the path name.
 * @param start_number The next number to be checked (0 by default).
 * @return The path name that is not yet existing.
 */
std::string AddUnusedNumberToPath(std::string const &path,
                                  unsigned int const start_number) {

  std::string target;
  if (start_number == 0) {
    target = path;
  } else {
    target = path + "-" + std::to_string(start_number);
  }

  if (CheckIfPathExists(target)) {
    target = AddUnusedNumberToPath(path, start_number + 1);
  }

  return target;
}

/**
 * @brief Writes a file containing human readable (ascii) text to disk.
 * @param filename The name of the file to be created.
 * @param content The text to be written into the file.
 */
void WriteTextBasedFile(std::string const &filename,
                        std::string const &content) {
  std::ofstream output_stream(filename, std::ios::trunc);
  output_stream << content;
  output_stream.flush();
  output_stream.close();
}

/**
 * @brief Appends content to an existing file containing human readable (ascii)
 * text to disk.
 * @param filename The name of the file to be appended.
 * @param content The text to be written into the file.
 */
void AppendToTextBasedFile(std::string const &filename,
                           std::string const &content) {
  std::ofstream output_stream(filename, std::ios::app);
  output_stream << content;
  output_stream.flush();
  output_stream.close();
}

} // namespace FileUtilities
