//===------------------------- file_utilities.h ---------------------------===//
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
#ifndef FILE_UTILITIES_H
#define FILE_UTILITIES_H

#include <filesystem>
#include <string>

namespace FileUtilities {

// operations for file and path modification
std::string RemoveFileExtension(std::string const &filename);
std::string RemoveFilePath(std::string const &filename);
std::string ChangeFileExtension(std::string const &filename_with_path,
                                std::string const extension);
std::string GetFileExtension(std::string const &filename);
bool CheckIfPathExists(std::string const &path);
bool CreateFolder(std::string const &path);
std::string AddUnusedNumberToPath(std::string const &path,
                                  unsigned int const start_number = 0);
void WriteTextBasedFile(std::string const &filename,
                        std::string const &content);
void AppendToTextBasedFile(std::string const &filename,
                           std::string const &content);
std::filesystem::path
CreateOutputBaseFolder(std::filesystem::path const &inputfile);

} // namespace FileUtilities

#endif // FILE_UTILITIES_H
