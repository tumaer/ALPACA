//===------------------------- xml_utilities.h ----------------------------===//
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
#ifndef XML_UTILITIES_H
#define XML_UTILITIES_H

#include <string>
#include <tinyxml2.h>
#include <vector>

/**
 * @brief In this name space all functions are defiend that simplify the reading
 * process of xml input files.
 */
namespace XmlUtilities {

tinyxml2::XMLElement const *GetChild(tinyxml2::XMLDocument const &xml_document,
                                     std::vector<std::string> child_names);
tinyxml2::XMLElement const *GetChild(tinyxml2::XMLElement const *parent_node,
                                     std::vector<std::string> child_names);
std::vector<tinyxml2::XMLElement const *>
GetChilds(tinyxml2::XMLElement const *parent_node,
          std::string const &child_name);
bool ChildExists(tinyxml2::XMLDocument const &xml_document,
                 std::vector<std::string> child_names);
bool ChildExists(tinyxml2::XMLElement const *parent_node,
                 std::string child_name);
double ReadDouble(tinyxml2::XMLElement const *node);
int ReadInt(tinyxml2::XMLElement const *node);
unsigned int ReadUnsignedInt(tinyxml2::XMLElement const *node);
std::int64_t ReadLongInt(tinyxml2::XMLElement const *node);
std::uint64_t ReadUnsignedLongInt(tinyxml2::XMLElement const *node);
std::string ReadString(tinyxml2::XMLElement const *node);
std::vector<double> ReadTimeStamps(tinyxml2::XMLElement const *parent_node);

} // namespace XmlUtilities

#endif // XML_UTILITIES_H
