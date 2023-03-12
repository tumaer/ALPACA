//===----------------- xml_multi_resolution_reader.cpp --------------------===//
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
#include "input_output/input_reader/multi_resolution_reader/xml_multi_resolution_reader.h"

#include "input_output/utilities/xml_utilities.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief Default constructor for the multiresolution reader for xml-type input
 * files.
 * @param inputfile The xml input file document holding all information of the
 * user inputs (shared pointer to provide document for different readers).
 */
XmlMultiResolutionReader::XmlMultiResolutionReader(
    std::shared_ptr<tinyxml2::XMLDocument> inputfile)
    : MultiResolutionReader(), xml_input_file_(std::move(inputfile)) {
  /** Empty besides initializer list and base class constructor call */
}

/**
 * @brief See base class definition.
 */
double XmlMultiResolutionReader::DoReadNodeSizeOnLevelZero() const {
  // Obtain correct node
  tinyxml2::XMLElement const *node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "domain", "nodeSize"});
  return XmlUtilities::ReadDouble(node);
}

/**
 * @brief See base class definition.
 */
int XmlMultiResolutionReader::DoReadNumberOfNodes(
    Direction const direction) const {
  // Get the correct component name
  std::string const component = direction == Direction::X   ? "x"
                                : direction == Direction::Y ? "y"
                                                            : "z";
  // obtain correct node
  tinyxml2::XMLElement const *node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "domain", "nodeRatio", component});
  return XmlUtilities::ReadInt(node);
}

/**
 * @brief See base class definition.
 */
int XmlMultiResolutionReader::DoReadMaximumLevel() const {
  // Obtain correct node
  tinyxml2::XMLElement const *level_node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "multiResolution", "maximumLevel"});
  return XmlUtilities::ReadInt(level_node);
}

/**
 * @brief See base class definition.
 */
int XmlMultiResolutionReader::DoReadEpsilonLevelReference() const {
  // Obtain correct nodes
  tinyxml2::XMLElement const *level_node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "multiResolution",
                         "refinementCriterion", "levelOfEpsilonReference"});
  return XmlUtilities::ReadInt(level_node);
}

/**
 * @brief See base class definition.
 */
double XmlMultiResolutionReader::DoReadEpsilonReference() const {
  // Obtain correct node
  tinyxml2::XMLElement const *epsilon_node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "multiResolution",
                         "refinementCriterion", "epsilonReference"});
  return XmlUtilities::ReadDouble(epsilon_node);
}
