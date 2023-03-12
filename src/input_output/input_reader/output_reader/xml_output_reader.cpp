//===----------------------- xml_output_reader.cpp ------------------------===//
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
#include "input_output/input_reader/output_reader/xml_output_reader.h"

#include "input_output/utilities/xml_utilities.h"

/**
 * @brief Default constructor for the output reader for xml-type input files.
 * @param inputfile The xml input file document holding all information of the
 * user inputs (shared pointer to provide document for different readers).
 */
XmlOutputReader::XmlOutputReader(
    std::shared_ptr<tinyxml2::XMLDocument> inputfile)
    : OutputReader(), xml_input_file_(std::move(inputfile)) {
  /** Empty besides initializer list and base class constructor call */
}

/**
 * @brief See base class definition.
 */
double XmlOutputReader::DoReadTimeNamingFactor() const {
  // Obtain correct node
  tinyxml2::XMLElement const *node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "output", "timeNamingFactor"});
  return XmlUtilities::ReadDouble(node);
}

/**
 * @brief See base class definition.
 */
std::string
XmlOutputReader::DoReadOutputTimesType(OutputType const output_type) const {
  // specify the correct tag for the given output type
  std::string const output_tag = output_type == OutputType::Interface
                                     ? "interfaceOutput"
                                     : "standardOutput";
  // Obtain correct nodes
  tinyxml2::XMLElement const *type_node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "output", output_tag, "type"});

  return XmlUtilities::ReadString(type_node);
}

/**
 * @brief See base class definition.
 */
double
XmlOutputReader::DoReadOutputInterval(OutputType const output_type) const {
  // specify the correct tag for the given output type
  std::string const output_tag = output_type == OutputType::Interface
                                     ? "interfaceOutput"
                                     : "standardOutput";
  // Obtain correct node
  tinyxml2::XMLElement const *interval_node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "output", output_tag, "interval"});

  return XmlUtilities::ReadDouble(interval_node);
}

/**
 * @brief See base class definition.
 */
std::vector<double>
XmlOutputReader::DoReadOutputTimeStamps(OutputType const output_type) const {
  // specify the correct tag for the given output type
  std::string const output_tag = output_type == OutputType::Interface
                                     ? "interfaceOutput"
                                     : "standardOutput";
  // Obtain correct node
  tinyxml2::XMLElement const *stamp_node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "output", output_tag, "stamps"});

  return XmlUtilities::ReadTimeStamps(stamp_node);
}
