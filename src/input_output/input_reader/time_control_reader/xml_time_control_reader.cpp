//===------------------- xml_time_control_reader.cpp ----------------------===//
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
#include "input_output/input_reader/time_control_reader/xml_time_control_reader.h"

#include "input_output/utilities/xml_utilities.h"

/**
 * @brief Default constructor for the time-control reader for xml-type input
 * files.
 * @param inputfile The xml input file document holding all information of the
 * user inputs (shared pointer to provide document for different readers).
 */
XmlTimeControlReader::XmlTimeControlReader(
    std::shared_ptr<tinyxml2::XMLDocument> inputfile)
    : TimeControlReader(), xml_input_file_(std::move(inputfile)) {
  /** Empty besides initializer list and base class constructor call */
}

/**
 * @brief See base class definition.
 */
double XmlTimeControlReader::DoReadStartTime() const {
  tinyxml2::XMLElement const *node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "timeControl", "startTime"});
  return XmlUtilities::ReadDouble(node);
}

/**
 * @brief See base class definition.
 */
double XmlTimeControlReader::DoReadEndTime() const {
  tinyxml2::XMLElement const *node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "timeControl", "endTime"});
  return XmlUtilities::ReadDouble(node);
}

/**
 * @brief See base class definition.
 */
double XmlTimeControlReader::DoReadCFLNumber() const {
  tinyxml2::XMLElement const *node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "timeControl", "CFLNumber"});
  return XmlUtilities::ReadDouble(node);
}
