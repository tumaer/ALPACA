//===--------------------- xml_restart_reader.cpp -------------------------===//
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
#include "input_output/input_reader/restart_reader/xml_restart_reader.h"

#include "input_output/utilities/xml_utilities.h"
#include "utilities/string_operations.h"

/**
 * @brief Default constructor for the restart reader for xml-type input files.
 * @param inputfile The xml input file document holding all information of the
 * user inputs (shared pointer to provide document for different readers).
 */
XmlRestartReader::XmlRestartReader(
    std::shared_ptr<tinyxml2::XMLDocument> inputfile)
    : RestartReader(), xml_input_file_(std::move(inputfile)) {
  /** Empty besides initializer list and base class constructor call */
}

/**
 * @brief See base class definition.
 */
std::string XmlRestartReader::DoReadRestoreMode() const {
  // Obtain correct node
  tinyxml2::XMLElement const *node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "restart", "restore", "mode"});
  return XmlUtilities::ReadString(node);
}

/**
 * @brief See base class definition.
 */
std::string XmlRestartReader::DoReadRestoreFilename() const {
  // Obtain correct node
  tinyxml2::XMLElement const *node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "restart", "restore", "fileName"});
  return XmlUtilities::ReadString(node);
}

/**
 * @brief See base class definition.
 */
std::string XmlRestartReader::DoReadSnapshotTimesType() const {
  // Obtain correct node
  tinyxml2::XMLElement const *node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "restart", "snapshots", "type"});
  return XmlUtilities::ReadString(node);
}

/**
 * @brief See base class definition.
 */
int XmlRestartReader::DoReadSnapshotInterval() const {
  // Obtain correct node
  tinyxml2::XMLElement const *node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "restart", "snapshots", "interval"});
  return XmlUtilities::ReadInt(node);
}

/**
 * @brief See base class definition.
 */
int XmlRestartReader::DoReadSnapshotIntervalsToKeep() const {
  // Obtain correct node
  tinyxml2::XMLElement const *node = XmlUtilities::GetChild(
      *xml_input_file_,
      {"configuration", "restart", "snapshots", "intervalsToKeep"});
  return XmlUtilities::ReadInt(node);
}

/**
 * @brief See base class definition.
 */
std::vector<double> XmlRestartReader::DoReadSnapshotTimeStamps() const {
  // Obtain correct node
  tinyxml2::XMLElement const *node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "restart", "snapshots", "stamps"});
  return XmlUtilities::ReadTimeStamps(node);
}
