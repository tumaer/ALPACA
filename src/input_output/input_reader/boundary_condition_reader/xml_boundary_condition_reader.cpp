//===---------------- xml_boundary_condition_reader.cpp -------------------===//
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
#include "input_output/input_reader/boundary_condition_reader/xml_boundary_condition_reader.h"

#include "input_output/utilities/xml_utilities.h"

/**
 * @brief Default constructor for the boundary condition reader for xml-type
 * input files.
 * @param inputfile The xml input file document holding all information of the
 * user inputs (shared pointer to provide document for different readers).
 */
XmlBoundaryConditionReader::XmlBoundaryConditionReader(
    std::shared_ptr<tinyxml2::XMLDocument> inputfile)
    : BoundaryConditionReader(), xml_input_file_(std::move(inputfile)) {
  /** Empty besides initializer list and base class constructor call */
}

/**
 * @brief See base class definition.
 */
std::string XmlBoundaryConditionReader::DoReadMaterialBoundaryType(
    BoundaryLocation const location) const {
  // get the correct name
  std::string const location_name(BoundaryLocationToString(location, false));
  // get the correct node and read name
  tinyxml2::XMLElement const *const node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "domain", "boundaryConditions",
                         "material", location_name});
  return XmlUtilities::ReadString(node);
}

/**
 * @brief See base class definition.
 */
std::string XmlBoundaryConditionReader::DoReadLevelSetBoundaryType(
    BoundaryLocation const location) const {
  // get the correct name
  std::string const location_name(BoundaryLocationToString(location, false));
  // get the correct node and read name
  tinyxml2::XMLElement const *const node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "domain", "boundaryConditions",
                         "levelSet", location_name});
  return XmlUtilities::ReadString(node);
}

/**
 * @brief See base class definition.
 */
double XmlBoundaryConditionReader::DoReadMaterialFixedValueBoundaryCondition(
    BoundaryLocation const location, std::string const &variable) const {

  // Define the string to be used for the given location
  std::string const fixed_value_name("values" +
                                     BoundaryLocationToString(location, true));

  // Get the correct nodes for the material
  tinyxml2::XMLElement const *fixed_value_node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "domain", "boundaryConditions",
                         "material", fixed_value_name, variable});
  return XmlUtilities::ReadDouble(fixed_value_node);
}
