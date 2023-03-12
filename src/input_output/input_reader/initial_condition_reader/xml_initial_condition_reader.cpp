//===--------------- xml_initial_condition_reader.cpp ---------------------===//
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
#include "input_output/input_reader/initial_condition_reader/xml_initial_condition_reader.h"
#include "input_output/utilities/xml_utilities.h"

#include "user_specifications/compile_time_constants.h"

/**
 * @brief Default constructor for the initial condition reader for xml-type
 * input files.
 * @param inputfile The xml input file document holding all information of the
 * user inputs (shared pointer to provide document for different readers).
 */
XmlInitialConditionReader::XmlInitialConditionReader(
    std::shared_ptr<tinyxml2::XMLDocument> inputfile)
    : InitialConditionReader(), xml_input_file_(std::move(inputfile)) {
  /** Empty besides initializer list and base class constructor call */
}

/**
 * @brief See base class definition.
 */
std::string XmlInitialConditionReader::DoReadMaterialInitialConditions(
    unsigned int const material_index) const {
  // Define the material name
  std::string const material_name("material" + std::to_string(material_index));
  // Get the correct node
  tinyxml2::XMLElement const *material_node = XmlUtilities::GetChild(
      *xml_input_file_,
      {"configuration", "domain", "initialConditions", material_name});

  return XmlUtilities::ReadString(material_node);
}

/**
 * @brief See base class definition.
 */
std::string XmlInitialConditionReader::DoReadLevelsetInitializerType(
    unsigned int const levelset_index) const {
  // Get the correct levelset node (which always must exist)
  std::string const levelset_name("levelSet" + std::to_string(levelset_index));
  tinyxml2::XMLElement const *levelset_node = XmlUtilities::GetChild(
      *xml_input_file_,
      {"configuration", "domain", "initialConditions", levelset_name});
  if (XmlUtilities::ChildExists(levelset_node, "type")) {
    return XmlUtilities::ReadString(
        XmlUtilities::GetChild(levelset_node, {"type"}));
  } else {
    return "";
  }
}

/**
 * @brief See base class definition.
 */
std::string XmlInitialConditionReader::DoReadLevelsetInitializerInput(
    unsigned int const levelset_index) const {
  // Get the correct levelset node (which always must exist)
  std::string const levelset_name("levelSet" + std::to_string(levelset_index));
  tinyxml2::XMLElement const *levelset_node = XmlUtilities::GetChild(
      *xml_input_file_,
      {"configuration", "domain", "initialConditions", levelset_name});
  // If the input tag exists. Read the data from it otherwise read from the
  // levelset node directly.
  if (XmlUtilities::ChildExists(levelset_node, "input")) {
    return XmlUtilities::ReadString(
        XmlUtilities::GetChild(levelset_node, {"input"}));
  } else {
    return XmlUtilities::ReadString(levelset_node);
  }
}

/**
 * @brief See base class definition.
 */
std::vector<std::tuple<std::string, double, double, std::uint64_t>>
XmlInitialConditionReader::DoReadParametricLevelsetInitializerVariables(
    unsigned int const levelset_index) const {
  /// Return array
  std::vector<std::tuple<std::string, double, double, std::uint64_t>>
      variables{};

  // Get the correct levelset node (which always must exist)
  std::string const levelset_name("levelSet" + std::to_string(levelset_index));
  tinyxml2::XMLElement const *levelset_node = XmlUtilities::GetChild(
      *xml_input_file_,
      {"configuration", "domain", "initialConditions", levelset_name});
  // Get all child nodes having the name parameter (if no are present, the
  // vector is empty)
  std::vector<tinyxml2::XMLElement const *> variable_nodes =
      XmlUtilities::GetChilds(levelset_node, "variable");
  // Loop through all elements to read the variables into a map
  for (auto variable_node : variable_nodes) {
    // Read all values
    std::string const name = XmlUtilities::ReadString(
        XmlUtilities::GetChild(variable_node, {"name"}));
    double const start = XmlUtilities::ReadDouble(
        XmlUtilities::GetChild(variable_node, {"start"}));
    double const end = XmlUtilities::ReadDouble(
        XmlUtilities::GetChild(variable_node, {"end"}));
    std::uint64_t const points = XmlUtilities::ReadUnsignedLongInt(
        XmlUtilities::GetChild(variable_node, {"points"}));
    // Add the variables to the tuple vector
    variables.push_back(std::make_tuple(name, start, end, points));
  }
  return variables;
}

/**
 * @brief See base class definition.
 */
std::array<double, 3>
XmlInitialConditionReader::DoReadParametricLevelsetInitializerReferencePoint(
    unsigned int const levelset_index) const {
  // Return array
  std::array<double, 3> reference_point = {0.0, 0.0, 0.0};

  // Define the name of the levelset tag to be read and read the levelset tag
  // (raises error if not existent)
  std ::string const levelset_name("levelSet" + std::to_string(levelset_index));
  tinyxml2::XMLElement const *ref_node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "domain", "initialConditions",
                         levelset_name, "referencePoint"});
  // Read the different coordinates depending on the dimension
  std::vector<std::string> const coordinates = {"x", "y", "z"};
  for (unsigned int dim = 0; dim < DTI(CC::DIM()); dim++) {
    reference_point[dim] = XmlUtilities::ReadDouble(
        XmlUtilities::GetChild(ref_node, {coordinates[dim]}));
  }
  return reference_point;
}

/**
 * @brief See base class definition.
 */
std::vector<std::array<double, 6>>
XmlInitialConditionReader::DoReadLevelsetInitializerBoundingBoxes(
    unsigned int const levelset_index) const {
  /// Return array and limit
  std::vector<std::array<double, 6>> bounding_boxes{};
  double const pos_limit = std::numeric_limits<double>::max();
  double const neg_limit = std::numeric_limits<double>::lowest();

  // Define the name of the levelset tag to be read and read the levelset tag
  // (raises error if not existent)
  std ::string const levelset_name("levelSet" + std::to_string(levelset_index));
  tinyxml2::XMLElement const *levelset_node = XmlUtilities::GetChild(
      *xml_input_file_,
      {"configuration", "domain", "initialConditions", levelset_name});
  // Get all child nodes having the name bounding box (if no are present, the
  // vector is empty)
  std::vector<tinyxml2::XMLElement const *> bounding_box_nodes =
      XmlUtilities::GetChilds(levelset_node, "boundingBox");

  // Loop through all elements to read the bounding box completely
  for (auto bounding_box_node : bounding_box_nodes) {
    // Default bounding box
    std::array<double, 6> bounding_box = {neg_limit, pos_limit, neg_limit,
                                          pos_limit, neg_limit, pos_limit};
    // Coordinates for each dimension
    std::vector<std::string> const coordinates = {"x", "y", "z"};
    for (unsigned int dim = 0; dim < DTI(CC::DIM()); dim++) {
      if (XmlUtilities::ChildExists(bounding_box_node,
                                    coordinates[dim] + "Min")) {
        bounding_box[2 * dim] = XmlUtilities::ReadDouble(XmlUtilities::GetChild(
            bounding_box_node, {coordinates[dim] + "Min"}));
      }
      if (XmlUtilities::ChildExists(bounding_box_node,
                                    coordinates[dim] + "Max")) {
        bounding_box[2 * dim + 1] =
            XmlUtilities::ReadDouble(XmlUtilities::GetChild(
                bounding_box_node, {coordinates[dim] + "Max"}));
      }
    }
    bounding_boxes.push_back(bounding_box);
  }
  return bounding_boxes;
}
