//===-------------------- xml_material_reader.cpp -------------------------===//
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
#include "input_output/input_reader/material_reader/xml_material_reader.h"

#include "input_output/utilities/xml_utilities.h"

/**
 * @brief Default constructor for the material reader for xml-type input files.
 * @param inputfile The xml input file document holding all information of the
 * user inputs (shared pointer to provide document for different readers).
 */
XmlMaterialReader::XmlMaterialReader(
    std::shared_ptr<tinyxml2::XMLDocument> inputfile)
    : MaterialReader(), xml_input_file_(std::move(inputfile)) {
  /** Empty besides initializer list and base class constructor call */
}

/**
 * @brief Reads the number of materials used for the simulation
 * @return number of materials
 */
int XmlMaterialReader::DoReadNumberOfMaterials() const {
  // Get correct node
  tinyxml2::XMLElement const *node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "materials", "numberOfMaterials"});
  return XmlUtilities::ReadInt(node);
}

/**
 * @brief See base class definition.
 */
std::string
XmlMaterialReader::DoReadMaterialType(unsigned int const material_index) const {
  // Define the material tag
  std::string const material_tag("material" + std::to_string(material_index));
  if (XmlUtilities::ChildExists(*xml_input_file_, {"configuration", "materials",
                                                   material_tag, "type"})) {
    // Get correct node (directly to the final node to get error back
    // propagation in case)
    tinyxml2::XMLElement const *type_node = XmlUtilities::GetChild(
        *xml_input_file_, {"configuration", "materials", material_tag, "type"});

    return XmlUtilities::ReadString(type_node);
  } else {
    return "";
  }
}

/**
 * @brief See base class definition.
 */
std::string XmlMaterialReader::DoReadEquationOfStateName(
    unsigned int const material_index) const {
  // Define the material tag
  std::string const material_tag("material" + std::to_string(material_index));
  // Get correct node (directly to the final node to get error back propagation
  // in case)
  tinyxml2::XMLElement const *type_node = XmlUtilities::GetChild(
      *xml_input_file_,
      {"configuration", "materials", material_tag, "equationOfState", "type"});

  return XmlUtilities::ReadString(type_node);
}

/**
 * @brief See base class definition.
 */
std::unordered_map<std::string, double>
XmlMaterialReader::DoReadEquationOfStateData(
    unsigned int const material_index) const {
  // Define the material tag and move to the eos node
  std::string const material_tag("material" + std::to_string(material_index));
  tinyxml2::XMLElement const *eos_node = XmlUtilities::GetChild(
      *xml_input_file_,
      {"configuration", "materials", material_tag, "equationOfState"});
  // Declare the map
  std::unordered_map<std::string, double> eos_properties;
  // Fill the map with all values that are present
  eos_node = eos_node->FirstChildElement();
  while (eos_node != nullptr) {
    std::string const name = eos_node->Name();
    // except the type of the eos ( read separately )
    if (name != "type") {
      eos_properties[name] = XmlUtilities::ReadDouble(eos_node);
    }
    eos_node = eos_node->NextSiblingElement();
  }

  return eos_properties;
}

/**
 * @brief See base class definition.
 */
double XmlMaterialReader::DoReadFixedValue(
    std::vector<unsigned int> const &material_indices,
    MaterialProperty const property) const {
  // Obtain the correct tags dependent on the input data
  std::string const material_tag(
      MaterialReader::MaterialInputTag(material_indices));
  std::string const property_tag(MaterialPropertyToString(property, false));

  // In the following two options are tried
  try {
    // First try to read the fixed value directly from the property node
    // Obtain the correct node for the given material tag and property. Multiple
    // material indices imply material pairings.
    tinyxml2::XMLElement const *property_node =
        material_indices.size() > 1
            ? XmlUtilities::GetChild(*xml_input_file_,
                                     {"configuration", "materialPairings",
                                      material_tag, property_tag})
            : XmlUtilities::GetChild(
                  *xml_input_file_, {"configuration", "materials", material_tag,
                                     "properties", property_tag});

    return XmlUtilities::ReadDouble(property_node);

  } catch (std::invalid_argument const &) {
    // Second try to read the fixed value from the additional tag <fixedValue>
    // in the property node Obtain the correct node for the given material tag
    // and property. Multiple material indices imply material pairings.
    tinyxml2::XMLElement const *property_node =
        material_indices.size() > 1
            ? XmlUtilities::GetChild(*xml_input_file_,
                                     {"configuration", "materialPairings",
                                      material_tag, property_tag, "fixedValue"})
            : XmlUtilities::GetChild(
                  *xml_input_file_, {"configuration", "materials", material_tag,
                                     "properties", property_tag, "fixedValue"});

    return XmlUtilities::ReadDouble(property_node);
  }
}

/**
 * @brief See base class definition.
 */
std::string XmlMaterialReader::DoReadModelName(
    std::vector<unsigned int> const &material_indices,
    MaterialProperty const property) const {
  // Obtain the correct tags dependent on the inout data
  std::string const material_tag(
      MaterialReader::MaterialInputTag(material_indices));
  std::string const property_tag(MaterialPropertyToString(property, false));
  // Obtain the correct node for the given material tag and property. Multiple
  // aterial indices imply interface property. Directly move to the model name
  // node to get error back propagation in case
  tinyxml2::XMLElement const *model_name_node =
      material_indices.size() > 1
          ? XmlUtilities::GetChild(
                *xml_input_file_, {"configuration", "materialPairings",
                                   material_tag, property_tag, "model", "name"})
          : XmlUtilities::GetChild(*xml_input_file_,
                                   {"configuration", "materials", material_tag,
                                    "properties", property_tag, "model",
                                    "name"});

  // return the correct string
  return XmlUtilities::ReadString(model_name_node);
}

/**
 * @brief See base class definition.
 */
std::unordered_map<std::string, double> XmlMaterialReader::DoReadModelData(
    std::vector<unsigned int> const &material_indices,
    MaterialProperty const property) const {
  // Obtain the correct tags dependent on the inout data
  std::string const material_tag(
      MaterialReader::MaterialInputTag(material_indices));
  std::string const property_tag(MaterialPropertyToString(property, false));
  // Obtain the correct node for the given material tag and property. Multiple
  // aterial indices imply interface property. Directly move to the model
  // parameter node to get error back propagation in case
  tinyxml2::XMLElement const *parameter_node =
      material_indices.size() > 1
          ? XmlUtilities::GetChild(*xml_input_file_,
                                   {"configuration", "materialPairings",
                                    material_tag, property_tag, "model",
                                    "parameter"})
          : XmlUtilities::GetChild(*xml_input_file_,
                                   {"configuration", "materials", material_tag,
                                    "properties", property_tag, "model",
                                    "parameter"});

  // Declare map to be filled and get the first element of the map
  std::unordered_map<std::string, double> model_properties;
  parameter_node = parameter_node->FirstChildElement();
  // Loop through the whole property block and add everything except model
  // related data
  while (parameter_node != nullptr) {
    std::string const name(parameter_node->Name());
    // read parameter
    model_properties[name] = XmlUtilities::ReadDouble(parameter_node);
    parameter_node = parameter_node->NextSiblingElement();
  }
  // Return the filled map
  return model_properties;
}
