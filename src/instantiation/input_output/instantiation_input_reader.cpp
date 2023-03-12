//===----------------- instantiation_input_reader.cpp ---------------------===//
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
#include "instantiation/input_output/instantiation_input_reader.h"

#include <memory>
#include <stdexcept>

#include "input_output/input_reader/input_definitions.h"
#include "input_output/log_writer/log_writer.h"
#include "input_output/utilities/file_utilities.h"
#include <tinyxml2.h>

#include "input_output/input_reader/boundary_condition_reader/xml_boundary_condition_reader.h"
#include "input_output/input_reader/dimensionalization_reader/xml_dimensionalization_reader.h"
#include "input_output/input_reader/initial_condition_reader/xml_initial_condition_reader.h"
#include "input_output/input_reader/material_reader/xml_material_reader.h"
#include "input_output/input_reader/multi_resolution_reader/xml_multi_resolution_reader.h"
#include "input_output/input_reader/output_reader/xml_output_reader.h"
#include "input_output/input_reader/restart_reader/xml_restart_reader.h"
#include "input_output/input_reader/source_term_reader/xml_source_term_reader.h"
#include "input_output/input_reader/time_control_reader/xml_time_control_reader.h"

namespace Instantiation {

/**
 * @brief Instantiates the full input reader class with the given input file.
 * @param input_filename Name of the file use for input.
 * @return The fully instantiated InputReader class.
 */
InputReader InstantiateInputReader(std::string const &input_filename) {
  // Check whether the file exists
  if (!FileUtilities::CheckIfPathExists(input_filename)) {
    throw std::logic_error("Input file " + input_filename + " does not exist!");
  }

  // Determine the input type
  InputType const input_type(
      StringToInputType(FileUtilities::GetFileExtension(input_filename)));
  // Instantiate correct reader
  switch (input_type) {
  case InputType::Xml: {
    // Open the file (here std::make_shared not possible)
    // shared pinter required to distribute the open input file on different
    // reader
    std::shared_ptr<tinyxml2::XMLDocument> input_file(
        new tinyxml2::XMLDocument);
    tinyxml2::XMLError error = input_file->LoadFile(input_filename.c_str());
    // Check if eversthing worked properly
    if (error != tinyxml2::XML_SUCCESS) {
      throw std::logic_error("Syntax error parsing the XML inputfile file, "
                             "check opening and closing tags!");
    }
    // Create the input reader properly
    LogWriter &logger = LogWriter::Instance();
    logger.LogMessage("Instantiating xml input reader");
    return InputReader(
        input_filename, input_type,
        std::make_unique<XmlMaterialReader const>(input_file),
        std::make_unique<XmlBoundaryConditionReader const>(input_file),
        std::make_unique<XmlInitialConditionReader const>(input_file),
        std::make_unique<XmlMultiResolutionReader const>(input_file),
        std::make_unique<XmlDimensionalizationReader const>(input_file),
        std::make_unique<XmlOutputReader const>(input_file),
        std::make_unique<XmlRestartReader const>(input_file),
        std::make_unique<XmlSourceTermReader const>(input_file),
        std::make_unique<XmlTimeControlReader const>(input_file));
  }

  default: {
    throw std::invalid_argument(
        "Input file extension not known! Cannot choose correct reader!");
  }
  }
}
} // namespace Instantiation
