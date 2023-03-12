//===------------------- xml_source_term_reader.cpp -----------------------===//
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
#include "input_output/input_reader/source_term_reader/xml_source_term_reader.h"

#include "input_output/utilities/xml_utilities.h"
#include <string>

/**
 * @brief Default constructor for the source term reader for xml-type input
 * files.
 * @param inputfile The xml input file document holding all information of the
 * user inputs (shared pointer to provide document for different readers).
 */
XmlSourceTermReader::XmlSourceTermReader(
    std::shared_ptr<tinyxml2::XMLDocument> inputfile)
    : SourceTermReader(), xml_input_file_(std::move(inputfile)) {
  /** Empty besides initializer list and base class constructor call */
}

/**
 * @brief See base class definition.
 */
double XmlSourceTermReader::DoReadGravity(Direction const direction) const {
  // Get the correct component name
  std::string const component = direction == Direction::X   ? "x"
                                : direction == Direction::Y ? "y"
                                                            : "z";
  // obtain correct node
  tinyxml2::XMLElement const *node = XmlUtilities::GetChild(
      *xml_input_file_, {"configuration", "sourceTerms", "gravity", component});
  return XmlUtilities::ReadDouble(node);
}
