//===---------------- xml_boundary_condition_reader.h ---------------------===//
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
#ifndef XML_BOUNDARY_CONDITION_READER_H
#define XML_BOUNDARY_CONDITION_READER_H

#include <memory>

#include "input_output/input_reader/boundary_condition_reader/boundary_condition_reader.h"
#include <tinyxml2.h>

/**
 * @brief Class that implements the actual reading procedure of boundary
 * condition data from input files of xml-type. Here, no consistency checks of
 * the read parameter are done. Only the validity of the correct variable type
 *        (double, int, string) is done.
 */
class XmlBoundaryConditionReader : public BoundaryConditionReader {

  // The already openend xml input file (must be shared pointer to distribute
  // input file on different readers)
  std::shared_ptr<tinyxml2::XMLDocument const> const xml_input_file_;

  // Functions that are required from base class
  std::string
  DoReadMaterialBoundaryType(BoundaryLocation const location) const override;
  std::string
  DoReadLevelSetBoundaryType(BoundaryLocation const location) const override;
  double DoReadMaterialFixedValueBoundaryCondition(
      BoundaryLocation const location,
      std::string const &variable) const override;

public:
  XmlBoundaryConditionReader() = delete;
  explicit XmlBoundaryConditionReader(
      std::shared_ptr<tinyxml2::XMLDocument> inputfile);
  ~XmlBoundaryConditionReader() = default;
  XmlBoundaryConditionReader(XmlBoundaryConditionReader const &) = delete;
  XmlBoundaryConditionReader &
  operator=(XmlBoundaryConditionReader const &) = delete;
  XmlBoundaryConditionReader(XmlBoundaryConditionReader &&) = delete;
  XmlBoundaryConditionReader &operator=(XmlBoundaryConditionReader &&) = delete;
};

#endif // XML_BOUNDARY_CONDITION_READER_H
