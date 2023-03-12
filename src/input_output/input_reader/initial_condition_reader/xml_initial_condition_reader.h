//===----------------- xml_initial_condition_reader.h ---------------------===//
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
#ifndef XML_INITIAL_CONDITION_READER_H
#define XML_INITIAL_CONDITION_READER_H

#include <memory>

#include "input_output/input_reader/initial_condition_reader/initial_condition_reader.h"
#include <string>
#include <tinyxml2.h>

/**
 * @brief Class that implements the actual reading procedure of initial
 * condition data from input files of xml-type. Here, no consistency checks of
 * the read parameter are done. Only the validity of the correct variable type
 *        (double, int, string) is done.
 */
class XmlInitialConditionReader : public InitialConditionReader {

  // The already openend xml input file (must be shared pointer to distribute
  // input file on different readers)
  std::shared_ptr<tinyxml2::XMLDocument const> const xml_input_file_;

  // Functions that are required from base class
  std::string DoReadMaterialInitialConditions(
      unsigned int const material_index) const override;
  std::string DoReadLevelsetInitializerType(
      unsigned int const material_index) const override;
  std::string DoReadLevelsetInitializerInput(
      unsigned int const levelset_index) const override;
  std::vector<std::tuple<std::string, double, double, std::uint64_t>>
  DoReadParametricLevelsetInitializerVariables(
      unsigned int const levelset_index) const override;
  std::array<double, 3> DoReadParametricLevelsetInitializerReferencePoint(
      unsigned int const levelset_index) const override;
  std::vector<std::array<double, 6>> DoReadLevelsetInitializerBoundingBoxes(
      unsigned int const material_index) const override;

public:
  XmlInitialConditionReader() = delete;
  explicit XmlInitialConditionReader(
      std::shared_ptr<tinyxml2::XMLDocument> inputfile);
  ~XmlInitialConditionReader() = default;
  XmlInitialConditionReader(XmlInitialConditionReader const &) = delete;
  XmlInitialConditionReader &
  operator=(XmlInitialConditionReader const &) = delete;
  XmlInitialConditionReader(XmlInitialConditionReader &&) = delete;
  XmlInitialConditionReader &operator=(XmlInitialConditionReader &&) = delete;
};

#endif // XML_INITIAL_CONDITION_READER_H
