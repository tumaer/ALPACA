//===---------------------- xml_material_reader.h -------------------------===//
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
#ifndef XML_MATERIAL_READER_H
#define XML_MATERIAL_READER_H

#include <memory>

#include "input_output/input_reader/material_reader/material_reader.h"
#include <tinyxml2.h>

/**
 * @brief Class that implements the actual reading procedure of material data
 * from input files of xml-type. Here, no consistency checks of the read
 * parameter are done. Only the validity of the correct variable type (double,
 * int, string) is done.
 */
class XmlMaterialReader : public MaterialReader {

  // The already openend xml input file (must be shared pointer to distribute
  // input file on different readers)
  std::shared_ptr<tinyxml2::XMLDocument const> const xml_input_file_;

  // Functions that are required from base class
  int DoReadNumberOfMaterials() const override;
  std::string
  DoReadMaterialType(unsigned int const material_index) const override;
  std::string
  DoReadEquationOfStateName(unsigned int const material_index) const override;
  std::unordered_map<std::string, double>
  DoReadEquationOfStateData(unsigned int const material_index) const override;
  double DoReadFixedValue(std::vector<unsigned int> const &material_indices,
                          MaterialProperty const property) const override;
  std::string DoReadModelName(std::vector<unsigned int> const &material_indices,
                              MaterialProperty const property) const override;
  std::unordered_map<std::string, double>
  DoReadModelData(std::vector<unsigned int> const &material_indices,
                  MaterialProperty const property) const override;

public:
  XmlMaterialReader() = delete;
  explicit XmlMaterialReader(std::shared_ptr<tinyxml2::XMLDocument> inputfile);
  ~XmlMaterialReader() = default;
  XmlMaterialReader(XmlMaterialReader const &) = delete;
  XmlMaterialReader &operator=(XmlMaterialReader const &) = delete;
  XmlMaterialReader(XmlMaterialReader &&) = delete;
  XmlMaterialReader &operator=(XmlMaterialReader &&) = delete;
};

#endif // XML_MATERIAL_READER_H
