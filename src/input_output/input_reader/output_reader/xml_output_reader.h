//===----------------------- xml_output_reader.h --------------------------===//
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
#ifndef XML_OUTPUT_READER_H
#define XML_OUTPUT_READER_H

#include <memory>

#include "input_output/input_reader/output_reader/output_reader.h"
#include <tinyxml2.h>

/**
 * @brief Class that implements the actual reading procedure of output data from
 * input files of xml-type. Here, no consistency checks of the read parameter
 * are done. Only the validity of the correct variable type (double, int,
 * string) is done.
 */
class XmlOutputReader : public OutputReader {

  // The already openend xml input file (must be shared pointer to distribute
  // input file on different readers)
  std::shared_ptr<tinyxml2::XMLDocument const> const xml_input_file_;

  // Functions that are required from base class
  double DoReadTimeNamingFactor() const override;
  std::string
  DoReadOutputTimesType(OutputType const output_type) const override;
  double DoReadOutputInterval(OutputType const output_type) const override;
  std::vector<double>
  DoReadOutputTimeStamps(OutputType const output_type) const override;

public:
  XmlOutputReader() = delete;
  explicit XmlOutputReader(std::shared_ptr<tinyxml2::XMLDocument> inputfile);
  ~XmlOutputReader() = default;
  XmlOutputReader(XmlOutputReader const &) = delete;
  XmlOutputReader &operator=(XmlOutputReader const &) = delete;
  XmlOutputReader(XmlOutputReader &&) = delete;
  XmlOutputReader &operator=(XmlOutputReader &&) = delete;
};

#endif // XML_OUTPUT_READER_H
