//===-------------------- xml_time_control_reader.h -----------------------===//
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
#ifndef XML_TIME_CONTROL_READER_H
#define XML_TIME_CONTROL_READER_H

#include <memory>

#include "input_output/input_reader/time_control_reader/time_control_reader.h"
#include <tinyxml2.h>

/**
 * @brief Class that implements the actual reading procedure of time control
 * data from input files of xml-type. Here, no consistency checks of the read
 * parameter are done. Only the validity of the correct variable type (double,
 * int, string) is done.
 */
class XmlTimeControlReader : public TimeControlReader {

  // The already openend xml input file (must be shared pointer to distribute
  // input file on different readers)
  std::shared_ptr<tinyxml2::XMLDocument const> const xml_input_file_;

  // Functions that are required from base class
  double DoReadStartTime() const override;
  double DoReadEndTime() const override;
  double DoReadCFLNumber() const override;

public:
  XmlTimeControlReader() = delete;
  explicit XmlTimeControlReader(
      std::shared_ptr<tinyxml2::XMLDocument> inputfile);
  ~XmlTimeControlReader() = default;
  XmlTimeControlReader(XmlTimeControlReader const &) = default;
  XmlTimeControlReader &operator=(XmlTimeControlReader const &) = delete;
  XmlTimeControlReader(XmlTimeControlReader &&) = default;
  XmlTimeControlReader &operator=(XmlTimeControlReader &&) = delete;
};

#endif // XML_TIME_CONTROL_READER_H
