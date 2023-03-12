//===------------------------- input_reader.h -----------------------------===//
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
#ifndef INPUT_READER_H
#define INPUT_READER_H

#include <filesystem>
#include <memory>

#include "input_output/input_reader/boundary_condition_reader/boundary_condition_reader.h"
#include "input_output/input_reader/dimensionalization_reader/dimensionalization_reader.h"
#include "input_output/input_reader/initial_condition_reader/initial_condition_reader.h"
#include "input_output/input_reader/input_definitions.h"
#include "input_output/input_reader/material_reader/material_reader.h"
#include "input_output/input_reader/multi_resolution_reader/multi_resolution_reader.h"
#include "input_output/input_reader/output_reader/output_reader.h"
#include "input_output/input_reader/restart_reader/restart_reader.h"
#include "input_output/input_reader/source_term_reader/source_term_reader.h"
#include "input_output/input_reader/time_control_reader/time_control_reader.h"

/**
 * @brief The input reader class serves as a holder for all base class readers
 * that can be used to read data from the input file. THe input reader does not
 * read any data but serves as a distributor of all different readers.
 */
class InputReader {

private:
  // Member to store the input type and filename
  std::string const input_filename_;
  InputType const input_type_;

  // All readere that are required from base class
  std::unique_ptr<MaterialReader const> const material_reader_;
  std::unique_ptr<BoundaryConditionReader const> const
      boundary_condition_reader_;
  std::unique_ptr<InitialConditionReader const> const initial_condition_reader_;
  std::unique_ptr<MultiResolutionReader const> const multi_resolution_reader_;
  std::unique_ptr<DimensionalizationReader const> const
      dimensionalization_reader_;
  std::unique_ptr<OutputReader const> const output_reader_;
  std::unique_ptr<RestartReader const> const restart_reader_;
  std::unique_ptr<SourceTermReader const> const source_term_reader_;
  std::unique_ptr<TimeControlReader const> const time_control_reader_;

public:
  explicit InputReader(
      std::string const &input_filename, InputType const input_type,
      std::unique_ptr<MaterialReader const> material_reader,
      std::unique_ptr<BoundaryConditionReader const> boundary_condition_reader,
      std::unique_ptr<InitialConditionReader const> initial_condition_reader,
      std::unique_ptr<MultiResolutionReader const> multi_resolution_reader,
      std::unique_ptr<DimensionalizationReader const> dimensionalization_reader,
      std::unique_ptr<OutputReader const> output_reader,
      std::unique_ptr<RestartReader const> restart_reader,
      std::unique_ptr<SourceTermReader const> source_term_reader,
      std::unique_ptr<TimeControlReader const> time_control_reader);
  InputReader() = delete;
  virtual ~InputReader() = default;
  InputReader(InputReader const &) = delete;
  InputReader &operator=(InputReader const &) = delete;
  InputReader(InputReader &&) = delete;
  InputReader &operator=(InputReader &&) = delete;

  // Return functions of the input reader
  TEST_VIRTUAL std::filesystem::path GetInputFile() const;
  TEST_VIRTUAL InputType GetInputType() const;
  TEST_VIRTUAL MaterialReader const &GetMaterialReader() const;
  TEST_VIRTUAL BoundaryConditionReader const &
  GetBoundaryConditionReader() const;
  TEST_VIRTUAL InitialConditionReader const &GetInitialConditionReader() const;
  TEST_VIRTUAL MultiResolutionReader const &GetMultiResolutionReader() const;
  DimensionalizationReader const &GetDimensionalizationReader() const;
  TEST_VIRTUAL OutputReader const &GetOutputReader() const;
  TEST_VIRTUAL RestartReader const &GetRestartReader() const;
  TEST_VIRTUAL SourceTermReader const &GetSourceTermReader() const;
  TEST_VIRTUAL TimeControlReader const &GetTimeControlReader() const;
};
#endif // INPUT_READER_H
