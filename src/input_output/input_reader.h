/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
*  \\                                                                                    *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA is a MPI-parallelized C++ code framework to simulate compressible multiphase    *
* flow physics. It allows for advanced high-resolution sharp-interface modeling          *
* empowered with efficient multiresolution compression. The modular code structure       *
* offers a broad flexibility to select among many most-recent numerical methods covering *
* WENO/T-ENO, Riemann solvers (complete/incomplete), strong-stability preserving Runge-  *
* Kutta time integration schemes, level set methods and many more.                       *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* LICENSE                                                                                *
*                                                                                        *
* ALPACA - Adaptive Level-set PArallel Code Alpaca                                       *
* Copyright (C) 2020 Nikolaus A. Adams and contributors (see AUTHORS list)               *
*                                                                                        *
* This program is free software: you can redistribute it and/or modify it under          *
* the terms of the GNU General Public License as published by the Free Software          *
* Foundation version 3.                                                                  *
*                                                                                        *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY        *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A        *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.               *
*                                                                                        *
* You should have received a copy of the GNU General Public License along with           *
* this program (gpl-3.0.txt).  If not, see <https://www.gnu.org/licenses/gpl-3.0.html>   *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* THIRD-PARTY tools                                                                      *
*                                                                                        *
* Please note, several third-party tools are used by ALPACA. These tools are not shipped *
* with ALPACA but available as git submodule (directing to their own repositories).      *
* All used third-party tools are released under open-source licences, see their own      *
* license agreement in 3rdParty/ for further details.                                    *
*                                                                                        *
* 1. tiny_xml           : See LICENSE_TINY_XML.txt for more information.                 *
* 2. expression_toolkit : See LICENSE_EXPRESSION_TOOLKIT.txt for more information.       *
* 3. FakeIt             : See LICENSE_FAKEIT.txt for more information                    *
* 4. Catch2             : See LICENSE_CATCH2.txt for more information                    *
* 5. ApprovalTests.cpp  : See LICENSE_APPROVAL_TESTS.txt for more information            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, February 10th, 2021                                                            *
*                                                                                        *
*****************************************************************************************/
#ifndef INPUT_READER_H
#define INPUT_READER_H

#include <filesystem>
#include <memory>

#include "input_output/input_reader/input_definitions.h"
#include "input_output/input_reader/boundary_condition_reader/boundary_condition_reader.h"
#include "input_output/input_reader/initial_condition_reader/initial_condition_reader.h"
#include "input_output/input_reader/material_reader/material_reader.h"
#include "input_output/input_reader/multi_resolution_reader/multi_resolution_reader.h"
#include "input_output/input_reader/dimensionalization_reader/dimensionalization_reader.h"
#include "input_output/input_reader/output_reader/output_reader.h"
#include "input_output/input_reader/restart_reader/restart_reader.h"
#include "input_output/input_reader/source_term_reader/source_term_reader.h"
#include "input_output/input_reader/time_control_reader/time_control_reader.h"

/**
 * @brief The input reader class serves as a holder for all base class readers that can be used to read data from the input file.
 *        THe input reader does not read any data but serves as a distributor of all different readers.
 */
class InputReader {

private:
   // Member to store the input type and filename
   std::string const input_filename_;
   InputType const input_type_;

   // All readere that are required from base class
   std::unique_ptr<MaterialReader const> const material_reader_;
   std::unique_ptr<BoundaryConditionReader const> const boundary_condition_reader_;
   std::unique_ptr<InitialConditionReader const> const initial_condition_reader_;
   std::unique_ptr<MultiResolutionReader const> const multi_resolution_reader_;
   std::unique_ptr<DimensionalizationReader const> const dimensionalization_reader_;
   std::unique_ptr<OutputReader const> const output_reader_;
   std::unique_ptr<RestartReader const> const restart_reader_;
   std::unique_ptr<SourceTermReader const> const source_term_reader_;
   std::unique_ptr<TimeControlReader const> const time_control_reader_;

public:
   explicit InputReader( std::string const& input_filename,
                         InputType const input_type,
                         std::unique_ptr<MaterialReader const> material_reader,
                         std::unique_ptr<BoundaryConditionReader const> boundary_condition_reader,
                         std::unique_ptr<InitialConditionReader const> initial_condition_reader,
                         std::unique_ptr<MultiResolutionReader const> multi_resolution_reader,
                         std::unique_ptr<DimensionalizationReader const> dimensionalization_reader,
                         std::unique_ptr<OutputReader const> output_reader,
                         std::unique_ptr<RestartReader const> restart_reader,
                         std::unique_ptr<SourceTermReader const> source_term_reader,
                         std::unique_ptr<TimeControlReader const> time_control_reader );
   InputReader()                     = delete;
   virtual ~InputReader()            = default;
   InputReader( InputReader const& ) = delete;
   InputReader& operator=( InputReader const& ) = delete;
   InputReader( InputReader&& )                 = delete;
   InputReader& operator=( InputReader&& ) = delete;

   // Return functions of the input reader
   TEST_VIRTUAL std::filesystem::path GetInputFile() const;
   TEST_VIRTUAL InputType GetInputType() const;
   TEST_VIRTUAL MaterialReader const& GetMaterialReader() const;
   TEST_VIRTUAL BoundaryConditionReader const& GetBoundaryConditionReader() const;
   TEST_VIRTUAL InitialConditionReader const& GetInitialConditionReader() const;
   TEST_VIRTUAL MultiResolutionReader const& GetMultiResolutionReader() const;
   DimensionalizationReader const& GetDimensionalizationReader() const;
   TEST_VIRTUAL OutputReader const& GetOutputReader() const;
   TEST_VIRTUAL RestartReader const& GetRestartReader() const;
   TEST_VIRTUAL SourceTermReader const& GetSourceTermReader() const;
   TEST_VIRTUAL TimeControlReader const& GetTimeControlReader() const;
};
#endif// INPUT_READER_H
