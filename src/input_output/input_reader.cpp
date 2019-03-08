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
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, July 1st, 2020                                                                 *
*                                                                                        *
*****************************************************************************************/
#include "input_output/input_reader.h"

/**
 * @brief Standard constructor for the input reader class 
 * @param input_filename Name of the input file that is used 
 * @param input_type Identifier of the input type 
 * @param material_reader Base class to the material reader 
 * @param boundary_condition_reader Base class to the boundary condition reader 
 * @param initial_condition_reader Base class to the initial condition reader 
 * @param multi_resolution_reader Base class to the multi resolution reader 
 * @param dimensionalization_reader Base class to the dimensionalization reader 
 * @param output_reader Base class to the output reader 
 * @param restart_reader Base class to the restart reader 
 * @param source_term_reader Base class to the source term reader 
 * @param time_control_reader Base class to the time control reader 
 */
InputReader::InputReader( std::string const& input_filename, 
                          InputType const input_type,
                          std::unique_ptr<MaterialReader const> material_reader, 
                          std::unique_ptr<BoundaryConditionReader const> boundary_condition_reader,
                          std::unique_ptr<InitialConditionReader const> initial_condition_reader, 
                          std::unique_ptr<MultiResolutionReader const> multi_resolution_reader,
                          std::unique_ptr<DimensionalizationReader const> dimensionalization_reader, 
                          std::unique_ptr<OutputReader const> output_reader,
                          std::unique_ptr<RestartReader const> restart_reader, 
                          std::unique_ptr<SourceTermReader const> source_term_reader, 
                          std::unique_ptr<TimeControlReader const> time_control_reader ) :
   // Start initializer list 
   input_filename_( input_filename ),
   input_type_( input_type ),
   material_reader_( std::move( material_reader ) ),
   boundary_condition_reader_( std::move( boundary_condition_reader ) ),
   initial_condition_reader_( std::move( initial_condition_reader ) ),
   multi_resolution_reader_( std::move( multi_resolution_reader ) ),
   dimensionalization_reader_( std::move( dimensionalization_reader ) ),
   output_reader_( std::move( output_reader ) ),
   restart_reader_( std::move( restart_reader ) ),
   source_term_reader_( std::move( source_term_reader ) ),
   time_control_reader_( std::move( time_control_reader ) ) {
   /** Empty besides initializer list */
}

/**
 * @brief Gives the filename used as an input 
 * @return input filename 
 */
std::string InputReader::GetInputFile() const {
   return input_filename_;
}

/**
 * @brief Gives the input type  
 * @return Input type identifier 
 */
InputType InputReader::GetInputType() const {
   return input_type_;
}

/**
 * @brief Gives the instance of the material reader
 * @return reference to the material reader
 */
MaterialReader const& InputReader::GetMaterialReader() const {
   return *material_reader_;
}

/**
 * @brief Gives the instance of the boundary condition reader
 * @return reference to the boundary condition reader
 */
BoundaryConditionReader const& InputReader::GetBoundaryConditionReader() const {
   return *boundary_condition_reader_;
}

/**
 * @brief Gives the instance of the initial condition reader
 * @return reference to the initial condition reader
 */
InitialConditionReader const& InputReader::GetInitialConditionReader() const {
   return *initial_condition_reader_;
}

/**
 * @brief Gives the instance of the multiresolution reader
 * @return reference to the multiresolution reader
 */
MultiResolutionReader const& InputReader::GetMultiResolutionReader() const {
   return *multi_resolution_reader_;
}

/**
 * @brief Gives the instance of the non-dimensionalizaion reader
 * @return reference to the non-dimensionalizaion reader
 */
DimensionalizationReader const& InputReader::GetDimensionalizationReader() const {
   return *dimensionalization_reader_;
}

/**
 * @brief Gives the instance of the output reader
 * @return reference to the output reader
 */
OutputReader const& InputReader::GetOutputReader() const {
   return *output_reader_;
}

/**
 * @brief Gives the instance of the restart reader
 * @return reference to the restart reader
 */
RestartReader const& InputReader::GetRestartReader() const {
   return *restart_reader_;
}

/**
 * @brief Gives the instance of the source term reader
 * @return reference to the source term reader
 */
SourceTermReader const& InputReader::GetSourceTermReader() const {
   return *source_term_reader_;
}

/**
 * @brief Gives the instance of the time control reader
 * @return reference to the time control reader
 */
TimeControlReader const& InputReader::GetTimeControlReader() const {
   return *time_control_reader_;
}