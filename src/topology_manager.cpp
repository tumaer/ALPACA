/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*  \\\\                                                                                  *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA                                                                                 *
* Copyright (c) 2017 Nikolaus A. Adams and contributors (see AUTHORS list)               *
* All rights reserved.                                                                   *
*                                                                                        *
* Chair of Aerodynamics and Fluid Mechanics                                              *
* Technical University of Munich                                                         *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
* This project has received funding from the European Reseach Council (ERC)              *
* under the European Union's Horizon 2020 research and innovation programme              *
* (grant agreement No 667483).                                                           *
*                                                                                        *
* ERC Advanced Grant No 667483, Prof. Dr. Nikolaus A. Adams:                             *
* "NANOSHOCK - Manufacturing Shock Interactions for Innovative Nanoscale Processes"      *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Redistribution and use in source and binary forms, with or without                     *
* modification, are permitted provided that the following conditions are met:            *
*                                                                                        *
* 1. Redistributions of source code must retain the above copyright notice,              *
*    this list of conditions and the following disclaimer.                               *
*                                                                                        *
* 2. Redistributions in binary form must reproduce the above copyright notice            *
*    this list of conditions and the following disclaimer in the documentation           *
*    and/or other materials provided with the distribution.                              *
*                                                                                        *
* 3. Neither the name of the copyright holder nor the names of its                       *
*    contributors may be used to endorse or promote products derived from this           *
*    software without specific prior written permission.                                 *
*                                                                                        *
* 4. Any redistribution of substantial fractions of the code as a                        *
*    different project should preserve the word ALPACA in the name                       *
*    of the code                                                                         *
*                                                                                        *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"            *
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE              *
* IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE            *
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE              *
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR                    *
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF                   *
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS               *
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN                *
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)                *
* ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE            *
* POSSIBILITY OF SUCH DAMAGE.                                                            *
*                                                                                        *
* Please note, several third-party tools are used within the ALPACA code under           *
* their own license agreement.                                                           *
*                                                                                        *
* 1. xdmf_writer        : Licensed by Technische Universitaet Muenchen                   *
*                         See 'COPYING_XDMF_WRITER' for more information.                *
*                                                                                        *
* 2. tiny_xml           : This software is provided 'as-is', without any express or      *
*                         implied warranty. In no event will the authors be held         *
*                         liable for any damages arising from the use of this software.  *
*                         See COPYING_TINY_XMLfor more information.                      *
*                                                                                        *
* 3. expression_toolkit : Free use of The C++ Mathematical Expression Toolkit Library is *
*                         permitted under the guidelines and in accordance with the most *
*                         current version of the Common Public License.                  *
*                         http://www.opensource.org/licenses/cpl1.0.php                  *
*                         See COPYING_EXPRESSION_TOOLKITfor more information.            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* AUTHORS                                                                                *
*                                                                                        *
*   Prof. Dr. Nikolaus A. Adams                                                          *
*                                                                                        *
*   Dr. Stefan Adami                                                                     *
*   Vladimir Bogdanov                                                                    *
*   Nico Fleischmann                                                                     *
*   Nils Hoppe                                                                           *
*   Naeimeh Hosseini                                                                     *
*   Jakob Kaiser                                                                         *
*   Aleksandr Lunkov                                                                     *
*   Thomas Paula                                                                         *
*   Josef Winter                                                                         *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
*   nanoshock@aer.mw.tum.de                                                              *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, December 15th 2017                                                             *
*                                                                                        *
*****************************************************************************************/

#include "topology_manager.h"

#include "node.h"

#include <algorithm>
#include <utility>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <iostream>
#include <numeric>

//Specialized Template Functions have to go first. Otherwise compiler cannot "see" them.

/**
 * @brief Datatypes for 3D Simulations. See Meta function.
 */
template<>
void TopologyManager::CreateDataTypes<Dimension::Three>() {

  const int block_size[3] = { CC::TCX(), CC::TCY(), CC::TCZ() };

  //Boundary Exchange
  const int x_boundary_size[3] = { CC::HS(), CC::TCY(), CC::TCZ() };
  const int y_boundary_size[3] = { CC::TCX(), CC::HS(), CC::TCZ() };
  const int z_boundary_size[3] = { CC::TCX(), CC::TCY(), CC::HS() };

  const int east_boundary_start_index[3] = { CC::FHH(), 0, 0 };
  const int west_boundary_start_index[3] = { 0, 0, 0 };
  const int north_boundary_start_index[3] = { 0, CC::FHH(), 0 };
  const int south_boundary_start_index[3] = { 0, 0, 0 };
  const int top_boundary_start_index[3] = { 0, 0, CC::FHH() };
  const int bottom_boundary_start_index[3] = { 0, 0, 0 };

  const int east_slice_start_index[3] = { CC::FHH() - CC::HS(), 0, 0 };
  const int west_slice_start_index[3] = { CC::HS(), 0, 0 };
  const int north_slice_start_index[3] = { 0, CC::FHH() - CC::HS(), 0 };
  const int south_slice_start_index[3] = { 0, CC::HS(), 0 };
  const int top_slice_start_index[3] = { 0, 0, CC::FHH() - CC::HS() };
  const int bottom_slice_start_index[3] = { 0, 0, CC::HS() };

  MPI_Type_create_subarray(3, block_size, x_boundary_size, east_boundary_start_index,   MPI_ORDER_C,MPI_DOUBLE,&east_boundary_array_);
  MPI_Type_create_subarray(3, block_size, x_boundary_size, west_boundary_start_index,   MPI_ORDER_C,MPI_DOUBLE,&west_boundary_array_);
  MPI_Type_create_subarray(3, block_size, y_boundary_size, north_boundary_start_index,  MPI_ORDER_C,MPI_DOUBLE,&north_boundary_array_);
  MPI_Type_create_subarray(3, block_size, y_boundary_size, south_boundary_start_index,  MPI_ORDER_C,MPI_DOUBLE,&south_boundary_array_);
  MPI_Type_create_subarray(3, block_size, z_boundary_size, top_boundary_start_index,    MPI_ORDER_C,MPI_DOUBLE,&top_boundary_array_);
  MPI_Type_create_subarray(3, block_size, z_boundary_size, bottom_boundary_start_index, MPI_ORDER_C,MPI_DOUBLE,&bottom_boundary_array_);

  MPI_Type_commit(&east_boundary_array_);
  MPI_Type_commit(&west_boundary_array_);
  MPI_Type_commit(&north_boundary_array_);
  MPI_Type_commit(&south_boundary_array_);
  MPI_Type_commit(&top_boundary_array_);
  MPI_Type_commit(&bottom_boundary_array_);

  MPI_Type_create_subarray(3, block_size, x_boundary_size, east_slice_start_index,   MPI_ORDER_C,MPI_DOUBLE,&east_boundary_slice_);
  MPI_Type_create_subarray(3, block_size, x_boundary_size, west_slice_start_index,   MPI_ORDER_C,MPI_DOUBLE,&west_boundary_slice_);
  MPI_Type_create_subarray(3, block_size, y_boundary_size, north_slice_start_index,  MPI_ORDER_C,MPI_DOUBLE,&north_boundary_slice_);
  MPI_Type_create_subarray(3, block_size, y_boundary_size, south_slice_start_index,  MPI_ORDER_C,MPI_DOUBLE,&south_boundary_slice_);
  MPI_Type_create_subarray(3, block_size, z_boundary_size, top_slice_start_index,    MPI_ORDER_C,MPI_DOUBLE,&top_boundary_slice_);
  MPI_Type_create_subarray(3, block_size, z_boundary_size, bottom_slice_start_index, MPI_ORDER_C,MPI_DOUBLE,&bottom_boundary_slice_);

  MPI_Type_commit(&east_boundary_slice_);
  MPI_Type_commit(&west_boundary_slice_);
  MPI_Type_commit(&north_boundary_slice_);
  MPI_Type_commit(&south_boundary_slice_);
  MPI_Type_commit(&top_boundary_slice_);
  MPI_Type_commit(&bottom_boundary_slice_);

  // Send Types for Jump Conditions -> We send the NON-Predicted values
  const int child_0_start_west_index[3] = { CC::PLCS(), CC::PLCS(), CC::PLCS() };
  const int child_0_start_south_index[3] = { CC::PLCS(), CC::PLCS(), CC::PLCS() };
  const int child_0_start_bottom_index[3] = { CC::PLCS(), CC::PLCS(), CC::PLCS() };

  const int child_1_start_east_index[3] = { CC::BPHCS(), CC::PLCS(), CC::PLCS() };
  const int child_1_start_south_index[3] = { CC::PHCS(), CC::PLCS(), CC::PLCS() };
  const int child_1_start_bottom_index[3] = { CC::PHCS(), CC::PLCS(), CC::PLCS() };

  const int child_2_start_west_index[3] = { CC::PLCS(), CC::PHCS(), CC::PLCS() };
  const int child_2_start_north_index[3] = { CC::PLCS(), CC::BPHCS(), CC::PLCS() };
  const int child_2_start_bottom_index[3] = { CC::PLCS(), CC::PHCS(), CC::PLCS() };

  const int child_3_start_east_index[3] = { CC::BPHCS(), CC::PHCS(), CC::PLCS() };
  const int child_3_start_north_index[3] = { CC::PHCS(), CC::BPHCS(), CC::PLCS() };
  const int child_3_start_bottom_index[3] = { CC::PHCS(), CC::PHCS(), CC::PLCS() };

  const int child_4_start_west_index[3] = { CC::PLCS(), CC::PLCS(), CC::PHCS() };
  const int child_4_start_south_index[3] = { CC::PLCS(), CC::PLCS(), CC::PHCS() };
  const int child_4_start_top_index[3] = { CC::PLCS(), CC::PLCS(), CC::BPHCS() };

  const int child_5_start_east_index[3] = { CC::BPHCS(), CC::PLCS(), CC::PHCS() };
  const int child_5_start_south_index[3] = { CC::PHCS(), CC::PLCS(), CC::PHCS() };
  const int child_5_start_top_index[3] = { CC::PHCS(), CC::PLCS(), CC::BPHCS() };

  const int child_6_start_west_index[3] = { CC::PLCS(), CC::PHCS(), CC::PHCS() };
  const int child_6_start_north_index[3] = { CC::PLCS(), CC::BPHCS(), CC::PHCS() };
  const int child_6_start_top_index[3] = { CC::PLCS(), CC::PHCS(), CC::BPHCS() };

  const int child_7_start_east_index[3] = { CC::BPHCS(), CC::PHCS(), CC::PHCS() };
  const int child_7_start_north_index[3] = { CC::PHCS(), CC::BPHCS(), CC::PHCS() };
  const int child_7_start_top_index[3] = { CC::PHCS(), CC::PHCS(), CC::BPHCS() };

  const int x_jump_size[3] = { CC::BPS(), CC::PHS(), CC::PHS() };
  const int y_jump_size[3] = { CC::PHS(), CC::BPS(), CC::PHS() };
  const int z_jump_size[3] = { CC::PHS(), CC::PHS(), CC::BPS() };

  MPI_Type_create_subarray(3, block_size, x_jump_size, child_1_start_east_index, MPI_ORDER_C,MPI_DOUBLE,&east_jump_child_1_);
  MPI_Type_create_subarray(3, block_size, x_jump_size, child_3_start_east_index, MPI_ORDER_C,MPI_DOUBLE,&east_jump_child_3_);
  MPI_Type_create_subarray(3, block_size, x_jump_size, child_5_start_east_index, MPI_ORDER_C,MPI_DOUBLE,&east_jump_child_5_);
  MPI_Type_create_subarray(3, block_size, x_jump_size, child_7_start_east_index, MPI_ORDER_C,MPI_DOUBLE,&east_jump_child_7_);

  MPI_Type_create_subarray(3, block_size, x_jump_size, child_0_start_west_index, MPI_ORDER_C,MPI_DOUBLE,&west_jump_child_0_);
  MPI_Type_create_subarray(3, block_size, x_jump_size, child_2_start_west_index, MPI_ORDER_C,MPI_DOUBLE,&west_jump_child_2_);
  MPI_Type_create_subarray(3, block_size, x_jump_size, child_4_start_west_index, MPI_ORDER_C,MPI_DOUBLE,&west_jump_child_4_);
  MPI_Type_create_subarray(3, block_size, x_jump_size, child_6_start_west_index, MPI_ORDER_C,MPI_DOUBLE,&west_jump_child_6_);

  MPI_Type_create_subarray(3, block_size, y_jump_size, child_2_start_north_index, MPI_ORDER_C,MPI_DOUBLE,&north_jump_child_2_);
  MPI_Type_create_subarray(3, block_size, y_jump_size, child_3_start_north_index, MPI_ORDER_C,MPI_DOUBLE,&north_jump_child_3_);
  MPI_Type_create_subarray(3, block_size, y_jump_size, child_6_start_north_index, MPI_ORDER_C,MPI_DOUBLE,&north_jump_child_6_);
  MPI_Type_create_subarray(3, block_size, y_jump_size, child_7_start_north_index, MPI_ORDER_C,MPI_DOUBLE,&north_jump_child_7_);

  MPI_Type_create_subarray(3, block_size, y_jump_size, child_0_start_south_index, MPI_ORDER_C,MPI_DOUBLE,&south_jump_child_0_);
  MPI_Type_create_subarray(3, block_size, y_jump_size, child_1_start_south_index, MPI_ORDER_C,MPI_DOUBLE,&south_jump_child_1_);
  MPI_Type_create_subarray(3, block_size, y_jump_size, child_4_start_south_index, MPI_ORDER_C,MPI_DOUBLE,&south_jump_child_4_);
  MPI_Type_create_subarray(3, block_size, y_jump_size, child_5_start_south_index, MPI_ORDER_C,MPI_DOUBLE,&south_jump_child_5_);

  MPI_Type_create_subarray(3, block_size, z_jump_size, child_4_start_top_index, MPI_ORDER_C,MPI_DOUBLE,&top_jump_child_4_);
  MPI_Type_create_subarray(3, block_size, z_jump_size, child_5_start_top_index, MPI_ORDER_C,MPI_DOUBLE,&top_jump_child_5_);
  MPI_Type_create_subarray(3, block_size, z_jump_size, child_6_start_top_index, MPI_ORDER_C,MPI_DOUBLE,&top_jump_child_6_);
  MPI_Type_create_subarray(3, block_size, z_jump_size, child_7_start_top_index, MPI_ORDER_C,MPI_DOUBLE,&top_jump_child_7_);

  MPI_Type_create_subarray(3, block_size, z_jump_size, child_0_start_bottom_index, MPI_ORDER_C,MPI_DOUBLE,&bottom_jump_child_0_);
  MPI_Type_create_subarray(3, block_size, z_jump_size, child_1_start_bottom_index, MPI_ORDER_C,MPI_DOUBLE,&bottom_jump_child_1_);
  MPI_Type_create_subarray(3, block_size, z_jump_size, child_2_start_bottom_index, MPI_ORDER_C,MPI_DOUBLE,&bottom_jump_child_2_);
  MPI_Type_create_subarray(3, block_size, z_jump_size, child_3_start_bottom_index, MPI_ORDER_C,MPI_DOUBLE,&bottom_jump_child_3_);

  MPI_Type_commit(&east_jump_child_1_);
  MPI_Type_commit(&east_jump_child_3_);
  MPI_Type_commit(&east_jump_child_5_);
  MPI_Type_commit(&east_jump_child_7_);
  MPI_Type_commit(&west_jump_child_0_);
  MPI_Type_commit(&west_jump_child_2_);
  MPI_Type_commit(&west_jump_child_4_);
  MPI_Type_commit(&west_jump_child_6_);
  MPI_Type_commit(&north_jump_child_2_);
  MPI_Type_commit(&north_jump_child_3_);
  MPI_Type_commit(&north_jump_child_6_);
  MPI_Type_commit(&north_jump_child_7_);
  MPI_Type_commit(&south_jump_child_0_);
  MPI_Type_commit(&south_jump_child_1_);
  MPI_Type_commit(&south_jump_child_4_);
  MPI_Type_commit(&south_jump_child_5_);
  MPI_Type_commit(&top_jump_child_4_);
  MPI_Type_commit(&top_jump_child_5_);
  MPI_Type_commit(&top_jump_child_6_);
  MPI_Type_commit(&top_jump_child_7_);
  MPI_Type_commit(&bottom_jump_child_0_);
  MPI_Type_commit(&bottom_jump_child_1_);
  MPI_Type_commit(&bottom_jump_child_2_);
  MPI_Type_commit(&bottom_jump_child_3_);
}

/**
 * @brief Frees MPI datatypes used in three dimensional simulations. See Meta function.
 */
template<>
void TopologyManager::FreeTypes<Dimension::Three>() {

  //No Jump Halo data types
  MPI_Type_free(&east_boundary_array_);
  MPI_Type_free(&east_boundary_slice_);
  MPI_Type_free(&west_boundary_array_);
  MPI_Type_free(&west_boundary_slice_);
  MPI_Type_free(&north_boundary_array_);
  MPI_Type_free(&north_boundary_slice_);
  MPI_Type_free(&south_boundary_array_);
  MPI_Type_free(&south_boundary_slice_);
  MPI_Type_free(&top_boundary_array_);
  MPI_Type_free(&top_boundary_slice_);
  MPI_Type_free(&bottom_boundary_array_);
  MPI_Type_free(&bottom_boundary_slice_);

  //Jump Halo data types
  MPI_Type_free(&east_jump_child_1_);
  MPI_Type_free(&east_jump_child_3_);
  MPI_Type_free(&east_jump_child_5_);
  MPI_Type_free(&east_jump_child_7_);
  MPI_Type_free(&west_jump_child_0_);
  MPI_Type_free(&west_jump_child_2_);
  MPI_Type_free(&west_jump_child_4_);
  MPI_Type_free(&west_jump_child_6_);
  MPI_Type_free(&north_jump_child_2_);
  MPI_Type_free(&north_jump_child_3_);
  MPI_Type_free(&north_jump_child_6_);
  MPI_Type_free(&north_jump_child_7_);
  MPI_Type_free(&south_jump_child_0_);
  MPI_Type_free(&south_jump_child_1_);
  MPI_Type_free(&south_jump_child_4_);
  MPI_Type_free(&south_jump_child_5_);
  MPI_Type_free(&top_jump_child_4_);
  MPI_Type_free(&top_jump_child_5_);
  MPI_Type_free(&top_jump_child_6_);
  MPI_Type_free(&top_jump_child_7_);
  MPI_Type_free(&bottom_jump_child_0_);
  MPI_Type_free(&bottom_jump_child_1_);
  MPI_Type_free(&bottom_jump_child_2_);
  MPI_Type_free(&bottom_jump_child_3_);
}

/**
 * @brief Datatypes for 2D Simulations. See Meta function.
 */
template<>
void TopologyManager::CreateDataTypes<Dimension::Two>() {

  const int block_size[3] = { CC::TCX(), CC::TCY(), CC::TCZ() };

  //Boundary Exchange
  const int x_boundary_size[3] = { CC::HS(), CC::TCY(), CC::TCZ() };
  const int y_boundary_size[3] = { CC::TCX(), CC::HS(), CC::TCZ() };

  const int east_boundary_start_index[3] = { CC::FHH(), 0, 0 };
  const int west_boundary_start_index[3] = { 0, 0, 0 };
  const int north_boundary_start_index[3] = { 0, CC::FHH(), 0 };
  const int south_boundary_start_index[3] = { 0, 0, 0 };

  const int east_slice_start_index[3] = { CC::FHH() - CC::HS(), 0, 0 };
  const int west_slice_start_index[3] = { CC::HS(), 0, 0 };
  const int north_slice_start_index[3] = { 0, CC::FHH() - CC::HS(), 0 };
  const int south_slice_start_index[3] = { 0, CC::HS(), 0 };

  MPI_Type_create_subarray(3, block_size, x_boundary_size, east_boundary_start_index,  MPI_ORDER_C,MPI_DOUBLE,&east_boundary_array_);
  MPI_Type_create_subarray(3, block_size, x_boundary_size, west_boundary_start_index,  MPI_ORDER_C,MPI_DOUBLE,&west_boundary_array_);
  MPI_Type_create_subarray(3, block_size, y_boundary_size, north_boundary_start_index, MPI_ORDER_C,MPI_DOUBLE,&north_boundary_array_);
  MPI_Type_create_subarray(3, block_size, y_boundary_size, south_boundary_start_index, MPI_ORDER_C,MPI_DOUBLE,&south_boundary_array_);

  MPI_Type_commit(&east_boundary_array_);
  MPI_Type_commit(&west_boundary_array_);
  MPI_Type_commit(&north_boundary_array_);
  MPI_Type_commit(&south_boundary_array_);

  MPI_Type_create_subarray(3, block_size, x_boundary_size, east_slice_start_index,  MPI_ORDER_C,MPI_DOUBLE,&east_boundary_slice_);
  MPI_Type_create_subarray(3, block_size, x_boundary_size, west_slice_start_index,  MPI_ORDER_C,MPI_DOUBLE,&west_boundary_slice_);
  MPI_Type_create_subarray(3, block_size, y_boundary_size, north_slice_start_index, MPI_ORDER_C,MPI_DOUBLE,&north_boundary_slice_);
  MPI_Type_create_subarray(3, block_size, y_boundary_size, south_slice_start_index, MPI_ORDER_C,MPI_DOUBLE,&south_boundary_slice_);

  MPI_Type_commit(&east_boundary_slice_);
  MPI_Type_commit(&west_boundary_slice_);
  MPI_Type_commit(&north_boundary_slice_);
  MPI_Type_commit(&south_boundary_slice_);

  // Send Types for Jump Conditions -> We send the NON-Predicted values
  const int child_0_start_west_index[3] = { CC::PLCS(), CC::PLCS(), 0 };
  const int child_0_start_south_index[3] = { CC::PLCS(), CC::PLCS(), 0 };

  const int child_1_start_east_index[3] = { CC::BPHCS(), CC::PLCS(), 0 };
  const int child_1_start_south_index[3] = { CC::PHCS(), CC::PLCS(), 0 };

  const int child_2_start_west_index[3] = { CC::PLCS(), CC::PHCS(), 0 };
  const int child_2_start_north_index[3] = { CC::PLCS(), CC::BPHCS(), 0 };

  const int child_3_start_east_index[3] = { CC::BPHCS(), CC::PHCS(), 0 };
  const int child_3_start_north_index[3] = { CC::PHCS(), CC::BPHCS(), 0 };

  const int x_jump_size[3] = { CC::BPS(), CC::PHS(), 1 };
  const int y_jump_size[3] = { CC::PHS(), CC::BPS(), 1 };

  MPI_Type_create_subarray(3, block_size, x_jump_size, child_1_start_east_index, MPI_ORDER_C,MPI_DOUBLE,&east_jump_child_1_);
  MPI_Type_create_subarray(3, block_size, x_jump_size, child_3_start_east_index, MPI_ORDER_C,MPI_DOUBLE,&east_jump_child_3_);

  MPI_Type_create_subarray(3, block_size, x_jump_size, child_0_start_west_index, MPI_ORDER_C,MPI_DOUBLE,&west_jump_child_0_);
  MPI_Type_create_subarray(3, block_size, x_jump_size, child_2_start_west_index, MPI_ORDER_C,MPI_DOUBLE,&west_jump_child_2_);

  MPI_Type_create_subarray(3, block_size, y_jump_size, child_2_start_north_index, MPI_ORDER_C,MPI_DOUBLE,&north_jump_child_2_);
  MPI_Type_create_subarray(3, block_size, y_jump_size, child_3_start_north_index, MPI_ORDER_C,MPI_DOUBLE,&north_jump_child_3_);

  MPI_Type_create_subarray(3, block_size, y_jump_size, child_0_start_south_index, MPI_ORDER_C,MPI_DOUBLE,&south_jump_child_0_);
  MPI_Type_create_subarray(3, block_size, y_jump_size, child_1_start_south_index, MPI_ORDER_C,MPI_DOUBLE,&south_jump_child_1_);

  MPI_Type_commit(&east_jump_child_1_);
  MPI_Type_commit(&east_jump_child_3_);

  MPI_Type_commit(&west_jump_child_0_);
  MPI_Type_commit(&west_jump_child_2_);

  MPI_Type_commit(&north_jump_child_2_);
  MPI_Type_commit(&north_jump_child_3_);

  MPI_Type_commit(&south_jump_child_0_);
  MPI_Type_commit(&south_jump_child_1_);

}

/**
 * @brief Frees MPI datatypes used in 2D simulations. See Meta function.
 */
template<>
void TopologyManager::FreeTypes<Dimension::Two>() {

  //No Jump Halo data types
  MPI_Type_free(&east_boundary_array_);
  MPI_Type_free(&east_boundary_slice_);
  MPI_Type_free(&west_boundary_array_);
  MPI_Type_free(&west_boundary_slice_);
  MPI_Type_free(&north_boundary_array_);
  MPI_Type_free(&north_boundary_slice_);
  MPI_Type_free(&south_boundary_array_);
  MPI_Type_free(&south_boundary_slice_);

  //Jump Halo data types
  MPI_Type_free(&east_jump_child_1_);
  MPI_Type_free(&east_jump_child_3_);
  MPI_Type_free(&west_jump_child_0_);
  MPI_Type_free(&west_jump_child_2_);
  MPI_Type_free(&north_jump_child_2_);
  MPI_Type_free(&north_jump_child_3_);
  MPI_Type_free(&south_jump_child_0_);
  MPI_Type_free(&south_jump_child_1_);
}

/**
 * @brief Datatypes for 1D Simulations. See Meta function.
 */
template<>
void TopologyManager::CreateDataTypes<Dimension::One>() {

  const int block_size[3] = { CC::TCX(), CC::TCY(), CC::TCZ() };

  //Boundary Exchange
  const int x_boundary_size[3] = { CC::HS(), CC::TCY(), CC::TCZ() };

  const int east_boundary_start_index[3] = { CC::FHH(), 0, 0 };
  const int west_boundary_start_index[3] = { 0, 0, 0 };

  const int east_slice_start_index[3] = { CC::FHH() - CC::HS(), 0, 0 };
  const int west_slice_start_index[3] = { CC::HS(), 0, 0 };

  MPI_Type_create_subarray(3, block_size, x_boundary_size, east_boundary_start_index, MPI_ORDER_C,MPI_DOUBLE,&east_boundary_array_);
  MPI_Type_create_subarray(3, block_size, x_boundary_size, west_boundary_start_index, MPI_ORDER_C,MPI_DOUBLE,&west_boundary_array_);

  MPI_Type_commit(&east_boundary_array_);
  MPI_Type_commit(&west_boundary_array_);

  MPI_Type_create_subarray(3, block_size, x_boundary_size, east_slice_start_index, MPI_ORDER_C,MPI_DOUBLE,&east_boundary_slice_);
  MPI_Type_create_subarray(3, block_size, x_boundary_size, west_slice_start_index, MPI_ORDER_C,MPI_DOUBLE,&west_boundary_slice_);

  MPI_Type_commit(&east_boundary_slice_);
  MPI_Type_commit(&west_boundary_slice_);

  // Send Types for Jump Conditions -> We send the NON-Predicted values
  const int child_0_start_west_index[3] = { CC::PLCS(), 0, 0 };

  const int child_1_start_east_index[3] = { CC::BPHCS(), 0, 0 };

  const int x_jump_size[3] = { CC::BPS(), 1, 1 };

  MPI_Type_create_subarray(3, block_size, x_jump_size, child_1_start_east_index, MPI_ORDER_C,MPI_DOUBLE,&east_jump_child_1_);

  MPI_Type_create_subarray(3, block_size, x_jump_size, child_0_start_west_index, MPI_ORDER_C,MPI_DOUBLE,&west_jump_child_0_);

  MPI_Type_commit(&east_jump_child_1_);

  MPI_Type_commit(&west_jump_child_0_);

}

/**
 * @brief Frees MPI datatypes used in 1D simulations. See Meta function.
 */
template<>
void TopologyManager::FreeTypes<Dimension::One>() {

  //No Jump Halo data types
  MPI_Type_free(&east_boundary_array_);
  MPI_Type_free(&east_boundary_slice_);
  MPI_Type_free(&west_boundary_array_);
  MPI_Type_free(&west_boundary_slice_);

  //Jump Halo data types
  MPI_Type_free(&east_jump_child_1_);
  MPI_Type_free(&west_jump_child_0_);
}

/**
 * @brief Clean-Up Constuctor
 */
TopologyManager::~TopologyManager() {
  FreeTypes<CC::DIM()>();
}

/**
 * @brief Gives the rank id from MPI directly as int. Avoids handle creation, e.g. for const members in initializer list.
 * @return Rank id
 */
int TopologyManager::ForwardRankId() const {
  int rank_id = -1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_id);
  return rank_id;
}

/**
 * @brief Gives the number of ranks in the MPI communicator. Avoids handle creation, e.g. for const members in initializer list.
 * @return Communicator Size which is the number of ranks.
 */
int TopologyManager::ForwardCommunicatorSize() const {
  int communicator_size = -1;
  MPI_Comm_size(MPI_COMM_WORLD,&communicator_size);
  return communicator_size;
}

/**
 * @brief Gives the MPI_TAG_UB. Avoids handle creation, e.g. for const members in initializer list.
 * @return MPI_TAG_UB
 */
int TopologyManager::ForwardTagUb() const {
  int *tag_ub;
  int flag;
  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &flag);
  return *tag_ub;
}

/**
 * @brief Default constructor. Creates the local and global Id list on level zero and initializes the MPI datatypes.
 * @param setup Reference used to incorporate user settings, in particular the geometry of the domain.
 */
TopologyManager::TopologyManager(const SimulationSetup& setup) : setup_(setup), rank_id_(ForwardRankId()), number_of_ranks_(ForwardCommunicatorSize()), mpi_tag_ub_(ForwardTagUb()), forest_{} {
  std::uint64_t id = Node::HeadbitOnLevel(0);

  std::vector<std::uint64_t> initialization_list;
  /*
   * The Nodes are created in a spatial fashion traversing through X, Y and finally Z.
   * The ids in this traversal are not continuous due to the implicit shadow levels in the tree
   * Therefore some non straightforward index magic needs to be applied
   * Initialization happens only on Level 0
   */
  for (unsigned short i = 0; i < setup_.GetLevelZeroBlocksZ(); ++i) {
      for (unsigned short j = 0; j < setup_.GetLevelZeroBlocksY(); ++j) {
          for (unsigned short k = 0; k < setup_.GetLevelZeroBlocksX(); ++k) {
              forest_.emplace_back(TopologyNode(id, 0));
              initialization_list.push_back(id);
              id = Node::EastNeigborOfNodeWithId(initialization_list.back()); //find eastern neighbor of the just created Node
          }
          // Index magic to create the correct node once the (outer) loop counters are reseted
          id = Node::NorthNeigborOfNodeWithId(initialization_list[setup_.GetLevelZeroBlocksX() * (j + setup_.GetLevelZeroBlocksY() * i)]);
      }
      id = Node::TopNeigborOfNodeWithId(initialization_list[(setup_.GetLevelZeroBlocksY() * setup_.GetLevelZeroBlocksX() * i)]);
  }

  //Topology Tree Node creation
  forest_.shrink_to_fit();

  // Now we distribute the global nodes accordingly:
  unsigned int number_of_nodes = forest_.size();
  int remainder = number_of_nodes % number_of_ranks_;
  int fraction = 0;

  std::vector<int> rank_map;

  for (int i = 0; i < number_of_ranks_; ++i) {
      fraction = 0;
      fraction += number_of_nodes / number_of_ranks_;
      if (i < remainder) {
          fraction++;
      }
      rank_map.push_back(fraction);
  }

  //Topology Tree Rank setup
  int sum = 0;
  for (int i = 0; i < number_of_ranks_; i++) {
      for (int j = 0; j < rank_map[i]; j++) {
          forest_[sum].SetCurrentRankOfLeaf(forest_[sum].Id(), i);
          sum++;
      }
  }

  CreateDataTypes<CC::DIM()>();
}

/**
 * @brief Updates the id list based on the recorded refinements, coarsening and weights.
 */
void TopologyManager::UpdateIdList() {

  int length;
  std::vector<int> all_lengths(number_of_ranks_);

  //Tree update
  //refine
  length = local_refine_list_.size();
  MPI_Allgather(&length, 1, MPI_INT, all_lengths.data(), 1, MPI_INT,MPI_COMM_WORLD);
  int sum = std::accumulate(all_lengths.begin(), all_lengths.end(), 0);
  std::vector<std::uint64_t> global_refine_list(sum);
  std::vector<int> displacements(number_of_ranks_);

  sum = 0;
  for (int i = 0; i < number_of_ranks_; i++) {
      displacements[i] = sum;
      sum += all_lengths[i];
  }
  MPI_Allgatherv(local_refine_list_.data(), length, MPI_LONG_LONG_INT, global_refine_list.data(), all_lengths.data(), displacements.data(), MPI_LONG_LONG_INT,MPI_COMM_WORLD);

  //global_refine_list.erase(global_refine_list.begin()+all_displs[rank_id_], global_refine_list.begin()+all_displs[rank_id_]+length);
  for (const auto& refine_id : global_refine_list) {
    forest_[PositionOfNodeInZeroTopology(refine_id)].Refine(refine_id);
  }
  local_refine_list_.clear();

  //weight exchange

  length=std::get<0>(local_weight_list_).size();
  MPI_Allgather(&length, 1, MPI_INT, all_lengths.data(), 1, MPI_INT,MPI_COMM_WORLD);
  sum = std::accumulate(all_lengths.begin(), all_lengths.end(), 0);
  std::tuple<std::vector<std::uint64_t>, std::vector<int>> globalWeigthList;
  std::get<0>(globalWeigthList).reserve(sum);
  std::get<1>(globalWeigthList).reserve(sum);

  sum=0;
  for (int i = 0; i < number_of_ranks_; i++) {
      displacements[i] = sum;
      sum += all_lengths[i];
  }

  MPI_Allgatherv(std::get<0>(local_weight_list_).data(), length, MPI_LONG_LONG_INT, std::get<0>(globalWeigthList).data(), all_lengths.data(), displacements.data(), MPI_LONG_LONG_INT,MPI_COMM_WORLD);
  MPI_Allgatherv(std::get<1>(local_weight_list_).data(), length, MPI_INT, std::get<1>(globalWeigthList).data(), all_lengths.data(), displacements.data(), MPI_INT,MPI_COMM_WORLD);

  for(int i=0;i<sum;i++){
    int position = PositionOfNodeInZeroTopology(std::get<0>(globalWeigthList)[i]);
                         // Forest size cannot exceed 128x128x128 = 10^6, fits in int (10^9 positive values), cast is safe.
    if(position >= 0 && position < static_cast<int>(forest_.size())){
        forest_[position].SetWeightOfId(std::get<0>(globalWeigthList)[i], std::get<1>(globalWeigthList)[i]);
    }else{
        throw std::logic_error("TopologyManager::UpdateIdsList root tree does not exist!");
    }
  }

  //coarse
  length = local_coarse_list_.size();
  MPI_Allgather(&length, 1, MPI_INT, all_lengths.data(), 1, MPI_INT,MPI_COMM_WORLD);
  sum = std::accumulate(all_lengths.begin(), all_lengths.end(), 0);
  std::vector<std::uint64_t> global_coarse_list(sum);

  sum = 0;
  for (int i = 0; i < number_of_ranks_; i++) {
      displacements[i] = sum;
      sum += all_lengths[i];
  }

  MPI_Allgatherv(local_coarse_list_.data(), length, MPI_LONG_LONG_INT, global_coarse_list.data(), all_lengths.data(), displacements.data(), MPI_LONG_LONG_INT,MPI_COMM_WORLD);

  for (const auto& coarse_id : global_coarse_list) {
      forest_[PositionOfNodeInZeroTopology(coarse_id)].Coarse(coarse_id);
  }
  local_coarse_list_.clear();
  GenerateTags();
}

/**
 * @brief Marks the node with the given id for refinement.
 * @param id The id of the leaf that is to be refined.
 * @note The actual refinement happens in a bundled fashion in another function. $No checks for leaves are performed caller is responsible$
 */
void TopologyManager::RefineNodeWithId(const std::uint64_t id) {
  local_refine_list_.push_back(id);
}

/**
 * @brief Adds the parent whose children may be coarsened to the coarsening list. The actual data deletion is than bundled using this list.
 * @param parent_id The id of the node that is to be made a leaf.
 */
void TopologyManager::CoarseNodeWithId(const std::uint64_t parent_id) {
  local_coarse_list_.push_back(parent_id);
}

/**
 * @brief Determines the rank which holds the node of given id.
 * @param id The unique id of the node.
 * @return Rank id for the requested node.
 */
int TopologyManager::GetRankOfNode(const std::uint64_t id) const {
  return forest_[PositionOfNodeInZeroTopology(id)].GetRank(id);
}

/**
 * @brief Gives a list which indicates which node should go from which mpi rank onto which mpi rank.
 * @param level The level on which balancing is to be executed.
 * @return A vector of all nodes and their current as well as their future mpi rank.
 */
std::vector<std::tuple<std::uint64_t, int, int>> TopologyManager::GetLoadBalancedTopology() {

  AssignTargetRankToLeaves();
  AssignBalancedLoad();
  std::vector<std::tuple<std::uint64_t, int, int>> ids_current_future_rank_map;

  ListNodeToBalance(ids_current_future_rank_map);

  return ids_current_future_rank_map;
}

/**
 * @brief Gives a unique tag for a specific node to be used in non-blocking MPI communication.
 * @param id The unique id of the node which should receive the information with this tag.
 * @return Tag.
 * @note The tag range is restricted by the MPI implementation to 2^(31) - 1. This is not enough to use the unique_id of the nodes as tag.
 */
int TopologyManager::GetTagForNode(const std::uint64_t id) const {
  int tag = forest_[PositionOfNodeInZeroTopology(id)].TagForId(id);
  int tag_test = tag;
  if( ((tag_test << 5)) > mpi_tag_ub_) {throw std::logic_error("You ran out of tags");}
  return tag;
}

/**
 * @brief Indicates whether a node exists in the global Tree, does not make implications about local tree
 * @param id id of the node one is looking for
 * @return true if node exists, false otherwise
 */
bool TopologyManager::NodeExists(const std::uint64_t id) const {

  int position_in_zero_topology = PositionOfNodeInZeroTopology(id);
                                                                  // Forest size cannot exceed 128x128x128 = 10^6, fits in int (10^9 positive values), cast is safe.
  if (position_in_zero_topology >= 0 && position_in_zero_topology < (static_cast<int>(forest_.size())) ) {
      return forest_[position_in_zero_topology].NodeExists(id);
  } else {
      return false;
  }
}

/**
 * @brief Gives the current maximum level of any global node.
 *        Can be less than the user set maximum level (if no interesting physics are present, or at initialization).
 * @return Globally Maximal Present Level.
 */
unsigned int TopologyManager::GetCurrentMaximumLevel() const {

  std::vector<unsigned int> tree_depths;
  tree_depths.reserve(forest_.size());
  for(const TopologyNode& tree : forest_) {
    tree_depths.emplace_back(tree.GetDepth());
  }

  return *std::max_element(tree_depths.begin(),tree_depths.end()) -1; //Tree Depth starts with 1 by definition, vs levels starts at 0.
}

/**
 * @brief Gives the ids of all globally existing nodes which descent from the given id.
 *        I.e. children, grand-children great-grand-children, ...
 * @param id Unique id of node whose descendants are searched for.
 * @return All ids of globally existing descendants.
 */
std::vector<std::uint64_t> TopologyManager::DescendantIdsOfNode(const std::uint64_t id) const {

  std::vector<std::uint64_t> descendants;
  std::vector<std::uint64_t> append_list;

  for(const auto& child_id : Node::IdsOfChildren(id)) {
    if(NodeExists(child_id)) {
        append_list = DescendantIdsOfNode(child_id);
        descendants.insert(descendants.end(),append_list.begin(), append_list.end());
        descendants.push_back(child_id);
    }
  }

  return descendants;
}

/**
 * @brief Indicates whether the invoking MPI rank holds the given node. $Throws exception if node does not exist!$
 * @param id Unique id of the node to be checked for self-ownership.
 * @return true if node is on the same rank as invoker, false otherwise.
 */
bool TopologyManager::NodeIsOnMyRank(const std::uint64_t id) const {
  if (NodeExists(id)) {
      return GetRankOfNode(id) == rank_id_;
  } else {
      throw std::logic_error("Node Ownership cannot be checked - Node does not exist");
  }
}

/**
 * @brief Determines if the specified node has a jump at the given direction.
 *        I.e. face does not have a global boundary and no neighbor on the same level exists.
 * @param id Unique id of the node under consideration.
 * @param location Direction of interest.
 * @return True if the Face is a Jump, false otherwise.
 */
bool TopologyManager::FaceIsJump(const std::uint64_t id, const BoundaryLocation location) const {

  bool is_external = Node::IsExternalBoundary(location, id, setup_);

  // If the neigbor does not exist and it is not an external BC we have a jump
      switch (location) {
      case BoundaryLocation::eEast:
          return (!NodeExists(Node::EastNeigborOfNodeWithId(id)) && !is_external);
          break;
      case BoundaryLocation::eWest:
          return (!NodeExists(Node::WestNeigborOfNodeWithId(id)) && !is_external);
          break;
      case BoundaryLocation::eNorth:
          return (!NodeExists(Node::NorthNeigborOfNodeWithId(id)) && !is_external);
          break;
      case BoundaryLocation::eSouth:
          return (!NodeExists(Node::SouthNeigborOfNodeWithId(id)) && !is_external);
          break;
      case BoundaryLocation::eTop:
          return (!NodeExists(Node::TopNeigborOfNodeWithId(id)) && !is_external);
          break;
      case BoundaryLocation::eBottom:
          return (!NodeExists(Node::BottomNeigborOfNodeWithId(id)) && !is_external);
          break;
      default:
          throw std::invalid_argument("Boundary Location not found - Impossible Error");
          break;
      }

}

/**
 * @brief Gives a list of ids of all globally present leaves.
 * @return Leaf Ids.
 */
std::vector<std::uint64_t> TopologyManager::GetLeafIds() const {

  std::vector<std::uint64_t> leaves;
  for (TopologyNode node : forest_) {
      node.GetLeafIds(leaves);
  }
  return leaves;
}

/**
 * @brief Gives the ids of all globally present leaves on the specified level.
 * @param level The level of interest.
 * @return The list of leaf ids.
 */
std::vector<std::uint64_t> TopologyManager::LeafIdsOnLevel(const unsigned int level) const {
  std::vector<std::uint64_t> leaves = GetLeafIds();
  std::vector<std::uint64_t> results;
  for (const auto& leaf : leaves) {
      if (Node::LevelOfNode(leaf) == level) {
          results.push_back(leaf);
      }
  }

  return results;
}

/**
 * @brief Give the datatype needed to send data large enough to use a prediction in order to
 *        fill this nodes halo at the given direction.
 * @param jump_boundary The node id and the direction that is to receive values.
 * @return MPI datatype needed for sending to this nodes jump from the parent.
 */
MPI_Datatype TopologyManager::GetDatatypeForJumpBoundary(const std::pair<std::uint64_t, BoundaryLocation> jump_boundary) const {

  /* TODO Make it templated to get rid of switch construct, but requries to lot of refactoring in MR because of conestexpr-ness in template parameter, postponed. */

  const unsigned short int position = Node::PositionOfNodeAmongSibilings(jump_boundary.first);
  const BoundaryLocation location = jump_boundary.second;

  switch (location) {
  case BoundaryLocation::eEast: {
      switch (position) {
      case 1:
          return east_jump_child_1_;
          break;
      case 3:
          return east_jump_child_3_;
          break;
      case 5:
          return east_jump_child_5_;
          break;
      case 7:
          return east_jump_child_7_;
          break;
      default:
          throw std::invalid_argument("Boundary Location not found in Tree FacIsJump - Impossible Error");
          break;
      }
  }
      break;
  case BoundaryLocation::eWest: {
      switch (position) {
      case 0:
          return west_jump_child_0_;
          break;
      case 2:
          return west_jump_child_2_;
          break;
      case 4:
          return west_jump_child_4_;
          break;
      case 6:
          return west_jump_child_6_;
          break;
      default:
          throw std::invalid_argument("Boundary Location not found in Tree FacIsJump - Impossible Error");
          break;
      }
  }
      break;
  case BoundaryLocation::eNorth: {
      switch (position) {
      case 2:
          return north_jump_child_2_;
          break;
      case 3:
          return north_jump_child_3_;
          break;
      case 6:
          return north_jump_child_6_;
          break;
      case 7:
          return north_jump_child_7_;
          break;
      default:
          throw std::invalid_argument("Boundary Location not found in Tree FacIsJump - Impossible Error");
          break;
      }
  }
      break;
  case BoundaryLocation::eSouth: {
      switch (position) {
      case 0:
          return south_jump_child_0_;
          break;
      case 1:
          return south_jump_child_1_;
          break;
      case 4:
          return south_jump_child_4_;
          break;
      case 5:
          return south_jump_child_5_;
          break;
      default:
          throw std::invalid_argument("Boundary Location not found in Tree FacIsJump - Impossible Error");
          break;
      }
  }
      break;
  case BoundaryLocation::eTop: {
      switch (position) {
      case 4:
          return top_jump_child_4_;
          break;
      case 5:
          return top_jump_child_5_;
          break;
      case 6:
          return top_jump_child_6_;
          break;
      case 7:
          return top_jump_child_7_;
          break;
      default:
          throw std::invalid_argument("Boundary Location not found in Tree FacIsJump - Impossible Error");
          break;
      }
  }
      break;
  case BoundaryLocation::eBottom: {
      switch (position) {
      case 0:
          return bottom_jump_child_0_;
          break;
      case 1:
          return bottom_jump_child_1_;
          break;
      case 2:
          return bottom_jump_child_2_;
          break;
      case 3:
          return bottom_jump_child_3_;
          break;
      default:
          throw std::invalid_argument("Boundary Location not found in Tree FacIsJump - Impossible Error");
          break;
      }
  }
      break;
  default:
      throw std::invalid_argument("Boundary Location not found in Tree FacIsJump - Impossible Error");
      break;
  }
}

/**
 * @brief Gives the MPI Datatype to send data into an internal no-jump exchange.
 * @param location The direction of the halo to be sent into.
 * @return The Datatype needed to send into the right sub_array.
 */
MPI_Datatype TopologyManager::NoJumpBoundaryDataType(const BoundaryLocation location) const {
// TODO Temaplte this function as soon as Boundary condition classes get templated.
  switch (location) {
  case BoundaryLocation::eEast:
      return east_boundary_array_;
      break;
  case BoundaryLocation::eWest:
      return west_boundary_array_;
      break;
  case BoundaryLocation::eNorth:
      return north_boundary_array_;
      break;
  case BoundaryLocation::eSouth:
      return south_boundary_array_;
      break;
  case BoundaryLocation::eTop:
      return top_boundary_array_;
      break;
  case BoundaryLocation::eBottom:
      return bottom_boundary_array_;
      break;
  default:
      throw std::invalid_argument("Boundary Location not found in Tree FacIsJump - Impossible Error");
      break;
  }

}

/**
 * @brief Gives the MPI Datatype to send data from an internal no-jump exchange.
 * @param location The direction of the halo to be sent from.
 * @return The Datatype needed to send from the right sub_array.
 */
MPI_Datatype TopologyManager::NoJumpSliceDataType(const BoundaryLocation location) const {
// TODO Template this function as soon as Boundary condition classes get templated.
  switch (location) {
  case BoundaryLocation::eEast:
      return east_boundary_slice_;
      break;
  case BoundaryLocation::eWest:
      return west_boundary_slice_;
      break;
  case BoundaryLocation::eNorth:
      return north_boundary_slice_;
      break;
  case BoundaryLocation::eSouth:
      return south_boundary_slice_;
      break;
  case BoundaryLocation::eTop:
      return top_boundary_slice_;
      break;
  case BoundaryLocation::eBottom:
      return bottom_boundary_slice_;
      break;
  default:
      throw std::invalid_argument("Boundary Location not found in Tree FacIsJump - Impossible Error");
      break;
  }
}

/**
 * @brief Assigns the target rank (rank on which the node SHOULD reside) to all leaf nodes.
 *        Uses either a linear or Hilbert Traversal to determine the target rank.
 */
void TopologyManager::AssignTargetRankToLeaves() {

  std::vector<std::vector<std::tuple<unsigned int, int>>> ids_current_future_rank_map; //[level][(count,rank)]
  ids_current_future_rank_map.reserve(number_of_ranks_);

  std::vector<int> weightList(WeightsOnLevels());

  for(unsigned int level=0; level <= setup_.GetMaximumLevel(); level++) {
    int weight_per_rank = weightList[level] / number_of_ranks_;
    int remainder       = weightList[level] % number_of_ranks_;
    ids_current_future_rank_map.push_back(std::vector<std::tuple<unsigned int, int>>());
    for(int rank=0; rank < number_of_ranks_; rank++) {
        if(rank < remainder) {
            ids_current_future_rank_map[level].push_back(std::make_tuple(weight_per_rank+1, rank));
        } else {
            ids_current_future_rank_map[level].push_back(std::make_tuple(weight_per_rank , rank));
        }
    }
    std::reverse(ids_current_future_rank_map[level].begin(), ids_current_future_rank_map[level].end());
  }

#ifndef HILBERT
  for(TopologyNode& node : forest_) {
      node.SetTargetRankForLeaf(ids_current_future_rank_map);
  }
#else

  //We travers the domain in a hilbert curve traversal, as we have shadow levels, however, some considerations on the forest_ have to be done so the curve is correct.

  bool x_inverse = false; //Inverse: reverse traversal when going back.
  bool y_inverse = false;
  for (unsigned short i = 0; i < setup_.GetLevelZeroBlocksZ(); ++i) {
    if(!y_inverse){
        for (short j = 0; j < setup_.GetLevelZeroBlocksY(); ++j) {
            if(!x_inverse){
                for (short k = 0; k < setup_.GetLevelZeroBlocksX(); ++k) {
                    int z = i*setup_.GetLevelZeroBlocksY()*setup_.GetLevelZeroBlocksX()+j*setup_.GetLevelZeroBlocksX()+k; //Implicit type conversions, but do no harm
                    forest_[z].SetTargetRankForLeaf(ids_current_future_rank_map, HilbertPosition::z_x_y);
                }
            }else{
                for (short k = setup_.GetLevelZeroBlocksX()-1; k >=0 ; --k) {
                    int z=i*setup_.GetLevelZeroBlocksY()*setup_.GetLevelZeroBlocksX()+j*setup_.GetLevelZeroBlocksX()+k; //Implicit type conversions, but do no harm
                    forest_[z].SetTargetRankForLeaf(ids_current_future_rank_map, HilbertPosition::_z_xy);
                }
            }
            x_inverse = !x_inverse;
        }
    } else{
        for (short j = setup_.GetLevelZeroBlocksY()-1 ; j >= 0; --j) {
            if(!x_inverse){
                for (short k = 0; k < setup_.GetLevelZeroBlocksX(); ++k) {
                    int z = i*setup_.GetLevelZeroBlocksY()*setup_.GetLevelZeroBlocksX()+j*setup_.GetLevelZeroBlocksX()+k; //Implicit type conversions, but do no harm
                    forest_[z].SetTargetRankForLeaf(ids_current_future_rank_map, HilbertPosition::z_x_y);
                }
            } else{
                for (short k = setup_.GetLevelZeroBlocksX()-1; k >=0 ; --k) {
                    int z=i*setup_.GetLevelZeroBlocksY()*setup_.GetLevelZeroBlocksX()+j*setup_.GetLevelZeroBlocksX()+k; //Implicit type conversions, but do no harm
                    forest_[z].SetTargetRankForLeaf(ids_current_future_rank_map, HilbertPosition::_z_xy);
                }
            }
            x_inverse = !x_inverse;
          }
      }
      y_inverse = !y_inverse;
  }

#endif

  // Sanity Check.
  for(unsigned int i=0; i <= setup_.GetMaximumLevel(); i++) {
      if(std::get<0>(ids_current_future_rank_map[i].back())!=0)
          throw std::logic_error("TopologyManager::AssignTargetRankToLeaves Some miscalculations in the distribution of leaves per rank");
  }
}

/**
 * @brief Assigns new tags to all nodes. Must be called after every refinement and coarsening in order to guarantee unique tags.
 */
void TopologyManager::GenerateTags() {
  int tag = 0;
  for(TopologyNode& tree : forest_) {
    tree.AssignTag(tag);
  }
}

/**
 * @brief Gives the position of the zero level ancestor node identified by the given id in the forest.
 * @param id The id of the node whose ancestor's position is to be determined
 * @return The index of the ancestor node in the zero topology. -1 If no such node could be found.
 */
int TopologyManager::PositionOfNodeInZeroTopology(const std::uint64_t id) const {

  std::uint64_t level_zero_id = id;
  while(Node::LevelOfNode(level_zero_id) != 0) {
    level_zero_id = Node::ParentIdOfNode(level_zero_id);
  }

  auto node_iterator = std::find_if(forest_.begin(),forest_.end(),[&level_zero_id](const TopologyNode& node){return node.Id() == level_zero_id;});

  if(node_iterator == forest_.end()) {
    return -1;
  }  else {
    return std::distance(forest_.begin(),node_iterator);
  }
}

/**
 * @brief Calculates a balanced distribution of nodes among the MPI ranks and assigns the determined rank to the nodes.
 *        Does not directly shift nodes among ranks! "Prepares for sending Load".
 */
void TopologyManager::AssignBalancedLoad() {
  for(TopologyNode& node : forest_) {
      node.BalanceTargetRanks();
  }
}

/**
 * @brief Gives a list of all nodes, that need to be balanced, i.e. shifted to another MPI rank.
 * @param ids_current_future_rank_map Indirect return parameter.
 * @note Lists the ranks to be balanced as tuple of their id, their current rank and the rank they are supposed to be shifted to
 */
void TopologyManager::ListNodeToBalance(std::vector<std::tuple<std::uint64_t, int, int>>& ids_current_future_rank_map) {
  for (TopologyNode& node : forest_) {
      node.ListUnbalancedNodes(rank_id_, ids_current_future_rank_map);
  }
}

/**
 * @brief Prints some statistics about how many messages were sent via MPI routines and the Distribution of Nodes (leaves) on the ranks etc.
 * @param logger The logger used for terminal output and log file creation.
 */
void TopologyManager::PrintStatistics(LogWriter& logger) {

  unsigned int maximum_level = setup_.GetMaximumLevel() + 1;
  if (MyRankId() == 0) {
    std::vector<std::vector<int>> levels(maximum_level,std::vector<int>(number_of_ranks_,0));

    for(unsigned int level = 0; level < maximum_level; level++) {
        std::vector<std::uint64_t> leaves = LeafIdsOnLevel(level);
        for (std::uint64_t id : leaves) {
            int rank = GetRankOfNode(id);
            levels[level][rank]++;
        }
    }

    //print
    logger.LogMessage("leave rank distribution: \n rank:\t",true,true);

    for (int rank = 0; rank < number_of_ranks_; rank++) {
        logger.LogMessage(std::to_string(rank) + "\t",true,true);
    }
    logger.LogMessage("\n",true,true);

    for(unsigned int level = 0; level < maximum_level; level++) {
        logger.LogMessage(std::to_string(level),true,true);
        for (int rank = 0; rank < number_of_ranks_; rank++) {
            logger.LogMessage("\t" + std::to_string(levels[level][rank]),true,true);
        }
        logger.LogMessage("\n",true,true);
    }
  }
}

/**
 * @brief Gives a list with the combined computation load (weight) on all levels.
 * @return List of summed weights on each level.
 */
std::vector<int> TopologyManager::WeightsOnLevels() const {
  std::vector<int> weights_on_level(setup_.GetMaximumLevel()+1); //If max level is = 2 we need 3 elements "0, 1, 2".
  for(const TopologyNode& node : forest_){
    node.ChildWeight(weights_on_level);
  }
  return weights_on_level;
}

/**
 * @brief Marks the given node with the given weight, which will be considered in the next load balancing.
 *        $CURRENTLY NOT USED, WILL BE IN MULITPHASE VERSION$
 * @param id The id of the node of interest.
 * @param weight The weight of the node of interest.
 * @note Does not directly assign the new node but buffers it in a list for one block update.
 */
void TopologyManager::MarkNodeWeight(const uint64_t id, const unsigned int weight){
  std::get<0>(local_weight_list_).push_back(id);
  std::get<1>(local_weight_list_).push_back(weight);
}

/**
 * @brief Gives the ids of all locally available nodes on level zero.
 * @return All ids of all nodes on level zero on the current mpi rank.
 * @note Function is used to initialize the Tree properly.
 */
std::vector<std::uint64_t> TopologyManager::LocalLevelZeroIds() const {
  std::vector<std::uint64_t> local_nodes;
  for(const TopologyNode& root_node : forest_) {
    if(root_node.GetRank() == rank_id_) {
        local_nodes.push_back(root_node.Id());
    }
  }
  return local_nodes;
}

/**
 * @brief Gives out the ids of all globally existent nodes on the specified level.
 * @param level Level of interest.
 * @return Ids of Nodes on level.
 */
std::vector<std::uint64_t> TopologyManager::GlobalIdsOnLevel(const unsigned int level) const {
  std::vector<std::uint64_t> ids;
  for(const TopologyNode& node : forest_) {
    node.IdsOnLevel(level,ids);
  }
  return ids;
}
