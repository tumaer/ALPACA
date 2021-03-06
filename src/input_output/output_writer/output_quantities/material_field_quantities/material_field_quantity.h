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
#ifndef MATERIAL_FIELD_QUANTITY_H
#define MATERIAL_FIELD_QUANTITY_H

#include "unit_handler.h"
#include "input_output/output_writer/output_quantity.h"
#include "input_output/output_writer/output_quantities/material_field_quantities/material_field_quantity_definitions.h"

/**
 * @brief The MaterialFieldOutputQuantity class handles the output to the filesystem in Xdmf+HDF5 file format for ParaView. MaterialFieldOutputQuantity must not change any data.
 *
 *        This quantity is used to write all material fields to the hdf5/xdmf file if desired. For activating different fields use the "output_constants.h" file. If a new
 *        material field needs to be added, refer to the "material_field_quantity_definitions.h". Here, nothing has to be changed.
 */
class MaterialFieldQuantity : public OutputQuantity {

private:
   // struct containing all data required for the output
   MaterialFieldQuantityData const quantity_data_;
   // Conservative buffer type
   ConservativeBufferType const buffer_type_;

   // Append functions required from base class
   void DoComputeCellData( Node const& node, std::vector<double>& cell_data, unsigned long long int& cell_data_counter ) const override;
   void DoComputeDebugCellData( Node const& node, std::vector<double>& cell_data, unsigned long long int& cell_data_counter, MaterialName const material ) const override;

public:
   MaterialFieldQuantity() = delete;
   explicit MaterialFieldQuantity( UnitHandler const& unit_handler,
                                   MaterialManager const& material_manager,
                                   std::string const& quantity_name,
                                   std::array<bool, 3> const output_flags,
                                   MaterialFieldQuantityData const& quantity_data,
                                   ConservativeBufferType const buffer_type = ConservativeBufferType::Average );
   virtual ~MaterialFieldQuantity()                      = default;
   MaterialFieldQuantity( MaterialFieldQuantity const& ) = delete;
   MaterialFieldQuantity& operator=( MaterialFieldQuantity const& ) = delete;
   MaterialFieldQuantity( MaterialFieldQuantity&& )                 = delete;
   MaterialFieldQuantity& operator=( MaterialFieldQuantity&& ) = delete;
};

#endif// MATERIAL_FIELD_QUANTITY_H
