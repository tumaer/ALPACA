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
#ifndef OUTPUT_QUANTITY_H
#define OUTPUT_QUANTITY_H

#include <hdf5.h>
#include "topology/node.h"
#include "unit_handler.h"
#include "materials/material_definitions.h"
#include "materials/material_manager.h"
#include "input_output/output_writer/output_definitions.h"

/**
 * @brief The OutputQuantity class handles the output to the filesystem in Xdmf+HDF5 file format for ParaView. OutputQuantity must not change any data.
 * @note This class works as a base class where all field output quantities must inherit from that write cell data to the appropriate file.
 */
class OutputQuantity {

protected:  
   // unit_handler for the dimensionalization of variables 
   UnitHandler const& unit_handler_;
   // material manager to allow computations for different materials and use properties/eos data or pairing data 
   MaterialManager const& material_manager_;
   // flags to specify where the quantities should be used (0: standard, 1: interface, 2:debug)
   std::array<bool, 3> const output_flags_;
   // name of the complete quantity used in the hdf5 + xdmf file for naming 
   std::string const quantity_name_;
   // dimensions of the quantity ( { 1 } : Scalar, { n } : vector, { n x m } : matrix )
   std::array<unsigned int,2> const dimensions_;

   /**
    * @brief Compute values for the data vector that is written to the hdf5 file (standard, interface mode)
    * @param node Node data that should be added
    * @param cell_data_counter Actual position of the index in the cell data vector
    * @param cell_data Vector in which the data is written 
    * 
    * @note pure virtual function that MUST be implemented by all derived classes
    */ 
   virtual void DoComputeCellData( Node const& node, std::vector<double>& cell_data, unsigned long long int & cell_data_counter ) const = 0;
   
   /**
    * @brief Compute values for the data vector that is written to the hdf5 file (debug mode)
    * @param node Node data that should be added
    * @param cell_data_counter Actual position of the index in the cell data vector
    * @param cell_data Vector in which the data is written 
    * @param material Material that should be considered (only used for material output quantities)
    * 
    * @note pure virtual function that MUST be implemented by all derived classes
    */ 
   virtual void DoComputeDebugCellData( Node const& node, std::vector<double>& cell_data, unsigned long long int & cell_data_counter, MaterialName const material ) const = 0;

   // protected constructor to prevent direct calls from outside
   explicit OutputQuantity( UnitHandler const& unit_handler, 
                            MaterialManager const& material_manager, 
                            std::string const& quantity_name,
                            std::array<bool, 3> const& output_flags, 
                            std::array<unsigned int,2> const& dimensions );

public:
   OutputQuantity() = delete;
   virtual ~OutputQuantity() = default;
   OutputQuantity( OutputQuantity const& ) = delete;
   OutputQuantity& operator=( OutputQuantity const& ) = delete;
   OutputQuantity( OutputQuantity&& ) = delete;
   OutputQuantity& operator=( OutputQuantity&& ) = delete;

   // Functions that fill the complete cell data vector with appropriate values
   void ComputeCellData( std::vector<std::reference_wrapper<Node const>> const& nodes, std::vector<double>& cell_data ) const;
   void ComputeDebugCellData( std::vector<std::reference_wrapper<Node const>> const& nodes, std::vector<double>& cell_data, MaterialName const material = MaterialName::MaterialOne ) const;

   // Creates and attribute string compliant with a xdmf reader
   std::string GetXdmfAttributeString( std::string const& filename, std::string const& group_name, hsize_t const number_of_values, std::string const prefix = "" ) const;

   // Additional return functions to provide data to outside 
   bool IsActive( OutputType const output_type ) const;
   std::string GetName() const;
   std::array<unsigned int, 2> GetDimensions() const;
};

#endif // OUTPUT_QUANTITY_H
