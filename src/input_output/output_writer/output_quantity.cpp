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
#include "output_quantity.h"
#include "input_output/utilities/xdmf_utilities.h"

/**
 * @brief Explicit constructor to be used to create an output quantity.
 * @param unit_handler Handler to allow dimensionalization.
 * @param material_manager Manager for handling of the materials.
 * @param output_flags Flags specifying for which output type the quantity should be written.
 * @param quantity_name Name of the given quantity used in the hdf5 and xdmf file.
 * @param dimensions Dimensions of the given quantity ( {1,1} : scalar, {3,1}: vector, {n,m}: matrix/tensor).
 */
OutputQuantity::OutputQuantity( UnitHandler const& unit_handler,
                                MaterialManager const& material_manager,
                                std::string const& quantity_name,
                                std::array<bool, 3> const& output_flags,
                                std::array<unsigned int, 2> const& dimensions ) :// Start initializer list
                                                                                  unit_handler_( unit_handler ),
                                                                                  material_manager_( material_manager ),
                                                                                  output_flags_( output_flags ),
                                                                                  quantity_name_( quantity_name ),
                                                                                  dimensions_( dimensions ) {
   /** Empty constructor besides initializer list */
}

/**
 * @brief Gives the name of the quantity used in the xdmf and hdf5 file.
 * @return Name of the quantity.
 */
std::string OutputQuantity::GetName() const {
   return quantity_name_;
}

/**
 * @brief Gives the dimensions of the quantity.
 * @return DImensions of the quantity.
 */
std::array<unsigned int, 2> OutputQuantity::GetDimensions() const {
   return dimensions_;
}

/**
 * @brief Checks whether the output quantity should be written for a given output type (0: standard, 1: interface, 2: debug).
 * @param output_type Output type to be checked.
 * @return true if output should be written, false if not.
 */
bool OutputQuantity::IsActive( OutputType const output_type ) const {
   return output_flags_[OTTI( output_type )];
}

/**
 * @brief Gives the appropriate attribute string for the xdmf file for this quantity (differentiation between multidimensional, vectorial and scalar quantities).
 * @param hdf5_filename HDF5 filename (without path) where the actual data has been written to.
 * @param group_name Name of the group the data was written.
 * @param number_of_global_cells Number of cells used for the complete quantity (globally on all ranks).
 * @param prefix Additional prefix that can be used for the attribute_name.
 * @return Compete Xdmf attribute string.
 */
std::string OutputQuantity::GetXdmfAttributeString( std::string const& hdf5_filename, std::string const& group_name, hsize_t const number_of_global_cells, std::string const prefix ) const {
   // Get the data item
   std::string const data_item( XdmfUtilities::DataItemString( hdf5_filename, group_name + "/" + prefix + quantity_name_, number_of_global_cells, dimensions_ ) );
   if( dimensions_.back() > 1 || dimensions_.front() > DTI( CC::DIM() ) ) {// multidimensional (second component dimension larger than one or first larger than the current dimension)
      // Both components are equal (nxn) and smaller or equal the dimension -> Assumption that it is a tensor
      // This is not relevant for the reading, but allows sometimes special computations in ParaView
      if( dimensions_.back() == dimensions_.front() && dimensions_.back() <= DTI( CC::DIM() ) ) {
         return XdmfUtilities::TensorAttributeString( prefix + quantity_name_, data_item );
      } else {
         return XdmfUtilities::MatrixAttributeString( prefix + quantity_name_, data_item );
      }
   } else if( dimensions_.front() > 1 ) {// vectorial (first component dimension larger than one and implicitly smaller than current dimension)
      return XdmfUtilities::VectorAttributeString( prefix + quantity_name_, data_item );
   } else {// scalar (both dimensions are one)
      return XdmfUtilities::ScalarAttributeString( prefix + quantity_name_, data_item );
   }
}

/**
 * @brief Computes the data of all internal cells and writes them into the output data vector.
 * @param nodes Vector of nodes for which the output should be done.
 * @param cell_data Vector in which the full set of cell data is written.
 *
 * @note The cell_data vector stores its values in a pre-defined order. For a n x m output quantity with and 3 cells the order is a follows:
 *                    cell 1                   |              cell 2                   |              cell 3                   |
 *       00 01 ... 0m 10 ... 1m ... n0 ... nm  | 00 01 ... 0m 10 ... 1m ... n0 ... nm  | 00 01 ... 0m 10 ... 1m ... n0 ... nm  |
 */
void OutputQuantity::ComputeCellData( std::vector<std::reference_wrapper<Node const>> const& nodes, std::vector<double>& cell_data ) const {

   // Counter for the cell data vector
   unsigned long long int cell_data_counter = 0;
   // node counter
   unsigned int node_counter = 0;
   // factor that specifies the number of cell data entries that must be/are written per node
   unsigned long long int cell_data_per_node = dimensions_[0] * dimensions_[1] * CC::ICX() * CC::ICY() * CC::ICZ();

   // Loop through all nodes to append correct number of data
   for( Node const& node : nodes ) {
      // Ensures that the data counter is at the correct position for this node
      cell_data_counter = cell_data_per_node * node_counter;
      // Call the function that is implemented by the derived class
      DoComputeCellData( node, cell_data, cell_data_counter );
      // Increment the node counter
      node_counter++;
   }
}

/**
 * @brief Computes the data of all cells (internal + halos) and writes them into the output data vector.
 * @param nodes Vector of nodes for which the output should be done.
 * @param cell_data Vector in which the full set of cell data is written.
 * @param material Material used for the output (not always used, only for material output quantities).
 *
 * @note The cell_data vector stores its values in a pre-defined order. For a n x m output quantity with and 3 cells the order is a follows:
 *                    cell 1                   |              cell 2                   |              cell 3                   |
 *       00 01 ... 0m 10 ... 1m ... n0 ... nm  | 00 01 ... 0m 10 ... 1m ... n0 ... nm  | 00 01 ... 0m 10 ... 1m ... n0 ... nm  |
 */
void OutputQuantity::ComputeDebugCellData( std::vector<std::reference_wrapper<Node const>> const& nodes, std::vector<double>& cell_data, MaterialName const material ) const {

   // Counter for the cell data vector
   unsigned long long int cell_data_counter = 0;
   // node counter
   unsigned int node_counter = 0;
   // factor that specifies the number of cell data entries that must be/are written per node
   unsigned long long int cell_data_per_node = dimensions_[0] * dimensions_[1] * CC::TCX() * CC::TCY() * CC::TCZ();

   // Loop through all nodes to append correct number of data
   for( Node const& node : nodes ) {
      // Ensures that the data counter is at the correct position for this node
      cell_data_counter = cell_data_per_node * node_counter;
      // Call the function that is implemented by the derived class
      DoComputeDebugCellData( node, cell_data, cell_data_counter, material );
      // Increment the node counter
      node_counter++;
   }
}
