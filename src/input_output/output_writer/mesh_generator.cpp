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
#include "mesh_generator.h"
#include "input_output/output_writer/xdmf_utilities.h"

/**
 * @brief Constructor for a generic mesh generator to be called from derived classes 
 * @param topology_manager Manager storing all topology information of all ranks 
 * @param flower The tree for the node informations on the current rank 
 * @param node_size_on_level_zero Already dimensionalized size of a single node on level zero 
 */ 
MeshGenerator::MeshGenerator( TopologyManager const& topology_manager, Tree const& flower, double const dimensionalized_node_size_on_level_zero ) :  
   topology_( topology_manager ),
   tree_( flower ),
   dimensionalized_node_size_on_level_zero_( dimensionalized_node_size_on_level_zero ) {
      /** Empty constructor besides initializer list */
}

/** 
 * @brief Returns the string for the topology (vertex ids of cells)
 * @param filename Name of the .hdf5 file where the actual data is found 
 * @param group_name Name of the group the topology is written into
 * @return Attribute string for the topology
 */ 
std::string MeshGenerator::GetXdmfTopologyString( std::string const& filename, std::string const& group_name ) const {
   hsize_t const global_number_cells = GetGlobalNumberOfCells();
   std::string const data_item("<DataItem NumberType=\"Int\" Format=\"HDF\" Dimensions=\"" + std::to_string( global_number_cells ) + " 8\"> " + filename + ":/" + group_name + "/" + vertex_ids_name_ + " </DataItem>\n");
   return XdmfUtilities::TopologyString( data_item, global_number_cells );
}

/** 
 * @brief Returns the string used for the geometry (vertex coordinates)
 * @param filename Name of the .hdf5 file where the actual data is found 
 * @param group_name Name of the group the geometry is written into
 * @return Attribute string for the geometry
 */ 
std::string MeshGenerator::GetXdmfGeometryString( std::string const& filename, std::string const& group_name ) const {
   hsize_t const global_number_vertices = GetGlobalDimensionsOfVertexCoordinates().front();
   std::string const data_item("<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"" + std::to_string( global_number_vertices ) + " 3\"> " + filename + ":/" + group_name + "/" + vertex_coordinates_name_ + " </DataItem>\n");
   return XdmfUtilities::GeometryString( data_item, global_number_vertices );
}

/**
 * @brief Appends the vertex IDs for the specific mesh generator to the hdf5 file. The vertex IDs represent 
 *        the cells, where one cell is built by 8 vertices
 * @param vertex_ids Vector where all vertex IDs are written into (indirect return)
 */  
void MeshGenerator::ComputeVertexIDs( std::vector<unsigned long long int> & vertex_ids ) const {
   // Compute the vertex IDs in the derived class
   DoComputeVertexIDs( vertex_ids );
}

/**
 * @brief Appends the vertex coordinates for the specific mesh generator to the hdf5 file 
 * @param vertex_coordinates Vector where all vertex coordinates are written into (indirect return)
 */     
void MeshGenerator::ComputeVertexCoordinates( std::vector<double> & vertex_coordinates ) const {
   // Compute the vertex coordinates in the derived class
   DoComputeVertexCoordinates( vertex_coordinates );
}