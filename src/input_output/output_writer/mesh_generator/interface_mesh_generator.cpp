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
#include "input_output/output_writer/mesh_generator/interface_mesh_generator.h"
#include "input_output/output_writer/mesh_generator/mesh_generator_utilities.h"
#include "communication/mpi_utilities.h"
#include "topology/id_information.h"

/**
 * @brief Constructor to create the interface mesh generator.
 * @param unit_handler Instance for dimensionalization of variables.
 * @param topology_manager Instance to provide node information on different ranks.
 * @param flower Instance to provide node information of current rank.
 * @param dimensionalized_node_size_on_level_zero Already dimensionalized size of a node on level zero.
 */
InterfaceMeshGenerator::InterfaceMeshGenerator( TopologyManager const& topology_manager,
                                                Tree const& flower,
                                                double const dimensionalized_node_size_on_level_zero ) : MeshGenerator( topology_manager, flower, dimensionalized_node_size_on_level_zero ) {
   /** Empty besides call of base class constructor */
}

/**
 * @brief See base class implementation.
 */
std::vector<std::reference_wrapper<Node const>> InterfaceMeshGenerator::DoGetLocalNodes() const {
   return tree_.InterfaceLeaves();
}

/**
 * @brief See base class implementation.
 */
hsize_t InterfaceMeshGenerator::DoGetGlobalNumberOfCells() const {
   return hsize_t( topology_.InterfaceLeafCount() ) * MeshGeneratorUtilities::NumberOfInternalCellsPerBlock();
}

/**
 * @brief See base class implementation.
 */
hsize_t InterfaceMeshGenerator::DoGetLocalNumberOfCells() const {
   return hsize_t( topology_.LocalInterfaceLeafIds().size() ) * MeshGeneratorUtilities::NumberOfInternalCellsPerBlock();
}

/**
 * @brief See base class implementation.
 */
hsize_t InterfaceMeshGenerator::DoGetLocalCellsStartIndex() const {
   return hsize_t( topology_.InterfaceLeafOffsetOfRank( MpiUtilities::MyRankId() ) ) * MeshGeneratorUtilities::NumberOfInternalCellsPerBlock();
}

/**
 * @brief See base class implementation.
 */
std::vector<hsize_t> InterfaceMeshGenerator::DoGetGlobalDimensionsOfVertexCoordinates() const {
   return { hsize_t( topology_.InterfaceLeafCount() ) * MeshGeneratorUtilities::NumberOfInternalVerticesPerBlock(), hsize_t( 3 ) };
}

/**
 * @brief See base class implementation.
 */
std::vector<hsize_t> InterfaceMeshGenerator::DoGetLocalDimensionsOfVertexCoordinates() const {
   return { hsize_t( topology_.LocalInterfaceLeafIds().size() ) * MeshGeneratorUtilities::NumberOfInternalVerticesPerBlock(), hsize_t( 3 ) };
}

/**
 * @brief See base class implementation.
 */
hsize_t InterfaceMeshGenerator::DoGetLocalVertexCoordinatesStartIndex() const {
   return hsize_t( topology_.InterfaceLeafOffsetOfRank( MpiUtilities::MyRankId() ) ) * MeshGeneratorUtilities::NumberOfInternalVerticesPerBlock();
}

/**
 * @brief See base class definition.
 */
void InterfaceMeshGenerator::DoComputeVertexCoordinates( std::vector<double>& vertex_coordinates ) const {

   // get the correct number of interface leaves for the rank
   std::vector<nid_t> local_interface_leaf_ids = topology_.LocalInterfaceLeafIds();
   // resize the vector to ensure enough memory for the cooridnates ( x,y,z coordinates for each vertex )
   vertex_coordinates.resize( local_interface_leaf_ids.size() * MeshGeneratorUtilities::NumberOfInternalVerticesPerBlock() * 3 );

   // Compute the correct coordinates of the vertices
   std::size_t vertex_counter = 0;
   for( auto const& id : local_interface_leaf_ids ) {
      double const block_size            = DomainSizeOfId( id, dimensionalized_node_size_on_level_zero_ );
      double const cell_size             = MeshGeneratorUtilities::CellSizeForBlockSize( block_size );
      std::array<double, 3> block_origin = DomainCoordinatesOfId( id, block_size );
      for( unsigned int k = 0; k <= CC::ICZ(); ++k ) {
         for( unsigned int j = 0; j <= CC::ICY(); ++j ) {
            for( unsigned int i = 0; i <= CC::ICX(); ++i ) {
               vertex_coordinates[vertex_counter]     = block_origin[0] + double( i ) * cell_size;
               vertex_coordinates[vertex_counter + 1] = block_origin[1] + double( j ) * cell_size;
               vertex_coordinates[vertex_counter + 2] = block_origin[2] + double( k ) * cell_size;
               vertex_counter += 3;
            }
         }
      }
   }
}

/**
 * @brief See base class definition.
 */
void InterfaceMeshGenerator::DoComputeVertexIDs( std::vector<unsigned long long int>& vertex_ids ) const {
   /************************************************************************/
   /** 1. Create full set of vertices */
   // Local leave definitions
   std::size_t number_of_local_interface_leaves = topology_.LocalInterfaceLeafIds().size();
   // Global offset between ranks (global vector is filled in the order (rank0, rank1, ..., rankN))
   unsigned int const offset = topology_.InterfaceLeafOffsetOfRank( MpiUtilities::MyRankId() );
   // Resize Vertex ID vector ( 8 vertices span one cell )
   vertex_ids.resize( number_of_local_interface_leaves * MeshGeneratorUtilities::NumberOfInternalCellsPerBlock() * 8 );
   unsigned long long int vertex_id_counter = 0;
   /**
    * Definition of offset parameters in y and z-direction to be specified for loop computation of vertices
    * Vertices are described in increasing order:
    * v( x,y,z ) = ( 0,0,0 )               -> ID: 0
    * v( x,y,z ) = ( total,0,0 )           -> ID: max
    * v( x,y,z ) = ( 0,1,0 )               -> ID: max + 1
    * v( x,y,z ) = ( 0,0,1 )               -> ID: ( max +1 ) * ( max + 1 )
    * v( x,y,z ) = ( total, total, total ) -> ID: ( max + 1 ) * ( max + 1 ) * ( max + 1 )
    */
   constexpr unsigned int j_ids_skew = ( CC::ICX() + 1 );
   constexpr unsigned int k_ids_skew = ( CC::ICX() + 1 ) * ( CC::ICY() + 1 );

   // loop through all number of leaves
   for( unsigned int leaves_counter = 0; leaves_counter < number_of_local_interface_leaves; leaves_counter++ ) {
      //store cell_vetrex offset for each node
      for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
         for( unsigned int j = 0; j < CC::ICY(); ++j ) {
            for( unsigned int i = 0; i < CC::ICX(); ++i ) {
               // Shift index for correct indexing in the global mesh (no duplicated IDs)
               unsigned long long int const shift = ( leaves_counter + offset ) * MeshGeneratorUtilities::NumberOfInternalVerticesPerBlock();
               // Add all vertices for one cell
               vertex_ids[vertex_id_counter]     = i + j * j_ids_skew + k * k_ids_skew + shift;
               vertex_ids[vertex_id_counter + 1] = ( i + 1 ) + j * j_ids_skew + k * k_ids_skew + shift;
               vertex_ids[vertex_id_counter + 2] = ( i + 1 ) + ( j + 1 ) * j_ids_skew + k * k_ids_skew + shift;
               vertex_ids[vertex_id_counter + 3] = i + ( j + 1 ) * j_ids_skew + k * k_ids_skew + shift;
               vertex_ids[vertex_id_counter + 4] = i + j * j_ids_skew + ( k + 1 ) * k_ids_skew + shift;
               vertex_ids[vertex_id_counter + 5] = ( i + 1 ) + j * j_ids_skew + ( k + 1 ) * k_ids_skew + shift;
               vertex_ids[vertex_id_counter + 6] = ( i + 1 ) + ( j + 1 ) * j_ids_skew + ( k + 1 ) * k_ids_skew + shift;
               vertex_ids[vertex_id_counter + 7] = i + ( j + 1 ) * j_ids_skew + ( k + 1 ) * k_ids_skew + shift;
               vertex_id_counter += 8;
            }
         }
      }
   }

   /************************************************************************/
   /** 2. Vertex filtering */
   // Not implemented yet -> There are doubled vertices in the output
}
