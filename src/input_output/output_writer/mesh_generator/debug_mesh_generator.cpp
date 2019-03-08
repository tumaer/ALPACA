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
#include "input_output/output_writer/mesh_generator/debug_mesh_generator.h"
#include "input_output/output_writer/mesh_generator/mesh_generator_utilities.h"
#include "communication/mpi_utilities.h"
#include "topology/id_information.h"

/**
 * @brief Constructor to create the standard finest level mesh generator
 * @param topology_manager Instance to provide node information on different ranks 
 * @param flower Instance to provide node information of current rank 
 * @param dimensionalized_node_size_on_level_zero Already dimensionalized size of a node on level zero 
 * @param number_of_z_nodes_on_level_zero Number of blocks on level zero in z-direction
 */
DebugMeshGenerator::DebugMeshGenerator( TopologyManager const& topology_manager, 
                                        Tree const& flower, 
                                        double const dimensionalized_node_size_on_level_zero,
                                        unsigned int const number_of_z_nodes_on_level_zero ) :
   MeshGenerator( topology_manager, flower, dimensionalized_node_size_on_level_zero ), 
   z_coordinates_offset_( dimensionalized_node_size_on_level_zero_ * ( number_of_z_nodes_on_level_zero + 1 ) ) {
   /** Empty besides call of base class constructor */
}

/**
 * @brief See base class implementation 
 */
std::vector<std::reference_wrapper<Node const>> DebugMeshGenerator::DoGetLocalNodes() const {
   return tree_.AllNodes();
}

/**
 * @brief See base class implementation
 */ 
hsize_t DebugMeshGenerator::DoGetGlobalNumberOfCells() const {
   return hsize_t( std::get<0>( topology_.NodeAndLeafCount() ) ) * MeshGeneratorUtilities::NumberOfTotalCellsPerBlock(); 
}

/**
 * @brief See base class implementation
 */ 
hsize_t DebugMeshGenerator::DoGetLocalNumberOfCells() const {
   return hsize_t( topology_.LocalNodeIds().size() ) * MeshGeneratorUtilities::NumberOfTotalCellsPerBlock();
}

/**
 * @brief See base class implementation
 */ 
hsize_t DebugMeshGenerator::DoGetLocalCellsStartIndex() const {
   return hsize_t( topology_.NodeOffsetOfRank( MpiUtilities::MyRankId() ) ) * MeshGeneratorUtilities::NumberOfTotalCellsPerBlock();
}

/**
 * @brief See base class implementation
 */ 
std::vector<hsize_t> DebugMeshGenerator::DoGetGlobalDimensionsOfVertexCoordinates() const {
   return { hsize_t( std::get<0>( topology_.NodeAndLeafCount() ) ) * MeshGeneratorUtilities::NumberOfTotalVerticesPerBlock(), hsize_t ( 3 ) };
}

/**
 * @brief See base class implementation
 */ 
std::vector<hsize_t> DebugMeshGenerator::DoGetLocalDimensionsOfVertexCoordinates() const {
   return { hsize_t( topology_.LocalNodeIds().size() ) * MeshGeneratorUtilities::NumberOfTotalVerticesPerBlock(), hsize_t ( 3 ) };
}

/**
 * @brief See base class implementation
 */ 
hsize_t DebugMeshGenerator::DoGetLocalVertexCoordinatesStartIndex() const {
   return hsize_t( topology_.NodeOffsetOfRank( MpiUtilities::MyRankId() ) ) * MeshGeneratorUtilities::NumberOfTotalVerticesPerBlock();
}

/**
 * @brief See base class definition
 */ 
void DebugMeshGenerator::DoComputeVertexCoordinates( std::vector<double> & vertex_coordinates ) const {

   // get the correct number of leaves fo the rank 
   std::vector<std::uint64_t> local_node_ids = topology_.LocalNodeIds();
   // resize the vector to ensure enough memory for the cooridnates ( x,y,z coordinates for each vertex )
   vertex_coordinates.resize( local_node_ids.size() * MeshGeneratorUtilities::NumberOfTotalVerticesPerBlock() * 3 );

   // Loop over all local leaves 
   std::size_t vertex_coordinates_counter = 0;
   for( auto const& id : local_node_ids ) {
      // Stores the current level of the node 
      unsigned int const level_of_node = LevelOfNode( id );
      /**
       * To obtain a gap between different levels add an additional block in z-direction for each level. This allows 
       * a representation of the full set of nodes and levels in a single ParaView file.
       */
      double const level_z_offset = level_of_node * z_coordinates_offset_;
      /**
       * Determine the cell size for the output. Since the debug output only serves as a helping tool, the real coordinates are 
       * not important. Therefore, the internal block size is divided in equal spaced cells. To obtain an empty gap between to blocks 
       * an additional cell ("+1") is considered. 
       */
      // Get all coordinate information for the current node ( cell size and origin )
      double const block_size = DomainSizeOfId( id, dimensionalized_node_size_on_level_zero_ );
      double const cell_size = block_size / double( CC::TCX() + 1 );
      std::array<double, 3> const block_origin = DomainCoordinatesOfId( id, block_size );
      // Loop through all total cells to append coordinates 
      for( unsigned int k = 0; k <= CC::TCZ(); ++k ) {
         for( unsigned int j = 0; j <= CC::TCY(); ++j ) {
            for( unsigned int i = 0; i <= CC::TCX(); ++i ) {               
               vertex_coordinates[vertex_coordinates_counter    ] = block_origin[0] + double( i ) * cell_size;
               vertex_coordinates[vertex_coordinates_counter + 1] = block_origin[1] + double( j ) * cell_size;
               vertex_coordinates[vertex_coordinates_counter + 2] = block_origin[2] + double( k ) * cell_size + level_z_offset;
               vertex_coordinates_counter += 3;
            }
         }
      }
   }
}


/**
 * @brief See base class definition
 */ 
void DebugMeshGenerator::DoComputeVertexIDs( std::vector<unsigned long long int> & vertex_ids ) const {
   // Local leave definitions 
   std::size_t const number_of_local_nodes = topology_.LocalNodeIds().size();
   // Global offset between ranks (global vector is filled in the order (rank0, rank1, ..., rankN))
   unsigned int const offset = topology_.NodeOffsetOfRank( MpiUtilities::MyRankId() );
   // Resize Vertex ID vector ( 8 vertices span one cell )
   vertex_ids.resize( number_of_local_nodes * MeshGeneratorUtilities::NumberOfTotalCellsPerBlock() * 8 );
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
   constexpr unsigned int j_ids_skew = ( CC::TCX() + 1 );
   constexpr unsigned int k_ids_skew = ( CC::TCX() + 1 ) * ( CC::TCY() + 1 );

   // loop through all number of nodes and cells  
   for( unsigned int nodes_counter = 0; nodes_counter < number_of_local_nodes; nodes_counter++ ) {
      for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int i = 0; i < CC::TCX(); ++i ) {
               // Shift index for correct indexing in the global mesh (no duplicated IDs)
               unsigned long long int const shift = ( nodes_counter + offset ) * MeshGeneratorUtilities::NumberOfTotalVerticesPerBlock();
               
               // Add all vertices for one cell 
               vertex_ids[vertex_id_counter    ] =    i    +    j    * j_ids_skew +    k    * k_ids_skew + shift;
               vertex_ids[vertex_id_counter + 1] = ( i+1 ) +    j    * j_ids_skew +    k    * k_ids_skew + shift;
               vertex_ids[vertex_id_counter + 2] = ( i+1 ) + ( j+1 ) * j_ids_skew +    k    * k_ids_skew + shift;
               vertex_ids[vertex_id_counter + 3] =    i    + ( j+1 ) * j_ids_skew +    k    * k_ids_skew + shift;
               vertex_ids[vertex_id_counter + 4] =    i    +    j    * j_ids_skew + ( k+1 ) * k_ids_skew + shift;
               vertex_ids[vertex_id_counter + 5] = ( i+1 ) +    j    * j_ids_skew + ( k+1 ) * k_ids_skew + shift;
               vertex_ids[vertex_id_counter + 6] = ( i+1 ) + ( j+1 ) * j_ids_skew + ( k+1 ) * k_ids_skew + shift;
               vertex_ids[vertex_id_counter + 7] =    i    + ( j+1 ) * j_ids_skew + ( k+1 ) * k_ids_skew + shift;
               vertex_id_counter += 8;
            }
         }
     }
   }

   /**
    * @note In other mesh_generators filtering operations are required to remove double placed vertices or to assign the 
    *       correct ID to double places vertices. In the debug Output this is not required, since an additional gap between 
    *       two blocks is used and the halo cells are also written. 
    */
}
