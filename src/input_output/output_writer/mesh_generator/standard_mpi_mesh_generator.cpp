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
#include "input_output/output_writer/mesh_generator/standard_mpi_mesh_generator.h"
#include "input_output/output_writer/mesh_generator/mpi_vertex_filter_utilities.h"
#include "input_output/output_writer/mesh_generator/mesh_generator_utilities.h"

/**
 * @brief Constructor to create the standard mpi mesh generator.
 * @param topology_manager Instance to provide node information on different ranks.
 * @param flower Instance to provide node information of current rank.
 * @param dimensionalized_node_size_on_level_zero Already dimensionalized size of a node on level zero.
 * @param is_mpi_filtering_active Flag whether mpi filtering should be used or not.
 */
StandardMpiMeshGenerator::StandardMpiMeshGenerator( TopologyManager const& topology_manager,
                                                    Tree const& flower,
                                                    double const dimensionalized_node_size_on_level_zero,
                                                    bool const mpi_filtering_active ) : MeshGenerator( topology_manager, flower, dimensionalized_node_size_on_level_zero ),
                                                                                        mpi_filtering_active_( mpi_filtering_active ) {
   /** Empty besides call of base class constructor */
}

/**
 * @brief See base class implementation.
 */
std::vector<std::reference_wrapper<Node const>> StandardMpiMeshGenerator::DoGetLocalNodes() const {
   return tree_.Leaves();
}

/**
 * @brief See base class implementation.
 */
hsize_t StandardMpiMeshGenerator::DoGetGlobalNumberOfCells() const {
   return hsize_t( std::get<1>( topology_.NodeAndLeafCount() ) ) * MeshGeneratorUtilities::NumberOfInternalCellsPerBlock();
}

/**
 * @brief See base class implementation.
 */
hsize_t StandardMpiMeshGenerator::DoGetLocalNumberOfCells() const {
   return hsize_t( topology_.LocalLeafIds().size() ) * MeshGeneratorUtilities::NumberOfInternalCellsPerBlock();
}

/**
 * @brief See base class implementation.
 */
hsize_t StandardMpiMeshGenerator::DoGetLocalCellsStartIndex() const {
   return hsize_t( topology_.LeafOffsetOfRank( MpiUtilities::MyRankId() ) ) * MeshGeneratorUtilities::NumberOfInternalCellsPerBlock();
}

/**
 * @brief See base class implementation.
 */
std::vector<hsize_t> StandardMpiMeshGenerator::DoGetGlobalDimensionsOfVertexCoordinates() const {
   return { hsize_t( std::get<1>( topology_.NodeAndLeafCount() ) ) * MeshGeneratorUtilities::NumberOfInternalVerticesPerBlock(), hsize_t( 3 ) };
}

/**
 * @brief See base class implementation.
 */
std::vector<hsize_t> StandardMpiMeshGenerator::DoGetLocalDimensionsOfVertexCoordinates() const {
   return { hsize_t( topology_.LocalLeafIds().size() ) * MeshGeneratorUtilities::NumberOfInternalVerticesPerBlock(), hsize_t( 3 ) };
}

/**
 * @brief See base class implementation.
 */
hsize_t StandardMpiMeshGenerator::DoGetLocalVertexCoordinatesStartIndex() const {
   return hsize_t( topology_.LeafOffsetOfRank( MpiUtilities::MyRankId() ) ) * MeshGeneratorUtilities::NumberOfInternalVerticesPerBlock();
}

/**
 * @brief See base class definition.
 *
 * @note The coordinates are created regardless they are needed at the end. If the mpi routine is used for filtering, some of the
 *       coordinates would not be needed anymore, but they still persist in the final vector to avoid complex reassignments.
 */
void StandardMpiMeshGenerator::DoComputeVertexCoordinates( std::vector<double>& vertex_coordinates ) const {
   // get the correct number of leaves fo the rank
   std::vector<nid_t> local_leaf_ids = topology_.LocalLeafIds();
   // resize the vector to ensure enough memory for the cooridnates ( x,y,z coordinates for each vertex )
   vertex_coordinates.resize( local_leaf_ids.size() * MeshGeneratorUtilities::NumberOfInternalVerticesPerBlock() * 3 );

   // Loop over all local leaves
   std::size_t vertex_coordinates_counter = 0;
   for( auto const& id : local_leaf_ids ) {
      // Get all information fo the current leave ( cell size and origin )
      double const block_size            = DomainSizeOfId( id, dimensionalized_node_size_on_level_zero_ );
      double const cell_size             = MeshGeneratorUtilities::CellSizeForBlockSize( block_size );
      std::array<double, 3> block_origin = DomainCoordinatesOfId( id, block_size );
      // Loop through all internal cells to append coordinates
      for( unsigned int k = 0; k <= CC::ICZ(); ++k ) {
         for( unsigned int j = 0; j <= CC::ICY(); ++j ) {
            for( unsigned int i = 0; i <= CC::ICX(); ++i ) {
               vertex_coordinates[vertex_coordinates_counter]     = block_origin[0] + double( i ) * cell_size;
               vertex_coordinates[vertex_coordinates_counter + 1] = block_origin[1] + double( j ) * cell_size;
               vertex_coordinates[vertex_coordinates_counter + 2] = block_origin[2] + double( k ) * cell_size;
               vertex_coordinates_counter += 3;
            }
         }
      }
   }
}

/**
 * @brief See base class definition.
 */
void StandardMpiMeshGenerator::DoComputeVertexIDs( std::vector<unsigned long long int>& vertex_ids ) const {
   /************************************************************************/
   /** 1. Create full set of vertices */
   // Local leave definitions
   std::size_t number_of_local_leaves = topology_.LocalLeafIds().size();
   // Global offset between ranks (global vector is filled in the order (rank0, rank1, ..., rankN))
   unsigned int const offset = topology_.LeafOffsetOfRank( MpiUtilities::MyRankId() );
   // Resize Vertex ID vector ( 8 vertices span one cell )
   vertex_ids.resize( number_of_local_leaves * MeshGeneratorUtilities::NumberOfInternalCellsPerBlock() * 8 );
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
   //Store offset in cell_vertex_ids for each leave; needed for vertex filtering
   std::vector<unsigned long long int> leave_offset;
   leave_offset.reserve( number_of_local_leaves );

   // loop through all number of leaves
   for( unsigned int leaves_counter = 0; leaves_counter < number_of_local_leaves; leaves_counter++ ) {
      //store cell_vetrex offset for each leave
      leave_offset.push_back( vertex_id_counter );
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
   // Carry out MPI sending and receiving operations to remove duplicated vertex ids of two neighboring cells
   // separation of x,y and z-axis
   // Only done if the vertex filter is set properly. Otherwise double created vertices are present in the grid
   if( mpi_filtering_active_ ) {
      FilterVertexIDs( vertex_ids, leave_offset );
   }
}

/**
 * @brief Filters double placed vertex IDs by an mpi routine
 * @param vertex_ids Vertex IDs to be filtered (indirect return)
 * @param leave_offset Offset index of each leaf node in the vertex ID vector (where they start)
 */
void StandardMpiMeshGenerator::FilterVertexIDs( std::vector<unsigned long long int>& vertex_ids, std::vector<unsigned long long int> const& leave_offset ) const {

   //create list of send vectors to each other rank; this is introduced to deal with the creation of tags
   std::vector<std::vector<std::pair<nid_t, nid_t>>> send_list_per_rank( MpiUtilities::NumberOfRanks() );
   std::vector<std::vector<std::pair<nid_t, nid_t>>> recv_list_per_rank( MpiUtilities::NumberOfRanks() );
   std::vector<MPI_Request> mpi_send_request_list;                      //vector to track all send requests
   std::vector<nid_t> const local_leaf_ids   = topology_.LocalLeafIds();//get local leaves
   unsigned int const number_of_local_leaves = local_leaf_ids.size();

   //x-Direction ( send West to East )
   MpiVertexFilter::InitializeSendRecvRankLists( send_list_per_rank, recv_list_per_rank, local_leaf_ids, MpiVertexFilter::Direction::X, topology_ );
   unsigned int const x_boundary_size = ( CC::ICZ() + 1 ) * ( CC::ICY() + 1 );
   //init send buffer for sends; exact size unknown -> size is upper bound
   std::vector<unsigned long long int> send_buffers( number_of_local_leaves * x_boundary_size );

   //send data
   for( unsigned int leaf_index = 0; leaf_index < number_of_local_leaves; leaf_index++ ) {
      nid_t node_id                              = local_leaf_ids[leaf_index];
      unsigned int send_buffer_offset            = leaf_index * x_boundary_size;
      unsigned long long int const vertex_offset = leave_offset[leaf_index];
      unsigned int const vertex_id_jump_y        = 8 * CC::ICX();
      unsigned int const vertex_id_jump_z        = 8 * CC::ICX() * CC::ICY();
      //load inner points
      MpiVertexFilter::LoadBoundaryBuffer( vertex_ids.data() + vertex_offset, send_buffers.data() + send_buffer_offset, vertex_id_jump_y, vertex_id_jump_z, CC::ICY(), CC::ICZ(), 0, 3, 4, 7 );
      //send data to all ranks necesarry
      MpiVertexFilter::SendFromNode( node_id, topology_.GetNeighboringLeaves( node_id, BoundaryLocation::West ), send_buffers.data() + send_buffer_offset, MpiVertexFilter::Direction::X, send_list_per_rank, mpi_send_request_list, topology_ );
   }

   //recv data
   for( unsigned int leaf_index = 0; leaf_index < local_leaf_ids.size(); leaf_index++ ) {
      nid_t const node_id = local_leaf_ids[leaf_index];
      std::vector<unsigned long long int> east_boundary_points( x_boundary_size );
      unsigned long long int const vertex_offset = leave_offset[leaf_index] + 8 * ( CC::ICX() - 1 );
      unsigned int const vertex_id_jump_y        = 8 * CC::ICX();
      unsigned int const vertex_id_jump_z        = 8 * CC::ICX() * CC::ICY();
      //init east vetrices
      MpiVertexFilter::LoadBoundaryBuffer( vertex_ids.data() + vertex_offset, east_boundary_points.data(), vertex_id_jump_y, vertex_id_jump_z, CC::ICY(), CC::ICZ(), 1, 2, 5, 6 );
      //receive data for node
      MpiVertexFilter::RecvToNode( node_id, topology_.GetNeighboringLeaves( node_id, BoundaryLocation::East ), east_boundary_points.data(), MpiVertexFilter::Direction::X, recv_list_per_rank, topology_ );

      //insert east boundary into cell vertices; iterate over every boundary block
      for( unsigned int k = 0; k < CC::ICZ(); k++ ) {
         for( unsigned int j = 0; j < CC::ICY(); j++ ) {
            const unsigned long long int vertex_offset = leave_offset[leaf_index];
            const unsigned long long int block_index   = 8 * ( CC::ICX() - 1 ) + 8 * j * CC::ICX() + 8 * k * CC::ICX() * CC::ICY() + vertex_offset;
            vertex_ids[block_index + 1]                = east_boundary_points[j + k * ( CC::ICY() + 1 )];
            vertex_ids[block_index + 2]                = east_boundary_points[( j + 1 ) + k * ( CC::ICY() + 1 )];
            vertex_ids[block_index + 5]                = east_boundary_points[j + ( k + 1 ) * ( CC::ICY() + 1 )];
            vertex_ids[block_index + 6]                = east_boundary_points[( j + 1 ) + ( k + 1 ) * ( CC::ICY() + 1 )];
         }
      }
   }

   //wait for all messages from x direction
   MPI_Waitall( mpi_send_request_list.size(), mpi_send_request_list.data(), MPI_STATUSES_IGNORE );

   //y direction ( send South to North )
   unsigned int const y_boundary_size = ( CC::ICZ() + 1 ) * ( CC::ICX() + 1 );

   if constexpr( CC::DIM() != Dimension::One ) {//in one-dimensional case no communication along the y-axis is needed

      MpiVertexFilter::InitializeSendRecvRankLists( send_list_per_rank, recv_list_per_rank, local_leaf_ids, MpiVertexFilter::Direction::Y, topology_ );
      mpi_send_request_list.clear();
      if( x_boundary_size != y_boundary_size )//only init new sendbuffer if size differs from existing one
         send_buffers.resize( number_of_local_leaves * y_boundary_size );
      for( unsigned int leaf_index = 0; leaf_index < number_of_local_leaves; leaf_index++ ) {
         nid_t const node_id                        = local_leaf_ids[leaf_index];
         unsigned int const send_buffer_offset      = leaf_index * y_boundary_size;
         unsigned long long int const vertex_offset = leave_offset[leaf_index];
         unsigned int const vertex_id_jump_x        = 8;
         unsigned int const vertex_id_jump_z        = 8 * CC::ICX() * CC::ICY();
         //load into buffer
         MpiVertexFilter::LoadBoundaryBuffer( vertex_ids.data() + vertex_offset, send_buffers.data() + send_buffer_offset, vertex_id_jump_x, vertex_id_jump_z, CC::ICX(), CC::ICZ(), 0, 1, 4, 5 );

         //send data to all ranks necessary
         MpiVertexFilter::SendFromNode( node_id, topology_.GetNeighboringLeaves( node_id, BoundaryLocation::South ), send_buffers.data() + send_buffer_offset, MpiVertexFilter::Direction::Y, send_list_per_rank, mpi_send_request_list, topology_ );
      }

      //recv all data
      for( unsigned int leaf_index = 0; leaf_index < local_leaf_ids.size(); leaf_index++ ) {
         nid_t const node_id = local_leaf_ids[leaf_index];
         //init north vetrices
         std::vector<unsigned long long int> north_boundary_points( y_boundary_size );
         unsigned long long int const vertex_offset = leave_offset[leaf_index] + 8 * CC::ICX() * ( CC::ICY() - 1 );
         unsigned int const vertex_id_jump_x        = 8;
         unsigned int const vertex_id_jump_z        = 8 * CC::ICX() * CC::ICY();
         MpiVertexFilter::LoadBoundaryBuffer( vertex_ids.data() + vertex_offset, north_boundary_points.data(), vertex_id_jump_x, vertex_id_jump_z, CC::ICX(), CC::ICZ(), 3, 2, 7, 6 );
         MpiVertexFilter::RecvToNode( node_id, topology_.GetNeighboringLeaves( node_id, BoundaryLocation::North ), north_boundary_points.data(), MpiVertexFilter::Direction::Y, recv_list_per_rank, topology_ );
         //insert north boundary into cell vertices; iterate over every boundary block
         for( unsigned int k = 0; k < CC::ICZ(); k++ ) {
            for( unsigned int i = 0; i < CC::ICX(); i++ ) {
               unsigned long long int const vertex_offset = leave_offset[leaf_index];
               unsigned long long int const block_index   = 8 * CC::ICX() * CC::ICY() * k + i * 8 + 8 * CC::ICX() * ( CC::ICY() - 1 ) + vertex_offset;
               vertex_ids[block_index + 3]                = north_boundary_points[i + k * ( CC::ICX() + 1 )];
               vertex_ids[block_index + 2]                = north_boundary_points[( i + 1 ) + k * ( CC::ICX() + 1 )];
               vertex_ids[block_index + 7]                = north_boundary_points[i + ( k + 1 ) * ( CC::ICX() + 1 )];
               vertex_ids[block_index + 6]                = north_boundary_points[( i + 1 ) + ( k + 1 ) * ( CC::ICX() + 1 )];
            }
         }
      }

      MPI_Waitall( mpi_send_request_list.size(), mpi_send_request_list.data(), MPI_STATUSES_IGNORE );
   }// if Dim != 1

   //z-direction ( send Bottom to Top )
   if constexpr( CC::DIM() == Dimension::Three ) {

      MpiVertexFilter::InitializeSendRecvRankLists( send_list_per_rank, recv_list_per_rank, local_leaf_ids, MpiVertexFilter::Direction::Z, topology_ );
      mpi_send_request_list.clear();
      unsigned int const z_boundary_size = ( CC::ICX() + 1 ) * ( CC::ICY() + 1 );
      if( y_boundary_size != z_boundary_size )
         send_buffers.resize( number_of_local_leaves * z_boundary_size );
      for( unsigned int leaf_index = 0; leaf_index < number_of_local_leaves; leaf_index++ ) {
         nid_t const node_id = local_leaf_ids[leaf_index];
         //Do MPI communication of boundaries
         unsigned int const send_buffer_offset = leaf_index * z_boundary_size;
         //load inner points
         unsigned long long int const vertex_offset = leave_offset[leaf_index];
         unsigned int const vertex_id_jump_x        = 8;
         unsigned int const vertex_id_jump_y        = 8 * CC::ICX();
         MpiVertexFilter::LoadBoundaryBuffer( vertex_ids.data() + vertex_offset, send_buffers.data() + send_buffer_offset, vertex_id_jump_x, vertex_id_jump_y, CC::ICX(), CC::ICY(), 0, 1, 3, 2 );
         //send data to all ranks necessary
         MpiVertexFilter::SendFromNode( node_id, topology_.GetNeighboringLeaves( node_id, BoundaryLocation::Bottom ), send_buffers.data() + send_buffer_offset, MpiVertexFilter::Direction::Z, send_list_per_rank, mpi_send_request_list, topology_ );
      }

      //recv all data
      for( unsigned int leaf_index = 0; leaf_index < local_leaf_ids.size(); leaf_index++ ) {
         nid_t const node_id = local_leaf_ids[leaf_index];
         //init boundary buffer
         std::vector<unsigned long long int> top_boundary_points( z_boundary_size );

         unsigned long long int const vertex_offset = leave_offset[leaf_index] + 8 * CC::ICX() * CC::ICY() * ( CC::ICZ() - 1 );
         unsigned int const vertex_id_jump_x        = 8;
         unsigned int const vertex_id_jump_y        = 8 * CC::ICX();
         MpiVertexFilter::LoadBoundaryBuffer( vertex_ids.data() + vertex_offset, top_boundary_points.data(), vertex_id_jump_x, vertex_id_jump_y, CC::ICX(), CC::ICY(), 4, 5, 7, 6 );
         MpiVertexFilter::RecvToNode( node_id, topology_.GetNeighboringLeaves( node_id, BoundaryLocation::Top ), top_boundary_points.data(), MpiVertexFilter::Direction::Z, recv_list_per_rank, topology_ );
         //insert boundary buffer into cell vertices; iterate over every boundary block
         for( unsigned int j = 0; j < CC::ICY(); j++ ) {
            for( unsigned int i = 0; i < CC::ICX(); i++ ) {
               unsigned long long int const vertex_offset = leave_offset[leaf_index];
               unsigned long long int const block_index   = 8 * ( i + j * CC::ICX() + CC::ICX() * CC::ICY() * ( CC::ICZ() - 1 ) ) + vertex_offset;
               vertex_ids[block_index + 4]                = top_boundary_points[i + j * ( CC::ICX() + 1 )];
               vertex_ids[block_index + 5]                = top_boundary_points[( i + 1 ) + j * ( CC::ICX() + 1 )];
               vertex_ids[block_index + 7]                = top_boundary_points[i + ( j + 1 ) * ( CC::ICX() + 1 )];
               vertex_ids[block_index + 6]                = top_boundary_points[( i + 1 ) + ( j + 1 ) * ( CC::ICX() + 1 )];
            }
         }
      }
      //wait for communication in z direction to be completed
      MPI_Waitall( mpi_send_request_list.size(), mpi_send_request_list.data(), MPI_STATUSES_IGNORE );
   }// if Dim == 3
}
