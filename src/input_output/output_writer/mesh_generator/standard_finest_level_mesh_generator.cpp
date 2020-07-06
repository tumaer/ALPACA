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
#include "input_output/output_writer/mesh_generator/standard_finest_level_mesh_generator.h"
#include "input_output/output_writer/mesh_generator/mesh_generator_utilities.h"
#include "communication/mpi_utilities.h"
#include "topology/id_information.h"

/**
 * @brief Constructor to create the standard finest level mesh generator.
 * @param topology_manager Instance to provide node information on different ranks.
 * @param flower Instance to provide node information of current rank.
 * @param dimensionalized_node_size_on_level_zero Already dimensionalized size of a node on level zero.
 * @param number_of_nodes_on_level_zero Number of nodes in all directions on level zero.
 */
StandardFinestLevelMeshGenerator::StandardFinestLevelMeshGenerator( TopologyManager const& topology_manager,
                                                                    Tree const& flower,
                                                                    double const dimensionalized_node_size_on_level_zero,
                                                                    std::array<unsigned int, 3> const number_of_nodes_on_level_zero ) :
   MeshGenerator( topology_manager, flower, dimensionalized_node_size_on_level_zero ),
   number_of_nodes_on_level_zero_( number_of_nodes_on_level_zero ) {
   /** Empty besides call of base class constructor and initializer list */
}

/**
 * @brief Counts the vertices per dimension for the current global maximum level of the simulation.
 * @return Number of vertices in each dimension on maximum level.
 */
std::array<unsigned long long int, 3> StandardFinestLevelMeshGenerator::GlobalNumberOfVerticesPerDimension() const {

   unsigned long long int&& resolution = 1 << topology_.GetCurrentMaximumLevel();
   unsigned long long int&& vertex_count_x = number_of_nodes_on_level_zero_[0] * resolution * CC::ICX() + 1;
   unsigned long long int&& vertex_count_y = number_of_nodes_on_level_zero_[1] * resolution * CC::ICY() + 1;
   unsigned long long int&& vertex_count_z = number_of_nodes_on_level_zero_[2] * resolution * CC::ICZ() + 1;
   return { vertex_count_x, vertex_count_y, vertex_count_z };
}

/**
 * @brief Returns the total global number of vertices on the finest level for the global maxmum level of the simulation.
 * @return Total number of vertices on finest level.
 */
unsigned long long int StandardFinestLevelMeshGenerator::TotalNumberOfVerticesOnFinestLevel() const {
   auto&& count_per_dimension = GlobalNumberOfVerticesPerDimension();
   return count_per_dimension[0] * count_per_dimension[1] * count_per_dimension[2];
}

/**
 * @brief Gives the local number of vertices present on this rank on finest level and start index.
 * @return Return the local number of vertices and start index.
 */
std::array<unsigned long long int, 2> StandardFinestLevelMeshGenerator::GetLocalNumberOfVerticesAndStartIndex() const {
   // Define the total number of vertices and current rank situation
   unsigned long long int total_number_vertices = TotalNumberOfVerticesOnFinestLevel();
   int const rank = MpiUtilities::MyRankId();
   int const number_of_ranks = MpiUtilities::NumberOfRanks();
   // Compute the number of local number of vertices and given start index
   int const remainder = total_number_vertices % number_of_ranks;
   unsigned long long int local_number_vertices = total_number_vertices / number_of_ranks;
   unsigned long long int const start_index = local_number_vertices * rank + ( rank < remainder ? rank : remainder );
   if(  rank  < remainder ) { local_number_vertices++; }
   return { local_number_vertices, start_index };
}

/**
 * @brief See base class implementation.
 */
std::vector<std::reference_wrapper<Node const>> StandardFinestLevelMeshGenerator::DoGetLocalNodes() const {
   return tree_.Leaves();
}

/**
 * @brief See base class implementation.
 */
hsize_t StandardFinestLevelMeshGenerator::DoGetGlobalNumberOfCells() const {
   return hsize_t( std::get<1>( topology_.NodeAndLeafCount() ) ) * MeshGeneratorUtilities::NumberOfInternalCellsPerBlock();
}

/**
 * @brief See base class implementation.
 */
hsize_t StandardFinestLevelMeshGenerator::DoGetLocalNumberOfCells() const {
   return hsize_t( topology_.LocalLeafIds().size() ) * MeshGeneratorUtilities::NumberOfInternalCellsPerBlock();
}

/**
 * @brief See base class implementation.
 */
hsize_t StandardFinestLevelMeshGenerator::DoGetLocalCellsStartIndex() const {
   return hsize_t( topology_.LeafOffsetOfRank( MpiUtilities::MyRankId() ) ) * MeshGeneratorUtilities::NumberOfInternalCellsPerBlock();
}

/**
 * @brief See base class implementation.
 */
std::vector<hsize_t> StandardFinestLevelMeshGenerator::DoGetGlobalDimensionsOfVertexCoordinates() const {
   return { hsize_t( TotalNumberOfVerticesOnFinestLevel() ), hsize_t ( 3 ) };
}

/**
 * @brief See base class implementation.
 */
std::vector<hsize_t> StandardFinestLevelMeshGenerator::DoGetLocalDimensionsOfVertexCoordinates() const {
   return { hsize_t( GetLocalNumberOfVerticesAndStartIndex()[0] ), hsize_t ( 3 ) };
}

/**
 * @brief See base class implementation.
 */
hsize_t StandardFinestLevelMeshGenerator::DoGetLocalVertexCoordinatesStartIndex() const {
   return hsize_t( GetLocalNumberOfVerticesAndStartIndex()[1] );
}

/**
 * @brief See base class definition.
 */
void StandardFinestLevelMeshGenerator::DoComputeVertexCoordinates( std::vector<double> & vertex_coordinates ) const {
   // Obtain correct local number of vertices and start index
   unsigned int const maximum_level = topology_.GetCurrentMaximumLevel();
   std::array<unsigned long long int, 2>&& local_number_vertices_and_startIndex = GetLocalNumberOfVerticesAndStartIndex();

   // resize the vector to ensure enough memory for the cooridnates ( x,y,z coordinates for each vertex )
   vertex_coordinates.resize( local_number_vertices_and_startIndex[0] * 3 );

   // Generate the vertex coordinates on the finest level
   // Remember: The mesh generator class already take the dimensionalized node size. Therfore, no additional dimensionalization is required.
   unsigned long long int const resolution = 1 << maximum_level;
   double const cell_size_x = dimensionalized_node_size_on_level_zero_ / resolution / CC::ICX();
   double const cell_size_y = dimensionalized_node_size_on_level_zero_ / resolution / CC::ICY();
   double const cell_size_z = dimensionalized_node_size_on_level_zero_ / resolution / CC::ICZ();

   // Generate vertex coordinates
   std::array<unsigned long long int, 3> const global_number_of_vertices = GlobalNumberOfVerticesPerDimension();
   for( unsigned long long int local_vertex_index = 0; local_vertex_index < local_number_vertices_and_startIndex[0]; local_vertex_index++ ) {
      // The index in in the global vector of all vertices is defined by the start index plus the local running index
      unsigned long long int  const global_vertex_index = local_vertex_index + local_number_vertices_and_startIndex[1];
      //invert 3D to 1D array indexing; global_vertex_index = x + X * y + X * Y * z
      unsigned long int const x =   global_vertex_index %   global_number_of_vertices[0];
      unsigned long int const y = ( global_vertex_index /   global_number_of_vertices[0] ) % global_number_of_vertices[1];
      unsigned long int const z =   global_vertex_index / ( global_number_of_vertices[0]   * global_number_of_vertices[1] );
      vertex_coordinates[local_vertex_index * 3    ] = x * cell_size_x;
      vertex_coordinates[local_vertex_index * 3 + 1] = y * cell_size_y;
      vertex_coordinates[local_vertex_index * 3 + 2] = z * cell_size_z;
   }
}

/**
 * @brief See base class definition.
 */
void StandardFinestLevelMeshGenerator::DoComputeVertexIDs( std::vector<unsigned long long int> & vertex_ids ) const {

   // Resize Vertex ID vector ( 8 vertices span one cell )
   vertex_ids.resize( topology_.LocalLeafIds().size() * MeshGeneratorUtilities::NumberOfInternalCellsPerBlock() * 8 );
   /**
    * Definition of offset parameters in x and y-direction to be specified for loop computation of vertices
    * Vertices are described in increasing order:
    * v( x,y,z ) = ( 0,0,0 )               -> ID: 0
    * v( x,y,z ) = ( total,0,0 )           -> ID: max
    * v( x,y,z ) = ( 0,1,0 )               -> ID: max + 1
    * v( x,y,z ) = ( 0,0,1 )               -> ID: ( max +1 ) * ( max + 1 )
    * v( x,y,z ) = ( total, total, total ) -> ID: ( max + 1 ) * ( max + 1 ) * ( max + 1 )
    */
   unsigned int maximum_level = topology_.GetCurrentMaximumLevel();
   std::array<unsigned long long int, 3>&& vertex_count_per_dimension = GlobalNumberOfVerticesPerDimension();
   unsigned long long int const vertex_count_x = vertex_count_per_dimension[0];
   unsigned long long int const vertex_count_y = vertex_count_per_dimension[1];

   //iterate over all leafs to generate vertex ids
   unsigned long long int vertex_id_counter = 0;
   unsigned long long int const resolution = 1 << maximum_level;
   for( std::uint64_t const& id : topology_.LocalLeafIds() ){
      // Additional distance to be considered between current level of node and maximum level where IDs are specified
      unsigned long long int const level_factor = resolution >> LevelOfNode( id );
      //Find index of nodes domain by coordinates function; Cast from double to long! -> for indices greater 2^53 this is incorrect
      std::array<double, 3> const block_origin = DomainCoordinatesOfId( id, 1 );
      for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
         for( unsigned int j = 0; j < CC::ICY(); ++j ) {
            for( unsigned int i = 0; i < CC::ICX(); ++i ) {
               unsigned long long int vertex_index = 0;
               unsigned long long int const index_x = ( block_origin[0] * CC::ICX() + i ) * level_factor;
               unsigned long long int const index_y = ( block_origin[1] * CC::ICY() + j ) * level_factor;
               unsigned long long int const index_z = ( block_origin[2] * CC::ICZ() + k ) * level_factor;

               vertex_index = index_x + index_y * vertex_count_x + index_z * vertex_count_x * vertex_count_y;
               vertex_ids[vertex_id_counter] = vertex_index;

               vertex_index = ( index_x + level_factor ) + index_y * vertex_count_x + index_z * vertex_count_x * vertex_count_y;
               vertex_ids[vertex_id_counter + 1] = vertex_index;

               vertex_index = ( index_x + level_factor ) + ( index_y + level_factor ) * vertex_count_x + index_z * vertex_count_x * vertex_count_y;
               vertex_ids[vertex_id_counter + 2] = vertex_index;

               vertex_index = index_x + ( index_y + level_factor ) * vertex_count_x + index_z * vertex_count_x * vertex_count_y;
               vertex_ids[vertex_id_counter + 3] = vertex_index;

               vertex_index = index_x + index_y * vertex_count_x + ( index_z + level_factor ) * vertex_count_x * vertex_count_y;
               vertex_ids[vertex_id_counter + 4] = vertex_index;

               vertex_index = ( index_x + level_factor ) + index_y * vertex_count_x + ( index_z + level_factor ) * vertex_count_x * vertex_count_y;
               vertex_ids[vertex_id_counter + 5] = vertex_index;

               vertex_index = ( index_x + level_factor ) + ( index_y + level_factor ) * vertex_count_x + ( index_z + level_factor ) * vertex_count_x * vertex_count_y;
               vertex_ids[vertex_id_counter + 6] = vertex_index;

               vertex_index = index_x + ( index_y + level_factor ) * vertex_count_x + ( index_z + level_factor ) * vertex_count_x * vertex_count_y;
               vertex_ids[vertex_id_counter + 7] = vertex_index;

               vertex_id_counter += 8;
            }
         }
      }
   }//for node id

   /**
    * @note: In other mesh_generators filtering operations are required to remove double placed vertices or to assign the
    *        correct ID to double placed vertices. In the finest level Output this is not required, since the IDs are mapped to the unique IDs on the finest
    *        level.
    */
}
