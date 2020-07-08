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
#include <catch.hpp>
#include "topology/topology_manager.h"
#include "test_mesh_generator_utilities.h"
#include "input_output/output_writer/mesh_generator/mesh_generator_utilities.h"
#include "input_output/output_writer/mesh_generator/standard_finest_level_mesh_generator.h"

/********************************************************************************************************************************************/
/*                                              TEST OF VERTEX COORDINATES                                                                  */
/********************************************************************************************************************************************/
/***********************************/
/* 1. Global dimensions check      */
/***********************************/
SCENARIO( "Standard mesh generator with finest level option: Dimension of global vertex coordinates vector is correctly computed", "[1rank],[2rank]" ) {

   GIVEN( "Underlying topology with Lmax being two" ) {

      // Parameter for the creation of the mesh geenrator
      std::array<unsigned int, 3> const number_of_blocks = { { 1, 1, 1 } };
      TopologyManager topology                           = TopologyManager( number_of_blocks, 1, 0 );
      Tree tree( topology, 1, 1.0 );
      std::unique_ptr<MeshGenerator const> mesh_generator = std::make_unique<StandardFinestLevelMeshGenerator const>( topology, tree, 1.0, number_of_blocks );

      WHEN( "The topology consists of just one node" ) {

         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );

         THEN( "The global vertex coordinates dimension vector has two entries with the first being equal the number of vertices (8 nodes on finest level) and the second being three" ) {
            std::vector<hsize_t> const global_vertex_coordinates_dimensions = mesh_generator->GetGlobalDimensionsOfVertexCoordinates();
            REQUIRE( global_vertex_coordinates_dimensions.size() == 2 );
            REQUIRE( global_vertex_coordinates_dimensions[0] == MeshGeneratorUtilities::NumberOfInternalVerticesPerBlock() );
            REQUIRE( global_vertex_coordinates_dimensions[1] == 3 );
         }
      }

      WHEN( "The node in the topology is refined ones (all on level 1)" ) {
         TestUtilities::RefineFirstNodeInTopology( topology );
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 9, 8 ) );

         THEN( "The global vertex coordinates dimension vector has two entries with the first being equal the number of vertices times the number of leafs (8) and the second being three" ) {
            std::vector<hsize_t> const global_vertex_coordinates_dimensions = mesh_generator->GetGlobalDimensionsOfVertexCoordinates();
            REQUIRE( global_vertex_coordinates_dimensions.size() == 2 );
            REQUIRE( global_vertex_coordinates_dimensions[0] == ( CC::ICX() * 2 + 1 ) * ( CC::ICY() * 2 + 1 ) * ( CC::ICZ() * 2 + 1 ) );
            REQUIRE( global_vertex_coordinates_dimensions[1] == 3 );
         }
      }
      WHEN( "The node in the topology is refined twice (all on level 1, and 8 on level 2)" ) {
         TestUtilities::RefineTwoNodesInTopology( topology );
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 17, 15 ) );

         THEN( "The global vertex coordinates dimension vector has two entries with the first being equal the number of vertices times the number of leafs (8) and the second being three" ) {
            std::vector<hsize_t> const global_vertex_coordinates_dimensions = mesh_generator->GetGlobalDimensionsOfVertexCoordinates();
            REQUIRE( global_vertex_coordinates_dimensions.size() == 2 );
            REQUIRE( global_vertex_coordinates_dimensions[0] == ( CC::ICX() * 4 + 1 ) * ( CC::ICY() * 4 + 1 ) * ( CC::ICZ() * 4 + 1 ) );
            REQUIRE( global_vertex_coordinates_dimensions[1] == 3 );
         }
      }
   }
}

/***********************************/
/* 2. Local dimensions check       */
/***********************************/
SCENARIO( "Standard mesh generator with finest level option: Vertex Coordinates local dimensions to write are correctly determined", "[1rank]" ) {

   GIVEN( "Underlying topology with Lmax being one" ) {
      // Parameter for the creation of the mesh geenrator
      std::array<unsigned int, 3> const number_of_blocks = { { 1, 1, 1 } };
      TopologyManager topology                           = TopologyManager( number_of_blocks, 1, 0 );
      Tree tree( topology, 1, 1.0 );
      std::unique_ptr<MeshGenerator const> mesh_generator = std::make_unique<StandardFinestLevelMeshGenerator const>( topology, tree, 1.0, number_of_blocks );

      WHEN( "The topology consists of just one node" ) {

         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );

         THEN( "The local vertex coordinates dimension vector has two entries with the first being equal the number of vertices and the second being three" ) {
            std::vector<hsize_t> const local_vertex_coordinates_dimensions = mesh_generator->GetLocalDimensionsOfVertexCoordinates();
            REQUIRE( local_vertex_coordinates_dimensions.size() == 2 );
            REQUIRE( local_vertex_coordinates_dimensions[0] == MeshGeneratorUtilities::NumberOfInternalVerticesPerBlock() );
            REQUIRE( local_vertex_coordinates_dimensions[1] == 3 );
         }
      }

      WHEN( "The node in the topology is refined" ) {
         TestUtilities::RefineFirstNodeInTopology( topology );
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 9, 8 ) );

         THEN( "The local vertex coordinates dimension vector has two entries with the first being equal the number of vertices time the number of leafs (8) and the second being three" ) {
            std::vector<hsize_t> const local_vertex_coordinates_dimensions = mesh_generator->GetLocalDimensionsOfVertexCoordinates();
            REQUIRE( local_vertex_coordinates_dimensions.size() == 2 );
            REQUIRE( local_vertex_coordinates_dimensions[0] == ( CC::ICX() * 2 + 1 ) * ( CC::ICY() * 2 + 1 ) * ( CC::ICZ() * 2 + 1 ) );
            REQUIRE( local_vertex_coordinates_dimensions[1] == 3 );
         }
      }

      WHEN( "The node in the topology is refined twice (all on level 1, and 8 on level 2)" ) {
         TestUtilities::RefineTwoNodesInTopology( topology );
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 17, 15 ) );

         THEN( "The global vertex coordinates dimension vector has two entries with the first being equal the number of vertices times the number of leafs (8) and the second being three" ) {
            std::vector<hsize_t> const global_vertex_coordinates_dimensions = mesh_generator->GetLocalDimensionsOfVertexCoordinates();
            REQUIRE( global_vertex_coordinates_dimensions.size() == 2 );
            REQUIRE( global_vertex_coordinates_dimensions[0] == ( CC::ICX() * 4 + 1 ) * ( CC::ICY() * 4 + 1 ) * ( CC::ICZ() * 4 + 1 ) );
            REQUIRE( global_vertex_coordinates_dimensions[1] == 3 );
         }
      }
   }
}

/***********************************/
/* 3. Start index check            */
/***********************************/
SCENARIO( "Standard mesh generator with finest level option: Vertex Coordinates start index to write is correctly determined", "[1rank]" ) {

   GIVEN( "Underlying topology with Lmax being one" ) {
      // Parameter for the creation of the mesh geenrator
      std::array<unsigned int, 3> const number_of_blocks = { { 1, 1, 1 } };
      TopologyManager topology                           = TopologyManager( number_of_blocks, 1, 0 );
      Tree tree( topology, 1, 1.0 );
      std::unique_ptr<MeshGenerator const> mesh_generator = std::make_unique<StandardFinestLevelMeshGenerator const>( topology, tree, 1.0, number_of_blocks );

      WHEN( "The topology consists of just one node" ) {

         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );

         THEN( "The local vertex coordinates start index is only value with entry 0" ) {
            hsize_t const local_vertex_coordinates_start_index = mesh_generator->GetLocalVertexCoordinatesStartIndex();
            REQUIRE( local_vertex_coordinates_start_index == 0 );
         }
      }

      WHEN( "The node in the topology is refined" ) {
         TestUtilities::RefineFirstNodeInTopology( topology );
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 9, 8 ) );

         THEN( "The local vertex coordinates start index is only value with entry 0" ) {
            hsize_t const local_vertex_coordinates_start_index = mesh_generator->GetLocalVertexCoordinatesStartIndex();
            REQUIRE( local_vertex_coordinates_start_index == 0 );
         }
      }
   }
}

/***********************************/
/* 4. Computation check            */
/***********************************/
SCENARIO( "Standard mesh generator with finest level option: Vertex coordinates node corner points are correctly written out", "[1rank]" ) {

   GIVEN( "Underlying topology with Lmax being one" ) {
      // Parameter for the creation of the mesh geenrator
      std::array<unsigned int, 3> const number_of_blocks = { { 1, 1, 1 } };
      TopologyManager topology                           = TopologyManager( number_of_blocks, 1, 0 );
      Tree tree( topology, 1, 1.0 );
      std::unique_ptr<MeshGenerator const> mesh_generator = std::make_unique<StandardFinestLevelMeshGenerator const>( topology, tree, 1.0, number_of_blocks );

      WHEN( "The topology consist of just one node" ) {

         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );

         std::vector<double> coordinates;
         mesh_generator->ComputeVertexCoordinates( coordinates );

         // Parameters to define the node counter for the three corners that are checked
         // 1: (x_max,0,0), 2: (xmax,ymax,0), 3: (xmax.ymax,zmax)
         std::array<unsigned long long int, 3> const corner_indices = { { ( CC::ICX() + 1 ) * 3 - 3, ( CC::ICX() + 1 ) * ( CC::ICY() + 1 ) * 3 - 3, ( CC::ICX() + 1 ) * ( CC::ICY() + 1 ) * ( CC::ICZ() + 1 ) * 3 - 3 } };

         THEN( "The coordinates vector has 3 * #of Leaves * #vertices per block entries and the first three all read zero and the last is one" ) {
            REQUIRE( coordinates.size() == MeshGeneratorUtilities::NumberOfInternalVerticesPerBlock() * 3 );

            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates[1] == Approx( 0.0 ) );
            REQUIRE( coordinates[2] == Approx( 0.0 ) );
            REQUIRE( coordinates[corner_indices[0]] == Approx( 1.0 ) );
            REQUIRE( coordinates[corner_indices[0] + 1] == Approx( 0.0 ) );
            REQUIRE( coordinates[corner_indices[0] + 2] == Approx( 0.0 ) );

            if constexpr( CC::DIM() != Dimension::One ) {
               REQUIRE( coordinates[corner_indices[1]] == Approx( 1.0 ) );
               REQUIRE( coordinates[corner_indices[1] + 1] == Approx( 1.0 ) );
               REQUIRE( coordinates[corner_indices[1] + 2] == Approx( 0.0 ) );
            }

            if constexpr( CC::DIM() == Dimension::Three ) {
               REQUIRE( coordinates[corner_indices[2]] == Approx( 1.0 ) );
               REQUIRE( coordinates[corner_indices[2] + 1] == Approx( 1.0 ) );
               REQUIRE( coordinates[corner_indices[2] + 2] == Approx( 1.0 ) );
            }
         }
      }

      WHEN( "The node in the topology is refined twice (all on level 1, and 8 on level 2)" ) {
         TestUtilities::RefineTwoNodesInTopology( topology );
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 17, 15 ) );

         std::vector<double> coordinates;
         mesh_generator->ComputeVertexCoordinates( coordinates );

         // Parameters to define the node counter for the three corners that are checked
         // 1: (x_max,0,0), 2: (xmax,ymax,0), 3: (xmax.ymax,zmax)
         std::array<unsigned long long int, 3> const corner_indices = { { ( CC::ICX() * 4 + 1 ) * 3 - 3, ( CC::ICX() * 4 + 1 ) * ( CC::ICY() * 4 + 1 ) * 3 - 3, ( CC::ICX() * 4 + 1 ) * ( CC::ICY() * 4 + 1 ) * ( CC::ICZ() * 4 + 1 ) * 3 - 3 } };

         THEN( "The coordinates vector has 3 * #of Leaves * #vertices * 8 #nodes and depending on the first coordinate of each node check corner points" ) {
            REQUIRE( coordinates.size() == ( CC::ICX() * 4 + 1 ) * ( CC::ICY() * 4 + 1 ) * ( CC::ICZ() * 4 + 1 ) * 3 );

            REQUIRE( coordinates[0] == Approx( 0.0 ) );
            REQUIRE( coordinates[1] == Approx( 0.0 ) );
            REQUIRE( coordinates[2] == Approx( 0.0 ) );
            REQUIRE( coordinates[corner_indices[0]] == Approx( 1.0 ) );
            REQUIRE( coordinates[corner_indices[0] + 1] == Approx( 0.0 ) );
            REQUIRE( coordinates[corner_indices[0] + 2] == Approx( 0.0 ) );

            if constexpr( CC::DIM() != Dimension::One ) {
               REQUIRE( coordinates[corner_indices[1]] == Approx( 1.0 ) );
               REQUIRE( coordinates[corner_indices[1] + 1] == Approx( 1.0 ) );
               REQUIRE( coordinates[corner_indices[1] + 2] == Approx( 0.0 ) );
            }

            if constexpr( CC::DIM() == Dimension::Three ) {
               REQUIRE( coordinates[corner_indices[2]] == Approx( 1.0 ) );
               REQUIRE( coordinates[corner_indices[2] + 1] == Approx( 1.0 ) );
               REQUIRE( coordinates[corner_indices[2] + 2] == Approx( 1.0 ) );
            }
         }
      }
   }
}

/********************************************************************************************************************************************/
/*                                                  TEST OF VERTEX IDs                                                                      */
/********************************************************************************************************************************************/
/***********************************/
/* 1. Global Dimension check       */
/***********************************/
SCENARIO( "Standard mesh generator with finest level option: Dimension of global vertex IDs vector is correctly computed", "[1rank],[2rank]" ) {

   GIVEN( "Underlying topology with Lmax being one" ) {

      // Parameter for the creation of the mesh geenrator
      std::array<unsigned int, 3> const number_of_blocks = { { 1, 1, 1 } };
      TopologyManager topology                           = TopologyManager( number_of_blocks, 1, 0 );
      Tree tree( topology, 1, 1.0 );
      std::unique_ptr<MeshGenerator const> mesh_generator = std::make_unique<StandardFinestLevelMeshGenerator const>( topology, tree, 1.0, number_of_blocks );

      WHEN( "The topology consists of just one node" ) {

         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );

         THEN( "The global vertex coordinates dimension vector has two entries with the first being equal the number of vertices and the second being three" ) {
            std::vector<hsize_t> const global_vertex_IDs_dimensions = mesh_generator->GetGlobalDimensionsOfVertexIDs();
            REQUIRE( global_vertex_IDs_dimensions.size() == 2 );
            REQUIRE( global_vertex_IDs_dimensions[0] == MeshGeneratorUtilities::NumberOfInternalCellsPerBlock() );
            REQUIRE( global_vertex_IDs_dimensions[1] == 8 );
         }
      }

      WHEN( "The node in the topology is refined" ) {
         TestUtilities::RefineFirstNodeInTopology( topology );
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 9, 8 ) );

         THEN( "The global vertex coordinates dimension vector has two entries with the first being equal the number of vertices time the number of leafs (8) and the second being three" ) {
            std::vector<hsize_t> const global_vertex_IDs_dimensions = mesh_generator->GetGlobalDimensionsOfVertexIDs();
            REQUIRE( global_vertex_IDs_dimensions.size() == 2 );
            REQUIRE( global_vertex_IDs_dimensions[0] == MeshGeneratorUtilities::NumberOfInternalCellsPerBlock() * 8 );
            REQUIRE( global_vertex_IDs_dimensions[1] == 8 );
         }
      }

      WHEN( "The node in the topology is refined twice (all on level 1, and 8 on level 2)" ) {
         TestUtilities::RefineTwoNodesInTopology( topology );
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 17, 15 ) );

         THEN( "The global vertex coordinates dimension vector has two entries with the first being equal the number of vertices time the number of leafs (8) and the second being three" ) {
            std::vector<hsize_t> const global_vertex_IDs_dimensions = mesh_generator->GetGlobalDimensionsOfVertexIDs();
            REQUIRE( global_vertex_IDs_dimensions.size() == 2 );
            REQUIRE( global_vertex_IDs_dimensions[0] == MeshGeneratorUtilities::NumberOfInternalCellsPerBlock() * 15 );
            REQUIRE( global_vertex_IDs_dimensions[1] == 8 );
         }
      }
   }
}

/***********************************/
/* 2. Local Dimension check        */
/***********************************/
SCENARIO( "Standard mesh generator with finest level option: Vertex IDs local dimensions to write are correctly determined", "[1rank]" ) {

   GIVEN( "Underlying topology with Lmax being one" ) {
      // Parameter for the creation of the mesh geenrator
      std::array<unsigned int, 3> const number_of_blocks = { { 1, 1, 1 } };
      TopologyManager topology                           = TopologyManager( number_of_blocks, 1, 0 );
      Tree tree( topology, 1, 1.0 );
      std::unique_ptr<MeshGenerator const> mesh_generator = std::make_unique<StandardFinestLevelMeshGenerator const>( topology, tree, 1.0, number_of_blocks );

      WHEN( "The topology consists of just one node" ) {

         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );

         THEN( "The local vertex IDs dimension vector has two entries with the first being equal the number of vertices and the second being three" ) {
            std::vector<hsize_t> const local_vertex_IDs_dimensions = mesh_generator->GetLocalDimensionsOfVertexIDs();
            REQUIRE( local_vertex_IDs_dimensions.size() == 2 );
            REQUIRE( local_vertex_IDs_dimensions[0] == MeshGeneratorUtilities::NumberOfInternalCellsPerBlock() );
            REQUIRE( local_vertex_IDs_dimensions[1] == 8 );
         }
      }

      WHEN( "The node in the topology is refined" ) {

         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );

         TestUtilities::RefineFirstNodeInTopology( topology );
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 9, 8 ) );

         THEN( "The local vertex IDs dimension vector has two entries with the first being equal the number of vertices time the number of leafs (8) and the second being three" ) {
            std::vector<hsize_t> const local_vertex_IDs_dimensions = mesh_generator->GetLocalDimensionsOfVertexIDs();
            REQUIRE( local_vertex_IDs_dimensions.size() == 2 );
            REQUIRE( local_vertex_IDs_dimensions[0] == MeshGeneratorUtilities::NumberOfInternalCellsPerBlock() * 8 );
            REQUIRE( local_vertex_IDs_dimensions[1] == 8 );
         }
      }

      WHEN( "The node in the topology is refined twice (all on level 1, and 8 on level 2)" ) {
         TestUtilities::RefineTwoNodesInTopology( topology );
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 17, 15 ) );

         THEN( "The local vertex IDs dimension vector has two entries with the first being equal the number of vertices time the number of leafs (8) and the second being three" ) {
            std::vector<hsize_t> const local_vertex_IDs_dimensions = mesh_generator->GetLocalDimensionsOfVertexIDs();
            REQUIRE( local_vertex_IDs_dimensions.size() == 2 );
            REQUIRE( local_vertex_IDs_dimensions[0] == MeshGeneratorUtilities::NumberOfInternalCellsPerBlock() * 15 );
            REQUIRE( local_vertex_IDs_dimensions[1] == 8 );
         }
      }
   }
}

/***********************************/
/* 3. Start index check            */
/***********************************/
SCENARIO( "Standard mesh generator with finest level option: Vertex IDs start index to write is correctly determined", "[1rank]" ) {

   GIVEN( "Underlying topology with Lmax being one" ) {
      // Parameter for the creation of the mesh geenrator
      std::array<unsigned int, 3> const number_of_blocks = { { 1, 1, 1 } };
      TopologyManager topology                           = TopologyManager( number_of_blocks, 1, 0 );
      Tree tree( topology, 1, 1.0 );
      std::unique_ptr<MeshGenerator const> mesh_generator = std::make_unique<StandardFinestLevelMeshGenerator const>( topology, tree, 1.0, number_of_blocks );

      WHEN( "The topology consists of just one node" ) {

         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );

         THEN( "The local vertex IDs start index is only value with entry 0" ) {
            hsize_t const local_vertex_IDs_start_index = mesh_generator->GetLocalVertexIDsStartIndex();
            REQUIRE( local_vertex_IDs_start_index == 0 );
         }
      }

      WHEN( "The node in the topology is refined" ) {
         TestUtilities::RefineFirstNodeInTopology( topology );
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 9, 8 ) );

         THEN( "The local vertex IDs start index is only value with entry 0" ) {
            hsize_t const local_vertex_IDs_start_index = mesh_generator->GetLocalVertexIDsStartIndex();
            REQUIRE( local_vertex_IDs_start_index == 0 );
         }
      }

      WHEN( "The node in the topology is refined twice (all on level 1, and 8 on level 2)" ) {
         TestUtilities::RefineTwoNodesInTopology( topology );
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 17, 15 ) );

         THEN( "The local vertex IDs start index is only value with entry 0" ) {
            hsize_t const local_vertex_IDs_start_index = mesh_generator->GetLocalVertexIDsStartIndex();
            REQUIRE( local_vertex_IDs_start_index == 0 );
         }
      }
   }
}

/********************************************************************************************************************************************/
/*                                             TEST COMBINATION VERTEX IDs AND COORDINATES                                               */
/********************************************************************************************************************************************/
SCENARIO( "Standard mesh generator with finest level option: Check that vertex IDs and coordinates are computed properly in combination", "[1rank]" ) {

   GIVEN( "Underlying topology with Lmax being one" ) {
      // Parameter for the creation of the mesh geenrator
      std::array<unsigned int, 3> const number_of_blocks = { { 1, 1, 1 } };
      TopologyManager topology                           = TopologyManager( number_of_blocks, 1, 0 );
      Tree tree( topology, 1, 1.0 );
      std::unique_ptr<MeshGenerator const> mesh_generator = std::make_unique<StandardFinestLevelMeshGenerator const>( topology, tree, 1.0, number_of_blocks );

      WHEN( "The topology consist of just one node" ) {

         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );

         std::vector<double> coordinates;
         mesh_generator->ComputeVertexCoordinates( coordinates );
         std::vector<unsigned long long int> vertex_ids;
         mesh_generator->ComputeVertexIDs( vertex_ids );

         THEN( "The coordinates of the IDs of each cell must be in an order where first x is increased, then y and last z" ) {
            // Check the total size of both vectors
            REQUIRE( coordinates.size() == MeshGeneratorUtilities::NumberOfInternalVerticesPerBlock() * 3 );
            REQUIRE( vertex_ids.size() == MeshGeneratorUtilities::NumberOfInternalCellsPerBlock() * 8 );
            // Check that the largest index in the vertex ids is not larger than the size of the coordinates vector
            REQUIRE( *std::max_element( vertex_ids.begin(), vertex_ids.end() ) < coordinates.size() );

            // Loop through all cells and check the coordinates of all vertex ids.
            /** 
             * A cell must be build in the following order (numbers correspond to vertex ID increment)
             *
             *                 7-------6
             *                /|      /|
             *               / |     / |
             *              3-------2  |
             * y   z        |  4----|--5
             *  ^ /         | /     | /
             *  |/          |/      |/
             *  --> x       0-------1
             */
            for( unsigned int cell = 0; cell < MeshGeneratorUtilities::NumberOfInternalCellsPerBlock(); cell++ ) {
               // Check that the x-coordinates of vertices are consistent
               REQUIRE( coordinates[vertex_ids[cell * 8 + 1] * 3] > coordinates[vertex_ids[cell * 8] * 3] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 2] * 3] == coordinates[vertex_ids[cell * 8 + 1] * 3] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 3] * 3] == coordinates[vertex_ids[cell * 8] * 3] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 4] * 3] == coordinates[vertex_ids[cell * 8] * 3] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 5] * 3] == coordinates[vertex_ids[cell * 8 + 1] * 3] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 6] * 3] == coordinates[vertex_ids[cell * 8 + 2] * 3] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 7] * 3] == coordinates[vertex_ids[cell * 8 + 3] * 3] );
               // Check that the y-coordinates of vertices are consistent
               REQUIRE( coordinates[vertex_ids[cell * 8 + 1] * 3 + 1] == coordinates[vertex_ids[cell * 8] * 3 + 1] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 2] * 3 + 1] > coordinates[vertex_ids[cell * 8 + 1] * 3 + 1] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 3] * 3 + 1] == coordinates[vertex_ids[cell * 8 + 2] * 3 + 1] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 4] * 3 + 1] == coordinates[vertex_ids[cell * 8] * 3 + 1] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 5] * 3 + 1] == coordinates[vertex_ids[cell * 8 + 4] * 3 + 1] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 6] * 3 + 1] == coordinates[vertex_ids[cell * 8 + 2] * 3 + 1] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 7] * 3 + 1] == coordinates[vertex_ids[cell * 8 + 3] * 3 + 1] );
               // Check that the z-coordinates of vertices are consistent
               REQUIRE( coordinates[vertex_ids[cell * 8 + 1] * 3 + 2] == coordinates[vertex_ids[cell * 8] * 3 + 2] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 2] * 3 + 2] == coordinates[vertex_ids[cell * 8 + 1] * 3 + 2] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 3] * 3 + 2] == coordinates[vertex_ids[cell * 8 + 2] * 3 + 2] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 4] * 3 + 2] > coordinates[vertex_ids[cell * 8] * 3 + 2] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 5] * 3 + 2] == coordinates[vertex_ids[cell * 8 + 4] * 3 + 2] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 6] * 3 + 2] == coordinates[vertex_ids[cell * 8 + 5] * 3 + 2] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 7] * 3 + 2] == coordinates[vertex_ids[cell * 8 + 6] * 3 + 2] );
            }
         }
      }

      WHEN( "The node in the topology is refined twice (all on level 1, and 8 on level 2)" ) {
         TestUtilities::RefineTwoNodesInTopology( topology );
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 17, 15 ) );

         std::vector<double> coordinates;
         mesh_generator->ComputeVertexCoordinates( coordinates );
         std::vector<unsigned long long int> vertex_ids;
         mesh_generator->ComputeVertexIDs( vertex_ids );

         THEN( "The coordinates of the IDs of each cell must be in an order where first x is increased, then y and last z" ) {
            // Check the total size of both vectors
            REQUIRE( coordinates.size() == ( CC::ICX() * 4 + 1 ) * ( CC::ICY() * 4 + 1 ) * ( CC::ICZ() * 4 + 1 ) * 3 );
            REQUIRE( vertex_ids.size() == MeshGeneratorUtilities::NumberOfInternalCellsPerBlock() * 8 * 15 );
            // Check that the largest index in the vertex ids is not larger than the size of the coordinates vector
            REQUIRE( *std::max_element( vertex_ids.begin(), vertex_ids.end() ) < coordinates.size() );

            // Loop through all cells and check the coordinates of all vertex ids.
            /** 
             * A cell must be build in the following order (numbers correspond to vertex ID increment)
             *
             *                 7-------6
             *                /|      /|
             *               / |     / |
             *              3-------2  |
             * y   z        |  4----|--5
             *  ^ /         | /     | /
             *  |/          |/      |/
             *  --> x       0-------1
             */
            for( unsigned int cell = 0; cell < 8 * MeshGeneratorUtilities::NumberOfInternalCellsPerBlock(); cell++ ) {
               // Check that the x-coordinates of vertices are consistent
               REQUIRE( coordinates[vertex_ids[cell * 8 + 1] * 3] > coordinates[vertex_ids[cell * 8] * 3] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 2] * 3] == coordinates[vertex_ids[cell * 8 + 1] * 3] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 3] * 3] == coordinates[vertex_ids[cell * 8] * 3] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 4] * 3] == coordinates[vertex_ids[cell * 8] * 3] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 5] * 3] == coordinates[vertex_ids[cell * 8 + 1] * 3] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 6] * 3] == coordinates[vertex_ids[cell * 8 + 2] * 3] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 7] * 3] == coordinates[vertex_ids[cell * 8 + 3] * 3] );
               // Check that the y-coordinates of vertices are consistent
               REQUIRE( coordinates[vertex_ids[cell * 8 + 1] * 3 + 1] == coordinates[vertex_ids[cell * 8] * 3 + 1] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 2] * 3 + 1] > coordinates[vertex_ids[cell * 8 + 1] * 3 + 1] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 3] * 3 + 1] == coordinates[vertex_ids[cell * 8 + 2] * 3 + 1] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 4] * 3 + 1] == coordinates[vertex_ids[cell * 8] * 3 + 1] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 5] * 3 + 1] == coordinates[vertex_ids[cell * 8 + 4] * 3 + 1] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 6] * 3 + 1] == coordinates[vertex_ids[cell * 8 + 2] * 3 + 1] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 7] * 3 + 1] == coordinates[vertex_ids[cell * 8 + 3] * 3 + 1] );
               // Check that the z-coordinates of vertices are consistent
               REQUIRE( coordinates[vertex_ids[cell * 8 + 1] * 3 + 2] == coordinates[vertex_ids[cell * 8] * 3 + 2] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 2] * 3 + 2] == coordinates[vertex_ids[cell * 8 + 1] * 3 + 2] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 3] * 3 + 2] == coordinates[vertex_ids[cell * 8 + 2] * 3 + 2] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 4] * 3 + 2] > coordinates[vertex_ids[cell * 8] * 3 + 2] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 5] * 3 + 2] == coordinates[vertex_ids[cell * 8 + 4] * 3 + 2] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 6] * 3 + 2] == coordinates[vertex_ids[cell * 8 + 5] * 3 + 2] );
               REQUIRE( coordinates[vertex_ids[cell * 8 + 7] * 3 + 2] == coordinates[vertex_ids[cell * 8 + 6] * 3 + 2] );
            }
         }
      }
   }
}