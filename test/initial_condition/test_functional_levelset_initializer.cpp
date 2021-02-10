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

#include <vector>
#include <string>
#include <cmath>

#include "initial_condition/functional_levelset_initializer.h"
#include "utilities/vector_utilities.h"
#include "topology/id_information.h"
#include "initial_condition/test_initial_condition_helper.h"

SCENARIO( "The functional initializer works properly for a quasi 1D interface.", "[1rank]" ) {
   GIVEN( "A topology with: node_size = 1.0, maximum level = 0, n_materials = 2." ) {
      // General data
      constexpr double node_size                        = 1.0;
      constexpr unsigned int max_level                  = 0;
      std::vector<std::vector<nid_t>> const level_nodes = TestInitialCondition::GetChildNodesForSingleRootNode( max_level );
      std::vector<MaterialName> const materials         = { MaterialName::MaterialOne, MaterialName::MaterialTwo };
      double const cell_size_finest_level               = CellSizeOfLevel( node_size, max_level );

      WHEN( "The interface lies at x=0.5 and no bounding box is given." ) {
         // Parametric expression data
         std::string const expression                            = "phi := 0.5 - x";
         std::vector<std::array<double, 6>> const bounding_boxes = {};
         // Create the initializer
         std::unique_ptr<LevelsetInitializer> initializer( std::make_unique<FunctionalLevelsetInitializer>( expression, bounding_boxes,
                                                                                                            materials, node_size, max_level ) );

         THEN( "All nodes must contain two materials." ) {
            for( nid_t const node_id : level_nodes[0] ) {
               std::vector<MaterialName> const found_materials = initializer->GetInitialMaterials( node_id );
               REQUIRE( found_materials.size() == 2 );
            }
         }
         THEN( "The levelset value must be same as the distance of the cell center x-coordinate to 0.5. The sign must be positive for cells below 0.5." ) {
            for( nid_t const node_id : level_nodes[max_level] ) {
               std::array<double, 3> const origin = DomainCoordinatesOfId( node_id, DomainSizeOfId( node_id, node_size ) );
               double levelset[CC::TCX()][CC::TCY()][CC::TCZ()];
               initializer->GetInitialLevelset( node_id, levelset );

               // Loop through all elements and check the data
               for( unsigned int i = 0; i < CC::TCX(); ++i ) {
                  double x = origin[0] + ( double( i ) - double( CC::FICX() ) + 0.5 ) * cell_size_finest_level;
                  for( unsigned int j = 0; j < CC::TCY(); ++j ) {
                     for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        REQUIRE( levelset[i][j][k] == Approx( ( 0.5 - x ) / cell_size_finest_level ) );
                     }
                  }
               }
            }
         }
      }
      WHEN( "The interface lies at y=0.5 and no bounding box is given." ) {
         // Parametric expression data
         std::string const expression                            = "phi := 0.5 - y";
         std::vector<std::array<double, 6>> const bounding_boxes = {};
         // Create the initializer
         std::unique_ptr<LevelsetInitializer> initializer( std::make_unique<FunctionalLevelsetInitializer>( expression, bounding_boxes,
                                                                                                            materials, node_size, max_level ) );

         THEN( "All nodes must contain two materials." ) {
            for( nid_t const node_id : level_nodes[0] ) {
               std::vector<MaterialName> const found_materials = initializer->GetInitialMaterials( node_id );
               REQUIRE( found_materials.size() == 2 );
            }
         }
         THEN( "The levelset value must be same as the distance of the cell center y-coordinate to 0.5. The sign must be positive for cells below 0.5." ) {
            for( nid_t const node_id : level_nodes[max_level] ) {
               std::array<double, 3> const origin = DomainCoordinatesOfId( node_id, DomainSizeOfId( node_id, node_size ) );
               double levelset[CC::TCX()][CC::TCY()][CC::TCZ()];
               initializer->GetInitialLevelset( node_id, levelset );

               // Loop through all elements and check the data
               for( unsigned int i = 0; i < CC::TCX(); ++i ) {
                  for( unsigned int j = 0; j < CC::TCY(); ++j ) {
                     double y = origin[1] + ( double( j ) - double( CC::FICY() ) + 0.5 ) * cell_size_finest_level;
                     for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        REQUIRE( levelset[i][j][k] == Approx( ( 0.5 - y ) / cell_size_finest_level ) );
                     }
                  }
               }
            }
         }
      }

      WHEN( "The interface lies at z=0.5 and no bounding box is given." ) {
         // Parametric expression data
         std::string const expression                            = "phi := -0.5 + z";
         std::vector<std::array<double, 6>> const bounding_boxes = {};
         // Create the initializer
         std::unique_ptr<LevelsetInitializer> initializer( std::make_unique<FunctionalLevelsetInitializer>( expression, bounding_boxes,
                                                                                                            materials, node_size, max_level ) );

         THEN( "All nodes must contain two materials." ) {
            for( nid_t const node_id : level_nodes[0] ) {
               std::vector<MaterialName> const found_materials = initializer->GetInitialMaterials( node_id );
               REQUIRE( found_materials.size() == 2 );
            }
         }
         THEN( "The levelset value must be same as the distance of the cell center z-coordinate to 0.5. The sign must be positive for cells above 0.5." ) {
            for( nid_t const node_id : level_nodes[max_level] ) {
               std::array<double, 3> const origin = DomainCoordinatesOfId( node_id, DomainSizeOfId( node_id, node_size ) );
               double levelset[CC::TCX()][CC::TCY()][CC::TCZ()];
               initializer->GetInitialLevelset( node_id, levelset );

               // Loop through all elements and check the data
               for( unsigned int i = 0; i < CC::TCX(); ++i ) {
                  for( unsigned int j = 0; j < CC::TCY(); ++j ) {
                     for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        double z = origin[2] + ( double( k ) - double( CC::FICZ() ) + 0.5 ) * cell_size_finest_level;
                        REQUIRE( levelset[i][j][k] == Approx( ( z - 0.5 ) / cell_size_finest_level ) );
                     }
                  }
               }
            }
         }
      }
   }
}

SCENARIO( "The functional initializer works properly for a 3D interface.", "[1rank]" ) {
   GIVEN( "A topology with: node_size = 1.0, maximum level = 1, n_materials = 2." ) {
      // General data
      constexpr double node_size                        = 1.0;
      constexpr unsigned int max_level                  = 1;
      std::vector<std::vector<nid_t>> const level_nodes = TestInitialCondition::GetChildNodesForSingleRootNode( max_level );
      std::vector<MaterialName> const materials         = { MaterialName::MaterialOne, MaterialName::MaterialTwo };
      double const cell_size_finest_level               = CellSizeOfLevel( node_size, max_level );

      WHEN( "The expression is a sphere with radius 0.125 and the bounding box is set to {0.0, 0.5, 0.0, 0.5, 0.0, 0.5}" ) {
         // Parametric expression data
         std::string const expression                            = "phi := -0.125 + sqrt( (x-0.25)*(x-0.25) + (y-0.25)*(y-0.25) + (z-0.25)*(z-0.25) );";
         std::vector<std::array<double, 6>> const bounding_boxes = { { 0.0, 0.5, 0.0, 0.5, 0.0, 0.5 } };
         // Create the initializer
         std::unique_ptr<LevelsetInitializer> initializer( std::make_unique<FunctionalLevelsetInitializer>( expression, bounding_boxes,
                                                                                                            materials, node_size, max_level ) );

         THEN( "Only the node in the west-south-bottom corner on level 1 should contain two materials. All others only the positive." ) {
            for( nid_t const node_id : level_nodes[max_level] ) {
               std::vector<MaterialName> const found_materials = initializer->GetInitialMaterials( node_id );
               if( IsExternalBoundary( BoundaryLocation::West, node_id, { 1, 1, 1 } ) &&
                   IsExternalBoundary( BoundaryLocation::South, node_id, { 1, 1, 1 } ) &&
                   IsExternalBoundary( BoundaryLocation::Bottom, node_id, { 1, 1, 1 } ) ) {
                  REQUIRE( found_materials.size() == 2 );
               } else {
                  REQUIRE( found_materials.size() == 1 );
                  REQUIRE( found_materials[0] == MaterialName::MaterialOne );
               }
            }
         }
         THEN( "The node on level zero should contain two materials." ) {
            for( nid_t const node_id : level_nodes[0] ) {
               std::vector<MaterialName> const found_materials = initializer->GetInitialMaterials( node_id );
               REQUIRE( found_materials.size() == 2 );
            }
         }
         THEN( "The levelset values of the nodes on level 1 should all be positive (except south-west-bottom). Outside bounding box the cut-off should be."
               "In the south-west-bottom the negative values are below the radius/cell_size." ) {
            for( nid_t const node_id : level_nodes[max_level] ) {
               std::array<double, 3> const origin = DomainCoordinatesOfId( node_id, DomainSizeOfId( node_id, node_size ) );
               double levelset[CC::TCX()][CC::TCY()][CC::TCZ()];
               initializer->GetInitialLevelset( node_id, levelset );

               // CHeck the node boundary
               bool const is_south_west_bottom = IsExternalBoundary( BoundaryLocation::West, node_id, { 1, 1, 1 } ) &&
                                                 IsExternalBoundary( BoundaryLocation::South, node_id, { 1, 1, 1 } ) &&
                                                 IsExternalBoundary( BoundaryLocation::Bottom, node_id, { 1, 1, 1 } );

               // Loop through all elements and check the data
               for( unsigned int i = 0; i < CC::TCX(); ++i ) {
                  double x = origin[0] + ( double( i ) - double( CC::FICX() ) + 0.5 ) * cell_size_finest_level;
                  for( unsigned int j = 0; j < CC::TCY(); ++j ) {
                     double y = origin[1] + ( double( j ) - double( CC::FICY() ) + 0.5 ) * cell_size_finest_level;
                     for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        double z = origin[2] + ( double( k ) - double( CC::FICZ() ) + 0.5 ) * cell_size_finest_level;
                        if( levelset[i][j][k] < 0.0 ) {
                           REQUIRE( is_south_west_bottom );
                           REQUIRE( std::abs( levelset[i][j][k] ) < 0.125 / cell_size_finest_level );
                        } else {
                           if( x > 0.5 or y > 0.5 or z > 0.5 ) {
                              REQUIRE( levelset[i][j][k] == CC::LSCOF() );
                           } else {
                              REQUIRE( levelset[i][j][k] > 0.0 );
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
}
