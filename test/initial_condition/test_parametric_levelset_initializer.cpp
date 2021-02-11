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
#include <catch2/catch.hpp>

#include <vector>
#include <string>
#include <cmath>

#include "initial_condition/parametric_levelset_initializer.h"
#include "utilities/vector_utilities.h"
#include "topology/id_information.h"
#include "initial_condition/test_initial_condition_helper.h"

namespace {
   constexpr double m_3_PI_2     = 3.0 * M_PI_2;
   constexpr double m_2_PI       = 2.0 * M_PI;
   std::string const string_pi   = std::to_string( M_PI );
   std::string const string_2_pi = std::to_string( m_2_PI );
}// namespace

SCENARIO( "The parametric interface coordinates are computed properly.", "[1rank]" ) {

   GIVEN( "A 2D parametric expression x=a*cos(phi), y=b*sin(phi) for elliptical corodinate computation.\n"
          "The spatial variable names {x, y} and the parametric variables theta phi.\n"
          "The semi axes are a=1.75 and b=0.5." ) {

      std::vector<std::string> const var_names = { "x", "y" };
      std::string const expression             = "x := 1.75*cos(phi); y := 0.5*sin(phi);";

      WHEN( "phi is in the range phi = [0, 2pi] with 101 points." ) {
         std::array<ParametricVariable, 2> variables = { ParametricVariable( "phi", 0.0, m_2_PI, 101 ),
                                                         ParametricVariable() };
         std::vector<std::array<double, 3>> const interface_coordinates( ComputeInterfaceCoordinates( expression, var_names, variables ) );

         THEN( "The number of interface coordinates should equal 121" ) {
            REQUIRE( interface_coordinates.size() == 101 );
         }
         THEN( "All x coordinates must be below equal 1.75 and y below 0.5. and z equal zero" ) {
            for( auto const& coord : interface_coordinates ) {
               REQUIRE( coord[0] <= 1.75 );
               REQUIRE( coord[1] <= 0.5 );
               REQUIRE( coord[2] == 0.0 );
            }
         }
      }
   }

   GIVEN( "A 3D parametric user expression for a 1D interface." ) {

      std::vector<std::string> const var_names    = { "x", "y", "z" };
      std::array<ParametricVariable, 2> variables = { ParametricVariable( "r", 0.0, 1.0, 21 ),
                                                      ParametricVariable( "s", 0.0, 1.0, 21 ) };

      WHEN( "The interface lies at x=0.475" ) {
         std::string const expression = "x := 0.475; y := r; z := s";
         std::vector<std::array<double, 3>> const interface_coordinates( ComputeInterfaceCoordinates( expression, var_names, variables ) );
         THEN( "The inteface coordinates should be 0.475 in x-direction. y and z should be multiple of 0.05." ) {
            REQUIRE( interface_coordinates.size() == 441 );
            for( auto const& coord : interface_coordinates ) {
               REQUIRE( coord[0] == 0.475 );
               REQUIRE( std::remainder( coord[1], 0.05 ) == Approx( 0.0 ).margin( 1e-16 ) );
               REQUIRE( std::remainder( coord[2], 0.05 ) == Approx( 0.0 ).margin( 1e-16 ) );
            }
         }
      }
      WHEN( "The interface lies at y=1.2345" ) {
         std::string const expression = "x := r; y := 1.2345; z := s";
         std::vector<std::array<double, 3>> const interface_coordinates( ComputeInterfaceCoordinates( expression, var_names, variables ) );
         THEN( "The inteface coordinates should be 1.2345 in y-direction. y and z should be multiple of 0.05." ) {
            REQUIRE( interface_coordinates.size() == 441 );
            for( auto const& coord : interface_coordinates ) {
               REQUIRE( std::remainder( coord[0], 0.05 ) == Approx( 0.0 ).margin( 1e-16 ) );
               REQUIRE( coord[1] == 1.2345 );
               REQUIRE( std::remainder( coord[2], 0.05 ) == Approx( 0.0 ).margin( 1e-16 ) );
            }
         }
      }
      WHEN( "The interface lies at z=-6.789" ) {
         std::string const expression = "x := r; y := s; z := -6.789";
         std::vector<std::array<double, 3>> const interface_coordinates( ComputeInterfaceCoordinates( expression, var_names, variables ) );
         THEN( "The inteface coordinates should be -6.789 in z-direction. y and z should be multiple of 0.05." ) {
            REQUIRE( interface_coordinates.size() == 441 );
            for( auto const& coord : interface_coordinates ) {
               REQUIRE( std::remainder( coord[0], 0.05 ) == Approx( 0.0 ).margin( 1e-16 ) );
               REQUIRE( std::remainder( coord[1], 0.05 ) == Approx( 0.0 ).margin( 1e-16 ) );
               REQUIRE( coord[2] == -6.789 );
            }
         }
      }
   }

   GIVEN( "A 3D parametric user expression of type x=r*sin(theta)*cos(phi), y=r*sin(theta)*sin(phi), z=r*cos(theta) for spherical corodinate computation.\n"
          "The spatial variable names {x, y, z} and two parametric variables theta and phi." ) {

      std::vector<std::string> const var_names = { "x", "y", "z" };

      WHEN( "r=1.75 and the angles have 11 points in the range theta = [-pi, pi) and phi = [0, 2pi)" ) {
         std::string const expression                = "x := 1.75*cos(theta)*cos(phi); y := 1.75*cos(theta)*sin(phi); z := 1.75*sin(theta);";
         std::array<ParametricVariable, 2> variables = { ParametricVariable( "theta", -M_PI, M_PI, 11 ),
                                                         ParametricVariable( "phi", 0.0, m_2_PI, 11 ) };
         std::vector<std::array<double, 3>> const interface_coordinates( ComputeInterfaceCoordinates( expression, var_names, variables ) );

         THEN( "The number of interface coordinates should equal 121" ) {
            REQUIRE( interface_coordinates.size() == 121 );
         }
         THEN( "The inteface coordinates should be 1.75 from the center {0,0,0} for all points." ) {
            for( auto const& coord : interface_coordinates ) {
               REQUIRE( VU::Distance( coord, { 0.0, 0.0, 0.0 } ) == Approx( 1.75 ) );
            }
         }
      }
      WHEN( "r=2.25, theta=0.5pi and phi has four points in the range [0, 2pi]" ) {
         std::string const expression                = "x := 2.25*sin(theta)*cos(phi); y := 2.25*sin(theta)*sin(phi); z := 2.25*cos(theta);";
         std::array<ParametricVariable, 2> variables = { ParametricVariable( "theta", M_PI_2, M_PI_2, 1 ),
                                                         ParametricVariable( "phi", 0.0, m_2_PI, 3 ) };
         std::vector<std::array<double, 3>> const interface_coordinates( ComputeInterfaceCoordinates( expression, var_names, variables ) );

         THEN( "The number of interface coordinates should equal 3" ) {
            REQUIRE( interface_coordinates.size() == 3 );
         }
         THEN( "The points should equal {+-r,0,0}" ) {
            for( auto const& coord : interface_coordinates ) {
               REQUIRE( std::abs( coord[0] ) == Approx( 2.25 ) );
               REQUIRE( coord[1] == Approx( 0.0 ).margin( 1e-15 ) );
               REQUIRE( coord[2] == Approx( 0.0 ).margin( 1e-15 ) );
            }
         }
      }
      WHEN( "r=2.25, theta=0.5pi and phi has two points in the range [0.5pi, 3/2pi]" ) {
         std::string const expression                = "x := 2.25*sin(theta)*cos(phi); y := 2.25*sin(theta)*sin(phi); z := 2.25*cos(theta);";
         std::array<ParametricVariable, 2> variables = { ParametricVariable( "theta", M_PI_2, M_PI_2, 1 ),
                                                         ParametricVariable( "phi", M_PI_2, m_3_PI_2, 2 ) };
         std::vector<std::array<double, 3>> const interface_coordinates( ComputeInterfaceCoordinates( expression, var_names, variables ) );

         THEN( "The number of interface coordinates should equal 2" ) {
            REQUIRE( interface_coordinates.size() == 2 );
         }
         THEN( "The points should equal {0,+-r,0}" ) {
            for( auto const& coord : interface_coordinates ) {
               REQUIRE( coord[0] == Approx( 0.0 ).margin( 1e-15 ) );
               REQUIRE( std::abs( coord[1] ) == Approx( 2.25 ) );
               REQUIRE( coord[2] == Approx( 0.0 ).margin( 1e-15 ) );
            }
         }
      }
   }
}

SCENARIO( "The parametric initializer works properly for a quasi 1D interface.", "[.slow1rank]" ) {
   GIVEN( "A topology with: node_size = 1.0, maximum level = 0, n_materials = 2. \n"
          "The expression parameters are in the range [0,1] and contain equal number of points to match the cell centers on level 0." ) {
      // General data
      constexpr double node_size                        = 1.0;
      constexpr unsigned int max_level                  = 0;
      std::vector<std::vector<nid_t>> const level_nodes = TestInitialCondition::GetChildNodesForSingleRootNode( max_level );
      std::vector<MaterialName> const materials         = { MaterialName::MaterialOne, MaterialName::MaterialTwo };
      double const cell_size_finest_level               = CellSizeOfLevel( node_size, max_level );
      // Parametric data
      std::uint64_t const number_of_points        = CC::TCX() + 1;
      double const start                          = -( double( CC::HS() ) - 0.5 ) * cell_size_finest_level;
      double const end                            = start + double( CC::TCX() ) * cell_size_finest_level;
      std::array<ParametricVariable, 2> variables = { ParametricVariable( "r", start, end, number_of_points ),
                                                      ParametricVariable( "s", start, end, number_of_points ) };
      WHEN( "The interface lies at x=0.5, the reference points is at {1.5, 0.5, 0.5} (outside of the domain) and no bounding box is given." ) {
         // Parametric expression data
         std::string const expression                            = "x := 0.5; y := r; z := s";
         std::array<double, 3> const ref_point                   = { 1.5, 0.5, 0.5 };
         std::vector<std::array<double, 6>> const bounding_boxes = {};
         // Create the initializer
         std::unique_ptr<LevelsetInitializer> initializer( std::make_unique<ParametricLevelsetInitializer>( expression, variables, ref_point, bounding_boxes,
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
      WHEN( "The interface lies at y=0.5, the reference points is at {0.5,1.5,0.5} (outside of the domain) and no bounding box is given." ) {
         // Parametric expression data
         std::string const expression                            = "x := r; y := 0.5; z := s";
         std::array<double, 3> const ref_point                   = { 0.5, 1.5, 0.5 };
         std::vector<std::array<double, 6>> const bounding_boxes = {};
         // Create the initializer
         std::unique_ptr<LevelsetInitializer> initializer( std::make_unique<ParametricLevelsetInitializer>( expression, variables, ref_point, bounding_boxes,
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

      WHEN( "The interface lies at z=0.5, the reference points is at {0.5,0.5,-1.5} (outside of the domain) and no bounding box is given." ) {
         // Parametric expression data
         std::string const expression                            = "x := r; y := s; z := 0.5";
         std::array<double, 3> const ref_point                   = { 0.5, 0.5, -1.5 };
         std::vector<std::array<double, 6>> const bounding_boxes = {};
         // Create the initializer
         std::unique_ptr<LevelsetInitializer> initializer( std::make_unique<ParametricLevelsetInitializer>( expression, variables, ref_point, bounding_boxes,
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

SCENARIO( "The parametric initializer works properly for a 3D interface.", "[.slow1rank]" ) {
   GIVEN( "A topology with: node_size = 1.0, maximum level = 1, n_materials = 2. \n"
          "The expression parameters are in the range [0,1] and the number of points are two times the number of cell corners on the finest level." ) {
      // General data
      constexpr double node_size                        = 1.0;
      constexpr unsigned int max_level                  = 1;
      std::vector<std::vector<nid_t>> const level_nodes = TestInitialCondition::GetChildNodesForSingleRootNode( max_level );
      std::vector<MaterialName> const materials         = { MaterialName::MaterialOne, MaterialName::MaterialTwo };
      double const cell_size_finest_level               = CellSizeOfLevel( node_size, max_level );
      // Parametric data
      std::uint64_t const number_of_points        = static_cast<std::uint64_t>( round( 1.0 / cell_size_finest_level ) ) * 2 + 1;
      std::array<ParametricVariable, 2> variables = { ParametricVariable( "r", 0, 1, number_of_points ),
                                                      ParametricVariable( "s", 0, 1, number_of_points ) };
      WHEN( "The expression is an ellipsoid with semi axis {r_x, r_y, r_z} = {0.125, 0.1, 0.05},\n"
            " the reference point and center are at {0.25,0.25,0.25}\n"
            " and the bounding box is set to {0.0, 0.5, 0.05, 0.45, 0.1, 0.4}" ) {
         // Parametric expression data
         std::string const expression = "x := 0.125 * sin(" + string_pi + " * s) * cos(" + string_2_pi + " * r) + 0.25;"
                                                                                                         "y := 0.1 * sin(" +
                                        string_pi + " * s) * sin(" + string_2_pi + " * r) + 0.25;"
                                                                                   "z := 0.05 * cos(" +
                                        string_pi + " * s) + 0.25";
         std::array<double, 3> const ref_point                   = { 0.25, 0.25, 0.25 };
         std::vector<std::array<double, 6>> const bounding_boxes = { { 0.0, 0.5, 0.0, 0.5, 0.0, 0.5 } };
         // Create the initializer
         std::unique_ptr<LevelsetInitializer> initializer( std::make_unique<ParametricLevelsetInitializer>( expression, variables, ref_point, bounding_boxes,
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
               "In the south-west-bottom the negative values are below the maximum radius/cell_size." ) {
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
