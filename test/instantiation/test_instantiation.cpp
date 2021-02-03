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
#include <fakeit.hpp>
#include <string>

#include "instantiation/materials/instantiation_material_manager.h"
#include "instantiation/topology/instantiation_topology_manager.h"
#include "instantiation/instantiation_initial_condition.h"
#include "instantiation/topology/instantiation_tree.h"
#include "materials/material_property_definitions.h"

namespace InstantiationConstants {
   constexpr unsigned int number_of_materials                          = 2;
   constexpr double first_material_bulk_viscosity                      = 4.0;
   constexpr double second_material_shear_viscosity                    = 5.0;
   constexpr double thermal_conductivity                               = 0.0;
   constexpr double specific_heat_capacity                             = 0.0;
   constexpr double reference_value                                    = 1.0;
   constexpr double gamma                                              = 1.4;
   constexpr double surface_tension_coefficient                        = 3.0;
   constexpr double node_size                                          = 1.0;
   constexpr double density_material_one_scenario_one                  = 0.4;
   constexpr double density_material_two_scenario_one                  = 0.9;
   constexpr double density_material_one_scenario_two                  = 0.7;
   constexpr double pressure_material_two_scenario_one                 = 0.5;
   constexpr double pressure_scenario_two                              = 1.2;
   std::array<unsigned int, 3> const number_of_nodes_on_level_zero_one = { 1, 1, 1 };
   std::array<unsigned int, 3> const number_of_nodes_on_level_zero_two = { 3, 2, 5 };
   nid_t const root_node_id                                            = IdSeed();
}// namespace InstantiationConstants

namespace {
   using namespace fakeit;
   std::string const scenario_one_material_one_density  = "density := " + std::to_string( InstantiationConstants::density_material_one_scenario_one ) + ";";
   std::string const scenario_one_material_one_velocity = "velocityX := 0.1;"
                                                          "velocityY := 0.0;"
                                                          "velocityZ := 0.0;";
   std::string const scenario_one_material_one_pressure = "pressure := pow(x,2);";
   std::string const scenario_one_material_two_density  = "density := " + std::to_string( InstantiationConstants::density_material_two_scenario_one ) + ";";
   std::string const scenario_one_material_two_velocity = "velocityX := 0.0;"
                                                          "velocityY := 0.0;"
                                                          "velocityZ := 0.0;";
   std::string const scenario_one_material_two_pressure = "pressure := " + std::to_string( InstantiationConstants::pressure_material_two_scenario_one ) + ";";

   std::string const scenario_two_material_one_density  = "density := " + std::to_string( InstantiationConstants::density_material_one_scenario_two ) + ";";
   std::string const scenario_two_material_one_velocity = "velocityX := 0.0;"
                                                          "velocityY := 0.0;"
                                                          "velocityZ := 0.0;";
   std::string const scenario_two_pressure              = "pressure := " + std::to_string( InstantiationConstants::pressure_scenario_two ) + ";";
   std::string const scenario_two_material_two_density  = "density := 1.0;";
   std::string const scenario_two_material_two_velocity = "velocityX := 0.0;"
                                                          "velocityY := 0.0;"
                                                          "velocityZ := 0.0;";

   std::string const scenario_one_levelset = "phi := 0.5-x;";
   std::string const scenario_two_levelset = "phi := 1.0;";

   constexpr MaterialName material_one = MaterialName::MaterialOne;
   constexpr MaterialName material_two = MaterialName::MaterialTwo;
   UnitHandler const unit_handler( InstantiationConstants::reference_value, InstantiationConstants::reference_value, InstantiationConstants::reference_value, InstantiationConstants::reference_value );

   Mock<MaterialReader> InstantiateMaterialReader( double const first_material_shear_viscosity = 0.0, double const second_material_bulk_viscosity = 0.0 ) {
      Mock<MaterialReader> material_reader;

      std::unordered_map<std::string, double> const eos_gas_data   = { { "gamma", InstantiationConstants::gamma }, { "backgroundPressure", 0.0 } };
      std::unordered_map<std::string, double> const eos_water_data = { { "gamma", InstantiationConstants::gamma }, { "A", 1.0 }, { "B", 1.0 }, { "rho0", 1.0 } };
      EquationOfStateName const eos_gas_name( EquationOfStateName::StiffenedGas );
      EquationOfStateName const eos_water_name( EquationOfStateName::WaterlikeFluid );

      When( Method( material_reader, ReadNumberOfMaterials ) ).AlwaysReturn( 2 );
      When( Method( material_reader, ReadMaterialType ) ).AlwaysReturn( MaterialType::Fluid );
      When( Method( material_reader, ReadEquationOfStateName ).Using( 2 ) ).AlwaysReturn( eos_water_name );
      When( Method( material_reader, ReadEquationOfStateData ).Using( 2 ) ).AlwaysReturn( eos_water_data );
      When( Method( material_reader, ReadEquationOfStateName ).Using( 1 ) ).AlwaysReturn( eos_gas_name );
      When( Method( material_reader, ReadEquationOfStateData ).Using( 1 ) ).AlwaysReturn( eos_gas_data );
      When( Method( material_reader, ReadFixedValue ).Using( { 1 }, MaterialProperty::BulkViscosity ) ).AlwaysReturn( InstantiationConstants::first_material_bulk_viscosity );
      When( Method( material_reader, ReadFixedValue ).Using( { 2 }, MaterialProperty::ShearViscosity ) ).AlwaysReturn( InstantiationConstants::second_material_shear_viscosity );
      When( Method( material_reader, ReadFixedValue ).Using( { 1, 2 }, MaterialProperty::SurfaceTensionCoefficient ) ).AlwaysReturn( InstantiationConstants::surface_tension_coefficient );
      When( Method( material_reader, ReadFixedValue ).Using( { 1 }, MaterialProperty::ShearViscosity ) ).AlwaysReturn( first_material_shear_viscosity );
      When( Method( material_reader, ReadFixedValue ).Using( { 2 }, MaterialProperty::BulkViscosity ) ).AlwaysReturn( second_material_bulk_viscosity );
      When( Method( material_reader, ReadFixedValue ).Using( { 1 }, MaterialProperty::ThermalConductivity ) ).AlwaysReturn( InstantiationConstants::thermal_conductivity );
      When( Method( material_reader, ReadFixedValue ).Using( { 2 }, MaterialProperty::ThermalConductivity ) ).AlwaysReturn( InstantiationConstants::thermal_conductivity );
      When( Method( material_reader, ReadFixedValue ).Using( { 1 }, MaterialProperty::SpecificHeatCapacity ) ).AlwaysReturn( InstantiationConstants::specific_heat_capacity );
      When( Method( material_reader, ReadFixedValue ).Using( { 2 }, MaterialProperty::SpecificHeatCapacity ) ).AlwaysReturn( InstantiationConstants::specific_heat_capacity );

      return material_reader;
   }

   Mock<MultiResolutionReader> InstantiateMultiResolutionReader( unsigned int const maximum_level, std::array<unsigned int, 3> const number_of_nodes, double const node_size = InstantiationConstants::node_size ) {
      Mock<MultiResolutionReader> multiresolution_reader;

      When( Method( multiresolution_reader, ReadNumberOfNodes ).Using( Direction::X ) ).AlwaysReturn( number_of_nodes[0] );
      When( Method( multiresolution_reader, ReadNumberOfNodes ).Using( Direction::Y ) ).AlwaysReturn( number_of_nodes[1] );
      When( Method( multiresolution_reader, ReadNumberOfNodes ).Using( Direction::Z ) ).AlwaysReturn( number_of_nodes[2] );
      When( Method( multiresolution_reader, ReadMaximumLevel ) ).Return( maximum_level );
      When( Method( multiresolution_reader, ReadNodeSizeOnLevelZero ) ).Return( node_size );

      return multiresolution_reader;
   }

   Mock<BoundaryConditionReader> InstantiateBoundaryConditionReader() {
      Mock<BoundaryConditionReader> boundary_condition_reader;

      When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::East ) ).AlwaysReturn( LevelSetBoundaryType::Symmetry );
      When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::West ) ).AlwaysReturn( LevelSetBoundaryType::Symmetry );
      When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::North ) ).AlwaysReturn( LevelSetBoundaryType::Symmetry );
      When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::South ) ).AlwaysReturn( LevelSetBoundaryType::Symmetry );
      When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::Bottom ) ).AlwaysReturn( LevelSetBoundaryType::Symmetry );
      When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::Top ) ).AlwaysReturn( LevelSetBoundaryType::Symmetry );

      When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::East ) ).AlwaysReturn( MaterialBoundaryType::Wall );
      When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::West ) ).AlwaysReturn( MaterialBoundaryType::ZeroGradient );
      When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::North ) ).AlwaysReturn( MaterialBoundaryType::ZeroGradient );
      When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::South ) ).AlwaysReturn( MaterialBoundaryType::ZeroGradient );
      When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::Bottom ) ).AlwaysReturn( MaterialBoundaryType::Symmetry );
      When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::Top ) ).AlwaysReturn( MaterialBoundaryType::Symmetry );

      return boundary_condition_reader;
   }

   Mock<InitialConditionReader> InstantiateInitialConditionReader( std::string const density_one, std::string const density_two,
                                                                   std::string const pressure_one, std::string const pressure_two,
                                                                   std::string const velocity_one, std::string const velocity_two, std::string const levelset ) {
      Mock<InitialConditionReader> initial_condition_reader;

      When( Method( initial_condition_reader, ReadMaterialInitialConditions ).Using( 1 ) ).Return( density_one + velocity_one + pressure_one );
      When( Method( initial_condition_reader, ReadMaterialInitialConditions ).Using( 2 ) ).Return( density_two + velocity_two + pressure_two );
      When( Method( initial_condition_reader, ReadLevelsetInitialConditions ).Using( 1 ) ).Return( levelset );

      return initial_condition_reader;
   }

   Mock<InputReader> InstantiateInputReaderForMaterialManager( Mock<MaterialReader>& material_reader ) {
      Mock<InputReader> input_reader;
      When( Method( input_reader, GetMaterialReader ) ).AlwaysReturn( material_reader.get() );

      return input_reader;
   }

   Mock<InputReader> InstantiateInputReaderForTopologyManager( Mock<MultiResolutionReader>& multiresolution_reader,
                                                               Mock<BoundaryConditionReader>& boundary_condition_reader ) {
      Mock<InputReader> input_reader;

      When( Method( input_reader, GetMultiResolutionReader ) ).AlwaysReturn( multiresolution_reader.get() );
      When( Method( input_reader, GetBoundaryConditionReader ) ).AlwaysReturn( boundary_condition_reader.get() );

      return input_reader;
   }

   Mock<InputReader> InstantiateInputReaderForInitialCondition( Mock<InitialConditionReader>& initial_condition_reader ) {
      Mock<InputReader> input_reader;

      When( Method( input_reader, GetInitialConditionReader ) ).AlwaysReturn( initial_condition_reader.get() );

      return input_reader;
   }

}// namespace

SCENARIO( "Instantiation of the material manager", "[1rank]" ) {
   GIVEN( "Mocks of material reader and input reader" ) {
      constexpr double shear_viscosity_reader_one_and_two_material_two = 1.1;
      constexpr double reader_three_shear_viscosity                    = 2.2;
      constexpr double bulk_viscosity_reader_one_and_two_material_one  = 1.7;
      constexpr double reader_three_bulk_viscosity                     = 3.3;
      Mock<MaterialReader> first_material_reader( InstantiateMaterialReader( shear_viscosity_reader_one_and_two_material_two, bulk_viscosity_reader_one_and_two_material_one ) );
      Mock<MaterialReader> second_material_reader( InstantiateMaterialReader( shear_viscosity_reader_one_and_two_material_two, bulk_viscosity_reader_one_and_two_material_one ) );
      Mock<MaterialReader> third_material_reader( InstantiateMaterialReader( reader_three_shear_viscosity, reader_three_bulk_viscosity ) );

      Mock<InputReader> first_input_reader  = InstantiateInputReaderForMaterialManager( first_material_reader );
      Mock<InputReader> second_input_reader = InstantiateInputReaderForMaterialManager( second_material_reader );
      Mock<InputReader> third_input_reader  = InstantiateInputReaderForMaterialManager( third_material_reader );

      WHEN( "A material manager with 2 materials is instantiated using mock readers" ) {
         MaterialManager const material_manager( Instantiation::InstantiateMaterialManager( first_input_reader.get(), unit_handler ) );
         THEN( "The values returned by the material manager are equal to expected values" ) {
            REQUIRE( material_manager.GetNumberOfMaterials() == InstantiationConstants::number_of_materials );
            REQUIRE( material_manager.GetMaterial( material_one ).GetShearViscosity() == shear_viscosity_reader_one_and_two_material_two );
            REQUIRE( material_manager.GetMaterial( material_one ).GetBulkViscosity() == InstantiationConstants::first_material_bulk_viscosity );
            REQUIRE( material_manager.GetMaterial( material_one ).GetEquationOfState().Gamma() == InstantiationConstants::gamma );
            REQUIRE( material_manager.GetMaterialPairing( material_one, material_two ).GetSurfaceTensionCoefficient() == InstantiationConstants::surface_tension_coefficient );
         }
      }

      WHEN( "The material values of the identical material managers are compared" ) {
         MaterialManager const first_material_manager( Instantiation::InstantiateMaterialManager( first_input_reader.get(), unit_handler ) );
         MaterialManager const second_material_manager( Instantiation::InstantiateMaterialManager( second_input_reader.get(), unit_handler ) );
         THEN( "The values of both material managers are equal" ) {
            REQUIRE( first_material_manager.GetNumberOfMaterials() == second_material_manager.GetNumberOfMaterials() );
            REQUIRE( first_material_manager.GetMaterial( material_one ).GetShearViscosity() == second_material_manager.GetMaterial( material_one ).GetShearViscosity() );
            REQUIRE( first_material_manager.GetMaterial( material_one ).GetBulkViscosity() == second_material_manager.GetMaterial( material_one ).GetBulkViscosity() );
            REQUIRE( first_material_manager.GetMaterial( material_one ).GetEquationOfState().Gamma() == second_material_manager.GetMaterial( material_one ).GetEquationOfState().Gamma() );
            REQUIRE( first_material_manager.GetMaterialPairing( material_one, material_two ).GetSurfaceTensionCoefficient() == second_material_manager.GetMaterialPairing( material_one, material_two ).GetSurfaceTensionCoefficient() );
         }
      }

      WHEN( "The shear and bulk viscosity values of first and third unidentical managers are compared" ) {
         MaterialManager const first_material_manager( Instantiation::InstantiateMaterialManager( first_input_reader.get(), unit_handler ) );
         MaterialManager const second_material_manager( Instantiation::InstantiateMaterialManager( third_input_reader.get(), unit_handler ) );
         THEN( "The values returned by both material managers are not equal" ) {
            REQUIRE_FALSE( first_material_manager.GetMaterial( material_one ).GetShearViscosity() == second_material_manager.GetMaterial( material_one ).GetShearViscosity() );
            REQUIRE_FALSE( first_material_manager.GetMaterial( material_two ).GetBulkViscosity() == second_material_manager.GetMaterial( material_two ).GetBulkViscosity() );
         }
      }
   }
}

SCENARIO( "Instantiation of topology manager", "[1rank]" ) {
   GIVEN( "Mocks of material, multiresolution, boundary condition and input readers and material manager" ) {
      Mock<MaterialReader> material_reader( InstantiateMaterialReader() );
      constexpr unsigned int maximum_level_reader_one_and_two = 3;
      constexpr unsigned int maximum_level_reader_three       = 1;
      constexpr double node_size_reader_three                 = 0.125;
      Mock<MultiResolutionReader> first_multiresolution( InstantiateMultiResolutionReader( maximum_level_reader_one_and_two, InstantiationConstants::number_of_nodes_on_level_zero_one ) );
      Mock<MultiResolutionReader> second_multiresolution( InstantiateMultiResolutionReader( maximum_level_reader_one_and_two, InstantiationConstants::number_of_nodes_on_level_zero_one ) );
      Mock<MultiResolutionReader> third_multiresolution( InstantiateMultiResolutionReader( maximum_level_reader_three, InstantiationConstants::number_of_nodes_on_level_zero_two, node_size_reader_three ) );
      Mock<BoundaryConditionReader> boundary_condition_reader( InstantiateBoundaryConditionReader() );

      Mock<InputReader> input_reader_for_material_manager( InstantiateInputReaderForMaterialManager( material_reader ) );
      Mock<InputReader> first_input_reader( InstantiateInputReaderForTopologyManager( first_multiresolution, boundary_condition_reader ) );
      Mock<InputReader> second_input_reader( InstantiateInputReaderForTopologyManager( second_multiresolution, boundary_condition_reader ) );
      Mock<InputReader> third_input_reader( InstantiateInputReaderForTopologyManager( third_multiresolution, boundary_condition_reader ) );

      MaterialManager const material_manager( Instantiation::InstantiateMaterialManager( input_reader_for_material_manager.get(), unit_handler ) );

      WHEN( "The topology manager is instantiated using the given mocks and material manager" ) {
         unsigned int const expected_maximum_level = 3;
         TopologyManager topology_manager( Instantiation::InstantiateTopologyManager( first_input_reader.get(), material_manager ) );
         THEN( "The returned values are equal to the expected values" ) {
            REQUIRE( topology_manager.GetMaximumLevel() == expected_maximum_level );
            REQUIRE( topology_manager.GetNumberOfNodesOnLevelZero() == InstantiationConstants::number_of_nodes_on_level_zero_one );
         }
      }

      WHEN( "Two topology managers are instantiated with identical initial conditions" ) {
         TopologyManager first_topology_manager( Instantiation::InstantiateTopologyManager( first_input_reader.get(), material_manager ) );
         TopologyManager second_topology_manager( Instantiation::InstantiateTopologyManager( second_input_reader.get(), material_manager ) );
         THEN( "The values returned by the both managers are equal" ) {
            REQUIRE( first_topology_manager.GetMaximumLevel() == second_topology_manager.GetMaximumLevel() );
            REQUIRE( first_topology_manager.GetNumberOfNodesOnLevelZero() == second_topology_manager.GetNumberOfNodesOnLevelZero() );
            REQUIRE( first_topology_manager.IsExternalTopologyBoundary( BoundaryLocation::East, InstantiationConstants::root_node_id ) == second_topology_manager.IsExternalTopologyBoundary( BoundaryLocation::East, InstantiationConstants::root_node_id ) );
         }
      }

      WHEN( "Two topology managers are instantiated with different initial conditions" ) {
         TopologyManager first_topology_manager( Instantiation::InstantiateTopologyManager( first_input_reader.get(), material_manager ) );
         TopologyManager second_topology_manager( Instantiation::InstantiateTopologyManager( third_input_reader.get(), material_manager ) );
         THEN( "The values returned by both managers are not equal" ) {
            REQUIRE_FALSE( first_topology_manager.GetMaximumLevel() == second_topology_manager.GetMaximumLevel() );
            REQUIRE( second_topology_manager.GetNumberOfNodesOnLevelZero() == InstantiationConstants::number_of_nodes_on_level_zero_two );
            REQUIRE( !second_topology_manager.IsExternalTopologyBoundary( BoundaryLocation::TopEast, InstantiationConstants::root_node_id ) );
         }
      }
   }
}

SCENARIO( "Instantiation of Tree", "[1rank]" ) {
   GIVEN( "Mocks of input reader, instance of material manager and topology manager" ) {
      Mock<MaterialReader> material_reader( InstantiateMaterialReader() );
      constexpr unsigned int maximum_level_reader_one_and_two = 3;
      constexpr unsigned int maximum_level_reader_three       = 1;
      constexpr double node_size_reader_three                 = 0.125;
      Mock<MultiResolutionReader> first_multiresolution( InstantiateMultiResolutionReader( maximum_level_reader_one_and_two, InstantiationConstants::number_of_nodes_on_level_zero_one ) );
      Mock<MultiResolutionReader> second_multiresolution( InstantiateMultiResolutionReader( maximum_level_reader_one_and_two, InstantiationConstants::number_of_nodes_on_level_zero_one ) );
      Mock<MultiResolutionReader> third_multiresolution( InstantiateMultiResolutionReader( maximum_level_reader_three, InstantiationConstants::number_of_nodes_on_level_zero_two, node_size_reader_three ) );
      Mock<BoundaryConditionReader> boundary_condition_reader( InstantiateBoundaryConditionReader() );

      Mock<InputReader> input_reader_for_material_manager( InstantiateInputReaderForMaterialManager( material_reader ) );
      Mock<InputReader> first_input_reader( InstantiateInputReaderForTopologyManager( first_multiresolution, boundary_condition_reader ) );
      Mock<InputReader> second_input_reader( InstantiateInputReaderForTopologyManager( second_multiresolution, boundary_condition_reader ) );
      Mock<InputReader> third_input_reader( InstantiateInputReaderForTopologyManager( third_multiresolution, boundary_condition_reader ) );

      MaterialManager const material_manager( Instantiation::InstantiateMaterialManager( input_reader_for_material_manager.get(), unit_handler ) );
      TopologyManager topology_manager( Instantiation::InstantiateTopologyManager( first_input_reader.get(), material_manager ) );

      WHEN( "The tree is instantiated using the input reader mock" ) {
         unsigned int const expected_nodes_on_level_zero = 0;
         Tree tree( Instantiation::InstantiateTree( first_input_reader.get(), topology_manager, unit_handler ) );
         THEN( "The obtained values are equal to the expected values" ) {
            REQUIRE( tree.GetNodeSizeOnLevelZero() == InstantiationConstants::node_size );
            REQUIRE( tree.NodesOnLevel( 0 ).size() == expected_nodes_on_level_zero );
         }
      }

      WHEN( "Two trees are instantiated using identical initial conditions" ) {
         Tree first_tree( Instantiation::InstantiateTree( first_input_reader.get(), topology_manager, unit_handler ) );
         Tree second_tree( Instantiation::InstantiateTree( second_input_reader.get(), topology_manager, unit_handler ) );
         THEN( "The two trees have same node size and number of nodes on level 0" ) {
            REQUIRE( first_tree.GetNodeSizeOnLevelZero() == second_tree.GetNodeSizeOnLevelZero() );
            REQUIRE( first_tree.NodesOnLevel( 0 ).size() == second_tree.NodesOnLevel( 0 ).size() );
         }
      }

      WHEN( "Two trees are instantiated using different node size on level zero" ) {
         Tree first_tree( Instantiation::InstantiateTree( first_input_reader.get(), topology_manager, unit_handler ) );
         Tree second_tree( Instantiation::InstantiateTree( third_input_reader.get(), topology_manager, unit_handler ) );
         THEN( "Node size returned by both trees are not equal but nodes on level 0 are 0, since no nodes are created yet" ) {
            REQUIRE_FALSE( first_tree.GetNodeSizeOnLevelZero() == second_tree.GetNodeSizeOnLevelZero() );
            REQUIRE( first_tree.NodesOnLevel( 0 ).size() == second_tree.NodesOnLevel( 0 ).size() );
         }
      }
   }
}

SCENARIO( "Instantiation of initial condition", "[1rank]" ) {
   //inital primestate are just internal cells (index shifting needed)
   constexpr std::size_t FIX = CC::FICX() - CC::HSSX();
   constexpr std::size_t FIY = CC::FICY() - CC::HSSY();
   constexpr std::size_t FIZ = CC::FICZ() - CC::HSSZ();
   constexpr std::size_t LIX = CC::LICX() - CC::HSSX();
   constexpr std::size_t LIY = CC::LICX() - CC::HSSY();
   constexpr std::size_t LIZ = CC::LICX() - CC::HSSZ();
   GIVEN( "Instances of material manager, topology manager and tree using mock readers" ) {
      Mock<MaterialReader> material_reader( InstantiateMaterialReader() );
      Mock<MultiResolutionReader> multiresolution( InstantiateMultiResolutionReader( 0, InstantiationConstants::number_of_nodes_on_level_zero_two ) );
      Mock<BoundaryConditionReader> boundary_condition_reader( InstantiateBoundaryConditionReader() );
      Mock<InitialConditionReader> first_initial_condition_reader( InstantiateInitialConditionReader( scenario_one_material_one_density, scenario_one_material_one_velocity,
                                                                                                      scenario_one_material_one_pressure, scenario_one_material_two_density,
                                                                                                      scenario_one_material_two_velocity, scenario_one_material_two_pressure, scenario_one_levelset ) );
      Mock<InitialConditionReader> second_initial_condition_reader( InstantiateInitialConditionReader( scenario_one_material_one_density, scenario_one_material_one_velocity,
                                                                                                       scenario_one_material_one_pressure, scenario_one_material_two_density,
                                                                                                       scenario_one_material_two_velocity, scenario_one_material_two_pressure, scenario_one_levelset ) );
      Mock<InitialConditionReader> third_initial_condition_reader( InstantiateInitialConditionReader( scenario_two_material_one_density, scenario_two_material_one_velocity,
                                                                                                      scenario_two_pressure, scenario_two_material_two_density,
                                                                                                      scenario_two_material_two_velocity, scenario_two_pressure, scenario_two_levelset ) );

      Mock<InputReader> input_reader_for_material_manager( InstantiateInputReaderForMaterialManager( material_reader ) );
      Mock<InputReader> input_reader_for_topology_manager( InstantiateInputReaderForTopologyManager( multiresolution, boundary_condition_reader ) );
      Mock<InputReader> first_input_reader( InstantiateInputReaderForInitialCondition( first_initial_condition_reader ) );
      Mock<InputReader> second_input_reader( InstantiateInputReaderForInitialCondition( second_initial_condition_reader ) );
      Mock<InputReader> third_input_reader( InstantiateInputReaderForInitialCondition( third_initial_condition_reader ) );

      MaterialManager const material_manager( Instantiation::InstantiateMaterialManager( input_reader_for_material_manager.get(), unit_handler ) );
      TopologyManager topology_manager( Instantiation::InstantiateTopologyManager( input_reader_for_topology_manager.get(), material_manager ) );
      Tree tree( Instantiation::InstantiateTree( input_reader_for_topology_manager.get(), topology_manager, unit_handler ) );

      WHEN( "An initial condition is instantiated using given mock input reader" ) {
         InitialCondition initial_condition( Instantiation::InstantiateInitialCondition( first_input_reader.get(), topology_manager, tree, material_manager, unit_handler ) );
         double initial_prime_states_material_one[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()];
         initial_condition.GetInitialPrimeStates( InstantiationConstants::root_node_id, material_one, initial_prime_states_material_one );
         double initial_prime_states_material_two[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()];
         initial_condition.GetInitialPrimeStates( InstantiationConstants::root_node_id, material_two, initial_prime_states_material_two );
         THEN( "The initial material and prime state values return by the initial condition are equal to the expected values" ) {
            REQUIRE( initial_condition.GetInitialMaterials( InstantiationConstants::root_node_id ).size() == InstantiationConstants::number_of_materials );
            constexpr double expected_pressure_material_one = ( 0.5 * ( 1.0 / CC::ICX() ) ) * ( 0.5 * ( 1.0 / CC::ICX() ) );
            REQUIRE( initial_prime_states_material_one[PTI( PrimeState::Density )][FIX][FIY][FIZ] == Approx( InstantiationConstants::density_material_one_scenario_one ) );
            REQUIRE( initial_prime_states_material_two[PTI( PrimeState::Density )][LIX][LIY][LIZ] == Approx( InstantiationConstants::density_material_two_scenario_one ) );
            REQUIRE( initial_prime_states_material_one[PTI( PrimeState::Pressure )][FIX][FIY][FIZ] == Approx( expected_pressure_material_one ) );
            REQUIRE( initial_prime_states_material_two[PTI( PrimeState::Pressure )][LIX][LIY][LIZ] == Approx( InstantiationConstants::pressure_material_two_scenario_one ) );
         }
      }

      WHEN( "Two initial conditions are instantiated with identical initial conditions" ) {
         InitialCondition first_initial_condition( Instantiation::InstantiateInitialCondition( first_input_reader.get(), topology_manager, tree, material_manager, unit_handler ) );
         InitialCondition second_initial_condition( Instantiation::InstantiateInitialCondition( second_input_reader.get(), topology_manager, tree, material_manager, unit_handler ) );
         double levelset_temp_one[CC::TCX()][CC::TCY()][CC::TCZ()];
         double levelset_temp_two[CC::TCX()][CC::TCY()][CC::TCZ()];
         double initial_prime_states_one[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()];
         double initial_prime_states_two[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()];

         first_initial_condition.GetInitialLevelset( InstantiationConstants::root_node_id, levelset_temp_one );
         second_initial_condition.GetInitialLevelset( InstantiationConstants::root_node_id, levelset_temp_two );

         first_initial_condition.GetInitialPrimeStates( InstantiationConstants::root_node_id, material_one, initial_prime_states_one );
         second_initial_condition.GetInitialPrimeStates( InstantiationConstants::root_node_id, material_one, initial_prime_states_two );
         THEN( "The number of materials, material and levelset values returned by both initial conditions are equal" ) {

            REQUIRE( first_initial_condition.GetInitialMaterials( InstantiationConstants::root_node_id ).size() == second_initial_condition.GetInitialMaterials( InstantiationConstants::root_node_id ).size() );

            for( unsigned int i = 0; i < CC::TCX(); i++ ) {
               for( unsigned int j = 0; j < CC::TCY(); j++ ) {
                  for( unsigned int k = 0; k < CC::TCZ(); k++ ) {
                     REQUIRE( levelset_temp_one[i][j][k] == levelset_temp_two[i][j][k] );
                  }
               }
            }
            REQUIRE( initial_prime_states_one[PTI( PrimeState::Density )][FIX][FIY][FIZ] == initial_prime_states_two[PTI( PrimeState::Density )][FIX][FIY][FIZ] );
            REQUIRE( initial_prime_states_one[PTI( PrimeState::Density )][LIX][LIY][LIZ] == initial_prime_states_two[PTI( PrimeState::Density )][LIX][LIY][LIZ] );
            REQUIRE( initial_prime_states_one[PTI( PrimeState::Pressure )][FIX][FIY][FIZ] == initial_prime_states_two[PTI( PrimeState::Pressure )][FIX][FIY][FIZ] );
            REQUIRE( initial_prime_states_one[PTI( PrimeState::Pressure )][LIX][LIY][LIZ] == initial_prime_states_two[PTI( PrimeState::Pressure )][LIX][LIY][LIZ] );
         }
      }

      WHEN( "Two initial conditions are instantiated with different initial conditions" ) {
         InitialCondition first_initial_condition( Instantiation::InstantiateInitialCondition( first_input_reader.get(), topology_manager, tree, material_manager, unit_handler ) );
         InitialCondition second_initial_condition( Instantiation::InstantiateInitialCondition( third_input_reader.get(), topology_manager, tree, material_manager, unit_handler ) );
         double levelset_temp_one[CC::TCX()][CC::TCY()][CC::TCZ()];
         double levelset_temp_two[CC::TCX()][CC::TCY()][CC::TCZ()];
         double initial_prime_states_one[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()];
         double initial_prime_states_two[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()];

         first_initial_condition.GetInitialLevelset( InstantiationConstants::root_node_id, levelset_temp_one );
         second_initial_condition.GetInitialLevelset( InstantiationConstants::root_node_id, levelset_temp_two );

         first_initial_condition.GetInitialPrimeStates( InstantiationConstants::root_node_id, material_one, initial_prime_states_one );
         second_initial_condition.GetInitialPrimeStates( InstantiationConstants::root_node_id, material_one, initial_prime_states_two );
         THEN( "The material and levelset values returned by both initial conditions are not equal" ) {
            for( unsigned int i = 0; i < CC::TCX(); i++ ) {
               for( unsigned int j = 0; j < CC::TCY(); j++ ) {
                  for( unsigned int k = 0; k < CC::TCZ(); k++ ) {
                     REQUIRE_FALSE( levelset_temp_one[i][j][k] == levelset_temp_two[i][j][k] );
                  }
               }
            }

            REQUIRE_FALSE( initial_prime_states_one[PTI( PrimeState::Density )][FIX][FIY][FIZ] == initial_prime_states_two[PTI( PrimeState::Density )][FIX][FIY][FIZ] );
            REQUIRE_FALSE( initial_prime_states_one[PTI( PrimeState::Pressure )][FIX][FIY][FIZ] == initial_prime_states_two[PTI( PrimeState::Pressure )][FIX][FIY][FIZ] );
            REQUIRE( initial_prime_states_two[PTI( PrimeState::Pressure )][FIX][FIY][FIZ] == Approx( InstantiationConstants::pressure_scenario_two ) );
            REQUIRE( initial_prime_states_two[PTI( PrimeState::Density )][FIX][FIY][FIZ] == Approx( InstantiationConstants::density_material_one_scenario_two ) );
         }
      }
   }
}
