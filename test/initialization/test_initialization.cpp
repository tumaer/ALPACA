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
* 5. ApprovalTests.cpp  : See LICENSE_APPROVAL_TESTS.txt for more information            *                                                                                        *
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
#include <ApprovalTests.hpp>
#include <algorithm>
#include <catch2/catch.hpp>
#include <fakeit.hpp>
#include <filesystem>

#include "instantiation/input_output/instantiation_restart_manager.h"
#include "instantiation/instantiation_communication_manager.h"
#include "instantiation/materials/instantiation_material_manager.h"
#include "instantiation/topology/instantiation_topology_manager.h"
#include "instantiation/halo_manager/instantiation_external_halo_manager.h"
#include "instantiation/halo_manager/instantiation_internal_halo_manager.h"
#include "instantiation/halo_manager/instantiation_halo_manager.h"
#include "instantiation/instantiation_multiresolution.h"
#include "instantiation/input_output/instantiation_output_writer.h"
#include "instantiation/input_output/instantiation_input_output_manager.h"
#include "instantiation/instantiation_modular_algorithm_assembler.h"
#include "instantiation/instantiation_initial_condition.h"
#include "instantiation/topology/instantiation_tree.h"
#include "unit_handler.h"
#include "input_output/hdf5/hdf5_manager.h"

namespace ApprovalBackend {
   /**
    * @brief Reads the content of a dataset of an HDF5 file into a vector and returns it.
    * @param file The location of the hdf file.
    * @param group_name The name of the group within the hdf file in which the dataset resides.
    * @param dataset_name The name of the dataset to be read.
    */
   std::vector<double> DatasetContent( std::string const file, std::string const group_name, std::string const dataset_name ) {
      Hdf5Manager& hdf5_manager( Hdf5Manager::Instance() );
      hdf5_manager.OpenFile( file, Hdf5Access::Read );
      hdf5_manager.OpenGroup( group_name );
      hsize_t const dataset_size = hdf5_manager.GetDatasetExtent( dataset_name );
      std::vector<double> dataset( dataset_size );
      hdf5_manager.ReadFullDataset( dataset_name, dataset.data() );
      hdf5_manager.CloseFile();
      return dataset;
   }

   /**
    * @brief Compares two datasets (with the same name in the same group) of two HDF5 files.
    * @param group_name The name of the group within the hdf file in which the dataset resides.
    * @param dataset_name The name of the dataset to be read.
    * @param file_a On of the HDF file to be used in the comparions.
    * @param file_b On of the HDF file to be used in the comparions.
    * @return True if the dataset differ. False otherwise.
    */
   bool DatasetContentDiffersInFiles( std::string const group_name, std::string const dataset_name, std::string const file_a, std::string const file_b ) {
      auto const dataset_of_file_a = DatasetContent( file_a, group_name, dataset_name );
      auto const dataset_of_file_b = DatasetContent( file_b, group_name, dataset_name );
      // Because we use trigonometric functions and such we need to include floating-point tolerances.
      return !std::equal( std::cbegin( dataset_of_file_a ), std::cend( dataset_of_file_a ), std::cbegin( dataset_of_file_b ), std::cend( dataset_of_file_b ), []( auto const& a, auto const& b ) { return a == Approx( b ).epsilon( ( std::numeric_limits<double>::epsilon() * 10 ) ); } );
   }

   /**
    * @brief Compares two datasets (with the same name in the same group) of two HDF5 files. If the dataset differ a message is logged.
    * @param group_name The name of the group within the hdf file in which the dataset resides.
    * @param dataset_name The name of the dataset to be read.
    * @param file_a On of the HDF file to be used in the comparions.
    * @param file_b On of the HDF file to be used in the comparions.
    * @return True if the dataset are identical. False otherwise.
    */
   bool FieldValuesAreEquivalentElseLog( std::string const value_name, std::string const group_name, std::string const dataset_name, std::string const file_a, std::string const file_b ) {
      if( DatasetContentDiffersInFiles( group_name, dataset_name, file_a, file_b ) ) {
         UNSCOPED_INFO( "-------> " + value_name + " values differ! <" + std::string( 38 - value_name.size(), '-' ) + "\n" );
         return false;
      }
      return true;
   }

   // We need a global switch for to stay competible with the ApprovalTests library in the override of the methods.
   static bool approving_restart = false;

   /**
    * @brief Custom comparator for HDF5 files to pass into the ApprovalTests library.
    */
   class H5FileComparator : public ApprovalTests::ApprovalComparator {
   public:
      /**
       * @brief Compares two (ALPACA output) HDF5 files for differences in the respective datasets.
       * @param received_file The received, ie freshly created file.
       * @param approved_file The approved, ie. ideal file.
       * @note Does not follow ALPACA naming convetion. Naming enforced by library.
       */
      bool contentsAreEquivalent( std::string const received_file, std::string const approved_file ) const override {
         bool result = true;

         if( approving_restart ) {
            result &= FieldValuesAreEquivalentElseLog( "CONSERVATIVES", "node_cell_data", "Conservatives", received_file, approved_file );
            result &= FieldValuesAreEquivalentElseLog( "LEVELSET", "node_cell_data", "Levelset", received_file, approved_file );
            result &= FieldValuesAreEquivalentElseLog( "PRIME STATES", "node_cell_data", "PrimeStates", received_file, approved_file );
         } else {
            result &= FieldValuesAreEquivalentElseLog( "DENSITY", "cell_data", "density", received_file, approved_file );
            result &= FieldValuesAreEquivalentElseLog( "INTERFACE_VELOCITY", "cell_data", "interface_velocity", received_file, approved_file );
            result &= FieldValuesAreEquivalentElseLog( "PRESSURE", "cell_data", "pressure", received_file, approved_file );
            result &= FieldValuesAreEquivalentElseLog( "VELOCITY", "cell_data", "velocity", received_file, approved_file );
            result &= FieldValuesAreEquivalentElseLog( "LEVELSET", "cell_data", "levelset", received_file, approved_file );
            result &= FieldValuesAreEquivalentElseLog( "CELL_VERTEX_IDS", "mesh_topology", "cell_vertex_IDs", received_file, approved_file );
            result &= FieldValuesAreEquivalentElseLog( "CELL_VERTEX_COORDINATES", "mesh_topology", "cell_vertex_coordinates", received_file, approved_file );
         }

         return result;
      }
   };
   // 'Linking' our new comparator into the library to be used on all HDF5 files.
   auto const disposer = ApprovalTests::FileApprover::registerComparatorForExtension( ".h5", std::make_shared<H5FileComparator>() );
}// namespace ApprovalBackend

namespace CaseSelection {
   enum class Case { SinglePhaseSingleNode,
                     ComplexTopology,
                     ComplexBoundaries,
                     ComplexInitialCondition };

   std::filesystem::path const base_folder = "./approval_output";

   /**
    * @brief Gives the parts of the directory or the HDF5-file path which change in the testcase.
    * @param test_case Indicator what test is to be run.
    * @param with_file_type_ending Indicator whether an HDF5-file or a directory is to be returned.
    * @param with_restart_ending Indicator whether the restart file or the output file is to be returned
    */
   std::filesystem::path NamedPartsOfCasePath( Case const test_case, bool const with_file_type_ending = false, bool const with_restart_ending = false ) {
      std::string const file_appendix = std::string( with_restart_ending ? "_restart" : "" ) + std::string( with_file_type_ending ? ".h5" : "" );
      switch( test_case ) {
         case Case::SinglePhaseSingleNode:
            return "single_phase_single_node" + file_appendix;
         case Case::ComplexInitialCondition:
            return "complex_initial_condition" + file_appendix;
         case Case::ComplexTopology:
            return "complex_topology" + file_appendix;
         default:// Case::ComplexBoundaries
            return "complex_boundaries" + file_appendix;
      }
   }

   /**
    * @brief Gives the path to directory of the respective testcase.
    * @param test_case Indicator what test is to be run.
    */
   std::filesystem::path CaseFolderPath( Case const test_case ) {
      return base_folder / NamedPartsOfCasePath( test_case );
   }

   /**
    * @brief Gives the path to approved output file for the respective testcase.
    * @param test_case Indicator what test is to be run.
    * @param restart Indicator whether the restart file or the output file is to be returned
    */
   std::filesystem::path ApprovedFile( Case const test_case, bool const restart = false ) {
      return "../test/initialization/" / NamedPartsOfCasePath( test_case, true, restart );
   }

   /**
    * @brief Gives the path to the generated HDF5-output file according to the tested case.
    * @param test_case Indicator what test is to be run.
    */
   std::filesystem::path OutputFile( Case const test_case ) {
      return CaseFolderPath( test_case ) / "domain/data_0.000000.h5";
   }

   /**
    * @brief Gives the path to the generated HDF5-restart file according to the tested case.
    * @param test_case Indicator what test is to be run.
    */
   std::filesystem::path RestartFile( Case const test_case ) {
      return CaseFolderPath( test_case ) / "restart/restart_0.000000.h5";
   }
}// namespace CaseSelection

namespace Mocking {

   namespace VaryingAttributes {
      /**
       * @brief Gives the inital condition string for the density of material one, based on the case under test.
       * @param scenario Indicator what test is to be run.
       */
      std::string DensityStringMaterialOne( CaseSelection::Case const scenario ) {
         switch( scenario ) {
            case CaseSelection::Case::ComplexTopology:
               return "density := z < 0.5 ? 1.0 : 2.0;\n";
            case CaseSelection::Case::ComplexInitialCondition:
               return "density := exp(-x);\n";
            default:
               return "density := 1.0;\n";
         }
      }

      /**
       * @brief Gives the inital x-velocity string for the density of material one, based on the case under test.
       * @param scenario Indicator what test is to be run.
       */
      std::string VelocityXStringMaterialOne( CaseSelection::Case const scenario ) {
         switch( scenario ) {
            case CaseSelection::Case::SinglePhaseSingleNode:
               return "velocityX := 2.0;\n";
            case CaseSelection::Case::ComplexInitialCondition:
               return "velocityX := tanh(x);\n";
            default:
               return "velocityX := 0.0;\n";
         }
      }

      /**
       * @brief Gives the inital y-velocity string for the density of material one, based on the case under test.
       * @param scenario Indicator what test is to be run.
       */
      std::string VelocityYStringMaterialOne( CaseSelection::Case const scenario ) {
         return scenario == CaseSelection::Case::SinglePhaseSingleNode ? "velocityY := 3.0;\n" : "velocityY := 0.0;\n";
      }

      /**
       * @brief Gives the inital condition string for the z-velocity of material one, based on the case under test.
       * @param scenario Indicator what test is to be run.
       */
      std::string VelocityZStringMaterialOne( CaseSelection::Case const scenario ) {
         switch( scenario ) {
            case CaseSelection::Case::SinglePhaseSingleNode:
               return "velocityZ := 4.0;\n";
            case CaseSelection::Case::ComplexInitialCondition:
               return "velocityZ := abs( y*y );\n";
            default:
               return "velocityZ := 0.0;\n";
         }
      }

      /**
       * @brief Gives the inital condition string for the pressure of material one, based on the case under test.
       * @param scenario Indicator what test is to be run.
       */
      std::string PressureStringMaterialOne( CaseSelection::Case const scenario ) {
         switch( scenario ) {
            case CaseSelection::Case::SinglePhaseSingleNode:
               return "pressure := 5.0;\n";
            case CaseSelection::Case::ComplexInitialCondition:
               return "pressure := 2.0 + sin( y * pi );\n";
            case CaseSelection::Case::ComplexTopology:
               return "pressure := z < 0.5 ? 1.0 : 2.0;\n";
            default:
               return "pressure := 1.0\n";
         }
      }

      /**
       * @brief Gives the inital condition string for the density of material two, based on the case under test.
       * @param scenario Indicator what test is to be run.
       */
      std::string DensityStringMaterialTwo( CaseSelection::Case const scenario ) {
         return scenario == CaseSelection::Case::ComplexInitialCondition ? "density :=  ( 1 + x ) + 0.25 * sin( y / z );" : "density := 0.125;";
      }

      /**
       * @brief Gives the inital condition string for the x-velocity of material two, based on the case under test.
       * @param scenario Indicator what test is to be run.
       */
      std::string VelocityXStringMaterialTwo( CaseSelection::Case const scenario ) {
         return scenario == CaseSelection::Case::ComplexInitialCondition ? "velocityX := y;" : "velocityX := 0.0;";
      }

      /**
       * @brief Gives the inital condition string for the y-velocity of material two, based on the case under test.
       * @param scenario Indicator what test is to be run.
       */
      std::string VelocityYStringMaterialTwo( CaseSelection::Case const scenario ) {
         return scenario == CaseSelection::Case::ComplexInitialCondition ? "velocityY := z;" : "velocityY := 0.0;";
      }

      /**
       * @brief Gives the inital condition string for the z-velocity of material two, based on the case under test.
       * @param scenario Indicator what test is to be run.
       */
      std::string VelocityZStringMaterialTwo( CaseSelection::Case const scenario ) {
         return scenario == CaseSelection::Case::ComplexInitialCondition ? "velocityZ := x;" : "velocityZ := 0.0;";
      }

      /**
       * @brief Gives the inital condition string for the pressure of material two, based on the case under test.
       * @param scenario Indicator what test is to be run.
       */
      std::string PressureStringMaterialTwo( CaseSelection::Case const scenario ) {
         return scenario == CaseSelection::Case::ComplexInitialCondition ? "pressure := ( 1.0 / x )*( y / z );" : "pressure := 1.0;";
      }

      /**
       * @brief Gives the inital condition string for the levelset, based on the case under test.
       * @param scenario Indicator what test is to be run.
       */
      std::string LevelsetString( CaseSelection::Case const scenario ) {
         switch( scenario ) {
            case CaseSelection::Case::ComplexInitialCondition:
               return "phi := 0.4 - sqrt( pow( x - 1, 2 ) + pow( y - 1, 2) );";
            case CaseSelection::Case::ComplexTopology:
               return "phi := 3.0 - x;";
            case CaseSelection::Case::ComplexBoundaries:
               return "phi := 0.5 - y;";
            default://SinglePhaseSingleNode
               return "phi := 1.0;";
         }
      }

      /**
       * @brief Gives the name the input file which matches the mocked setup, based on the case under test.
       * @param scenario Indicator what test is to be run.
       */
      std::string InputfileString( CaseSelection::Case const scenario ) {
         switch( scenario ) {
            case CaseSelection::Case::SinglePhaseSingleNode:
               return "single_phase_single_node.xml";
            case CaseSelection::Case::ComplexInitialCondition:
               return "complex_initial_conditions.xml";
            case CaseSelection::Case::ComplexTopology:
               return "complex_topology.xml";
            default://CaseSelection::Case::ComplexBoundaries:
               return "complex_boundaries.xml";
         }
      }
   }// namespace VaryingAttributes

   using namespace fakeit;

   /**
    * @brief Gives a material reader that acts according to the given test scenario.
    * @param scenario Indicator what test is to be run.
    */
   Mock<MaterialReader> InstantiateMaterialReader( CaseSelection::Case const scenario ) {
      std::unordered_map<std::string, double> const eos_gas_data   = { { "gamma", 1.4 }, { "backgroundPressure", 0.0 } };
      std::unordered_map<std::string, double> const eos_water_data = { { "gamma", 4.4 }, { "A", 1.0 }, { "B", 1.0 }, { "rho0", 1.0 } };
      EquationOfStateName const eos_gas_name( EquationOfStateName::StiffenedGas );
      EquationOfStateName const eos_water_name( EquationOfStateName::WaterlikeFluid );

      Mock<MaterialReader> material_reader;
      switch( scenario ) {
         case CaseSelection::Case::SinglePhaseSingleNode:
         case CaseSelection::Case::ComplexTopology:
            When( Method( material_reader, ReadNumberOfMaterials ) ).AlwaysReturn( 1 );
            break;
         default://ComplexInitialConditions, ComplexBoundaries
            When( Method( material_reader, ReadNumberOfMaterials ) ).AlwaysReturn( 2 );
            When( Method( material_reader, ReadEquationOfStateName ).Using( 2 ) ).AlwaysReturn( eos_water_name );
            When( Method( material_reader, ReadEquationOfStateData ).Using( 2 ) ).AlwaysReturn( eos_water_data );
            When( Method( material_reader, ReadFixedValue ).Using( { 2 }, MaterialProperty::BulkViscosity ) ).AlwaysReturn( 0.0 );
            When( Method( material_reader, ReadFixedValue ).Using( { 2 }, MaterialProperty::ShearViscosity ) ).AlwaysReturn( 0.0 );
            When( Method( material_reader, ReadFixedValue ).Using( { 1, 2 }, MaterialProperty::SurfaceTensionCoefficient ) ).AlwaysReturn( 0.0 );
            break;
      }

      When( Method( material_reader, ReadFixedValue ).Using( { 1 }, MaterialProperty::ShearViscosity ) ).AlwaysReturn( 0.0 );
      When( Method( material_reader, ReadEquationOfStateName ).Using( 1 ) ).AlwaysReturn( eos_gas_name );
      When( Method( material_reader, ReadEquationOfStateData ).Using( 1 ) ).AlwaysReturn( eos_gas_data );
      When( Method( material_reader, ReadFixedValue ).Using( { 1 }, MaterialProperty::BulkViscosity ) ).AlwaysReturn( 0.0 );
      When( Method( material_reader, ReadFixedValue ) ).AlwaysReturn( 0.0 );
      When( Method( material_reader, ReadMaterialType ) ).AlwaysReturn( MaterialType::Fluid );

      return material_reader;
   }

   /**
    * @brief Gives a multiresolution reader that acts according to the given test scenario.
    * @param scenario Indicator what test is to be run.
    */
   Mock<MultiResolutionReader> InstantiateMultiResolutionReader( CaseSelection::Case const scenario ) {
      Mock<MultiResolutionReader> multiresolution_reader;
      switch( scenario ) {
         case CaseSelection::Case::ComplexTopology:
            When( Method( multiresolution_reader, ReadNumberOfNodes ).Using( Direction::X ) ).AlwaysReturn( 1 );
            When( Method( multiresolution_reader, ReadNumberOfNodes ).Using( Direction::Y ) ).AlwaysReturn( 1 );
            When( Method( multiresolution_reader, ReadNumberOfNodes ).Using( Direction::Z ) ).AlwaysReturn( 2 );
            When( Method( multiresolution_reader, ReadMaximumLevel ) ).AlwaysReturn( 2 );
            break;
         default:
            When( Method( multiresolution_reader, ReadNumberOfNodes ).Using( Direction::X ) ).AlwaysReturn( 1 );
            When( Method( multiresolution_reader, ReadNumberOfNodes ).Using( Direction::Y ) ).AlwaysReturn( 1 );
            When( Method( multiresolution_reader, ReadNumberOfNodes ).Using( Direction::Z ) ).AlwaysReturn( 1 );
            When( Method( multiresolution_reader, ReadMaximumLevel ) ).AlwaysReturn( 0 );
            break;
      }
      When( Method( multiresolution_reader, ReadNodeSizeOnLevelZero ) ).AlwaysReturn( 1 );
      When( Method( multiresolution_reader, ReadEpsilonLevelReference ) ).AlwaysReturn( 1 );
      When( Method( multiresolution_reader, ReadEpsilonReference ) ).AlwaysReturn( 0.01 );

      return multiresolution_reader;
   }

   /**
    * @brief Gives a boundary condition reader that acts according to the given test scenario.
    * @param scenario Indicator what test is to be run.
    */
   Mock<BoundaryConditionReader> InstantiateBoundaryConditionReader( CaseSelection::Case const scenario ) {
      Mock<BoundaryConditionReader> boundary_condition_reader;
      switch( scenario ) {
         case CaseSelection::Case::ComplexBoundaries:
            When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::West ) ).AlwaysReturn( MaterialBoundaryType::ZeroGradient );
            When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::East ) ).AlwaysReturn( MaterialBoundaryType::Wall );
            When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::South ) ).AlwaysReturn( MaterialBoundaryType::Periodic );
            When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::North ) ).AlwaysReturn( MaterialBoundaryType::Periodic );
            When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::Top ) ).AlwaysReturn( MaterialBoundaryType::Symmetry );
            When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::Bottom ) ).AlwaysReturn( MaterialBoundaryType::FixedValue );
            When( Method( boundary_condition_reader, ReadMaterialFixedValueBoundaryConditions ).Using( BoundaryLocation::Bottom ) ).AlwaysReturn( { 1, 1, 0, 0, 0 } );

            When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::East ) ).AlwaysReturn( LevelSetBoundaryType::Symmetry );
            When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::West ) ).AlwaysReturn( LevelSetBoundaryType::Symmetry );
            When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::South ) ).AlwaysReturn( LevelSetBoundaryType::Periodic );
            When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::North ) ).AlwaysReturn( LevelSetBoundaryType::Periodic );
            When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::Top ) ).AlwaysReturn( LevelSetBoundaryType::ZeroGradient );
            When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::Bottom ) ).AlwaysReturn( LevelSetBoundaryType::ZeroGradient );
            break;
         default:
            When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::East ) ).AlwaysReturn( MaterialBoundaryType::Symmetry );
            When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::West ) ).AlwaysReturn( MaterialBoundaryType::Symmetry );
            When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::South ) ).AlwaysReturn( MaterialBoundaryType::Symmetry );
            When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::North ) ).AlwaysReturn( MaterialBoundaryType::Symmetry );
            When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::Top ) ).AlwaysReturn( MaterialBoundaryType::Symmetry );
            When( Method( boundary_condition_reader, ReadMaterialBoundaryType ).Using( BoundaryLocation::Bottom ) ).AlwaysReturn( MaterialBoundaryType::Symmetry );

            When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::East ) ).AlwaysReturn( LevelSetBoundaryType::Symmetry );
            When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::West ) ).AlwaysReturn( LevelSetBoundaryType::Symmetry );
            When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::South ) ).AlwaysReturn( LevelSetBoundaryType::Symmetry );
            When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::North ) ).AlwaysReturn( LevelSetBoundaryType::Symmetry );
            When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::Top ) ).AlwaysReturn( LevelSetBoundaryType::Symmetry );
            When( Method( boundary_condition_reader, ReadLevelsetBoundaryType ).Using( BoundaryLocation::Bottom ) ).AlwaysReturn( LevelSetBoundaryType::Symmetry );
            break;
      }

      return boundary_condition_reader;
   }

   /**
    * @brief Gives a source term reader that acts according to the given test scenario.
    * @param scenario Indicator what test is to be run.
    */
   Mock<SourceTermReader> InstantiateSourceTermReader() {
      Mock<SourceTermReader> source_term_reader;
      When( Method( source_term_reader, ReadGravity ).Using( Direction::X ) ).Return( 0 );
      When( Method( source_term_reader, ReadGravity ).Using( Direction::Y ) ).Return( 0 );
      When( Method( source_term_reader, ReadGravity ).Using( Direction::Z ) ).Return( 0 );
      return source_term_reader;
   }

   /**
    * @brief Gives a time control reader that acts according to the given test scenario.
    * @param scenario Indicator what test is to be run.
    */
   Mock<TimeControlReader> InstantiateTimeControlReader() {
      Mock<TimeControlReader> time_control_reader;
      When( Method( time_control_reader, ReadStartTime ) ).AlwaysReturn( 0.0 );
      When( Method( time_control_reader, ReadEndTime ) ).AlwaysReturn( 0.0 );
      When( Method( time_control_reader, ReadCFLNumber ) ).AlwaysReturn( 0.6 );
      return time_control_reader;
   }

   /**
    * @brief Gives a resart reader that acts according to the given test scenario.
    * @param scenario Indicator what test is to be run.
    */
   Mock<RestartReader> InstantiateRestartReader() {
      Mock<RestartReader> restart_reader;
      When( Method( restart_reader, ReadRestoreMode ) ).Return( RestoreMode::Off );
      When( Method( restart_reader, ReadSnapshotTimesType ) ).AlwaysReturn( SnapshotTimesType::Stamps );
      When( Method( restart_reader, ReadSnapshotTimeStamps ) ).AlwaysReturn( { 0.0 } );
      return restart_reader;
   }

   /**
    * @brief Gives an output reader that acts according to the given test scenario.
    * @param scenario Indicator what test is to be run.
    */
   Mock<OutputReader> InstantiateOutputReader() {
      Mock<OutputReader> output_reader;
      When( Method( output_reader, ReadOutputTimesType ).Using( OutputType::Standard ) ).Return( OutputTimesType::Interval );
      When( Method( output_reader, ReadOutputTimesType ).Using( OutputType::Interface ) ).Return( OutputTimesType::Off );
      When( Method( output_reader, ReadOutputInterval ).Using( OutputType::Standard ) ).Return( 0.000001 );
      When( Method( output_reader, ReadTimeNamingFactor ) ).AlwaysReturn( 1.e0 );
      return output_reader;
   }

   /**
    * @brief Gives an inital condition reader that acts according to the given test scenario.
    * @param scenario Indicator what test is to be run.
    */
   Mock<InitialConditionReader> InstantiateInitialConditionReader( CaseSelection::Case const scenario ) {
      using namespace VaryingAttributes;
      Mock<InitialConditionReader> initial_condition_reader;
      std::string material_one = DensityStringMaterialOne( scenario ) + VelocityXStringMaterialOne( scenario ) + VelocityYStringMaterialOne( scenario ) + VelocityZStringMaterialOne( scenario ) + PressureStringMaterialOne( scenario );
      std::string material_two = DensityStringMaterialTwo( scenario ) + VelocityXStringMaterialTwo( scenario ) + VelocityYStringMaterialTwo( scenario ) + VelocityZStringMaterialTwo( scenario ) + PressureStringMaterialTwo( scenario );
      std::string levelset     = LevelsetString( scenario );
      When( Method( initial_condition_reader, ReadMaterialInitialConditions ).Using( 1 ) ).AlwaysReturn( material_one );
      When( Method( initial_condition_reader, ReadMaterialInitialConditions ).Using( 2 ) ).AlwaysReturn( material_two );
      When( Method( initial_condition_reader, ReadLevelsetInitializerInput ).Using( 1 ) ).Return( levelset );
      When( Method( initial_condition_reader, ReadLevelsetInitializerType ) ).AlwaysReturn( LevelsetInitializerType::Functional );
      When( Method( initial_condition_reader, ReadLevelsetInitializerBoundingBoxes ) ).AlwaysReturn( std::vector<std::array<double, 6>>() );
      return initial_condition_reader;
   }

   /**
    * @brief Gives an input reader reader that acts according to the given test scenario.
    * @param scenario Indicator what test is to be run.
    */
   Mock<InputReader> InstantiateInputReader( Mock<MaterialReader>& material_reader, Mock<MultiResolutionReader>& multiresolution_reader,
                                             Mock<BoundaryConditionReader>& boundary_condition_reader, Mock<TimeControlReader>& time_control_reader,
                                             Mock<RestartReader>& restart_reader, Mock<OutputReader>& output_reader,
                                             Mock<InitialConditionReader>& initial_condition_reader, Mock<SourceTermReader>& source_term_reader, CaseSelection::Case const scenario ) {
      Mock<InputReader> input_reader;
      When( Method( input_reader, GetMaterialReader ) ).AlwaysReturn( material_reader.get() );
      When( Method( input_reader, GetBoundaryConditionReader ) ).AlwaysReturn( boundary_condition_reader.get() );
      When( Method( input_reader, GetMultiResolutionReader ) ).AlwaysReturn( multiresolution_reader.get() );
      When( Method( input_reader, GetTimeControlReader ) ).Return( time_control_reader.get() );
      When( Method( input_reader, GetRestartReader ) ).Return( restart_reader.get() );
      When( Method( input_reader, GetOutputReader ) ).Return( output_reader.get() );
      When( Method( input_reader, GetInputType ) ).Return( InputType::Xml );
      When( Method( input_reader, GetInputFile ) ).Return( VaryingAttributes::InputfileString( scenario ) );
      When( Method( input_reader, GetInitialConditionReader ) ).AlwaysReturn( initial_condition_reader.get() );
      When( Method( input_reader, GetTimeControlReader ) ).AlwaysReturn( time_control_reader.get() );
      When( Method( input_reader, GetSourceTermReader ) ).AlwaysReturn( source_term_reader.get() );
      return input_reader;
   }
}// namespace Mocking

namespace ScenarioTest {

   /**
    * @brief Executes a full unit test (GIVEN-WHEN-THEN) to test the initialization using HDF5-file approval.
    * @param test_case Indicator what test is to be run.
    */
   void RunTestsForCase( CaseSelection::Case const test_case ) {
      auto const case_base_folder = CaseFolderPath( test_case );
      {// Separate scope and long info message needed as approvals + catch library combine their message passing in strange ways ...
         INFO( "If you see only a single failure message the approval-output directory could not be created. Maybe a folder from a previous runs still exists?. Otherwise refer the the other failure messages first." );
         REQUIRE( std::filesystem::create_directories( case_base_folder ) );
      }
      GIVEN( "Input reader mocks providing initial conditions for this case" ) {
         fakeit::Mock<MaterialReader> material_reader( Mocking::InstantiateMaterialReader( test_case ) );
         fakeit::Mock<MultiResolutionReader> multiresolution_reader( Mocking::InstantiateMultiResolutionReader( test_case ) );
         fakeit::Mock<BoundaryConditionReader> boundary_condition_reader( Mocking::InstantiateBoundaryConditionReader( test_case ) );
         fakeit::Mock<TimeControlReader> time_control_reader( Mocking::InstantiateTimeControlReader() );
         fakeit::Mock<RestartReader> restart_reader( Mocking::InstantiateRestartReader() );
         fakeit::Mock<OutputReader> output_reader( Mocking::InstantiateOutputReader() );
         fakeit::Mock<InitialConditionReader> initial_condition_reader( Mocking::InstantiateInitialConditionReader( test_case ) );
         fakeit::Mock<SourceTermReader> source_term_reader( Mocking::InstantiateSourceTermReader() );

         fakeit::Mock<InputReader> input_reader( Mocking::InstantiateInputReader( material_reader, multiresolution_reader, boundary_condition_reader,
                                                                                  time_control_reader, restart_reader, output_reader, initial_condition_reader,
                                                                                  source_term_reader, test_case ) );

         WHEN( "The simulation is initialized" ) {
            UnitHandler const unit_handler( 1.0, 1.0, 1.0, 1.0 );
            //UnitHandler const unit_handler( 0.75 , 0.75, 0.75, 0.75 );
            MaterialManager const material_manager( Instantiation::InstantiateMaterialManager( input_reader.get(), unit_handler ) );
            TopologyManager topology_manager( Instantiation::InstantiateTopologyManager( input_reader.get(), material_manager ) );
            Multiresolution const multiresolution( Instantiation::InstantiateMultiresolution( input_reader.get(), topology_manager ) );
            Tree tree( Instantiation::InstantiateTree( input_reader.get(), topology_manager, unit_handler ) );
            CommunicationManager communication_manager( Instantiation::InstantiateCommunicationManager( topology_manager ) );
            ExternalHaloManager const external_halo_manager( Instantiation::InstantiateExternalHaloManager( input_reader.get(), unit_handler, material_manager ) );
            InternalHaloManager internal_halo_manager( Instantiation::InstantiateInternalHaloManager( topology_manager, tree, communication_manager, material_manager ) );
            HaloManager halo_manager( Instantiation::InstantiateHaloManager( topology_manager, tree, external_halo_manager, internal_halo_manager, communication_manager ) );
            OutputWriter const output_writer( Instantiation::InstantiateOutputWriter( topology_manager, tree, material_manager, unit_handler ) );
            RestartManager const restart_manager( Instantiation::InstantiateRestartManager( topology_manager, tree, unit_handler ) );
            InputOutputManager input_output_manager( Instantiation::InstantiateInputOutputManager( input_reader.get(), output_writer, restart_manager, unit_handler, case_base_folder ) );
            ModularAlgorithmAssembler modular_assembler( Instantiation::InstantiateModularAlgorithmAssembler( input_reader.get(), topology_manager, tree, communication_manager, halo_manager, multiresolution,
                                                                                                              material_manager, input_output_manager, unit_handler ) );

            auto const initial_condition( Instantiation::InstantiateInitialCondition( input_reader.get(), topology_manager, tree, material_manager, unit_handler ) );
            modular_assembler.Initialization( *initial_condition );

            THEN( "The inital output and restart matches the expected one" ) {
               {// separate scope for logging purposes
                  ApprovalBackend::approving_restart = false;
                  auto const output_file             = CaseSelection::OutputFile( test_case );
                  auto const approved_output_file    = CaseSelection::ApprovedFile( test_case );
                  INFO( "Approving output failed. Output directory is kept for manual inspection" );
                  ApprovalTests::FileApprover::verify( output_file, approved_output_file );
               }
               {// separate scope for logging purposes
                  ApprovalBackend::approving_restart = true;
                  auto const restart_file            = CaseSelection::RestartFile( test_case );
                  auto const approved_restart_file   = CaseSelection::ApprovedFile( test_case, true );
                  INFO( "Approving restart failed. Output directory is kept for manual inspection" );
                  ApprovalTests::FileApprover::verify( restart_file, approved_restart_file );
               }
               constexpr std::size_t expected_number_of_files_to_clean = 10;
               REQUIRE( std::filesystem::remove_all( case_base_folder ) == expected_number_of_files_to_clean );
            }
         }
      }
   }
}// namespace ScenarioTest

SCENARIO( "Initialization of a single node with one material", "[.slow1rank]" ) {
   constexpr CaseSelection::Case test_case = CaseSelection::Case::SinglePhaseSingleNode;
   ScenarioTest::RunTestsForCase( test_case );
}

SCENARIO( "Initialization of a complex topology", "[.slow1rank]" ) {
   constexpr CaseSelection::Case test_case = CaseSelection::Case::ComplexTopology;
   ScenarioTest::RunTestsForCase( test_case );
}

SCENARIO( "Initalization of complex initial condtions", "[.slow1rank]" ) {
   constexpr CaseSelection::Case test_case = CaseSelection::Case::ComplexInitialCondition;
   ScenarioTest::RunTestsForCase( test_case );
}

SCENARIO( "Initalization of complex boundary conditions", "[.slow1rank]" ) {
   constexpr CaseSelection::Case test_case = CaseSelection::Case::ComplexBoundaries;
   ScenarioTest::RunTestsForCase( test_case );
}
