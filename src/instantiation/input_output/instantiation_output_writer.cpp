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
#include "instantiation/input_output/instantiation_output_writer.h"

#include "user_specifications/compile_time_constants.h"
#include "user_specifications/output_constants.h"
// mesh generators for the topology generations
#include "input_output/output_writer/mesh_generator/standard_mpi_mesh_generator.h"
#include "input_output/output_writer/mesh_generator/standard_finest_level_mesh_generator.h"
#include "input_output/output_writer/mesh_generator/interface_mesh_generator.h"
#include "input_output/output_writer/mesh_generator/debug_mesh_generator.h"
// output quantities for the cell data generation
#include "input_output/output_writer/output_quantities/material_field_quantities/material_field_quantity.h"
#include "input_output/output_writer/output_quantities/interface_field_quantities/interface_field_quantity.h"
#include "input_output/output_writer/output_quantities/custom_material_quantities/partition_output.h"
#include "input_output/output_writer/output_quantities/custom_material_quantities/interface_tags_output.h"
#include "input_output/output_writer/output_quantities/custom_material_quantities/mach_number_output.h"
#include "input_output/output_writer/output_quantities/custom_material_quantities/numerical_schlieren_output.h"
#include "input_output/output_writer/output_quantities/custom_material_quantities/vorticity_absolute_output.h"
#include "input_output/output_writer/output_quantities/custom_material_quantities/helicity_output.h"
#include "input_output/output_writer/output_quantities/custom_material_quantities/baroclinicity_output.h"
#include "input_output/output_writer/output_quantities/custom_material_quantities/vortex_dilatation_output.h"
#include "input_output/output_writer/output_quantities/custom_material_quantities/vortex_stretching_output.h"

namespace {
   /**
    * @brief Checks whether any of the given array flags is active.
    * @param setting Array holding all setting flags.
    * @return True if any is active, False otherwise.
    */
   inline bool IsAnyActive( std::array<bool, 3> const& setting ) {
      return std::any_of( setting.begin(), setting.end(), []( bool const output ) { return output; } );
   }
}// namespace

namespace Instantiation {

   /**
   * @brief Add conservative and primitive state quantities to material output vector.
   * @param quantity the given MaterialFieldQuantity.
   * @param quantity_flags quantity output decision flags.
   * @param quantity_name quantity name.
   * @param unit_handler Instance for dimensionalization of quantities.
   * @param material_manager Instance for handling material properties and pairing properties.
   * @param debug_flag Enable special debug treatment.
   */
   void AddStateQuantity( MaterialFieldQuantityName const quantity,
                          std::array<bool, 3> const quantity_flags,
                          std::string const quantity_name,
                          std::vector<std::unique_ptr<OutputQuantity const>>& quantity_vector,
                          UnitHandler const& unit_handler,
                          MaterialManager const& material_manager,
                          bool const debug_flag ) {

      // Obtain data and verify it
      MaterialFieldQuantityData quantity_data( DataOfMaterialFieldQuantity( quantity ) );
      if( quantity_data.VerifyQuantity() ) {
         // Add quantities to the vector
         quantity_vector.emplace_back( std::make_unique<MaterialFieldQuantity const>( unit_handler, material_manager, quantity_name, quantity_flags, quantity_data ) );
         // Check if the debug output is activated and add the right-hand side buffer
         if( debug_flag ) {
            if( quantity_flags.back() && MaterialFieldOutputSettings::use_all_buffers_in_debug ) {
               quantity_vector.emplace_back(
                     std::make_unique<MaterialFieldQuantity const>( unit_handler, material_manager, quantity_name, std::array<bool, 3>( { false, false, true } ),
                                                                    quantity_data, ConservativeBufferType::RightHandSide ) );
            }
         }
      }
   }

   /**
   * @brief Add interface quantities to material output vector.
   * @param quantity the given InterfaceQuantity.
   * @param quantity_flags quantity output decision flags.
   * @param quantity_name quantity name.
   * @param unit_handler Instance for dimensionalization of quantities.
   * @param material_manager Instance for handling material properties and pairing properties.
   * @param debug_flag Enable special debug treatment.
   */
   void AddInterfaceQuantity( InterfaceFieldQuantityName const quantity,
                              std::array<bool, 3> const quantity_flags,
                              std::string const quantity_name,
                              std::vector<std::unique_ptr<OutputQuantity const>>& quantity_vector,
                              UnitHandler const& unit_handler,
                              MaterialManager const& material_manager,
                              bool const debug_flag ) {

      // Obtain data and verify it
      InterfaceFieldQuantityData quantity_data( DataOfInterfaceFieldQuantity( quantity ) );
      if( quantity_data.VerifyQuantity() ) {
         // Add quantity to vector
         quantity_vector.emplace_back(
               std::make_unique<InterfaceFieldQuantity const>( unit_handler, material_manager, quantity_name, quantity_flags, quantity_data ) );
         // Check if the debug output is activated and add the right-hand side and reinitialized buffer
         if( debug_flag ) {
            if( quantity_flags.back() && InterfaceFieldOutputSettings::use_all_buffers_in_debug ) {
               quantity_vector.emplace_back(
                     std::make_unique<InterfaceFieldQuantity const>( unit_handler, material_manager, quantity_name, std::array<bool, 3>( { false, false, true } ),
                                                                     quantity_data, InterfaceDescriptionBufferType::RightHandSide ) );
               quantity_vector.emplace_back(
                     std::make_unique<InterfaceFieldQuantity const>( unit_handler, material_manager, quantity_name, std::array<bool, 3>( { false, false, true } ),
                                                                     quantity_data, InterfaceDescriptionBufferType::Reinitialized ) );
               quantity_vector.emplace_back(
                     std::make_unique<InterfaceFieldQuantity const>( unit_handler, material_manager, quantity_name, std::array<bool, 3>( { false, false, true } ),
                                                                     quantity_data, InterfaceDescriptionBufferType::Integrated ) );
            }
         }
      }
   }

   /**
    * @brief Gives pointers of all material output quantities.
    * @param unit_handler Instance for dimensionalization of quantities.
    * @param material_manager Instance for handling material properties and pairing properties.
    * @return Vector of pointer to all material output quantities.
    */
   std::vector<std::unique_ptr<OutputQuantity const>> GetMaterialOutputQuantities( UnitHandler const& unit_handler, MaterialManager const& material_manager ) {
      // short name using
      namespace MFOS = MaterialFieldOutputSettings;

      // Declare vector to be returned
      std::vector<std::unique_ptr<OutputQuantity const>> output_quantities;

      /**************************************************************************************************/
      /*                                 MATERIAL FIELD QUANTITIES                                      */
      /**************************************************************************************************/
      // conservatives
      if( IsAnyActive( MFOS::Mass ) ) {
         AddStateQuantity( MaterialFieldQuantityName::Mass, MFOS::Mass, MFOS::MassName, output_quantities, unit_handler, material_manager, true );
      }
      if( IsAnyActive( MFOS::Momentum ) ) {
         AddStateQuantity( MaterialFieldQuantityName::Momentum, MFOS::Momentum, MFOS::MomentumName, output_quantities, unit_handler, material_manager, true );
      }
      if( IsAnyActive( MFOS::Energy ) ) {
         AddStateQuantity( MaterialFieldQuantityName::Energy, MFOS::Energy, MFOS::EnergyName, output_quantities, unit_handler, material_manager, true );
      }
      if( IsAnyActive( MFOS::GammaConservative ) ) {
         AddStateQuantity( MaterialFieldQuantityName::GammaConservative, MFOS::GammaConservative, MFOS::GammaConservativeName, output_quantities, unit_handler, material_manager, true );
      }
      if( IsAnyActive( MFOS::PiConservative ) ) {
         AddStateQuantity( MaterialFieldQuantityName::PiConservative, MFOS::PiConservative, MFOS::PiConservativeName, output_quantities, unit_handler, material_manager, true );
      }

      // prime states
      if( IsAnyActive( MFOS::Velocity ) ) {
         AddStateQuantity( MaterialFieldQuantityName::Velocity, MFOS::Velocity, MFOS::VelocityName, output_quantities, unit_handler, material_manager, false );
      }
      if( IsAnyActive( MFOS::Pressure ) ) {
         AddStateQuantity( MaterialFieldQuantityName::Pressure, MFOS::Velocity, MFOS::PressureName, output_quantities, unit_handler, material_manager, false );
      }
      if( IsAnyActive( MFOS::Temperature ) ) {
         AddStateQuantity( MaterialFieldQuantityName::Temperature, MFOS::Temperature, MFOS::TemperatureName, output_quantities, unit_handler, material_manager, false );
      }
      if( IsAnyActive( MFOS::Density ) ) {
         AddStateQuantity( MaterialFieldQuantityName::Density, MFOS::Density, MFOS::DensityName, output_quantities, unit_handler, material_manager, false );
      }
      if( IsAnyActive( MFOS::GammaPrimitive ) ) {
         AddStateQuantity( MaterialFieldQuantityName::GammaPrimitive, MFOS::GammaPrimitive, MFOS::GammaPrimitiveName, output_quantities, unit_handler, material_manager, false );
      }
      if( IsAnyActive( MFOS::PiPrimitive ) ) {
         AddStateQuantity( MaterialFieldQuantityName::PiPrimitive, MFOS::PiPrimitive, MFOS::PiPrimitiveName, output_quantities, unit_handler, material_manager, false );
      }

      // parameter output ( can only be used if models are used, otherwise the required buffers do not exist)
      if constexpr( CC::ParameterModelActive() ) {
         if( IsAnyActive( MFOS::ShearViscosity ) ) {
            AddStateQuantity( MaterialFieldQuantityName::ShearViscosity, MFOS::ShearViscosity, MFOS::ShearViscosityName, output_quantities, unit_handler, material_manager, false );
         }
         if( IsAnyActive( MFOS::ThermalConductivity ) ) {
            AddStateQuantity( MaterialFieldQuantityName::ThermalConductivity, MFOS::ThermalConductivity, MFOS::ThermalConductivityName, output_quantities, unit_handler, material_manager, false );
         }
      }
      /**************************************************************************************************/
      /*                                 CUSTOM QUANTITIES                                              */
      /**************************************************************************************************/
      // short name using
      namespace COS = CustomOutputSettings;

      if( IsAnyActive( COS::Partition ) ) {
         output_quantities.push_back( std::make_unique<PartitionOutput const>( unit_handler, material_manager, COS::PartitionName, COS::Partition ) );
      }
      if( IsAnyActive( COS::InterfaceTags ) ) {
         output_quantities.push_back( std::make_unique<InterfaceTagsOutput const>( unit_handler, material_manager, COS::InterfaceTagsName, COS::InterfaceTags ) );
      }
      if( IsAnyActive( COS::MachNumber ) ) {
         output_quantities.push_back( std::make_unique<MachNumberOutput const>( unit_handler, material_manager, COS::MachNumberName, COS::MachNumber ) );
      }
      if( IsAnyActive( COS::NumericalSchlieren ) ) {
         output_quantities.push_back( std::make_unique<NumericalSchlierenOutput const>( unit_handler, material_manager, COS::NumericalSchlierenName, COS::NumericalSchlieren ) );
      }
      if( IsAnyActive( COS::VorticityAbsolute ) ) {
         output_quantities.push_back( std::make_unique<VorticityAbsoluteOutput const>( unit_handler, material_manager, COS::VorticityAbsoluteName, COS::VorticityAbsolute ) );
      }
      if( IsAnyActive( COS::Helicity ) ) {
         output_quantities.push_back( std::make_unique<HelicityOutput const>( unit_handler, material_manager, COS::HelicityName, COS::Helicity ) );
      }
      if( IsAnyActive( COS::Baroclinicity ) ) {
         output_quantities.push_back( std::make_unique<BaroclinicityOutput const>( unit_handler, material_manager, COS::BaroclinicityName, COS::Baroclinicity ) );
      }
      if( IsAnyActive( COS::VortexDilatation ) ) {
         output_quantities.push_back( std::make_unique<VortexDilatationOutput const>( unit_handler, material_manager, COS::VortexDilatationName, COS::VortexDilatation ) );
      }
      if( IsAnyActive( COS::VortexStretching ) ) {
         output_quantities.push_back( std::make_unique<VortexStretchingOutput const>( unit_handler, material_manager, COS::VortexStretchingName, COS::VortexStretching ) );
      }

      return output_quantities;
   }

   /**
    * @brief Gives pointers of all interface output quantities.
    * @param unit_handler Instance for dimensionalization of quantities.
    * @param material_manager Instance for handling material properties and pairing properties.
    * @return Vector of pointer to all interface output quantities.
    */
   std::vector<std::unique_ptr<OutputQuantity const>> GetInterfaceOutputQuantities( UnitHandler const& unit_handler, MaterialManager const& material_manager ) {
      // short name using
      namespace IFOS = InterfaceFieldOutputSettings;

      // Declare vector to be returned
      std::vector<std::unique_ptr<OutputQuantity const>> output_quantities;
      /**************************************************************************************************/
      /*                                 INTERFACE FIELD QUANTITIES                                     */
      /**************************************************************************************************/
      // interface descriptions
      if( IsAnyActive( IFOS::Levelset ) ) {
         AddInterfaceQuantity( InterfaceFieldQuantityName::Levelset, IFOS::Levelset, IFOS::LevelsetName, output_quantities, unit_handler, material_manager, true );
      }
      if( IsAnyActive( IFOS::VolumeFraction ) ) {
         AddInterfaceQuantity( InterfaceFieldQuantityName::VolumeFraction, IFOS::VolumeFraction, IFOS::VolumeFractionName, output_quantities, unit_handler, material_manager, true );
      }

      // interface state output
      if( IsAnyActive( IFOS::InterfaceVelocity ) ) {
         AddInterfaceQuantity( InterfaceFieldQuantityName::InterfaceVelocity, IFOS::InterfaceVelocity, IFOS::InterfaceVelocityName, output_quantities, unit_handler, material_manager, false );
      }
      if( IsAnyActive( IFOS::PressurePositive ) ) {
         AddInterfaceQuantity( InterfaceFieldQuantityName::PressurePositive, IFOS::PressurePositive, IFOS::PressurePositiveName, output_quantities, unit_handler, material_manager, false );
      }
      if( IsAnyActive( IFOS::PressureNegative ) ) {
         AddInterfaceQuantity( InterfaceFieldQuantityName::PressureNegative, IFOS::PressureNegative, IFOS::PressureNegativeName, output_quantities, unit_handler, material_manager, false );
      }

      // interface parameter output ( can only be used if models are used, otherwise the required buffers do not exist)
      if constexpr( CC::InterfaceParameterModelActive() ) {
         if( IsAnyActive( IFOS::SurfaceTensionCoefficient ) ) {
            AddInterfaceQuantity( InterfaceFieldQuantityName::SurfaceTensionCoefficient, IFOS::SurfaceTensionCoefficient, IFOS::SurfaceTensionCoefficientName, output_quantities, unit_handler, material_manager, false );
         }
      }
      /**************************************************************************************************/
      /*                                 CUSTOM QUANTITIES                                              */
      /**************************************************************************************************/
      // No custom interface quantities are implemented yet

      return output_quantities;
   }

   /**
    * @brief Gives the standard mesh generator depending on the Vertex filter.
    * @param topology_manager Class providing global (on all ranks) node information.
    * @param tree Tree class providing local (on current rank) node information.
    * @param node_size_on_level_zero Size of one block on level zero.
    * @return Pointer to the base class of all mesh generators.
    */
   std::unique_ptr<MeshGenerator const> GetStandardMeshGenerator( TopologyManager const& topology_manager,
                                                                  Tree const& tree,
                                                                  double const node_size_on_level_zero,
                                                                  std::array<unsigned int, 3> const number_of_nodes_on_level_zero ) {

      if constexpr( CC::VERTEX_FILTER() == VertexFilterType::FinestLevel ) {
         return std::make_unique<StandardFinestLevelMeshGenerator const>( topology_manager, tree, node_size_on_level_zero, number_of_nodes_on_level_zero );
      } else {
         return std::make_unique<StandardMpiMeshGenerator const>( topology_manager, tree, node_size_on_level_zero, CC::VERTEX_FILTER() == VertexFilterType::Mpi );
      }
   }

   /**
    * @brief Instantiates the complete output writer class with the given input classes.
    * @param topology_manager Class providing global (on all ranks) node information.
    * @param tree Tree class providing local (on current rank) node information.
    * @param material_manager Instance providing initialized material data.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return The fully instantiated OutputWriter class as pointer (allows movements of it).
    */
   OutputWriter InstantiateOutputWriter( TopologyManager& topology_manager,
                                         Tree& tree,
                                         MaterialManager const& material_manager,
                                         UnitHandler const& unit_handler ) {

      // get the block size on level zero
      double const node_size_on_level_zero( unit_handler.DimensionalizeValue( tree.GetNodeSizeOnLevelZero(), UnitType::Length ) );
      std::array<unsigned int, 3> const number_of_nodes_on_level_zero( topology_manager.GetNumberOfNodesOnLevelZero() );

      return OutputWriter( GetStandardMeshGenerator( topology_manager, tree, node_size_on_level_zero, number_of_nodes_on_level_zero ),
                           std::make_unique<DebugMeshGenerator const>( topology_manager, tree, node_size_on_level_zero, number_of_nodes_on_level_zero[2] ),
                           std::make_unique<InterfaceMeshGenerator const>( topology_manager, tree, node_size_on_level_zero ),
                           GetMaterialOutputQuantities( unit_handler, material_manager ),
                           GetInterfaceOutputQuantities( unit_handler, material_manager ),
                           material_manager.GetNumberOfMaterials() );
   }
}// namespace Instantiation
