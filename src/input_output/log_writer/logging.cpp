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

#include "input_output/log_writer/log_writer.h"

#include "integrator/time_integrator_setup.h"
#include "levelset/levelset_advector/levelset_advector_setup.h"
#include "levelset/multi_phase_manager/multi_phase_manager_setup.h"
#include "stencils/spatial_derivative_stencils/derivative_stencil_setup.h"
#include "stencils/spatial_reconstruction_stencils/reconstruction_stencil_setup.h"
#include "user_specifications/compile_time_constants.h"
#include "user_specifications/riemann_solver_settings.h"
#include "user_specifications/equation_settings.h"
#include "user_specifications/state_reconstruction_settings.h"
#include "solvers/convective_term_contributions/riemann_solvers/riemann_solver_setup.h"

#define TOSTRING1( str ) #str
#define TOSTRING( str ) TOSTRING1( str )

namespace {
   /**
    * @brief Gives a pretty-print string telling whether or not a parameter is active.
    * @brief active Indicates wether or not the variable is active.
    * @brief parameter_message The message to preceed the (de-)active indicator string.
    * @note Once compiler support for constexpr strings is a thing, this function can be constexpr.
    */
   std::string ActiveStatusString( bool const active, std::string const& parameter_message ) {
      return parameter_message + ( active ? "ACTIVE" : "NOT active" );
   }
}// namespace

namespace Logging {
   /**
    * @brief Logs the compile-time decisions used in this compilation of the program.
    */
   void LogCompiledSettings() {
      LogWriter& logger = LogWriter::Instance();

      logger.LogMessage( std::string( "Based on git commit: " ) + TOSTRING( GITHASH ) );

      char const* val = getenv( "HOSTNAME" );
      if( val != NULL ) {
         logger.LogMessage( "Running on              : " + std::string( val ) );
      } else {
         logger.LogMessage( "Running on              : Hostname not known" );
      }
      logger.LogMessage( "Dimensions              : " + std::to_string( static_cast<int>( CC::DIM() ) ) );
      if constexpr( CC::DIM() == Dimension::Two ) {
         std::string message = "Axis symmetry";
         if constexpr( CC::Axisymmetric() ) {
            message.append( "           : ON" );
         } else {
            message.append( "           : OFF" );
         }
         logger.LogMessage( message );
      }
      logger.LogBreakLine();

      logger.LogMessage( ActiveStatusString( CC::InviscidExchangeActive(), "Inviscid exchange       : " ) );
      logger.LogMessage( ActiveStatusString( CC::GravityIsActive(), "Gravity                 : " ) );
      logger.LogMessage( ActiveStatusString( CC::ViscosityIsActive(), "Viscosity               : " ) );
      logger.LogMessage( ActiveStatusString( CC::CapillaryForcesActive(), "Surface tension         : " ) );
      logger.LogMessage( ActiveStatusString( CC::HeatConductionActive(), "Heat conduction         : " ) );
      logger.LogMessage( ActiveStatusString( CC::ScaleSeparationActive(), "Scale seperation        : " ) );
      logger.LogMessage( ActiveStatusString( CC::FUSY(), "FP control for symmetry : " ) );
      logger.LogMessage( ActiveStatusString( CC::SolidBoundaryActive(), "Solid boundaries        : " ) );
      if constexpr( active_equations == EquationSet::GammaModel ) {
         logger.LogMessage( ActiveStatusString( GammaModelSettings::EosSafeGuarding, "Safe equation of state  : " ) );
      }
      logger.LogBreakLine();

      logger.LogMessage( "Equation set                                  : " + SetToString( active_equations ) );
      logger.LogMessage( "Interface set                                 : " + SetToString( active_interface_quantities ) );
      logger.LogMessage( "Parameter set                                 : " + SetToString( active_parameters ) );
      logger.LogMessage( "Time integrator                               : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( TimeIntegratorSetup::Concretize<time_integrator>::type ).name() ) ) );

      if constexpr( convective_term_solver == ConvectiveTermSolvers::FluxSplitting ) {
         logger.LogMessage( "Flux Splitting Scheme                         : " + FluxSplittingToString( FluxSplittingSettings::flux_splitting_scheme ) );
         if constexpr( FluxSplittingSettings::flux_splitting_scheme == FluxSplitting::Roe_M || FluxSplittingSettings::flux_splitting_scheme == FluxSplitting::LocalLaxFriedrichs_M ) {
            logger.LogMessage( "Low-Mach-number limit factor                  : " + std::to_string( FluxSplittingSettings::low_mach_number_limit_factor ) );
         }
      }
      if constexpr( convective_term_solver == ConvectiveTermSolvers::FiniteVolume ) {
         logger.LogMessage( "Riemann solver                                : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( RiemannSolverSetup::Concretize<FiniteVolumeSettings::riemann_solver>::type ).name() ) ) );
         logger.LogMessage( "Signal Speed Selection                        : " + SignalSpeedToString( FiniteVolumeSettings::signal_speed_selection ) );
      }

      logger.LogMessage( "Reconstruction stencil                        : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( ReconstructionStencilSetup::Concretize<reconstruction_stencil>::type ).name() ) ) );
      logger.LogMessage( "Reconstruction Type                           : " + SetToString( state_reconstruction_type ) );
      logger.LogMessage( "Viscous fluxes reconstruction stencil         : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( ReconstructionStencilSetup::Concretize<viscous_fluxes_reconstruction_stencil>::type ).name() ) ) );
      logger.LogMessage( "Heat fluxes reconstruction stencil            : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( ReconstructionStencilSetup::Concretize<heat_fluxes_reconstruction_stencil>::type ).name() ) ) );
      logger.LogMessage( "Derivative stencil                            : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( DerivativeStencilSetup::Concretize<derivative_stencil>::type ).name() ) ) );
      logger.LogMessage( "Viscous fluxes derivative stencil cell center : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( DerivativeStencilSetup::Concretize<viscous_fluxes_derivative_stencil_cell_center>::type ).name() ) ) );
      logger.LogMessage( "Viscous fluxes derivative stencil cell face   : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( DerivativeStencilSetup::Concretize<viscous_fluxes_derivative_stencil_cell_face>::type ).name() ) ) );
      logger.LogMessage( "Heat fluxes derivative stencil cell center    : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( DerivativeStencilSetup::Concretize<heat_fluxes_derivative_stencil_cell_center>::type ).name() ) ) );
      logger.LogMessage( "Heat fluxes derivative stencil cell face      : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( DerivativeStencilSetup::Concretize<heat_fluxes_derivative_stencil_cell_face>::type ).name() ) ) );
      logger.LogMessage( "Curvature calculation derivative stencil      : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( DerivativeStencilSetup::Concretize<curvature_calculation_derivative_stencil>::type ).name() ) ) );
      logger.LogMessage( "Normal calculation derivative stencil         : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( DerivativeStencilSetup::Concretize<normal_calculation_derivative_stencil>::type ).name() ) ) );
      logger.LogBreakLine();
      logger.LogMessage( "Multi-phase manager        : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( MultiPhaseManagerSetup::Concretize<phase_manager>::type ).name() ) ) );
      logger.LogMessage( "LS Advection               : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( LevelsetAdvectorSetup::Concretize<levelset_advector>::type ).name() ) ) );
      logger.LogMessage( "LS Reinitialization        : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( LevelsetReinitializerSetup::Concretize<levelset_reinitializer>::type ).name() ) ) );
      logger.LogMessage( "Geometry calculator        : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( GeometryCalculatorSetup::Concretize<geometry_calculator>::type ).name() ) ) );
      logger.LogMessage( "Mixing method              : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( CutCellMixerSetup::Concretize<cut_cell_mixer>::type ).name() ) ) );
      logger.LogMessage( "Extension method           : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( GhostFluidExtenderSetup::Concretize<extender>::type_primestates ).name() ) ) );
      logger.LogMessage( "Interface Extension method : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( InterfaceExtenderSetup::Concretize<interface_extender>::type_states ).name() ) ) );
      logger.LogMessage( "Interface Riemann solver   : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( InterfaceRiemannSolverSetup::Concretize<interface_riemann_solver>::type ).name() ) ) );
      logger.LogMessage( "Scale Separation method    : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( ScaleSeparatorSetup::Concretize<scale_separator>::type ).name() ) ) );
      logger.LogMessage( "Buffer handler             : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( BufferHandlerSetup::Concretize<buffer_handler>::type ).name() ) ) );
      logger.LogBreakLine();
      logger.LogMessage( "Output Vertex Filer: " + VertexFilterTypeToString( CC::VERTEX_FILTER() ) );
      logger.LogBreakLine();
      logger.Flush();
   }
}// namespace Logging
