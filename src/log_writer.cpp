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
#include "log_writer.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <stdlib.h> /* getenv */
#include <typeinfo>
#include <functional>

#include "user_specifications/equation_settings.h"
#include "user_specifications/compile_time_constants.h"
#include "user_specifications/riemann_solver_settings.h"
#include "utilities/helper_functions.h"
#include "utilities/string_operations.h"
#include "enums/flux_splitting.h"

#include "solvers/riemann_solvers/riemann_solver_setup.h"
#include "prime_states/prime_state_handler.h"
#include "integrator/time_integrator_setup.h"
#include "levelset/multi_phase_manager/multi_phase_manager_setup.h"
#include "levelset/levelset_advector/levelset_advector_setup.h"
#include "levelset/multi_phase_manager/levelset_reinitializer/levelset_reinitializer_setup.h"
#include "levelset/multi_phase_manager/cut_cell_mixer/cut_cell_mixer_setup.h"
#include "interface_interaction/interface_riemann_solver/interface_riemann_solver_setup.h"

#include "user_specifications/stencil_setup.h"
#include "stencils/spatial_reconstruction_stencils/reconstruction_stencil_setup.h"
#include "stencils/spatial_derivative_stencils/derivative_stencil_setup.h"
#include "levelset/multi_phase_manager/ghost_fluid_extender/ghost_fluid_extender_setup.h"
#include "levelset/geometry/geometry_calculator_setup.h"
#include "levelset/multi_phase_manager/scale_separator/scale_separator_setup.h"
#include "levelset/multi_phase_manager/interface_extender/interface_extender_setup.h"
#include "levelset/multi_phase_manager/buffer_handler_setup.h"

// double expansion to get string from compiler variable:
#define TOSTRING1( str ) #str
#define TOSTRING( str ) TOSTRING1( str )

/**
 * @brief Default constructor, allows to set if all or just the root rank are to be logged.
 * @param save_all_ranks Decider if all ranks are to be logged. If false only root is logged.
 */
LogWriter::LogWriter( bool const save_all_ranks ) : logfile_name_( "Unnamed_Simulation.log" ),
                                                    rank_( ForwardRankId() ),
                                                    save_all_ranks_( save_all_ranks ) {

   if( rank_ == 0 ) {
      std::cout << "|******************************************************************************|\n"
                << "|*  \\\\                                                                        *|\n"
                << "|*  l '>                                                                      *|\n"
                << "|*  | |                                                                       *|\n"
                << "|*  | |                                                                       *|\n"
                << "|*  | alpaca~                                                                 *|\n"
                << "|*  ||    ||                                                                  *|\n"
                << "|*  ''    ''                                                                  *|\n"
                << "|*                                                                            *|\n"
                << "|*  THE AGE OF ALPACA HAS COME                                                *|\n"
                << "|*                                                                            *|\n"
                << "|******************************************************************************|\n"
                << "|* based on git commit:                                                       *|\n"
                << "|* " << std::setw( 74 ) << std::left << TOSTRING( GITHASH ) << " *|\n"
                << "|******************************************************************************|" << std::endl;
   }

   char const* val = getenv( "HOSTNAME" );
   if( val != NULL ) {
      LogMessage( "Running on              : " + std::string( val ), true, true );
   } else {
      LogMessage( "Running on              : Hostname not known", true, true );
   }
   LogMessage( "Dimensions              : " + std::to_string( static_cast<int>( CC::DIM() ) ) );
   if constexpr( CC::DIM() == Dimension::Two ) {
      std::string message = "Axis symmetry";
      if constexpr( CC::Axisymmetric() ) {
         message.append( "           : ON" );
      } else {
         message.append( "           : OFF" );
      }
      LogMessage( message, true, true );
   }
   AddBreakLine( true );

   // Logging the input parameters - Does not claim completeness.
   if constexpr( CC::InviscidExchangeActive() ) {
      LogMessage( "Inviscid exchange       : ACTIVE" );
   } else {
      LogMessage( "Inviscid exchange       : NOT active" );
   }
   if constexpr( CC::GravityIsActive() ) {
      LogMessage( "Gravity                 : ACTIVE" );
   } else {
      LogMessage( "Gravity                 : NOT active" );
   }
   if constexpr( CC::ViscosityIsActive() ) {
      LogMessage( "Viscosity               : ACTIVE" );
   } else {
      LogMessage( "Viscosity               : NOT active" );
   }
   if constexpr( CC::CapillaryForcesActive() ) {
      LogMessage( "Surface tension         : ACTIVE" );
   } else {
      LogMessage( "Surface tension         : NOT active" );
   }
   if constexpr( CC::HeatConductionActive() ) {
      LogMessage( "Heat conduction         : ACTIVE" );
   } else {
      LogMessage( "Heat conduction         : NOT active" );
   }
   if constexpr( CC::ScaleSeparationActive() ) {
      LogMessage( "Scale seperation        : ACTIVE" );
   } else {
      LogMessage( "Scale seperation        : NOT active" );
   }
   if constexpr( CC::FUSY() ) {
      LogMessage( "FP control for symmetry : ACTIVE" );
   } else {
      LogMessage( "FP control for symmetry : NOT active" );
   }
   if constexpr( CC::SolidBoundaryActive() ) {
      LogMessage( "Solid boundaries        : ACTIVE" );
   } else {
      LogMessage( "Solid boundaries        : NOT active" );
   }
   AddBreakLine( true );

   LogMessage( "Equation set                                  : " + SetToString( active_equations ) );
   LogMessage( "Interface set                                 : " + SetToString( active_interface_quantities ) );
   LogMessage( "Parameter set                                 : " + SetToString( active_parameters ) );
   LogMessage( "Time integrator                               : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( TimeIntegratorSetup::Concretize<time_integrator>::type ).name() ) ) );
   LogMessage( "Riemann solver                                : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( RiemannSolverSetup::Concretize<riemann_solver>::type ).name() ) ) );
   if constexpr( riemann_solver == RiemannSolvers::Roe ) {
      LogMessage( "Flux Splitting Scheme                         : " + FluxSplittingToString( RoeSolverSettings::flux_splitting_scheme ) );
   }
   if constexpr( riemann_solver == RiemannSolvers::Hllc || riemann_solver == RiemannSolvers::Hll ) {
      LogMessage( "Signal Speed Selection                        : " + SignalSpeedToString( HllSolverSettings::signal_speed_selection ) );
   }
   if constexpr( RoeSolverSettings::flux_splitting_scheme == FluxSplitting::Roe_M || RoeSolverSettings::flux_splitting_scheme == FluxSplitting::LocalLaxFriedrichs_M ) {
      LogMessage( "Low-Mach-number limit factor                  : " + std::to_string( RoeSolverSettings::low_mach_number_limit_factor ) );
   }
   LogMessage( "Reconstruction stencil                        : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( ReconstructionStencilSetup::Concretize<reconstruction_stencil>::type ).name() ) ) );
   LogMessage( "Viscous fluxes reconstruction stencil         : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( ReconstructionStencilSetup::Concretize<viscous_fluxes_reconstruction_stencil>::type ).name() ) ) );
   LogMessage( "Heat fluxes reconstruction stencil            : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( ReconstructionStencilSetup::Concretize<heat_fluxes_reconstruction_stencil>::type ).name() ) ) );
   LogMessage( "Derivative stencil                            : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( DerivativeStencilSetup::Concretize<derivative_stencil>::type ).name() ) ) );
   LogMessage( "Viscous fluxes derivative stencil cell center : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( DerivativeStencilSetup::Concretize<viscous_fluxes_derivative_stencil_cell_center>::type ).name() ) ) );
   LogMessage( "Viscous fluxes derivative stencil cell face   : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( DerivativeStencilSetup::Concretize<viscous_fluxes_derivative_stencil_cell_face>::type ).name() ) ) );
   LogMessage( "Heat fluxes derivative stencil cell center    : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( DerivativeStencilSetup::Concretize<heat_fluxes_derivative_stencil_cell_center>::type ).name() ) ) );
   LogMessage( "Heat fluxes derivative stencil cell face      : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( DerivativeStencilSetup::Concretize<heat_fluxes_derivative_stencil_cell_face>::type ).name() ) ) );
   LogMessage( "Curvature calculation derivative stencil      : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( DerivativeStencilSetup::Concretize<curvature_calculation_derivative_stencil>::type ).name() ) ) );
   LogMessage( "Normal calculation derivative stencil         : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( DerivativeStencilSetup::Concretize<normal_calculation_derivative_stencil>::type ).name() ) ) );
   AddBreakLine( true );
   LogMessage( "Multi-phase manager        : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( MultiPhaseManagerSetup::Concretize<phase_manager>::type ).name() ) ) );
   LogMessage( "LS Advection               : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( LevelsetAdvectorSetup::Concretize<levelset_advector>::type ).name() ) ) );
   LogMessage( "LS Reinitialization        : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( LevelsetReinitializerSetup::Concretize<levelset_reinitializer>::type ).name() ) ) );
   LogMessage( "Geometry calculator        : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( GeometryCalculatorSetup::Concretize<geometry_calculator>::type ).name() ) ) );
   LogMessage( "Mixing method              : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( CutCellMixerSetup::Concretize<cut_cell_mixer>::type ).name() ) ) );
   LogMessage( "Extension method           : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( GhostFluidExtenderSetup::Concretize<extender>::type_primestates ).name() ) ) );
   LogMessage( "Interface Extension method : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( InterfaceExtenderSetup::Concretize<interface_extender>::type_states ).name() ) ) );
   LogMessage( "Interface Riemann solver   : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( InterfaceRiemannSolverSetup::Concretize<interface_riemann_solver>::type ).name() ) ) );
   LogMessage( "Scale Separation method    : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( ScaleSeparatorSetup::Concretize<scale_separator>::type ).name() ) ) );
   LogMessage( "Buffer handler             : " + StringOperations::RemoveLeadingNumbers( std::string( typeid( BufferHandlerSetup::Concretize<buffer_handler>::type ).name() ) ) );
   AddBreakLine( true );
   LogMessage( "Output Vertex Filer: " + VertexFilterTypeToString( CC::VERTEX_FILTER() ) );
   AddBreakLine( true );
}

/**
 * @brief Writes the welcome message to the output file.
 */
void LogWriter::FlushWelcomeMessage() {

   int number_of_ranks;
   MPI_Comm_size( MPI_COMM_WORLD, &number_of_ranks );
   std::string const filename      = logfile_name_;
   std::string const rank_filename = filename + "_" + std::to_string( rank_ );
   auto rank_messages( FormatMessage( "Data from Rank " + std::to_string( rank_ ) + ":" ) );
   auto rank_message = rank_messages[0] + "\n";// Dirty Hack, we know the line is to short.

   if( rank_ == 0 ) {
      std::ofstream output_stream( filename, std::ios::app );
      output_stream << std::scientific << std::setprecision( 5 );
      output_stream << "|******************************************************************************|" << std::endl;
      output_stream << "|*  \\\\                                                                        *|" << std::endl;
      output_stream << "|*  l '>                                                                      *|" << std::endl;
      output_stream << "|*  | |                                                                       *|" << std::endl;
      output_stream << "|*  | |                                                                       *|" << std::endl;
      output_stream << "|*  | alpaca~                                                                 *|" << std::endl;
      output_stream << "|*  ||    ||                                                                  *|" << std::endl;
      output_stream << "|*  ''    ''                                                                  *|" << std::endl;
      output_stream << "|*                                                                            *|" << std::endl;
      output_stream << "|* TERMINAL OUTPUT HAS GONE - THE ALPACA LOGBOOK HAS COME                     *|" << std::endl;
      output_stream << "|*                                                                            *|" << std::endl;
      output_stream << "|******************************************************************************|" << std::endl;
      output_stream << "|* based on git commit:                                                       *|" << std::endl;
      output_stream << "|* " << std::setw( 74 ) << std::left << TOSTRING( GITHASH ) << " *|" << std::endl;
      output_stream << "|******************************************************************************|" << std::endl;
      output_stream.flush();
      output_stream.close();
   }

   if( save_all_ranks_ && number_of_ranks > 1 ) {
      std::ofstream output_stream( rank_filename, std::ios::app );
      output_stream << std::scientific << std::setprecision( 5 );
      output_stream << "|******************************************************************************|" << std::endl;
      output_stream << "|*  \\\\                                                                        *|" << std::endl;
      output_stream << "|*  l '>                                                                      *|" << std::endl;
      output_stream << "|*  | |                                                                       *|" << std::endl;
      output_stream << "|*  | |                                                                       *|" << std::endl;
      output_stream << "|*  | alpaca~                                                                 *|" << std::endl;
      output_stream << "|*  ||    ||                                                                  *|" << std::endl;
      output_stream << "|*  ''    ''                                                                  *|" << std::endl;
      output_stream << "|*                                                                            *|" << std::endl;
      output_stream << "|* TERMINAL OUTPUT HAS GONE - THE ALPACA LOGBOOK HAS COME                     *|" << std::endl;
      output_stream << "|*                                                                            *|" << std::endl;
      output_stream << "|******************************************************************************|" << std::endl;
      output_stream << "|* based on git commit:                                                       *|" << std::endl;
      output_stream << "|* " << std::setw( 74 ) << std::left << TOSTRING( GITHASH ) << " *|" << std::endl;
      output_stream << "|******************************************************************************|" << std::endl;
      output_stream << rank_message;
      output_stream.flush();
      output_stream.close();
   }
}

/**
 * @brief Writes the output file, and collects output from all ranks if necessary
 */
void LogWriter::Flush() {

   int number_of_ranks;
   MPI_Comm_size( MPI_COMM_WORLD, &number_of_ranks );
   std::string const filename      = logfile_name_;
   std::string const rank_filename = filename + "_" + std::to_string( rank_ );

   if( rank_ == 0 ) {
      std::cout << "|******************************************************************************|" << std::endl;
      std::ofstream output_stream( filename, std::ios::app );
      AddBreakLine();
      output_stream << log_;
      output_stream.flush();
      output_stream.close();
   }

   if( save_all_ranks_ && number_of_ranks > 1 ) {
      std::ofstream output_stream( rank_filename, std::ios::app );
      AddBreakLine();
      output_stream << log_;
      output_stream.flush();
      output_stream.close();
   }

   log_.clear();
}

/**
 * @brief      Plot an Alpaca crossing the terminal.
 *
 * @param[in]  percentage  Indicates which percentage of the final simulation time is already finished. Based on that the position of the
 *                         Alpaca is deternined.
 * @param[in]  fast_forward  Whether the Alpaca should rush to its current position to indicate that the simulation was restarted.
 */
void LogWriter::FlushAlpaca( double const percentage, bool const fast_forward ) {

   double cut_percentage              = std::max( 0.0, percentage );
   cut_percentage                     = std::min( cut_percentage, 1.0 );
   unsigned int const plot_percentage = (unsigned int)std::floor( cut_percentage * 64 + 9 );

   std::stringstream line_stream;
   AddBreakLine( true );
   line_stream << "  ";
   LogMessage( line_stream.str() );
   line_stream.str( "" );
   line_stream << std::setw( plot_percentage - 1 ) << " \\\\";
   LogMessage( line_stream.str() );
   line_stream.str( "" );
   if( fast_forward ) line_stream << std::setfill( '>' );
   line_stream << std::setw( plot_percentage + 1 ) << " l '>";
   if( fast_forward ) line_stream << std::setfill( ' ' );
   LogMessage( line_stream.str() );
   line_stream.str( "" );
   line_stream << std::setw( plot_percentage ) << " | |";
   LogMessage( line_stream.str() );
   line_stream.str( "" );
   if( fast_forward ) line_stream << std::setfill( '>' );
   line_stream << std::setw( plot_percentage ) << " | |";
   if( fast_forward ) line_stream << std::setfill( ' ' );
   LogMessage( line_stream.str() );
   line_stream.str( "" );
   line_stream << std::setw( plot_percentage ) << "~alpaca |";
   LogMessage( line_stream.str() );
   line_stream.str( "" );
   if( fast_forward ) line_stream << std::setfill( '>' );
   line_stream << std::setw( plot_percentage ) << " ||    ||";
   if( fast_forward ) line_stream << std::setfill( ' ' );
   LogMessage( line_stream.str() );
   line_stream.str( "" );
   line_stream << std::setw( plot_percentage ) << " ''    ''";
   LogMessage( line_stream.str() );
   line_stream.str( "" );
   line_stream << "  ";
   LogMessage( line_stream.str() );
   line_stream.str( "" );
   AddBreakLine( true );
}

/**
 * @brief      Appends a string to a delayed log message.
 *
 * @param[in]  delayed_log  The string that should be appended.
 */
void LogWriter::AppendDelayedLog( std::string const delayed_log ) {
   delayed_log_ += delayed_log;
}

/**
 * @brief Gives the rank id from MPI directly as int. Avoids handle creation, e.g. for const members in initializer list.
 * @return Rank id
 */
int LogWriter::ForwardRankId() const {
   int rank_id = -1;
   MPI_Comm_rank( MPI_COMM_WORLD, &rank_id );
   return rank_id;
}

/**
 * @brief Brings a message into the same ASCII-Layout.
 * @param message An arbitrary string message $MUST NOT CONTAIN TABS '\t' or NEWLINES '\n'$
 * @return Gives a string in homogeneous ASCII-Layout.
 */
std::vector<std::string> LogWriter::FormatMessage( std::string const& message ) const {
   unsigned int message_size    = message.size();
   unsigned int number_of_lines = std::ceil( double( message_size ) / double( 76.0 ) );

   std::vector<std::string> final_lines;

   for( unsigned int i = 0; i < number_of_lines; ++i ) {
      if( message.size() <= i * 75 + 75 ) {
         final_lines.emplace_back( message.begin() + ( i * 75 ), message.end() );
      } else {
         final_lines.emplace_back( message.begin() + ( i * 75 ), message.begin() + ( i * 75 ) + 75 );
      }
   }

   for( auto& line : final_lines ) {
      if( line.size() < 75 ) {
         int difference = 75 - line.size();
         line.insert( line.end(), difference, ' ' );
      }
      line.insert( line.end(), '*' );
      line.insert( line.end(), '|' );
      line.insert( 0, "|* " );
   }

   return final_lines;
}

/**
 * @brief Writes a message to the terminal (std::cout) and/or save the message in order to include it in the log file.
 * @param message String containing the message to be printed/logged.
 * @param print_to_terminal Decider if message is to be printed to std::cout.
 * @param save_in_logfile Decider if message is to be saved in the log file.
 */
void LogWriter::LogMessage( std::string const& message, bool const print_to_terminal, bool const save_in_logfile ) {

   std::vector<std::string> formatted( FormatMessage( message ) );

   if( print_to_terminal ) {
      if( rank_ == 0 ) {
         for( auto const& line : formatted ) {
            std::cout << std::scientific << std::setprecision( 5 ) << line << std::endl;
         }
      }
   }
   if( save_in_logfile ) {
      if( save_all_ranks_ ) {
         for( auto const& line : formatted ) {
            log_.append( line );
            log_.append( "\n" );
         }
      } else {
         if( rank_ == 0 ) {
            for( auto const& line : formatted ) {
               log_.append( line );
               log_.append( "\n" );
            }
         }
      }
   }
}

/**
 * @brief Writes a message to the terminal (std::cout) and/or save the message in order to include it in the log file.
 * @param message String containing the message to be printed/logged.
 * @param print_to_terminal Decider if message is to be printed to std::cout.
 * @param save_in_logfile Decider if message is to be saved in the log file.
 */
void LogWriter::LogLinebreakMessage( std::string const& message, bool const print_to_terminal, bool const save_in_logfile ) {
   // Instantiate vector containing split strings if needed
   std::vector<std::string> all_lines;
   std::string line;
   std::stringstream stream( message );
   while( std::getline( stream, line, '\n' ) ) {
      all_lines.push_back( line );
   }

   // Loop through all single_messages and call single send message function
   for( auto const& single_message : all_lines ) {
      LogMessage( single_message, print_to_terminal, save_in_logfile );
   }
}

/**
 * @brief Writes a delayed message to the terminal (std::cout) and/or save the message in order to include it in the log file.
 * @param print_to_terminal Decider if message is to be printed to std::cout.
 * @param save_in_logfile Decider if message is to be saved in the log file.
 */
void LogWriter::DelayedLogMessage( bool const print_to_terminal, bool const save_in_logfile ) {
   LogMessage( delayed_log_, print_to_terminal, save_in_logfile );
   delayed_log_.clear();
}

/**
 * @brief Introduces a breaking line to the log message to enhance clearness of the logging output.
 * @param print_to_terminal Decision wether break line is also written to the terminal cout. Default: false
 */
void LogWriter::AddBreakLine( bool const print_to_terminal ) {
   LogMessage( "**************************************************************************", print_to_terminal );
}

/**
 * @brief Sets the name of the log file.
 * @param name Log file name.
 */
void LogWriter::SetLogfileName( std::string const name ) {
   logfile_name_ = name;
}

/**
 * @brief This function is used to get the Logger. If no Logger exists yet it is created, otherwise the existing logger is passed back. "Singleton Constructor"
 * @param save_all_ranks Decider whether or not all ranks are to be logged. (Only relevant at first call!)
 * @return The logger instance.
 */
LogWriter& LogWriter::Instance( bool save_all_ranks ) {
   static LogWriter instance_( save_all_ranks );
   return instance_;
}
