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
#ifndef COMPILE_TIME_CONSTANTS_H
#define COMPILE_TIME_CONSTANTS_H

#include "boundary_condition/boundary_specifications.h"
#include <array>
#include <limits>
#include "enums/dimension_definition.h"
#include "enums/norms.h"
#include "enums/vertex_filter_type.h"

/**
 * @brief The CompileTimeConstants class holds all constants which are determined at compile time. Also the central place to change settings.
 */
class CompileTimeConstants {

   /*** TO BE SET BY EXPERIENCED USERS ***/
   static constexpr unsigned int internal_cells_per_block_and_dimension_ = 16; //Referred to as "IC"
   static constexpr unsigned int halo_width_ = 4; // Referred to as "HS"
   static constexpr unsigned int cells_per_dimension_with_halo_ = 2 * halo_width_ + internal_cells_per_block_and_dimension_; // Referred to as "TC"

   /*
    * The number of neighbors considered (in each direction) during the prediction. E.g. prediction_stencil = 2 means 5 cells  per dimension are considered:
    * The base the two left and two right neighbors
    */
   static constexpr unsigned int prediction_stencil_size_ = 2;

   /*
    * Alpha in \cite Roussel2013.
    * Needs to be adapted if time and/or space discretization schemes are changed.
    */
   static constexpr unsigned int space_time_discretization_order_ = 3;

   //Flags to enable/disable contributions for the right-hand side/advection terms of the Navier-Stokes equations
   static constexpr bool inviscid_exchange_active_ = true; // Enables advection term
   static constexpr bool gravitation_active_       = true; // Enables gravitational term
   static constexpr bool viscosity_active_         = true; // Enables viscosity term
   static constexpr bool heat_conduction_active_   = false; // Enables heat conduction term
   static constexpr bool capillary_forces_active_  = false;  // Enables surface tension term

   // Flags to activate material property models
   static constexpr bool viscosity_model_active_                     = false;
   static constexpr bool thermal_conductivity_model_active_          = false;
   static constexpr bool surface_tension_coefficient_model_active_   = false;

   // Flag to enforce symmetry for sums of more than two values
   static constexpr bool full_symmetry_active_     = true;
   // Flag to enable scale separation during levelset reinitialization
   static constexpr bool scale_separation_active_  = false;

   //flag for running axisymmetric simulations. The axial direction is always the y-direction, the radial
   //direction the x-direction
   static constexpr bool axisymmetric_ = false;
   // Flag to activate density dependency of Gruneisen coefficient (e.g., for NobleAbel EoS)
   static constexpr bool gruneisen_density_dependent_ = false;

   //Time step control
   static constexpr bool write_timestep_list_ = false; // All micro time steps are written into a .txt file
   static constexpr bool limit_end_time_ = true; // Last macro time step is adapted to exactly match user end time
   static constexpr double minimum_time_step_size_ = std::numeric_limits<double>::epsilon(); // Minimum micro time step size

   /* This factor defines the width of the levelset narrow band in terms of cellsizes
    * Outside of this range the levelset is fixed to this distance (levelset_cutoff_factor_ * cellsize)
    */
   static constexpr double levelset_cutoff_factor_ = 8.0;

   /*
    * Defines the threshold value for the extension into cut-cells for ghost fluid extension. If the volume fraction is lower than the given threshold, the extension
    * procedure is enabled. The mixing threshold is only relevant for cut-cells. All other cells are not affected.
    */
#if DIMENSION == 1
   static constexpr double mixing_threshold_ = 0.5; // 0.5 as mixing threshold for 1D is on purpose.
#elif DIMENSION == 2
   static constexpr double mixing_threshold_ = 0.5;
#else
   static constexpr double mixing_threshold_ = 0.6;
#endif

   // Width of the narrow band around cut-cells
   static constexpr unsigned int extension_band_ = 3; // Band for the ghost-fluid extension
   static constexpr unsigned int reinitialization_band_ = 4; // Band for the reinitialization

   // Norm used for the wavelet analysis triggering the refinement/coarsening of cells
   static constexpr Norm norm_for_wavelet_analysis_ = Norm::Linfinity;

   // Specification of the vertex filter for the standard output to remove doubled placed vertices
   static constexpr VertexFilterType output_vertex_filter_ = VertexFilterType::Mpi;

   /*** DEDUCED OR FIXED VALUES - MUST NOT BE CHANGED ***/

   // Macro "PERFORMANCE" set through makefile (only).

   // Macro is set through Makefile. Determines numbers of spatial dimensions considered in the simulation.
#if DIMENSION == 1
   static constexpr Dimension dimension_of_simulation_ = Dimension::One;
#elif DIMENSION == 2
   static constexpr Dimension dimension_of_simulation_ = Dimension::Two;
#else
#define HILBERT //Hilbert Load Balancing is only meaningful in 3D Simulations
   static constexpr Dimension dimension_of_simulation_ = Dimension::Three;
#endif

   // Multiresolution specifications
   static constexpr unsigned int number_of_children_ = dimension_of_simulation_ == Dimension::One ? 2 : (dimension_of_simulation_ == Dimension::Two ? 4 : 8);
   static constexpr unsigned int algorithmic_maximum_number_of_levels_ = 14; //MUST NOT BE CHANGED!

   // Number of natural plane domain sides (1D: East-West, 2D: East-West North-South, 3D: East-West North-South Top-Bottom)
   static constexpr std::underlying_type<Dimension>::type domain_sides_ = 2 * DTI(dimension_of_simulation_);

   //MUST BE IN THE ORDER AS THEIR UNDERLYING INDEX!
   static constexpr std::array<BoundaryLocation, 6> natural_boundary_sides_ = {BoundaryLocation::East, BoundaryLocation::West, BoundaryLocation::North,
                                                                               BoundaryLocation::South, BoundaryLocation::Top, BoundaryLocation::Bottom};
#if DIMENSION == 1
   static constexpr std::array<BoundaryLocation, 2> halo_boundary_sides_= { {
                                                                             BoundaryLocation::East, BoundaryLocation::West
                                                                          } };
#elif DIMENSION == 2
   static constexpr std::array<BoundaryLocation, 8> halo_boundary_sides_= { {
                                                                             BoundaryLocation::East, BoundaryLocation::West, BoundaryLocation::North, BoundaryLocation::South,
                                                                             // Diagonals in 2 Dimensions for Sticks in 3D and cubes in 2D
                                                                             BoundaryLocation::NorthEast, BoundaryLocation::NorthWest,
                                                                             BoundaryLocation::SouthEast,  BoundaryLocation::SouthWest,
                                                                          } };
#else
   static constexpr std::array<BoundaryLocation, 26> halo_boundary_sides_ = { { // #26
                                                                               BoundaryLocation::East, BoundaryLocation::West, BoundaryLocation::North, BoundaryLocation::South, BoundaryLocation::Top, BoundaryLocation::Bottom,
                                                                               // Diagonals in 2 Dimensions for Sticks in 3D and cubes in 2D
                                                                               BoundaryLocation::BottomNorth, BoundaryLocation::BottomSouth, BoundaryLocation::TopNorth, BoundaryLocation::TopSouth, // x-Axis Sticks
                                                                               BoundaryLocation::BottomEast, BoundaryLocation::BottomWest, BoundaryLocation::TopEast, BoundaryLocation::TopWest,  // y-Axis Sticks
                                                                               BoundaryLocation::NorthEast, BoundaryLocation::NorthWest, BoundaryLocation::SouthEast, BoundaryLocation::SouthWest,// z-Axis Sticks
                                                                               // Diagonals in 3 Dimensions for cubes in 3D
                                                                               BoundaryLocation::EastNorthTop, BoundaryLocation::EastNorthBottom, BoundaryLocation::EastSouthTop, BoundaryLocation::EastSouthBottom, //East-Side
                                                                               BoundaryLocation::WestNorthTop, BoundaryLocation::WestNorthBottom, BoundaryLocation::WestSouthTop, BoundaryLocation::WestSouthBottom, //West-Side
                                                                            } };
#endif

   static constexpr std::array <BoundaryLocation, domain_sides_> active_natural_boundary_sides_ = {BoundaryLocation::East, BoundaryLocation::West
#if DIMENSION != 1
      , BoundaryLocation::North, BoundaryLocation::South
#endif
#if DIMENSION == 3
      ,BoundaryLocation::Top, BoundaryLocation::Bottom
#endif
   };

   /*
    * Values for influencing the load distribution between ranks
    * Gives the number of coarsenings or refinements that are allowed on each rank before the topology is load balanced (chosen on experience)
    */
   static constexpr unsigned int topology_changes_until_load_balancing_ = 8;

   // Assertions for consistency checks
   static_assert((internal_cells_per_block_and_dimension_ % 4) == 0, "IC must be a multiple of four, stupid!");
   static_assert(internal_cells_per_block_and_dimension_ > 0, "IC must be greater zero, stupid!");
   static_assert((internal_cells_per_block_and_dimension_ + halo_width_ + halo_width_) < 32768, "IC must be smaller than 2^15 to fit signed int");
   static_assert((halo_width_ % 2) == 0, "Halo Width should be divisable by two");
   static_assert(halo_width_ < internal_cells_per_block_and_dimension_, "IC must be larger or equal than halo width, stupid!");
   static_assert(halo_width_ >= extension_band_, "Extension width must not be larger than halo size. With this setup, extension algorithm and tagging system do not work!");
   static_assert(halo_width_ >= reinitialization_band_, "Reinitialization width must not be larger than halo size. With this setup, reinitialization algorithm and tagging system do not work!");
   static_assert(reinitialization_band_ > extension_band_, "Reinitialization band must be larger than extension band to allow correct normal computation in last cell which gets extended values!");
   static_assert(((axisymmetric_ == true && dimension_of_simulation_ == Dimension::Two) || axisymmetric_ == false), "Axisymmetric case can only be run with DIM=2");

public:
   CompileTimeConstants() = delete;
   ~CompileTimeConstants() = default;
   CompileTimeConstants( CompileTimeConstants const& ) = delete;
   CompileTimeConstants& operator=( CompileTimeConstants const& ) = delete;
   CompileTimeConstants( CompileTimeConstants&& ) = delete;
   CompileTimeConstants& operator=( CompileTimeConstants&& ) = delete;

   /**
    * @brief Gives the dimension of the Simulation, i.e. 1D, 2D or 3D
    * @return Dimension.
    */
   static constexpr Dimension DIM() { return dimension_of_simulation_; }

   /**
    * @brief Gives the number of Borders the simulation domain has.
    * @return 1D: 2, 2D: 4, 3D:6
    */
   static constexpr unsigned int SIDES() { return domain_sides_; }

   /**@{
    * @brief Gives the number of Internal Cells in a block per dimension.
    * @return Number of internal cells. 1 if dimension does not exist
    */
   static constexpr unsigned int ICX() { return internal_cells_per_block_and_dimension_; }
   static constexpr unsigned int ICY() { return DIM()!=Dimension::One   ? internal_cells_per_block_and_dimension_ : 1; }
   static constexpr unsigned int ICZ() { return DIM()==Dimension::Three ? internal_cells_per_block_and_dimension_ : 1; }
   /**@}*/

   /**@{
    * @brief Gives the number of Total Cells in a block per dimension, i.e. number of internal cells per dimension + 2 * number of halo cells per dimension.
    * @return Number of total cells. 1 if dimension does not exist
    */
   static constexpr unsigned int TCX() { return cells_per_dimension_with_halo_; }
   static constexpr unsigned int TCY() { return DIM()!=Dimension::One   ? cells_per_dimension_with_halo_ : 1; }
   static constexpr unsigned int TCZ() { return DIM()==Dimension::Three ? cells_per_dimension_with_halo_ : 1; }
   /**@}*/

   /**@{
    * @brief Gives the index of the First Internal Cell in a block per dimension.
    * @return Index of first internal cell in block. 0 if dimension does not exist
    */
   static constexpr unsigned int FICX() { return halo_width_; }
   static constexpr unsigned int FICY() { return DIM()!=Dimension::One   ? halo_width_ : 0; }
   static constexpr unsigned int FICZ() { return DIM()==Dimension::Three ? halo_width_ : 0; }
   /**@}*/

   /**@{
    * @brief Gives the index of the Last Internal Cell in a block per dimension. I.e. the returned index must be included if the internal cells are of interest.
    * @return Index of the last internal cell in a block. 0 if dimension does not exist.
    */
   static constexpr unsigned int LICX() { return halo_width_ + internal_cells_per_block_and_dimension_ - 1; }
   static constexpr unsigned int LICY() { return DIM()!=Dimension::One   ? halo_width_ + internal_cells_per_block_and_dimension_ - 1 : 0; }
   static constexpr unsigned int LICZ() { return DIM()==Dimension::Three ? halo_width_ + internal_cells_per_block_and_dimension_ - 1 : 0; }
   /**@}*/

   /**
    * @brief Gives the size of the halo "HS = Halo Size" in its shortest direction.
    * @return Number of cells in the halo in shortest halo dimension.
    */
   static constexpr unsigned int HS() { return halo_width_; }

   /**@{
    * @brief Gives the Halo-Slice Size per dimension.
    * @return Index of first internal cell in block. 1 if dimension does not exist
    */
   static constexpr unsigned int HSSX() { return halo_width_; }
   static constexpr unsigned int HSSY() { return DIM()!=Dimension::One   ? halo_width_ : 1; }
   static constexpr unsigned int HSSZ() { return DIM()==Dimension::Three ? halo_width_ : 1; }
   /**@}*/

   /**@{
    * @brief Gives the index of the first cell in the high-index halo "FHH = First High-index Halo",
    *        i.e. the eastern halo in X-Direction, northern in Y-Direction and top in Z-Direction.
    * @return Start index of the high-index halo.
    */
   static constexpr unsigned int FHHX() { return halo_width_ + internal_cells_per_block_and_dimension_; }
   static constexpr unsigned int FHHY() { return DIM()!=Dimension::One   ? halo_width_ + internal_cells_per_block_and_dimension_ : 0; }
   static constexpr unsigned int FHHZ() { return DIM()==Dimension::Three ? halo_width_ + internal_cells_per_block_and_dimension_ : 0; }
   /**@}*/

   /**
    * @brief Indicates whether inviscid exchange processes are considered.
    * @return True if Euler equations are solved. False otherwise
    */
   static constexpr inline bool InviscidExchangeActive() { return inviscid_exchange_active_; }

   /**
    * @brief Indicates if gravity should be considered as source term.
    * @return True if gravity source terms are to be calculated. False otherwise.
    */
   static constexpr bool GravityIsActive() { return gravitation_active_; }

   /**
    * @brief Indicates if viscosity should be considered as source term.
    * @return True if viscosity source terms are to be calculated. False otherwise.
    */
   static constexpr bool ViscosityIsActive() { return viscosity_active_; }

   /**
    * @brief Indicates if heat conduction should be considered as source term.
    * @return True if heat conduction is to be calculated. False otherwise.
    */
   static constexpr bool HeatConductionActive() { return heat_conduction_active_; }

   /**
    * @brief Indicates if capillary forces should be considered as interface source term.
    * @return True if capillary forces are to be calculated. False otherwise.
    */
   static constexpr bool CapillaryForcesActive() { return capillary_forces_active_; }

   /**
    * @brief Indicates whether scale separation is to be done for multi-phase simulations or not.
    * @return True if scale separation is to be done. False otherwise.
    */
   static constexpr bool ScaleSeparationActive() { return scale_separation_active_; }

   /**
    * @brief Give a bool to indicate whether or not to run an axisymmetric simulation.
    * @return True if axisymmetric simulation is performed. False otherwise.
    */
   static constexpr bool Axisymmetric() { return axisymmetric_; }

   /**
    * @brief Indicates whether the Gruneisen parameter is density dependent.
    * @return True if the Gruneisen parameter is density dependent.
    */
   static constexpr bool GruneisenDensityDependent() { return gruneisen_density_dependent_; }

   /**
    * @brief Indicates whether the viscosity is calculated by a model or not.
    * @return True if a model is used.
    */
   static constexpr bool ShearViscosityModelActive() { return viscosity_model_active_; }

   /**
    * @brief Indicates whether the thermal conductivity is calculated by a model or not.
    * @return True if a model is used.
    */
   static constexpr bool ThermalConductivityModelActive() { return thermal_conductivity_model_active_; }

   /**
    * @brief Indicates whether the thermal conductivity is calculated by a model or not.
    * @return True if a model is used.
    */
   static constexpr bool SurfaceTensionCoefficientModelActive() { return surface_tension_coefficient_model_active_; }

   /**
    * @brief Indicates whether the thermal conductivity is calculated by a model or not.
    * @return True if a model is used.
    */
   static constexpr bool ParameterModelActive() { return viscosity_model_active_ || thermal_conductivity_model_active_; }

   /**
    * @brief Indicates whether the thermal conductivity is calculated by a model or not.
    * @return True if a model is used.
    */
   static constexpr bool InterfaceParameterModelActive() { return surface_tension_coefficient_model_active_; }

   /**
    * @brief Gives the number of the cells considered in the parent for a prediction into the child's halo (long side). "PHS = Prediction Halo Size".
    * @return Number of cells for prediction per dimension.
    */
   static constexpr unsigned int PHS() { return (cells_per_dimension_with_halo_/2) + 2 * prediction_stencil_size_; }

   /**
    * @brief Gives the number of cells considered in the parent during a prediction into the childs halo,
    *        e.g. at internal jumps. "BPS = Boundary Prediction Size".
    * @return Number of cells for boundary prediction per dimension.
    */
   static constexpr unsigned int BPS() { return (halo_width_/2) + (prediction_stencil_size_*2); }

   /**
    * @brief Gives the start index in the parent used during prediction to the high-index child, i.e. the eastern child in X-Direction.
    *        "PHCS = Prediction High-index Child Start index".
    * @return Starting index in Parent for prediction to high-index child.
    */
   static constexpr unsigned int PHCS() { return ( (cells_per_dimension_with_halo_/2) - prediction_stencil_size_ - (halo_width_/2) ); }

   /**
    * @brief Gives the start index in the parent used during prediction to the high-index child's halo.
    *        "BPHCS = Boundary Prediction High-index Child Start index".
    * @return Starting index in parent for prediction to high-index child's halo.
    */
   static constexpr unsigned int BPHCS() { return halo_width_ + internal_cells_per_block_and_dimension_ - (halo_width_/2); }

   /**
    * @brief Gives the start index in the parent used during prediction to the low-index child.
    *        "PLCS = Prediction Low-index Child Start index"
    * @return Start index in the parent for prediction to low-index child.
    */
   static constexpr unsigned int PLCS() { return (halo_width_/2) - prediction_stencil_size_; }

   /**
    * @brief Gives the index in the parent that 'lies over' the first internal cell of the low-index child.
    *        "PIOLCFIC = Parent Index Overlaying Low-index Child First Internal Cell".
    * @return Index in parent of first internal cell of low-index child.
    */
   static constexpr unsigned int PIOLCFICX() { return FICX(); }
   static constexpr unsigned int PIOLCFICY() { return FICY(); }
   static constexpr unsigned int PIOLCFICZ() { return FICZ(); }

   /**
    * @brief Gives the index in the parent which 'lies over' the first internal cell of the high-index child.
    *        "PIOHCFIC = Parent Index Overlaying High-Index Child's First Internal Cell".
    * @return Index in parent of first internal cell of high-index child.
    */
   static constexpr unsigned int PIOHCFICX() { return cells_per_dimension_with_halo_/2; }
   static constexpr unsigned int PIOHCFICY() { return DIM()!=Dimension::One   ? cells_per_dimension_with_halo_/2 : 0; }
   static constexpr unsigned int PIOHCFICZ() { return DIM()==Dimension::Three ? cells_per_dimension_with_halo_/2 : 0; }

   /**
    * @brief Gives the number of parent cells which 'lie over' the internal cells of a child.
    *        "PSOCIC = Parent Size Overlaying Children's Internal Cells".
    * @return Number of cells in parent which span over all internal cells of a child.
    */
   static constexpr unsigned int PSOCICX() { return (internal_cells_per_block_and_dimension_/2); }
   static constexpr unsigned int PSOCICY() { return DIM()!=Dimension::One   ? (internal_cells_per_block_and_dimension_/2) : 1; }
   static constexpr unsigned int PSOCICZ() { return DIM()==Dimension::Three ? (internal_cells_per_block_and_dimension_/2) : 1; }

   /**
    * @brief Gives the index in the parent which 'lies over' the low-index childrens halo's first cell.
    *        "PIOLCH = Parent Index Overlaying Low-Index Child Halo".
    * @return Index in parent of first halo cell of low-index child.
    */
   static constexpr unsigned int PIOLCH() { return (halo_width_/2); }

   /**
    * @brief Gives the index in the parent which 'lies over' the high-index childrens halo's first cell.
    *        "PIOHCH = Parent Index Overlaying High-Index Child Halo".
    * @return Index in parent of first halo cell of high-index child.
    */
   static constexpr unsigned int PIOHCH() { return ((cells_per_dimension_with_halo_/2) - (halo_width_/2)); }

   /**
    * @brief Gives the order coefficient, i.e. order of the Finite Volume scheme and the time integration.
    *        The coefficient is used to adjust the coarsening/refinement threshold. "MROC = Multi-Resolution Order Coefficient".
    * @return Order coefficient.
    */
   static constexpr unsigned int STDO() { return space_time_discretization_order_; }

   /**
    * @brief The maximum number of levels that can be used in the current implementation of the algorithm. Bottleneck is the geometry
    *        feature of the node indexing. Counts level zero as level! $Convenience function in order to make code changes only once and locally$
    *        AMNL = Algorithm Maximum Number of Level".
    * @return Number of Levels possible to simulate.
    */
   static constexpr unsigned int AMNL() {return algorithmic_maximum_number_of_levels_;}

   /**
    * @brief Gives a bool to decide if the end time is limited to the exact value specified in the inputfile.
    * @return End time decision.
    */
   static constexpr bool LET() { return limit_end_time_; }

   /**
    * @brief Gives the minimum time-step size as specified by the user.
    * @return Minimum time-step size.
    */
   static constexpr double MTS() { return minimum_time_step_size_; }

   /**
    * @brief Gives the number of children when the block is refined
    * @return Number of children = 2 (1D), 4 (2D), or 8 (3D)
    */
   static constexpr unsigned int NOC() { return number_of_children_; }

   /**
    * @brief Gives a bool to decide whether a list of the timestep sizes is written to a file or not.
    * @return Plot pressure decision for the current build.
    */
   static constexpr bool WTL() { return write_timestep_list_; }

   /**
    * @brief Gives a bool to indicate, whether or not changes to make the calculation fully symmetric are active.
    *        These changes might affect performance.
    * @return Symmetry decision for the current build.
    */
   static constexpr bool FUSY() { return full_symmetry_active_; }

   /**
    * @brief Gives the Cut-off factor for the levelset function. "LSCOF = Levelset Cut Off Factor".
    * @return Levelset cut-off factor.
    */
   static constexpr double LSCOF() { return levelset_cutoff_factor_; }

   /**
    * @brief Gives the Mixing Threshold for mixing and extending into cut cells. "MITH = MIxing ThresHold".
    * @return mixing threshold value.
    */
   static constexpr double MITH() { return mixing_threshold_; }

   /**
    * @brief Gives the size of the extension band. "EBW = Extension band width".
    * @return Extension band width.
    */
   static constexpr unsigned int EBW() { return extension_band_; }

   /**
    * @brief Gives the size of the reinitialization band. "RBW = Reinitialization band width".
    * @return Extension band width.
    */
   static constexpr unsigned int RBW() { return reinitialization_band_; }

   /**
    * @brief Gives a list of all boundary sides. "NBS = Natural Boundary Sides".
    * @return Array contains East, West, North, South, Top and Bottom.
    */
   static constexpr std::array<BoundaryLocation,6> NBS() { return natural_boundary_sides_; }

   /**
    * @brief Gives a list of all active boundary sides. "ANBS = Active Natural Boundary Sides"
    * @return East, West in 1D.
    * @return East, West, North, South in 2D.
    * @return East, West, North, South, Top, Bottom in 3D.
    */
   static constexpr std::array<BoundaryLocation, active_natural_boundary_sides_.size()> ANBS() { return active_natural_boundary_sides_; }

   /**
    * @brief Gives a list of all active boundary sides. "HBS = Halo Boundary Sides"
    * @return East, West in 1D.
    * @return East, West, North, South + 4 corners in 2D.
    * @return East, West, North, South, Top, Bottom + 12 Sticks + 8 Corners in 3D.
    */
   static constexpr std::array<BoundaryLocation,halo_boundary_sides_.size()> HBS() { return halo_boundary_sides_; }

   /**
    * @brief Gives the norm to be used for the wavelet analysis. "NFWA = Norm For Wavelet Analysis".
    * @return Norm type.
    */
   static constexpr Norm NFWA() { return norm_for_wavelet_analysis_; }

   /**
    * @brief Indicates whether or not duplicated vertices in the output files should be filtered. Can affect performance.
    * @return Output vertex filter decision.
    */
   static constexpr VertexFilterType VERTEX_FILTER() { return output_vertex_filter_; }

   /**
    * @brief Gives the number of topology changes that are allowed on each rank (refinements, coarsenings) before load load balancing
    * @return Number of topology changes
    */
   static constexpr unsigned int TCULB() { return topology_changes_until_load_balancing_; }
};

using CC = CompileTimeConstants;


#endif // COMPILE_TIME_CONSTANTS_H
