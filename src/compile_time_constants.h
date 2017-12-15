/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*  \\\\                                                                                  *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA                                                                                 *
* Copyright (c) 2017 Nikolaus A. Adams and contributors (see AUTHORS list)               *
* All rights reserved.                                                                   *
*                                                                                        *
* Chair of Aerodynamics and Fluid Mechanics                                              *
* Technical University of Munich                                                         *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
* This project has received funding from the European Reseach Council (ERC)              *
* under the European Union's Horizon 2020 research and innovation programme              *
* (grant agreement No 667483).                                                           *
*                                                                                        *
* ERC Advanced Grant No 667483, Prof. Dr. Nikolaus A. Adams:                             *
* "NANOSHOCK - Manufacturing Shock Interactions for Innovative Nanoscale Processes"      *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Redistribution and use in source and binary forms, with or without                     *
* modification, are permitted provided that the following conditions are met:            *
*                                                                                        *
* 1. Redistributions of source code must retain the above copyright notice,              *
*    this list of conditions and the following disclaimer.                               *
*                                                                                        *
* 2. Redistributions in binary form must reproduce the above copyright notice            *
*    this list of conditions and the following disclaimer in the documentation           *
*    and/or other materials provided with the distribution.                              *
*                                                                                        *
* 3. Neither the name of the copyright holder nor the names of its                       *
*    contributors may be used to endorse or promote products derived from this           *
*    software without specific prior written permission.                                 *
*                                                                                        *
* 4. Any redistribution of substantial fractions of the code as a                        *
*    different project should preserve the word ALPACA in the name                       *
*    of the code                                                                         *
*                                                                                        *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"            *
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE              *
* IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE            *
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE              *
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR                    *
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF                   *
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS               *
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN                *
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)                *
* ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE            *
* POSSIBILITY OF SUCH DAMAGE.                                                            *
*                                                                                        *
* Please note, several third-party tools are used within the ALPACA code under           *
* their own license agreement.                                                           *
*                                                                                        *
* 1. xdmf_writer        : Licensed by Technische Universitaet Muenchen                   *
*                         See 'COPYING_XDMF_WRITER' for more information.                *
*                                                                                        *
* 2. tiny_xml           : This software is provided 'as-is', without any express or      *
*                         implied warranty. In no event will the authors be held         *
*                         liable for any damages arising from the use of this software.  *
*                         See COPYING_TINY_XMLfor more information.                      *
*                                                                                        *
* 3. expression_toolkit : Free use of The C++ Mathematical Expression Toolkit Library is *
*                         permitted under the guidelines and in accordance with the most *
*                         current version of the Common Public License.                  *
*                         http://www.opensource.org/licenses/cpl1.0.php                  *
*                         See COPYING_EXPRESSION_TOOLKITfor more information.            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* AUTHORS                                                                                *
*                                                                                        *
*   Prof. Dr. Nikolaus A. Adams                                                          *
*                                                                                        *
*   Dr. Stefan Adami                                                                     *
*   Vladimir Bogdanov                                                                    *
*   Nico Fleischmann                                                                     *
*   Nils Hoppe                                                                           *
*   Naeimeh Hosseini                                                                     *
*   Jakob Kaiser                                                                         *
*   Aleksandr Lunkov                                                                     *
*   Thomas Paula                                                                         *
*   Josef Winter                                                                         *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
*   nanoshock@aer.mw.tum.de                                                              *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, December 15th 2017                                                             *
*                                                                                        *
*****************************************************************************************/

#ifndef COMPILE_TIME_CONSTANTS_H
#define COMPILE_TIME_CONSTANTS_H

#include "enums/dimension_definition.h"
#include "enums/flux_splitting.h"

// double expansion to get string from compiler variable:
#define TOSTRING1(str) #str
#define TOSTRING(str) TOSTRING1(str)

/**
 * @brief The CompileTimeConstants class holds all constants which are determined at compile time. Also the central place to change setting.
 */
class CompileTimeConstants {

    // Macro is set through Makefile. Determines numbers of spatial dimensions considered in the simulation.
    #if DIMENSION == 1
        static constexpr Dimension dimension_of_simulation = Dimension::One;
    #elif DIMENSION == 2
        static constexpr Dimension dimension_of_simulation = Dimension::Two;
    #else
        #define HILBERT //Hilbert Load Balancing is only meanigfull in 3D Simulations
        static constexpr Dimension dimension_of_simulation = Dimension::Three;
    #endif


    static constexpr FluxSplitting flux_splitting_scheme_ = FluxSplitting::Roe;

    static constexpr unsigned short start_ = 0;
    static constexpr unsigned int internal_cells_per_block_and_dimension_ = INTERNALCELLS; //Referred to as "IC"
    static constexpr unsigned int halo_width_ = 4;
    static constexpr unsigned int cells_per_dimension_with_halo_ = 2*halo_width_+internal_cells_per_block_and_dimension_; // Refered to as "TC"
    /* The number of neighbors considered (in each direction) during the predicton. E.g. prediction_stencil = 2 means 5 cells  per dimesnion are condisidered:
     * The base the two left and two right neibors
     */
    static constexpr unsigned int prediction_stencil_size_ = 2;

    static_assert((internal_cells_per_block_and_dimension_ % 4) == 0, "IC must be a multiple of four, stupid!");
    static_assert(internal_cells_per_block_and_dimension_ > 0, "IC must be a greater zero, stupid!");
    static_assert(internal_cells_per_block_and_dimension_ < 32768, "IC must be smaller than 2^15 to fit signed int");
    static_assert((halo_width_ % 2) == 0, "Halo Width should be divisable by two");

    static constexpr unsigned short multi_resolution_order_coefficient_ = 2;

    static constexpr unsigned short algorithmic_maximum_number_of_levels_ = 14; //MUST NOT BE CHANGED!

    static constexpr bool gravitation_active_ = true;
    static constexpr bool viscosity_active_ = true;

    static constexpr bool debug_printing_ = false;
    static constexpr bool debug_plotting_ = false;

    static constexpr bool testing_ = true;
    static constexpr bool profiling_ = false;

    static constexpr bool use_kaiser_adaptive_local_timestepping_ = true;

    static constexpr unsigned short number_of_children_  = dimension_of_simulation == Dimension::One ? 2 : (dimension_of_simulation == Dimension::Two ? 4 : 8);
    static constexpr unsigned short number_of_equations_ = dimension_of_simulation == Dimension::One ? 3 : (dimension_of_simulation == Dimension::Two ? 4 : 5);

public:

    /**
     * @brief Gives the hash of the git commit this build is based on
     * @return Git commit hash
     */
    static constexpr inline const char* GIT_COMMIT() {return TOSTRING(GITHASH);}

    /**
     * @brief Gives the entry of density in the block-buffer-data
     * @return Index of density
     */
    static constexpr inline unsigned short ID_RHO() {return 0;}

    /**
       * @brief Gives the entry of energy in the block-buffer-data
       * @return Index of energy
       */
    static constexpr inline unsigned short ID_ENERGY() {return 1;}

    /**
       * @brief Gives the entry of the x-momentum in the block-buffer-data
       * @return Index of x-momentum
       */
    static constexpr inline unsigned short ID_XMOM() {return 2;}

    /**
       * @brief Gives the entry of the y-momentum in the block-buffer-data
       * @return Index of y-momentum
       */
    static constexpr inline unsigned short ID_YMOM() {return 3;}

    /**
       * @brief Gives the entry of the z-momentum in the block-buffer-data
       * @return Index of z-momentum
       */
    static constexpr inline unsigned short ID_ZMOM() {return 4;}

    /**
     * @brief Gives the number of equations considered in the simulation, i.e. Euler Equations.
     * @return Number of Equations = 5 (rho, E, X-,Y-,Z-Momentum) for 3D
     * @return Number of Equations = 4 (rho, E, X-,Y-Momentum) for 2D
     * @return Number of Equations = 3 (rho, E, X-Momentum) for 1D
     */
    static constexpr inline unsigned short NoEq() {return number_of_equations_;}

    /**
     * @brief Gives the dimension of the Simulation, i.e. 1D, 2D or 3D
     * @return Dimension.
     */
    static constexpr inline Dimension DIM() {return dimension_of_simulation;}

    /**@{
     * @brief Gives the number of Internal Cells in a block per dimension.
     * @return Number of internal cells. 1 if dimension does not exist
     */
    static constexpr inline unsigned int ICX() {return internal_cells_per_block_and_dimension_;}
    static constexpr inline unsigned int ICY() {return dimension_of_simulation!=Dimension::One   ? internal_cells_per_block_and_dimension_ : 1;}
    static constexpr inline unsigned int ICZ() {return dimension_of_simulation==Dimension::Three ? internal_cells_per_block_and_dimension_ : 1;}
    /**@}*/

    /**@{
     * @brief Gives the number of Total Cells in a block per dimension, i.e. number of internal cells per dimension + 2 * number of halo cells per dimension.
     * @return Number of total cells. 1 if dimension does not exist
     */
    static constexpr inline unsigned int TCX() {return cells_per_dimension_with_halo_;}
    static constexpr inline unsigned int TCY() {return dimension_of_simulation!=Dimension::One   ? cells_per_dimension_with_halo_ : 1;}
    static constexpr inline unsigned int TCZ() {return dimension_of_simulation==Dimension::Three ? cells_per_dimension_with_halo_ : 1;}
    /**@}*/

    /**@{
     * @brief Gives the index of the First Internal Cell in a block per dimension.
     * @return Index of first internal cell in block. 0 if dimension does not exist
     */
    static constexpr inline unsigned int FICX() {return halo_width_;}
    static constexpr inline unsigned int FICY() {return dimension_of_simulation!=Dimension::One   ? halo_width_ : 0;}
    static constexpr inline unsigned int FICZ() {return dimension_of_simulation==Dimension::Three ? halo_width_ : 0;}
    /**@}*/

    /**@{
     * @brief Gives the index of the Last Internal Cell in a block per dimension. I.e. the returned index must be included if the internal cells are of interest.
     * @return Index of the last internal cell in a block. 0 if dimension does not exist
     */
    static constexpr inline unsigned int LICX() {return halo_width_ + internal_cells_per_block_and_dimension_ - 1;}
    static constexpr inline unsigned int LICY() {return dimension_of_simulation!=Dimension::One   ? (halo_width_ + internal_cells_per_block_and_dimension_ - 1) : 0;}
    static constexpr inline unsigned int LICZ() {return dimension_of_simulation==Dimension::Three ? (halo_width_ + internal_cells_per_block_and_dimension_ - 1) : 0;}
    /**@}*/

    /**
     * @brief Gives the size of the halo "HS = Halo Size" in its shortest direction.
     * @return Number of cells in the halo in shortest halo dimension.
     */
    static constexpr inline unsigned int HS() {return halo_width_;}

    /**
     * @brief Gives the index of the first cell in the high-index halo "FHH = First High-index Halo",
     *        i.e. the eastern halo in X-Direction, northern in Y-Direction and top in Z-Direction.
     * @return Start index of the high-index halo.
     */
    static constexpr inline unsigned int FHH() {return halo_width_ + internal_cells_per_block_and_dimension_;}

    /**
     * @brief Indicates if gravity should be considered as source term.
     * @return True if gravity source terms are to be calculated. False otherwise.
     */
    static constexpr inline bool GravityIsActive() {return gravitation_active_;}

    /**
     * @brief Indicates if viscosity should be considered as source term.
     * @return True if viscosity source terms are to be calculated. False otherwise.
     */
    static constexpr inline bool ViscosityIsActive() {return viscosity_active_;}

    /**
     * @brief Gives the number of the cells considered in the parent for a prediction into the childs halo (long side). "PHS = Prediction Halo Size".
     * @return Number of cells for prediction per dimension.
     */
    static constexpr inline unsigned int PHS() {return (cells_per_dimension_with_halo_/2) + 2 * prediction_stencil_size_;}

    /**
     * @brief Gives the number of cells considered in the parent during a prediction into the childs halo,
     *        e.g. at internal jumps. "BPS = Boundary Prediction Size".
     * @return Number of cells for boundary prediction per dimension.
     */
    static constexpr inline unsigned int BPS() {return (halo_width_/2) + (prediction_stencil_size_*2);}

    /**
     * @brief Gives the start index in the parent used during prediction to the high-index child, i.e. the eastern child in X-Direction.
     *        "PHCS = Prediction High-index Child Start index".
     * @return Starting index in Parent for prediction to high-index child.
     */
    static constexpr inline unsigned int PHCS() {return ( (cells_per_dimension_with_halo_/2) - prediction_stencil_size_ - (halo_width_/2) );}

    /**
     * @brief Gives the start index in the parent used during prediction to the high-index child's halo.
     *        "BPHCS = Boundary Prediction High-index Child Start index".
     * @return Starting index in parent for prediction to high-index child's halo.
     */
    static constexpr inline unsigned int BPHCS() {return halo_width_ + internal_cells_per_block_and_dimension_ - (halo_width_/2);}

    /**
     * @brief Gives the start index in the parent used during prediction to the low-index child.
     *        "PLCS = Prediction Low-index Child Start index"
     * @return Start index in the parent for prediction to low-index child.
     */
    static constexpr inline unsigned int PLCS() {return (halo_width_/2) - prediction_stencil_size_;}

    /**
     * @brief Gives the index in the parent that 'lies over' the first internal cell of the low-index child.
     *        "PIOLCFIC = Parent Index Overlaying Low-index Child First Internal Cell".
     * @return Index in parent of first internal cell of low-index child.
     */
    static constexpr inline unsigned int PIOLCFIC() {return halo_width_;}

    /**
     * @brief Gives the index in the parent which 'lies over' the first internal cell of the high-index child.
     *        "PIOHCFIC = Parent Index Overlaying High-Index Childs First Internal Cell".
     * @return Index in parent of first internal cell of high-index child.
     */
    static constexpr inline unsigned int PIOHCFIC() {return cells_per_dimension_with_halo_/2;}

    /**
     * @brief Gives the number of parent cells which 'lies over' the internal cells of a child.
     *        "PSOCIC = Parent Size Overlaying Childrens Internal Cells".
     * @return Number of cells in parent which span over all internal cells of a child.
     */
    static constexpr inline unsigned int PSOCIC() {return (internal_cells_per_block_and_dimension_/2);}

    /**
     * @brief Gives the index in the parent which 'lies over' the low-index childrens halo's first cell.
     *        "PIOLCH = Parent Index Overlaying Low-Index Child Halo".
     * @return Index in parent of first halo cell of low-index child.
     */
    static constexpr inline unsigned int PIOLCH() {return (halo_width_/2);}

    /**
     * @brief Gives the index in the parent which 'lies over' the high-index childrens halo's first cell.
     *        "PIOHCH = Parent Index Overlaying High-Index Child Halo".
     * @return Index in parent of first halo cell of high-index child.
     */
    static constexpr inline unsigned int PIOHCH() {return ((cells_per_dimension_with_halo_/2) - (halo_width_/2));}

    /**
     * @brief Gives the order coefficient, i.e. order of the Finite Volume scheme and the time integration.
     *        The coefficient is used to adjust the coarsening/refinement threshold. "MROC = Multi-Resolution Order Coefficient".
     * @return Order coefficient.
     */
    static constexpr inline unsigned short MROC() {return multi_resolution_order_coefficient_;}

    /**
     * @brief The maximum number of levels that can be used in the current implementation of the algorithm. Bottleneck is the geometry
     *        feature of the node indexing. Counts level zero as level! $Convenience function in order to make code changes only once and locally$
     *        ALNM = Algorithm Maximum Number of Level".
     * @return Number of Levels possible to simulate.
     */
    static constexpr inline unsigned short AMNL() {return algorithmic_maximum_number_of_levels_;}

    /**
     * @brief Gives a bool to decide if Kaiser Adaptive Local Time Stepping (KALTS) or the standard local time stepping scheme is to be used.
     * @return Time stepping decision
     */
    static constexpr inline bool KALTS() {return use_kaiser_adaptive_local_timestepping_;}

    /**
     * @brief Gives the number of children when the block is refined
     * @return Number of children = 2 (1D), 4 (2D), or 8 (3D)
     */
    static constexpr inline unsigned short NOC() {return number_of_children_;}

    /**
     * @brief Gives a bool to decide if additional debugging information should be provided. This function is not about debugging information for e.g.
     *        GNU Debugger, but to print additional output during a debug run or to create additional plots after separate functions. $Convenience function.
     *        The compiler should recognize unset debug if clauses and take them out in release-builds$
     * @return Debug logging decision for the current build.
     */
    static constexpr inline bool DEBUG_PRINT(){return debug_printing_;}

    /**
     * @brief Gives a bool to decide if additional debug plots are to be written out.
     * @return Debug plotting decision for the current build.
     */
    static constexpr inline bool DEBUG_PLOT(){return debug_plotting_;}

    /**
     * @brief Gives a bool to indicate, whether or not test (or verification) runs are being run. $Convenience function.
     *        The compiler should recognize unset debug if clauses and take them out in release-builds$
     * @return Test run decision for the current build.
     */
    static constexpr inline bool TEST() {return testing_;}

    /**
     * @brief Give a bool to indicate whether or not profiling (timing) runs are beeing run.
     * @return Profiling decision.
     */
    static constexpr inline bool PROFILE() {return profiling_;}


    /**
     * @brief Gives the flux splitting scheme of the Simulation, i.e. Roe, global Lax-Friedrichs, or local Lax-Friedrichs
     * @return FluxSplitting scheme.
     */
    static constexpr inline FluxSplitting FSS() {return flux_splitting_scheme_;}

};

typedef CompileTimeConstants CC;

#endif // COMPILE_TIME_CONSTANTS_H
