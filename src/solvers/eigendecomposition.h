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
#ifndef EIGENVALUE_CALCULATOR_H
#define EIGENVALUE_CALCULATOR_H

#include <cmath>
#include "utilities/mathematical_functions.h"

#include "materials/material_manager.h"
#include "block_definitions/block.h"
#include "user_specifications/compile_time_constants.h"
#include "user_specifications/riemann_solver_settings.h"
#include "enums/direction_definition.h"

/**
 * @brief The RoeEigenvalues class computes Roe eigenvalues and eigenvectors within one block. For further information consult
 * \cite Roe1981.
 * For characteristic decomposition, especially the eigenvectors used here, consult \cite Fedkiw1999a.
 */
class EigenDecomposition {

   const MaterialManager& material_manager_;

   static constexpr std::array<std::array<unsigned int, 3>, 3> direction_momentum_indices_ = {{
      {0, 1, 2},
      {1, 0, 2},
      {2, 0, 1}
   }};
   static constexpr std::array<std::array<unsigned int, 3>, 3> direction_eigenvector_indices_ = {{
      {1, 2, 3},
      {2, 1, 3},
      {3, 1, 2}
   }};

   // Using static to cheat constness (for this global variable okay - don't do this at home (we are what you a call Experts) ;)
   static double global_eigenvalues_[DTI(CC::DIM())][MF::ANOE()];

public:
   EigenDecomposition() = delete;
   explicit EigenDecomposition( MaterialManager const& material_manager);
   ~EigenDecomposition() = default;
   EigenDecomposition( EigenDecomposition const& ) = delete;
   EigenDecomposition& operator=( EigenDecomposition const& ) = delete;
   EigenDecomposition( EigenDecomposition&& ) = delete;
   EigenDecomposition& operator=( EigenDecomposition&& ) = delete;

   template<Direction DIR>
   void ComputeRoeEigendecomposition( std::pair<MaterialName const, Block> const& mat_block,
      double (&roe_eigenvectors_left)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][MF::ANOE()][MF::ANOE()],
      double (&roe_eigenvectors_right)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][MF::ANOE()][MF::ANOE()],
      double (&fluxfunction_eigenvalues)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][MF::ANOE()] ) const;

   void ComputeMaxEigenvaluesOnBlock(std::pair<MaterialName const, Block> const& mat_block, double (&eigenvalues)[DTI(CC::DIM())][MF::ANOE()]) const;
   void SetGlobalEigenvalues(double (&eigenvalues)[DTI(CC::DIM())][MF::ANOE()]) const;
   auto GetGlobalEigenvalues() const -> double const (&)[DTI(CC::DIM())][MF::ANOE()];
};

namespace {
   /**
    * @brief Save the wavespeeds in the fluxfunction_eigenvalues buffer.
    * Taking care of num-of-eq and multiple material waves for higher dimensions.
    * @param eigenvalues (indirect return parameter) .
    * @param u_minus_c .
    * @param u .
    * @param u_plus_c .
    */
   inline void SaveForAllFields(double (&eigenvalues)[MF::ANOE()], double const u_minus_c, double const u, double const u_plus_c){
      eigenvalues[0]            = u_minus_c;
      eigenvalues[MF::ANOE()-1] = u_plus_c;
      for(unsigned int l = 1; l < MF::ANOE()-1; ++l) {
         eigenvalues[l] = u;
      }
   }
}

/**
 * @brief Computes the Roe left and right eigenvectors and the Roe eigenvalues in x-direction according to \cite Fedkiw1999a.
 * @param mat_block The block and material information of the phase under consideration.
 * @param roe_eigenvectors_left Reference to an array which is filled with the computed eigenvectors (indirect return parameter).
 * @param roe_eigenvectors_right Reference to an array which is filled with the computed eigenvectors (indirect return parameter).
 * @param fluxfunction_eigenvalues Reference to an array which is filled with the computed eigenvalues (indirect return parameter).
 * @note Hotpath function.
 */
template<Direction DIR>
void EigenDecomposition::ComputeRoeEigendecomposition( std::pair<MaterialName const, Block> const& mat_block,
   double (&roe_eigenvectors_left)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][MF::ANOE()][MF::ANOE()],
   double (&roe_eigenvectors_right)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][MF::ANOE()][MF::ANOE()],
   double (&fluxfunction_eigenvalues)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][MF::ANOE()] ) const {

   constexpr int total_to_internal_offset_x = CC::FICX() - 1;
   constexpr int total_to_internal_offset_y = CC::DIM() != Dimension::One   ? int(CC::FICY()) - 1 : -1;
   constexpr int total_to_internal_offset_z = CC::DIM() == Dimension::Three ? int(CC::FICZ()) - 1 : -1;

   constexpr unsigned int start_x = DIR == Direction::X ? CC::FICX() - 1 : CC::FICX();
   constexpr unsigned int start_y = DIR == Direction::Y ? CC::FICY() - 1 : CC::FICY();
   constexpr unsigned int start_z = DIR == Direction::Z ? CC::FICZ() - 1 : CC::FICZ();

   constexpr unsigned int x_varying = DIR == Direction::X ? 1 : 0;
   constexpr unsigned int y_varying = DIR == Direction::Y ? 1 : 0;
   constexpr unsigned int z_varying = DIR == Direction::Z ? 1 : 0;

   // indices of principal and first/second minor momentum/velocity within the three-packs MF::AME and MF::AV
   constexpr unsigned int principal = direction_momentum_indices_[DTI(DIR)][0];
   constexpr unsigned int minor1    = direction_momentum_indices_[DTI(DIR)][1];
   constexpr unsigned int minor2    = direction_momentum_indices_[DTI(DIR)][2];

   // Indices of eigenvectors for characteristic fields due to principal and first/second minor momentum
   constexpr unsigned int ev_principal = direction_eigenvector_indices_[DTI(DIR)][0];
   constexpr unsigned int ev_minor1    = direction_eigenvector_indices_[DTI(DIR)][1];
   constexpr unsigned int ev_minor2    = direction_eigenvector_indices_[DTI(DIR)][2];

   // Access the pair's elements directly.
   auto const& [material, block] = mat_block;

   // We need to use the conservative buffer for density as the prime state is only consistent (if zero) after last RK stage
   Conservatives const& conservatives = block.GetAverageBuffer();
   double const    (&density)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(Equation::Mass);
   double const     (&energy)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(Equation::Energy);

   PrimeStates const& prime_states = block.GetPrimeStateBuffer();
   double const   (&pressure)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetPrimeStateBuffer(PrimeState::Pressure);

   double const gruneisen_coefficient_material = CC::GruneisenDensityDependent() ? 0.0 : material_manager_.GetMaterial(material).GetEquationOfState().GetGruneisen();

   for( unsigned int i = start_x; i <= CC::LICX(); ++i ) {
      for( unsigned int j = start_y; j <= CC::LICY(); ++j ) {
         for( unsigned int k = start_z; k <= CC::LICZ(); ++k ) {
            // Indices of neighbor cell
            unsigned int const in = i + x_varying;
            unsigned int const jn = j + y_varying;
            unsigned int const kn = k + z_varying;
            /**
             * This if statement is necessary due to the ghost-fluid method. In ghost-material cells which do not lie on the extension band, i.e. therein
             * we do not have extended or integrated values, the density is zero. Therefore, we cannot compute Roe eigenvalues in those cells.
             */
            if(density[i][j][k] <= 0.0 || density[in][jn][kn] <= 0.0) continue;

            // References to write directly to the buffers of the currently considered cell face
            auto& left_eigenvector  = roe_eigenvectors_left   [i - total_to_internal_offset_x][j - total_to_internal_offset_y][k - total_to_internal_offset_z];
            auto& right_eigenvector = roe_eigenvectors_right  [i - total_to_internal_offset_x][j - total_to_internal_offset_y][k - total_to_internal_offset_z];
            auto& eigenvalues       = fluxfunction_eigenvalues[i - total_to_internal_offset_x][j - total_to_internal_offset_y][k - total_to_internal_offset_z];

            // extract required conservatives and primes
            double const rho_target                 =       density[i][j][k];
            double const one_rho_target             = 1.0 / density[i][j][k];
            double const energy_target              =     energy[i][j][k];
            double const x_momentum_target          = conservatives[Equation::MomentumX][i][j][k];
            double const y_momentum_target          = CC::DIM() != Dimension::One ? conservatives[Equation::MomentumY][i][j][k] : 0.0;
            double const z_momentum_target          = CC::DIM() == Dimension::Three ? conservatives[Equation::MomentumZ][i][j][k] : 0.0;
            double const principal_velocity_target  = prime_states[MF::AV()[principal]][i][j][k];
            double const minor1_velocity_target     = CC::DIM() != Dimension::One   ? prime_states[MF::AV()[minor1]][i][j][k] : 0.0;
            double const minor2_velocity_target     = CC::DIM() == Dimension::Three ? prime_states[MF::AV()[minor2]][i][j][k] : 0.0;
            double const pressure_target            = pressure[i][j][k];
            double const generalized_psi_target     = material_manager_.GetMaterial(material).GetEquationOfState().GetPsi(pressure_target, one_rho_target);
            double const gruneisen_target           = CC::GruneisenDensityDependent() ? material_manager_.GetMaterial(material).GetEquationOfState().GetGruneisen(rho_target) : 0.0;

            double const rho_neighbor               =       density[in][jn][kn];
            double const one_rho_neighbor           = 1.0 / density[in][jn][kn];
            double const energy_neighbor            =     energy[in][jn][kn];
            double const x_momentum_neighbor        = conservatives[Equation::MomentumX][in][jn][kn];
            double const y_momentum_neighbor        = CC::DIM() != Dimension::One ? conservatives[Equation::MomentumY][in][jn][kn] : 0.0;
            double const z_momentum_neighbor        = CC::DIM() == Dimension::Three ? conservatives[Equation::MomentumZ][in][jn][kn] : 0.0;
            double const principal_velocity_neighbor= prime_states[MF::AV()[principal]][in][jn][kn];
            double const minor1_velocity_neighbor   = CC::DIM() != Dimension::One   ? prime_states[MF::AV()[minor1]][in][jn][kn] : 0.0;
            double const minor2_velocity_neighbor   = CC::DIM() == Dimension::Three ? prime_states[MF::AV()[minor2]][in][jn][kn] : 0.0;
            double const pressure_neighbor          =   pressure[in][jn][kn];
            double const generalized_psi_neighbor   = material_manager_.GetMaterial(material).GetEquationOfState().GetPsi(pressure_neighbor, one_rho_neighbor);
            double const gruneisen_neighbor         = CC::GruneisenDensityDependent() ? material_manager_.GetMaterial(material).GetEquationOfState().GetGruneisen(rho_neighbor) : 0.0;

            // compute expensive and frequently used temporaries
            double const sqrt_rho_target   = std::sqrt(rho_target);
            double const sqrt_rho_neighbor = std::sqrt(rho_neighbor);
            double const rho_div = 1.0 / ( sqrt_rho_target + sqrt_rho_neighbor );
            double const density_roe_ave = sqrt_rho_target * sqrt_rho_neighbor;
            double const one_density_roe_ave = 1.0 / density_roe_ave;

            // compute Roe averages
            double const principal_velocity_roe_ave = ((principal_velocity_target * sqrt_rho_target) + (principal_velocity_neighbor * sqrt_rho_neighbor)) * rho_div;
            double const minor1_velocity_roe_ave = ((minor1_velocity_target * sqrt_rho_target) + (minor1_velocity_neighbor * sqrt_rho_neighbor)) * rho_div;
            double const minor2_velocity_roe_ave = ((minor2_velocity_target * sqrt_rho_target) + (minor2_velocity_neighbor * sqrt_rho_neighbor)) * rho_div;
            double const generalized_psi_roe_ave = ((generalized_psi_target * sqrt_rho_target) + (generalized_psi_neighbor * sqrt_rho_neighbor)) * rho_div;
            double const gruneisen_roe_ave = CC::GruneisenDensityDependent() ? ((gruneisen_target * sqrt_rho_target) + (gruneisen_neighbor * sqrt_rho_neighbor)) * rho_div : gruneisen_coefficient_material;
            double const enthalpy_roe_ave = ( material_manager_.GetMaterial(material).GetEquationOfState().GetEnthalpy(rho_target,   x_momentum_target,   y_momentum_target,   z_momentum_target,   energy_target )   * sqrt_rho_target
                                            + material_manager_.GetMaterial(material).GetEquationOfState().GetEnthalpy(rho_neighbor, x_momentum_neighbor, y_momentum_neighbor, z_momentum_neighbor, energy_neighbor ) * sqrt_rho_neighbor)
                                            * rho_div;

            // Absolute roe averaged velocity, speed of sound
            double const q_squared  = DimensionAwareConsistencyManagedSum( principal_velocity_roe_ave * principal_velocity_roe_ave, minor1_velocity_roe_ave * minor1_velocity_roe_ave, minor2_velocity_roe_ave * minor2_velocity_roe_ave);
            double const tmp = DimensionAwareConsistencyManagedSum( (principal_velocity_target - principal_velocity_neighbor) * (principal_velocity_target - principal_velocity_neighbor),
                                                         (minor1_velocity_target    - minor1_velocity_neighbor)    * (minor1_velocity_target    - minor1_velocity_neighbor),
                                                         (minor2_velocity_target    - minor2_velocity_neighbor)    * (minor2_velocity_target    - minor2_velocity_neighbor) );
            double const pressure_over_density_roe_ave = ((pressure_target*one_rho_target) * sqrt_rho_target + (pressure_neighbor*one_rho_neighbor) * sqrt_rho_neighbor) * rho_div
                                                           + 0.5 * density_roe_ave * (rho_div * rho_div) * tmp;
            double const cc = generalized_psi_roe_ave + gruneisen_roe_ave * pressure_over_density_roe_ave;
            double const one_cc = 1.0 / cc;
            double const c = std::sqrt(cc);

            // LEFT EIGENVECTORS **********************************************
            // ****************************************************************

            // Eigenvector for lambda = u-c
            left_eigenvector[         0     ][ETI(Equation::Mass)       ] = 0.5 * one_cc * (gruneisen_roe_ave*q_squared - gruneisen_roe_ave*enthalpy_roe_ave + (principal_velocity_roe_ave+c) * c);
            left_eigenvector[         0     ][ETI(MF::AME()[principal]) ] = 0.5 * one_cc * (-principal_velocity_roe_ave * gruneisen_roe_ave - c);
            if constexpr( CC::DIM() != Dimension::One )
               left_eigenvector[      0     ][ETI(MF::AME()[minor1])    ] = 0.5 * one_cc * (-minor1_velocity_roe_ave*gruneisen_roe_ave);
            if constexpr( CC::DIM() == Dimension::Three )
               left_eigenvector[      0     ][ETI(MF::AME()[minor2])    ] = 0.5 * one_cc * (-minor2_velocity_roe_ave*gruneisen_roe_ave);
            left_eigenvector[         0     ][ETI(Equation::Energy)     ] = 0.5 * one_cc * gruneisen_roe_ave;

            // Eigenvector for lambda = u related to principal momentum direction
            left_eigenvector[   ev_principal][ETI(Equation::Mass)       ] =  one_cc * (enthalpy_roe_ave - q_squared);
            left_eigenvector[   ev_principal][ETI(MF::AME()[principal]) ] =  one_cc * principal_velocity_roe_ave;
            if constexpr( CC::DIM() != Dimension::One )
               left_eigenvector[ev_principal][ETI(MF::AME()[minor1])    ] =  one_cc * minor1_velocity_roe_ave;
            if constexpr( CC::DIM() == Dimension::Three )
               left_eigenvector[ev_principal][ETI(MF::AME()[minor2])    ] =  one_cc * minor2_velocity_roe_ave;
            left_eigenvector[   ev_principal][ETI(Equation::Energy)     ] = -one_cc;

            if constexpr( CC::DIM() != Dimension::One ) {
               // Additional eigenvector for lambda = u related to second momentum equation (= first minor momentum)
               left_eigenvector[   ev_minor1][ETI(Equation::Mass)       ] = minor1_velocity_roe_ave * one_density_roe_ave;
               left_eigenvector[   ev_minor1][ETI(MF::AME()[principal]) ] = 0.0;
               left_eigenvector[   ev_minor1][ETI(MF::AME()[minor1])    ] = -one_density_roe_ave;
               if constexpr( CC::DIM() == Dimension::Three )
                  left_eigenvector[ev_minor1][ETI(MF::AME()[minor2])    ] = 0.0;
               left_eigenvector[   ev_minor1][ETI(Equation::Energy)     ] = 0.0;
            }

            if constexpr( CC::DIM() == Dimension::Three ) {
               // Additional eigenvector for lambda = u related to third momentum equation (= second minor momentum)
               left_eigenvector[   ev_minor2][ETI(Equation::Mass)       ] = -minor2_velocity_roe_ave * one_density_roe_ave;
               left_eigenvector[   ev_minor2][ETI(MF::AME()[principal]) ] = 0.0;
               left_eigenvector[   ev_minor2][ETI(MF::AME()[minor1])    ] = 0.0;
               left_eigenvector[   ev_minor2][ETI(MF::AME()[minor2])    ] = one_density_roe_ave;
               left_eigenvector[   ev_minor2][ETI(Equation::Energy)     ] = 0.0;
            }

            // eigenvector for lambda = u+c
            left_eigenvector[   MF::ANOE()-1][ETI(Equation::Mass)       ] = 0.5 * one_cc * (gruneisen_roe_ave*q_squared - gruneisen_roe_ave*enthalpy_roe_ave - (principal_velocity_roe_ave-c) * c);
            left_eigenvector[   MF::ANOE()-1][ETI(MF::AME()[principal]) ] = 0.5 * one_cc * (-principal_velocity_roe_ave*gruneisen_roe_ave + c);
            if constexpr( CC::DIM() != Dimension::One )
               left_eigenvector[MF::ANOE()-1][ETI(MF::AME()[minor1])    ] = 0.5 * one_cc * (-minor1_velocity_roe_ave*gruneisen_roe_ave);
            if constexpr( CC::DIM() == Dimension::Three )
               left_eigenvector[MF::ANOE()-1][ETI(MF::AME()[minor2])    ] = 0.5 * one_cc * (-minor2_velocity_roe_ave*gruneisen_roe_ave);
            left_eigenvector[   MF::ANOE()-1][ETI(Equation::Energy)     ] = 0.5 * one_cc * gruneisen_roe_ave;


            // RIGHT EIGENVECTORS *********************************************
            // ****************************************************************

            // Mass equation entries of right eigenvectors
            right_eigenvector[   ETI(Equation::Mass)      ][      0     ] = 1.0;
            right_eigenvector[   ETI(Equation::Mass)      ][ev_principal] = gruneisen_roe_ave;
            if constexpr( CC::DIM() != Dimension::One )
               right_eigenvector[ETI(Equation::Mass)      ][ev_minor1   ] = 0.0;
            if constexpr( CC::DIM() == Dimension::Three )
               right_eigenvector[ETI(Equation::Mass)      ][ev_minor2   ] = 0.0;
            right_eigenvector[   ETI(Equation::Mass)      ][MF::ANOE()-1] = 1.0;

            // principal momentum equation entries
            right_eigenvector[   ETI(MF::AME()[principal])][      0     ] = principal_velocity_roe_ave - c;
            right_eigenvector[   ETI(MF::AME()[principal])][ev_principal] = gruneisen_roe_ave * principal_velocity_roe_ave;
            if constexpr( CC::DIM() != Dimension::One )
               right_eigenvector[ETI(MF::AME()[principal])][ev_minor1   ] = 0.0;
            if constexpr( CC::DIM() == Dimension::Three )
               right_eigenvector[ETI(MF::AME()[principal])][ev_minor2   ] = 0.0;
            right_eigenvector[   ETI(MF::AME()[principal])][MF::ANOE()-1] = principal_velocity_roe_ave + c;

            if constexpr( CC::DIM() != Dimension::One ) {
               // first minor momentum equation entries
               right_eigenvector[   ETI(MF::AME()[minor1])][      0     ] = minor1_velocity_roe_ave;
               right_eigenvector[   ETI(MF::AME()[minor1])][ev_principal] = gruneisen_roe_ave * minor1_velocity_roe_ave;
               right_eigenvector[   ETI(MF::AME()[minor1])][ev_minor1   ] = -density_roe_ave;
               if constexpr( CC::DIM() == Dimension::Three )
                  right_eigenvector[ETI(MF::AME()[minor1])][ev_minor2   ] = 0.0;
               right_eigenvector[   ETI(MF::AME()[minor1])][MF::ANOE()-1] = minor1_velocity_roe_ave;
            }

            if constexpr( CC::DIM() == Dimension::Three ) {
               // second minor momentum equation entries
               right_eigenvector[ETI(MF::AME()[minor2])   ][      0     ] = minor2_velocity_roe_ave;
               right_eigenvector[ETI(MF::AME()[minor2])   ][ev_principal] = gruneisen_roe_ave * minor2_velocity_roe_ave;
               right_eigenvector[ETI(MF::AME()[minor2])   ][ev_minor1   ] = 0.0;
               right_eigenvector[ETI(MF::AME()[minor2])   ][ev_minor2   ] = density_roe_ave;
               right_eigenvector[ETI(MF::AME()[minor2])   ][MF::ANOE()-1] = minor2_velocity_roe_ave;
            }

            // Energy equation entries
            right_eigenvector[   ETI(Equation::Energy)    ][      0     ] = enthalpy_roe_ave - principal_velocity_roe_ave * c;
            right_eigenvector[   ETI(Equation::Energy)    ][ev_principal] = gruneisen_roe_ave * enthalpy_roe_ave - c*c;
            if constexpr( CC::DIM() != Dimension::One )
               right_eigenvector[ETI(Equation::Energy)    ][ev_minor1   ] = -minor1_velocity_roe_ave * density_roe_ave;
            if constexpr( CC::DIM() == Dimension::Three )
               right_eigenvector[ETI(Equation::Energy)    ][ev_minor2   ] = density_roe_ave * minor2_velocity_roe_ave;
            right_eigenvector[   ETI(Equation::Energy)    ][MF::ANOE()-1] = enthalpy_roe_ave + principal_velocity_roe_ave * c;


            // EIGENVALUES ****************************************************
            // ****************************************************************

            // Flux-function artificial viscosity -> Roe eigenvalues
            if constexpr( RoeSolverSettings::flux_splitting_scheme == FluxSplitting::Roe ) {
               SaveForAllFields( eigenvalues,
                  std::abs(principal_velocity_roe_ave - c),
                  std::abs(principal_velocity_roe_ave),
                  std::abs(principal_velocity_roe_ave + c) );
            }

            // Flux-function artificial viscosity -> Roe-M eigenvalues
            if constexpr( RoeSolverSettings::flux_splitting_scheme == FluxSplitting::Roe_M ) {
               SaveForAllFields( eigenvalues,
                  std::abs( principal_velocity_roe_ave - std::min( c, RoeSolverSettings::low_mach_number_limit_factor * std::abs( principal_velocity_roe_ave ) ) ),
                  std::abs( principal_velocity_roe_ave ),
                  std::abs( principal_velocity_roe_ave + std::min( c, RoeSolverSettings::low_mach_number_limit_factor * std::abs( principal_velocity_roe_ave ) ) ) );
            }

            // Flux-function artificial viscosity -> Local-Lax-Friedrichs eigenvalues
            if constexpr( RoeSolverSettings::flux_splitting_scheme == FluxSplitting::LocalLaxFriedrichs ) {
               double const c_target   = material_manager_.GetMaterial(material).GetEquationOfState().GetSpeedOfSound(rho_target,   pressure_target);
               double const c_neighbor = material_manager_.GetMaterial(material).GetEquationOfState().GetSpeedOfSound(rho_neighbor, pressure_neighbor);

               SaveForAllFields( eigenvalues,
                  std::max(std::abs(principal_velocity_target - c_target),std::abs(principal_velocity_neighbor - c_neighbor)),
                  std::max(std::abs(principal_velocity_target), std::abs(principal_velocity_neighbor)),
                  std::max(std::abs(principal_velocity_target + c_target),std::abs(principal_velocity_neighbor + c_neighbor)) );
            }

            // Flux-function artificial viscosity -> LLF-M eigenvalues
            if constexpr( RoeSolverSettings::flux_splitting_scheme == FluxSplitting::LocalLaxFriedrichs_M ) {
               double const c_target   = material_manager_.GetMaterial(material).GetEquationOfState().GetSpeedOfSound(rho_target,   pressure_target);
               double const c_neighbor = material_manager_.GetMaterial(material).GetEquationOfState().GetSpeedOfSound(rho_neighbor, pressure_neighbor);

               double const c_target_m   = std::min( RoeSolverSettings::low_mach_number_limit_factor * std::abs( principal_velocity_target )  , c_target );
               double const c_neighbor_m = std::min( RoeSolverSettings::low_mach_number_limit_factor * std::abs( principal_velocity_neighbor ), c_neighbor );

               SaveForAllFields( eigenvalues,
                  std::max(std::abs(principal_velocity_target - c_target_m),std::abs(principal_velocity_neighbor - c_neighbor_m)),
                  std::max(std::abs(principal_velocity_target), std::abs(principal_velocity_neighbor)),
                  std::max(std::abs(principal_velocity_target + c_target_m),std::abs(principal_velocity_neighbor + c_neighbor_m)) );
            }

            // Flux-function artificial viscosity -> Global-Lax-Friedrichs or LaxFriedrichs scheme
            if constexpr( RoeSolverSettings::flux_splitting_scheme == FluxSplitting::GlobalLaxFriedrichs ) {
               for(unsigned int l = 0; l < MF::ANOE(); ++l) {
                  eigenvalues[l] = global_eigenvalues_[0][l];
               }
            }
         }
      }
   }
}

#endif // EIGENVALUE_CALCULATOR_H
