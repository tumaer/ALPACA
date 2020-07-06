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
#include "eigendecomposition.h"
#include "utilities/mathematical_functions.h"
#include "user_specifications/numerical_setup.h"

#include <cmath>

double EigenDecomposition::global_eigenvalues_[DTI(CC::DIM())][MF::ANOE()];

/**
 * @brief Standard constructor using an already existing MaterialManager.
 * @param material_manager The MaterialManager provides the correct equation of state for a given Material.
 */
EigenDecomposition::EigenDecomposition( MaterialManager const& material_manager) : material_manager_(material_manager)
{
   /*Empty besides initializer list*/
}

/**
 * @brief Computes the global Lax-Friedrichs eigenvalues (maxima) within the given block.
 * @param mat_block The block in which the eigenvalues are to be computed.
 * @param eigenvalues Indirect return parameter.
 */
void EigenDecomposition::ComputeMaxEigenvaluesOnBlock( std::pair<MaterialName const, Block> const& mat_block, double (&eigenvalues)[DTI(CC::DIM())][MF::ANOE()]) const {

   // We save u-c,u,u+c, i.e. three values, for the eigenvalues per direction
   std::array<double,3> max_eigenvalue_x = {0.0};
   std::array<double,3> max_eigenvalue_y = {0.0};
   std::array<double,3> max_eigenvalue_z = {0.0};

   // Access the pair's elements directly.
   auto const& [material, block] = mat_block;

   double const (&density)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetPrimeStateBuffer(PrimeState::Density);

   for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
         for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
            /**
              * This if statement is necessary due to the ghost-fluid method. In ghost-material cells which do not lie on the extension band, i.e. therein
              * we do not have extended or integrated values, the density is zero. Therefore, we cannot compute Lax-Friedrichs eigenvalues in those cells.
              */
            if(density[i][j][k] <= 0.0) continue;

            double const c = material_manager_.GetMaterial(material).GetEquationOfState().GetSpeedOfSound(density[i][j][k], block.GetPrimeStateBuffer(PrimeState::Pressure)[i][j][k]);
            double const u = block.GetPrimeStateBuffer(PrimeState::VelocityX)[i][j][k];

            max_eigenvalue_x[0] = std::max( max_eigenvalue_x[0], std::abs(u - c) );
            max_eigenvalue_x[1] = std::max( max_eigenvalue_x[1], std::abs(u) );
            max_eigenvalue_x[2] = std::max( max_eigenvalue_x[2], std::abs(u + c) );

            if constexpr(CC::DIM() != Dimension::One) {
               double const v = block.GetPrimeStateBuffer(PrimeState::VelocityY)[i][j][k];

               max_eigenvalue_y[0] = std::max( max_eigenvalue_y[0], std::abs(v - c) );
               max_eigenvalue_y[1] = std::max( max_eigenvalue_y[1], std::abs(v) );
               max_eigenvalue_y[2] = std::max( max_eigenvalue_y[2], std::abs(v + c) );
            }

            if constexpr(CC::DIM() == Dimension::Three) {
               double const w = block.GetPrimeStateBuffer(PrimeState::VelocityZ)[i][j][k];

               max_eigenvalue_z[0] = std::max( max_eigenvalue_z[0], std::abs(w - c) );
               max_eigenvalue_z[1] = std::max( max_eigenvalue_z[1], std::abs(w) );
               max_eigenvalue_z[2] = std::max( max_eigenvalue_z[2], std::abs(w + c) );
            }
         }
      }
   }

   SaveForAllFields( eigenvalues[0], max_eigenvalue_x[0], max_eigenvalue_x[1], max_eigenvalue_x[2] );

   if constexpr(CC::DIM() != Dimension::One) {
      SaveForAllFields( eigenvalues[1], max_eigenvalue_y[0], max_eigenvalue_y[1], max_eigenvalue_y[2] );
   }

   if constexpr(CC::DIM() == Dimension::Three) {
      SaveForAllFields( eigenvalues[2], max_eigenvalue_z[0], max_eigenvalue_z[1], max_eigenvalue_z[2] );
   }
}

/**
 * @brief Stores the global Lax-Friedrichs eigenvalues for later usage.
 * @param eigenvalues The eigenvalues to be set.
 */
void EigenDecomposition::SetGlobalEigenvalues( double (&eigenvalues)[DTI(CC::DIM())][MF::ANOE()] ) const {
   for( unsigned int d = 0; d < DTI(CC::DIM()); ++d ) {
      for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
         global_eigenvalues_[d][e] = eigenvalues[d][e];
      }
   }
}

/**
 * @brief Gives the stored global Lax-Friedrichs eigenvalues.
 * @return The buffer holding the global eigenvalues.
 */
auto EigenDecomposition::GetGlobalEigenvalues() const -> double const (&)[DTI(CC::DIM())][MF::ANOE()] {
   return global_eigenvalues_;
}
