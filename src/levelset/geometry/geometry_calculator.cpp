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
#include "geometry_calculator.h"
#include "stencils/stencil_utilities.h"
#include "enums/interface_tag_definition.h"
#include "utilities/mathematical_functions.h"
#include "user_specifications/two_phase_constants.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

/**
 * @brief Calculates level-set values at the corners of cell (i,j,k). Calculations are based on linear interpolations.
 * @param cell_corner_levelset The levelset value at cell corners. Indirect return parameter.
 * @param levelset The level-set field.
 * @param i Index i of the cell of interest.
 * @param j Index j of the cell of interest.
 * @param k Index k of the cell of interest.
 */
void GetLevelsetAtCellCorners(double (&cell_corner_levelset)[cell_box_size_x][cell_box_size_y][cell_box_size_z],
   double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
   unsigned int const i, unsigned int const j, unsigned int const k) {

   //indices of bottom left are chosen for the basic box that gets shifted later to determine all boxes in an iterative way
   unsigned int const index_i = i-1;
   unsigned int const index_j = CC::DIM() != Dimension::One   ? j-1 : 0;
   unsigned int const index_k = CC::DIM() == Dimension::Three ? k-1 : 0;
   unsigned int const index_ip1 = i;
   unsigned int const index_jp1 = CC::DIM() != Dimension::One   ? j : 0;
   unsigned int const index_kp1 = CC::DIM() == Dimension::Three ? k : 0;

   //number of boxes in each direction
   unsigned int const rmax = 2;
   unsigned int const smax = CC::DIM() != Dimension::One   ? 2 : 1;
   unsigned int const tmax = CC::DIM() == Dimension::Three ? 2 : 1;

   // degenerates automatically in 2D and 1D, however always uses 8 corner values which are mirrored then into the not considered dimension
   std::vector<double> levelset_at_corners;

   //get levelset values at the corners of the cell
   for(unsigned int r = 0; r < rmax; ++r) {
      for(unsigned int s = 0; s < smax; ++s) {
         for(unsigned int t = 0; t < tmax; ++t) {
            unsigned int const index_i_r   = index_i   + r;
            unsigned int const index_ip1_r = index_ip1 + r;
            unsigned int const index_j_s   = index_j   + s;
            unsigned int const index_jp1_s = index_jp1 + s;
            unsigned int const index_k_t   = index_k   + t;
            unsigned int const index_kp1_t = index_kp1 + t;

            levelset_at_corners = {levelset[index_i_r][index_j_s]  [index_k_t]   , levelset[index_ip1_r][index_j_s]  [index_k_t],
                                   levelset[index_i_r][index_jp1_s][index_k_t]   , levelset[index_ip1_r][index_jp1_s][index_k_t],
                                   levelset[index_i_r][index_j_s]  [index_kp1_t] , levelset[index_ip1_r][index_j_s]  [index_kp1_t],
                                   levelset[index_i_r][index_jp1_s][index_kp1_t] , levelset[index_ip1_r][index_jp1_s][index_kp1_t]};
            std::sort( levelset_at_corners.begin(), levelset_at_corners.end() );
            cell_corner_levelset[r][s][t] = 0.125 * (std::accumulate(levelset_at_corners.begin(), levelset_at_corners.end(), 0.0));
         }
      }
   }
}

/**
* @brief Calculates the level-set values at the corners of all subcells which appear when the parent cell gets refined once, i.e.
*        at the center of this cell, its corners, and the the center of each patch and edge. Calculations are based on linear interpolations.
* @param corner_levelset The array in which to store the level-set values at the corners of the subcells. Indirect return parameter.
* @param levelset The level-set field.
* @param i Index i of the cell of interest.
* @param j Index j of the cell of interest.
* @param k Index k of the cell of interest.
*/
void GetLevelsetAtSubcellCorners(double (&subcell_corner_levelset)[subcell_box_size_x][subcell_box_size_y][subcell_box_size_z],
   double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
   unsigned int const i, unsigned int const j, unsigned int const k) {


   //indices of bottom left are chosen for the basic box that gets shifted later to determine all boxes in an iterative way
   unsigned int const index_i = i-1;
   unsigned int const index_j = CC::DIM() != Dimension::One   ? j-1 : 0;
   unsigned int const index_k = CC::DIM() == Dimension::Three ? k-1 : 0;
   unsigned int const index_ip1 = i;
   unsigned int const index_jp1 = CC::DIM() != Dimension::One   ? j : 0;
   unsigned int const index_kp1 = CC::DIM() == Dimension::Three ? k : 0;

   //number of boxes in each direction
   unsigned int const rmax = 2;
   unsigned int const smax = CC::DIM() != Dimension::One   ? 2 : 1;
   unsigned int const tmax = CC::DIM() == Dimension::Three ? 2 : 1;

   //indices that are used while iterating through all 8 boxes
   unsigned int index_i_r;
   unsigned int index_ip1_r;
   unsigned int index_j_s;
   unsigned int index_jp1_s;
   unsigned int index_k_t;
   unsigned int index_kp1_t;

   //levelset at the center of the box
   subcell_corner_levelset[rmax-1][smax-1][tmax-1] = levelset[index_ip1][index_jp1][index_kp1];

   // degenerates automatically in 2D and 1D, however always uses 8 corner values which are mirrored then into the not considered dimension
   std::vector<double> levelset_at_corners;

   //get levelset values at the corners of the cell
   for(unsigned int r = 0; r < rmax; ++r) {
      for(unsigned int s = 0; s < smax; ++s) {
         for(unsigned int t = 0; t < tmax; ++t) {
            index_i_r   = index_i   + r;
            index_ip1_r = index_ip1 + r;
            index_j_s   = index_j   + s;
            index_jp1_s = index_jp1 + s;
            index_k_t   = index_k   + t;
            index_kp1_t = index_kp1 + t;

            levelset_at_corners = {levelset[index_i_r][index_j_s]  [index_k_t]   , levelset[index_ip1_r][index_j_s]  [index_k_t],
                                   levelset[index_i_r][index_jp1_s][index_k_t]   , levelset[index_ip1_r][index_jp1_s][index_k_t],
                                   levelset[index_i_r][index_j_s]  [index_kp1_t] , levelset[index_ip1_r][index_j_s]  [index_kp1_t],
                                   levelset[index_i_r][index_jp1_s][index_kp1_t] , levelset[index_ip1_r][index_jp1_s][index_kp1_t]};
            std::sort( levelset_at_corners.begin(), levelset_at_corners.end() );
            subcell_corner_levelset[2*r][2*s][2*t] = 0.125 * (std::accumulate(levelset_at_corners.begin(), levelset_at_corners.end(), 0.0));
            levelset_at_corners.clear();
         }
      }
   }

   std::vector<double> levelset_edge_vector = {0.0, 0.0, 0.0, 0.0};

   //edges for 2D and 3D case
   if constexpr(CC::DIM() != Dimension::One) {

      for(unsigned int s = 0; s < smax; ++s) {
         for(unsigned int t = 0; t < tmax; ++t) {
            index_j_s   = index_j   + s;
            index_jp1_s = index_jp1 + s;
            index_k_t   = index_k   + t;
            index_kp1_t = index_kp1 + t;

            levelset_edge_vector = {levelset[i][index_j_s]  [index_k_t], levelset[i][index_j_s] [index_kp1_t],
                                    levelset[i][index_jp1_s][index_k_t], levelset[i][index_jp1_s][index_kp1_t]};
            std::sort( levelset_edge_vector.begin(), levelset_edge_vector.end() );
            subcell_corner_levelset[1][2*s][2*t] = 0.25 * (std::accumulate(levelset_edge_vector.begin(), levelset_edge_vector.end(), 0.0));
            levelset_edge_vector.clear();
         }
      }

      for(unsigned int r = 0; r < rmax; ++r) {
         for(unsigned int t = 0; t < tmax; ++t) {
            index_i_r   = index_i   + r;
            index_ip1_r = index_ip1 + r;
            index_k_t   = index_k   + t;
            index_kp1_t = index_kp1 + t;

            levelset_edge_vector = {levelset[index_i_r][j][index_k_t],   levelset[index_ip1_r][j][index_k_t],
                                    levelset[index_i_r][j][index_kp1_t], levelset[index_ip1_r][j][index_kp1_t]};
            std::sort( levelset_edge_vector.begin(), levelset_edge_vector.end() );
            subcell_corner_levelset[2*r][1][2*t] = 0.25 * (std::accumulate(levelset_edge_vector.begin(), levelset_edge_vector.end(), 0.0));
            levelset_edge_vector.clear();
         }
      }
   }

   //edges and patches for 3D cases
   if constexpr(CC::DIM() == Dimension::Three) {

      // Calculation of edges
      for(unsigned int r = 0; r < rmax; ++r) {
         for(unsigned int s = 0; s < smax; ++s) {
            index_i_r   = index_i   + r;
            index_ip1_r = index_ip1 + r;
            index_j_s   = index_j   + s;
            index_jp1_s = index_jp1 + s;

            levelset_edge_vector = {levelset[index_i_r]  [index_j_s]  [k], levelset[index_ip1_r][index_j_s]  [k],
                                    levelset[index_i_r]  [index_jp1_s][k], levelset[index_ip1_r][index_jp1_s][k]};
            std::sort( levelset_edge_vector.begin(), levelset_edge_vector.end() );
            subcell_corner_levelset[2*r][2*s][1] = 0.25 * (std::accumulate(levelset_edge_vector.begin(), levelset_edge_vector.end(), 0.0));
            levelset_edge_vector.clear();
         }
      }

      // Calculation of patches
      std::vector<double> levelset_patch_vector;

      for(unsigned int r = 0; r < rmax; ++r) {

         levelset_patch_vector = {levelset[i-1+r][j][k], levelset[i+r]  [j][k]};
         std::sort( levelset_patch_vector.begin(), levelset_patch_vector.end() );
         subcell_corner_levelset[2*r][1][1] = 0.5 * (std::accumulate(levelset_patch_vector.begin(), levelset_patch_vector.end(), 0.0));
         levelset_patch_vector.clear();

         levelset_patch_vector = {levelset[i][j-1+r][k], levelset[i][j+r][k]};
         std::sort( levelset_patch_vector.begin(), levelset_patch_vector.end() );
         subcell_corner_levelset[1][2*r][1] = 0.5 * (std::accumulate(levelset_patch_vector.begin(), levelset_patch_vector.end(), 0.0));
         levelset_patch_vector.clear();

         levelset_patch_vector = {levelset[i][j][k-1+r], levelset[i][j][k+r]};
         std::sort( levelset_patch_vector.begin(), levelset_patch_vector.end() );
         subcell_corner_levelset[1][1][2*r] = 0.5 * (std::accumulate(levelset_patch_vector.begin(), levelset_patch_vector.end(), 0.0));
         levelset_patch_vector.clear();
      }
   }
}

/**
 * @brief Tests whether cell (i,j,k) is a cut cell. Implementation of a sign-change based criterion.
 * @param levelset The level-set field.
 * @param i Index i of the cell of interest.
 * @param j Index j of the cell of interest.
 * @param k Index k of the cell of interest.
 * @return Bool whether the cell (i,j,k) is a cut cell (true) or not (false).
 */
template<>
bool IsCutCell<CutCellCriteria::SignChangeBased>(double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const k) {

#ifndef PERFORMANCE
   //check if one runs this function in a correct cell
   bool cell_is_outmost = (i==0) || (i >= CC::TCX()-1) || (j >= CC::TCY()) || (k >= CC::TCZ());
   if constexpr(CC::DIM() != Dimension::One) {
      cell_is_outmost = cell_is_outmost || (j==0) || (j == CC::TCY()-1);
   }
   if constexpr(CC::DIM() == Dimension::Three) {
      cell_is_outmost = cell_is_outmost || (k==0) || (k == CC::TCZ()-1);
   }

   if(cell_is_outmost){
      throw std::invalid_argument("You cannot run IsCutCell function for this cell.");
   }
#endif

   double levelset_cube[subcell_box_size_x][subcell_box_size_y][subcell_box_size_z];
   for(unsigned int r = 0; r < subcell_box_size_x; ++r) {
      for(unsigned int s = 0; s < subcell_box_size_y; ++s) {
         for(unsigned int t = 0; t < subcell_box_size_z; ++t) {
            levelset_cube[r][s][t] = 0.0;
         }
      }
   }

   GetLevelsetAtSubcellCorners(levelset_cube, levelset, i, j, k);

   for(unsigned int r = 0; r < subcell_box_size_x; ++r) {
      for(unsigned int s = 0; s < subcell_box_size_y; ++s) {
         for(unsigned int t = 0; t < subcell_box_size_z; ++t) {
            /* With this if statement we check whether all elements in levelset_cube have the same sign. If this is not the case,
               the cell is a cut cell and we return true */
            if(Signum(levelset_cube[r][s][t]) != Signum(levelset_cube[0][0][0])) {return true;}
         }
      }
   }

   return false;

}

/**
 * @brief Tests whether cell (i,j,k) is a cut cell. Implementation of a value-based criterion.
 * @param levelset The level-set field.
 * @param i Index i of the cell of interest.
 * @param j Index j of the cell of interest.
 * @param k Index k of the cell of interest.
 * @return Bool whether the cell (i,j,k) is a cut cell (true) or not (false).
 */
template<>
bool IsCutCell<CutCellCriteria::ValueBased>(double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const k) {

#if DIMENSION == 1
   constexpr double largest_cut_cell_levelset_value = 0.51; // Must be larger than 1 / 2
#elif DIMENSION == 2
   constexpr double largest_cut_cell_levelset_value = 0.75; // Must be larger than sqrt(2) / 2
#else
   constexpr double largest_cut_cell_levelset_value = 1.0; // Must be larger than sqrt(3) / 2
#endif

   return std::abs(levelset[i][j][k]) < largest_cut_cell_levelset_value;

}

/**
 * @brief Calculates the cell-normal vector, which is the gradient of the level-set field. It is computed with a second-order
 *        central scheme.
 * @param levelset The level-set field.
 * @param i Index i of the cell of interest.
 * @param j Index j of the cell of interest.
 * @param k Index k of the cell of interest.
 * @param material_sign Indicates for which phase the inward pointing normal is calculated ( > 0 -> positive phase, < 0 -> negative phase). Default: positive material.
 * @return Normal unit vector to the interface in cell (i,j,k).
 */
std::array<double, 3> GetNormal(double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
   unsigned int const i, unsigned int const j, unsigned int const k, std::int8_t const material_sign) {

   std::array<double, 3> normal = {0.0, 0.0, 0.0};
   if constexpr( GeometryCalculationSettings::GeometryStencilType == GeometryStencilType::Reconstruction ) {
      using ReconstructionStencil = ReconstructionStencilSetup::Concretize<geometry_reconstruction_stencil>::type;
      normal = SU::GradientVector<ReconstructionStencil>(levelset, i, j, k, 1.0);
   } else if constexpr( GeometryCalculationSettings::GeometryStencilType == GeometryStencilType::Derivative ) {
      using DerivativeStencil = DerivativeStencilSetup::Concretize<normal_calculation_derivative_stencil>::type;
      normal = SU::GradientVector<DerivativeStencil>(levelset, i, j, k, 1.0);
   }
   double const one_length = 1.0 / std::sqrt(DimensionAwareConsistencyManagedSum(normal[0] * normal[0] , normal[1] * normal[1] , normal[2] * normal[2]));
   double const multiplicator = material_sign * one_length;
   std::transform(normal.begin(), normal.end(), normal.begin(), [&multiplicator](double const& normal_component) {return normal_component * multiplicator;});
   return normal;

}

/**
 * @brief Computes the interface curvature and stores it in the curvature buffer.
 * @param node The node for which the interface curvature is computed.
 * @param curvature The curvature buffer. Indirect return parameter.
 */
void ComputeInterfaceCurvature(Node const& node, double (&curvature)[CC::TCX()][CC::TCY()][CC::TCZ()]) {

   /**
    * @brief A constexpr to have the number of dimensions as a double.
    */
   constexpr double dimension_double = double(DTI(CC::DIM()));

   using ReconstructionStencil = ReconstructionStencilSetup::Concretize<geometry_reconstruction_stencil>::type;
   using Stencil = DerivativeStencilSetup::Concretize<curvature_calculation_derivative_stencil>::type;

   static_assert(2 * (ReconstructionStencil::DownstreamStencilSize() + 1) <= CC::HS(), "Halo size not enough to calculate curvature with the chosen stencil. Increase the halo size in compile_time_constants.h or choose another stencil to calculate curvature!");

   /**
    * Offsets in order to also calculate first derivatives in halo cells. This is necessary for the calculation of
    * second derivatives.
    */
   constexpr unsigned int offset_x = ReconstructionStencil::DownstreamStencilSize() + 1;
   constexpr unsigned int offset_y = CC::DIM() != Dimension::One ? ReconstructionStencil::DownstreamStencilSize() + 1: 0;
   constexpr unsigned int offset_z = CC::DIM() == Dimension::Three ? ReconstructionStencil::DownstreamStencilSize() + 1: 0;

   double const (&levelset_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetReinitializedBuffer(InterfaceDescription::Levelset);
   std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
   double const cell_size = node.GetCellSize();
   double const one_cell_size = 1.0 / cell_size;

   double levelset_x[CC::TCX()][CC::TCY()][CC::TCZ()];
   double levelset_y[CC::TCX()][CC::TCY()][CC::TCZ()];
   double levelset_z[CC::TCX()][CC::TCY()][CC::TCZ()];
   for(unsigned int i = 0; i < CC::TCX(); ++i){
      for(unsigned int j = 0; j < CC::TCY(); ++j){
         for(unsigned int k = 0; k < CC::TCZ(); ++k){
            levelset_x[i][j][k] = 0.0;
            levelset_y[i][j][k] = 0.0;
            levelset_z[i][j][k] = 0.0;
         }
      }
   }

   std::array<std::array<double, 3>, 3> second_derivatives;
   std::array<double, 3> gradient;

   for(unsigned int i = 0 + offset_x; i < CC::TCX() - offset_x; ++i){
      for(unsigned int j = 0 + offset_y; j < CC::TCY() - offset_y; ++j){
         for(unsigned int k = 0 + offset_z; k < CC::TCZ() - offset_z; ++k){
            if(std::abs(interface_tags[i][j][k]) <= ITTI(IT::ReinitializationBand)){
               if constexpr( GeometryCalculationSettings::GeometryStencilType == GeometryStencilType::Reconstruction ) {
                  gradient = SU::GradientVector<ReconstructionStencil>(levelset_reinitialized, i, j, k, 1.0);
               } else if constexpr( GeometryCalculationSettings::GeometryStencilType == GeometryStencilType::Derivative ) {
                  gradient = SU::GradientVector<Stencil>(levelset_reinitialized, i, j, k, 1.0);
               }
               levelset_x[i][j][k] = gradient[0];
               levelset_y[i][j][k] = gradient[1];
               levelset_z[i][j][k] = gradient[2];
            }
         }
      }
   }

   for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
         for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
            if(std::abs(interface_tags[i][j][k]) <= ITTI(IT::NewCutCell)){
               if constexpr( GeometryCalculationSettings::GeometryStencilType == GeometryStencilType::Reconstruction ) {
                  second_derivatives = SU::JacobianMatrix<ReconstructionStencil>(levelset_x, levelset_y, levelset_z, i, j, k, cell_size);
               } else if constexpr( GeometryCalculationSettings::GeometryStencilType == GeometryStencilType::Derivative ) {
                  second_derivatives = SU::JacobianMatrix<Stencil>(levelset_x, levelset_y, levelset_z, i, j, k, cell_size);
               }

               double const first_summand = ConsistencyManagedSum(
                  levelset_x[i][j][k] * levelset_x[i][j][k] * ( second_derivatives[1][1] + second_derivatives[2][2] )
                  , levelset_y[i][j][k] * levelset_y[i][j][k] * ( second_derivatives[0][0] + second_derivatives[2][2] )
                  , levelset_z[i][j][k] * levelset_z[i][j][k] * ( second_derivatives[0][0] + second_derivatives[1][1] ) );

               double const second_summand = - 2.0 * ConsistencyManagedSum(
                  ( levelset_x[i][j][k] * levelset_y[i][j][k] ) * second_derivatives[0][1]
                  , ( levelset_x[i][j][k] * levelset_z[i][j][k] ) * second_derivatives[0][2]
                  , ( levelset_y[i][j][k] * levelset_z[i][j][k] ) * second_derivatives[1][2] );

               double const length_of_normal = std::sqrt( DimensionAwareConsistencyManagedSum( levelset_x[i][j][k] * levelset_x[i][j][k], levelset_y[i][j][k] * levelset_y[i][j][k], levelset_z[i][j][k] * levelset_z[i][j][k] ) );
               double const denominator = length_of_normal * length_of_normal * length_of_normal;

               double single_cell_curvature = (first_summand + second_summand);

               if( CC::Axisymmetric() ) {
                  double const radius = node.GetBlockCoordinateX() + (double(i) - double(CC::FICX()) + 0.5) * cell_size;
                  single_cell_curvature += levelset_x[i][j][k] * (levelset_x[i][j][k] * levelset_x[i][j][k] + levelset_y[i][j][k] * levelset_y[i][j][k])/radius;
               }

               single_cell_curvature /= denominator;

               double const corrected_curvature = (dimension_double - 1.0) * single_cell_curvature / (dimension_double - 1.0 - levelset_reinitialized[i][j][k]*cell_size*single_cell_curvature);

               /**
                * Limit the curvature as described in \cite kang2000boundary.
                */
               double const limited_corrected_curvature = Sign(corrected_curvature) * std::min( std::abs(corrected_curvature), one_cell_size );

               curvature[i][j][k] = limited_corrected_curvature;
            }
         }
      }
   }

}
