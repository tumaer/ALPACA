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
#include "two_phase_scale_separator.h"

#include "enums/interface_tag_definition.h"
#include "utilities/mathematical_functions.h"
#include "user_specifications/two_phase_constants.h"

/**
 * @brief Default constructor of the TwoPhaseScaleSeparator class.
 * @param material_manager Instance of a MaterialManager, which already has been initialized according to the user input.
 * @param halo_manager Instance of a HaloManager which provides MPI-related methods.
 */
TwoPhaseScaleSeparator::TwoPhaseScaleSeparator( MaterialManager const& material_manager, HaloManager & halo_manager ) :
   ScaleSeparator( material_manager, halo_manager )
{
   // Empty Constructor, besides call of base class constructor.
}

namespace {

/**
 * Loop start indices: set to 0 in case dimension is not considered; otherwise set to 1
 */
unsigned int i_lower = 1;
unsigned int j_lower = CC::DIM() != Dimension::One   ? 1 : 0;
unsigned int k_lower = CC::DIM() == Dimension::Three ? 1 : 0;

/**
 * Loop end indices: set to 1 in case dimensions is not considered, otherwise to the end of the block minus one
 */
unsigned int i_upper = CC::TCX() - 1;
unsigned int j_upper = CC::DIM() != Dimension::One   ? CC::TCY()-1 : 1;
unsigned int k_upper = CC::DIM() == Dimension::Three ? CC::TCZ()-1 : 1;

/**
 * Offset for current cell: 0 if dimension is not considered, 1 otherwise
 */
unsigned int i_offset = 1;
unsigned int j_offset = CC::DIM() != Dimension::One   ? 1 : 0;
unsigned int k_offset = CC::DIM() == Dimension::Three ? 1 : 0;

/**
 * Stencil size for level-set recomputation
 */
int stencil_width = CC::HS() - 1;
int stencil_width_i = stencil_width;
int stencil_width_j = CC::DIM() != Dimension::One ?   stencil_width : 0;
int stencil_width_k = CC::DIM() == Dimension::Three ? stencil_width : 0;

/**
 * Box size
 */
constexpr unsigned int BSX = 3;
constexpr unsigned int BSY = CC::DIM() != Dimension::One   ? 3 : 1;
constexpr unsigned int BSZ = CC::DIM() == Dimension::Three ? 3 : 1;

// According to \cite Han2015
#if DIMENSION == 1
constexpr double stimulus_response_constant = 0.51; // Must be larger than 1 / 2
#elif DIMENSION == 2
constexpr double stimulus_response_constant = 0.75; // Must be larger than sqrt(2) / 2
#else
constexpr double stimulus_response_constant = 0.9; // Must be larger than sqrt(3) / 2
#endif

/**
 * @brief Calculate the resolved levelset value for cell (i, j, k) with a forward tracing method starting from cut cell (i+l, j+m, k+n). This procedure is
 * described in \cite Fu2016b.
 * @param levelset The level-set field for which the reinitialized value has to be calculated.
 * @param i The x-index of the cell where the resolved value should be calculated.
 * @param j The y-index of the cell where the resolved value should be calculated.
 * @param k The z-index of the cell where the resolved value should be calculated.
 * @param l The relative index shift in x-direction where the cut cell is located.
 * @param m The relative index shift in y-direction where the cut cell is located.
 * @param n The relative index shift in z-direction where the cut cell is located.
 * @return The reinitialized value.
 */
double GetResolvedValue(double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()], int const i, int const j, int const k, int const l, int const m, int const n){

   std::array<double, 3> const normal = GetNormal(levelset, i + l, j + m, k + n);

   double const r_x =                                 -l + std::min( levelset[i+l][j+m][k+n], std::sqrt(0.5 * (double)DTI(CC::DIM()))) * normal[0];
   double const r_y = CC::DIM() != Dimension::One   ? -m + std::min( levelset[i+l][j+m][k+n], std::sqrt(0.5 * (double)DTI(CC::DIM()))) * normal[1] : 0.0;
   double const r_z = CC::DIM() == Dimension::Three ? -n + std::min( levelset[i+l][j+m][k+n], std::sqrt(0.5 * (double)DTI(CC::DIM()))) * normal[2] : 0.0;

   constexpr double one_third = 1.0 / 3.0;
   constexpr double thirtyfive_thirds = 35.0 / 3.0;

   //d, r and h as shown in figure 3 in \cite Fu2016b.
   double const r = std::sqrt(r_x*r_x + r_y*r_y + r_z*r_z);
   double cosinus_of_angle_between_r_and_d = 0.0;
   cosinus_of_angle_between_r_and_d += r_x * normal[0];
   if constexpr(CC::DIM() != Dimension::One)   cosinus_of_angle_between_r_and_d += r_y * normal[1];
   if constexpr(CC::DIM() == Dimension::Three) cosinus_of_angle_between_r_and_d += r_z * normal[2];
   cosinus_of_angle_between_r_and_d = std::abs(cosinus_of_angle_between_r_and_d) / r;

   double const h = r * std::sqrt(std::abs(1.0 - cosinus_of_angle_between_r_and_d * cosinus_of_angle_between_r_and_d));
   double const d = r * cosinus_of_angle_between_r_and_d;

   //calculate reinitialized value as described in equations 15 and 16 in \cite Fu2016b. Limiters to avoid unrealistic values - resolved value is usually still a cut-cell or close to one.
   double const lambda = std::max( std::min(h * one_third, 1.0), 0.0); //equation 16 in \cite Fu2016b.
   double const factor_squared = (1.0 - lambda) * (1.0 - lambda);
   double const beta = factor_squared * factor_squared * factor_squared * (1.0 + 6.0 * lambda + thirtyfive_thirds * lambda * lambda); //equation 16 in \cite Fu2016b.
   double const resolved_value = beta * d + (1.0 - beta) * r; //equation 15 in \cite Fu2016b.

   return resolved_value;
}

/**
 * @brief Performs a scale separation procedure according to \cite Luo2016.
 * @param node The node for which scale separation is done.
 */
void ScaleSeparationProcedure(Node& node){

   double (&levelset_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetReinitializedBuffer(InterfaceDescription::Levelset);

   std::int8_t (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();

   std::int8_t interface_tags_positive_shift[CC::TCX()][CC::TCY()][CC::TCZ()];
   std::int8_t interface_tags_negative_shift[CC::TCX()][CC::TCY()][CC::TCZ()];

   double levelset_positive_shift[CC::TCX()][CC::TCY()][CC::TCZ()];
   double levelset_negative_shift[CC::TCX()][CC::TCY()][CC::TCZ()];

   //compute shifted level-set fields
   for(unsigned int i = 0; i < CC::TCX(); ++i) {
      for(unsigned int j = 0; j < CC::TCY(); ++j) {
         for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            levelset_positive_shift[i][j][k] = levelset_reinitialized[i][j][k] + stimulus_response_constant;
            levelset_negative_shift[i][j][k] = levelset_reinitialized[i][j][k] - stimulus_response_constant;
            interface_tags_positive_shift[i][j][k] = Signum(levelset_positive_shift[i][j][k]) * ITTI(IT::BulkPhase);
            interface_tags_negative_shift[i][j][k] = Signum(levelset_negative_shift[i][j][k]) * ITTI(IT::BulkPhase);
         } //k
      } //j
   } //i

   /**
    * Calculate cut-cells for positive and negative shifted levelset field. This can be done in
    * 1. inner cells
    * 2. halo cells except for the last layer of halo cells
    */
   for(unsigned int i = i_lower; i < i_upper; ++i) {
      for(unsigned int j = j_lower; j < j_upper; ++j) {
         for(unsigned int k = k_lower; k < k_upper; ++k) {
            if(IsCutCell<GeometryCalculationSettings::CutCellCriteria>(levelset_positive_shift,i,j,k)) interface_tags_positive_shift[i][j][k] = ITTI(IT::OldCutCell);
            if(IsCutCell<GeometryCalculationSettings::CutCellCriteria>(levelset_negative_shift,i,j,k)) interface_tags_negative_shift[i][j][k] = ITTI(IT::OldCutCell);
         }
      }
   }

   /**
    * Determine non-resolved cells.
    */
   for(unsigned int i = CC::FICX() - i_offset; i < CC::LICX() + i_offset + 1; ++i) {
      for(unsigned int j = CC::FICY() - j_offset; j < CC::LICY() + j_offset + 1; ++j) {
         for(unsigned int k = CC::FICZ() - k_offset; k < CC::LICZ() + k_offset + 1; ++k) {

            std::int8_t flag_sign_change = 0;

            //check for cells with level-set value smaller than scale separation cut-off and neighbor with different sign
            if(std::abs(levelset_reinitialized[i][j][k]) < stimulus_response_constant){
               for(unsigned int r = 0; r < BSX; r++) {
                  for(unsigned int s = 0; s < BSY; s++) {
                     for(unsigned int t = 0; t < BSZ; t++) {
                        if(levelset_reinitialized[i][j][k] * levelset_reinitialized[i + r - i_offset][j + s - j_offset][k + t - k_offset] < 0.0){
                           flag_sign_change += 1;
                        }
                     } //k neighbourhood
                  } //j neighbourhood
               } //i neighbourhood
            }

            if(std::abs(interface_tags[i][j][k]) <= ITTI(IT::NewCutCell) || flag_sign_change != 0) { //Non-resolved cells can either be cut-cells or non-cut-cells with a level-set value smaller than the threshold defined in CC
               std::int8_t flag_positive = 0;
               std::int8_t flag_negative = 0;

               for(unsigned int r = 0; r < BSX; r++) {
                  for(unsigned int s = 0; s < BSY; s++) {
                     for(unsigned int t = 0; t < BSZ; t++) {
                        if(interface_tags_positive_shift[i + r - i_offset][j + s - j_offset][k + t - k_offset] == ITTI(IT::OldCutCell)) {
                           flag_positive += 1;
                        }
                        if(interface_tags_negative_shift[i + r - i_offset][j + s - j_offset][ k + t -k_offset] == ITTI(IT::OldCutCell)) {
                           flag_negative += 1;
                        }
                     } //k neighbourhood
                  } //j neighbourhood
               } //i neighbourhood

               if(flag_negative == 0 && flag_positive != 0) {
                  interface_tags[i][j][k] = -ITTI(IT::ScaleSeparatedCell);
               }
               if(flag_positive == 0 && flag_negative != 0) {
                  interface_tags[i][j][k] = +ITTI(IT::ScaleSeparatedCell);
               }

            }
         } //k
      } //j
   } //i

   /**
    * Calculate scale separated levelset field.
    */
   int flag = 0;
   for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
         for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
            if(std::abs(interface_tags[i][j][k]) == ITTI(IT::ScaleSeparatedCell)) {

               flag = 0;

               levelset_reinitialized[i][j][k] = CC::LSCOF();
               for(int l = -stencil_width_i; l < stencil_width_i+1; l++) {
                  for(int m = -stencil_width_j; m < stencil_width_j+1; m++) {
                     for(int n = -stencil_width_k; n < stencil_width_k+1; n++) {

                        //compute re-constructed level-set. Cells with missing negative-shifted interface as neighbor are reconstructed from the closest negatively shifted interface, and analogously for the positive shift.
                        if(interface_tags[i][j][k] == +ITTI(IT::ScaleSeparatedCell) && interface_tags_positive_shift[i + l][j + m][k + n] == ITTI(IT::OldCutCell)){
                           levelset_reinitialized[i][j][k] = std::min(levelset_reinitialized[i][j][k], GetResolvedValue(levelset_positive_shift, i, j, k, l, m, n));
                           flag = 1;
                        } else if(interface_tags[i][j][k] == -ITTI(IT::ScaleSeparatedCell) && interface_tags_negative_shift[i + l][j + m][k + n] == ITTI(IT::OldCutCell)){
                           levelset_reinitialized[i][j][k] = std::min(levelset_reinitialized[i][j][k], GetResolvedValue(levelset_negative_shift, i, j, k, l, m, n));
                           flag = 1;
                        }
                     } //k neighbourhood
                  } //j neighbourhood
               } //i neighbourhood

               if(flag == 0) {
                  if(interface_tags[i][j][k] == -ITTI(IT::ScaleSeparatedCell)) {
                     levelset_reinitialized[i][j][k] = -CC::LSCOF() + stimulus_response_constant;
                  } else {
                     levelset_reinitialized[i][j][k] = CC::LSCOF() - stimulus_response_constant;
                  }
               } else {
                  levelset_reinitialized[i][j][k] = double( Signum(interface_tags[i][j][k]) ) * (levelset_reinitialized[i][j][k] - stimulus_response_constant);
               }
            }
         } //k
      } //j
   } //i
}
}

/**
 * @brief Performs a scale separation procedure for a given vector of nodes.
 * @param nodes The nodes for which scale separation is done.
 */
void TwoPhaseScaleSeparator::SeparateScalesImplementation(std::vector<std::reference_wrapper<Node>> const& nodes) const {
   for(Node& node : nodes) {
      ScaleSeparationProcedure(node);
   }
   halo_manager_.InterfaceHaloUpdateOnLmax( InterfaceBlockBufferType::LevelsetReinitialized );
}
