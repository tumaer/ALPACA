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
#include "geometry_calculator_marching_cubes.h"
#include "utilities/mathematical_functions.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <bitset>
#include <functional>

namespace {

/**
 * @brief      Function to calculate the part of a square which is covered by the positive side of the level-set field.
 *             The letters M respectively P behind the "Aperture" in the function name indicate which sign is assumed for
 *             the respective level-set value of the corner.
 *
 * @param[in]  face  The level-set values of the corners of the square.
 *
 * @return     The aperture.
 */
constexpr double ApertureMMMP(std::array<double, 4> const face) { // 8
   return 0.5 * ((face[3] / (face[3] - face[0]) ) * (face[3] / (face[3] - face[2]) ));
}
constexpr double ApertureMMMM(std::array<double, 4> const face) { // 0
#ifndef PERFORMANCE
   (void) face; //Avoid compiler warning
#endif
   return 0.0;
}
constexpr double AperturePMMM(std::array<double, 4> const face) { // 1
   return 0.5 * ((face[0] / (face[0] - face[1]) ) * (face[0] / (face[0] - face[3]) ));
}
constexpr double ApertureMPMM(std::array<double, 4> const face) { // 2
   return 0.5 * ((face[1] / (face[1] - face[0]) ) * (face[1] / (face[1] - face[2]) ));
}
constexpr double AperturePPMM(std::array<double, 4> const face) { // 3
   double const tmp0 = face[0] / (face[0] - face[3]);
   double const tmp1 = face[1] / (face[1] - face[2]);
   return 0.5 * (tmp0 + tmp1);
}
constexpr double ApertureMMPM(std::array<double, 4> const face) { // 4
   return 0.5 * ((face[2] / (face[2] - face[1]) ) * (face[2] / (face[2] - face[3]) ));
}
constexpr double AperturePMPM(std::array<double, 4> const face) { // 5
   double const center = ConsistencyManagedSum(face);
   return (center > 0.0) * (1.0 - (ApertureMMMP(face) + ApertureMPMM(face))) + (center <= 0.0) * (AperturePMMM(face) + ApertureMMPM(face));
}
constexpr double ApertureMPPM(std::array<double, 4> const face) { // 6
   double const tmp1 = face[1] / (face[1] - face[0]);
   double const tmp2 = face[2] / (face[2] - face[3]);
   return 0.5 * (tmp1 + tmp2);
}
constexpr double AperturePPPM(std::array<double, 4> const face) { // 7
   return 1.0 - ApertureMMMP(face);
}
constexpr double AperturePMMP(std::array<double, 4> const face) { // 9
   double const tmp0 = face[0] / (face[0] - face[1]);
   double const tmp3 = face[3] / (face[3] - face[2]);
   return 0.5 * (tmp0 + tmp3);
}
constexpr double ApertureMPMP(std::array<double, 4> const face) { // 10
   double const center = ConsistencyManagedSum(face);
   return (center > 0.0) * (1.0 - (AperturePMMM(face) + ApertureMMPM(face))) + (center <= 0.0) * (ApertureMPMM(face) + ApertureMMMP(face));
}
constexpr double AperturePPMP(std::array<double, 4> const face) { // 11
   return 1.0 - ApertureMMPM(face);
}
constexpr double ApertureMMPP(std::array<double, 4> const face) { // 12
   double const tmp2 = face[2] / (face[2] - face[1]);
   double const tmp3 = face[3] / (face[3] - face[0]);
   return 0.5 * (tmp2 + tmp3);
}
constexpr double AperturePMPP(std::array<double, 4> const face) { // 13
   return 1.0 - ApertureMPMM(face);
}
constexpr double ApertureMPPP(std::array<double, 4> const face) { // 14
   return 1.0 - AperturePMMM(face);
}
constexpr double AperturePPPP(std::array<double, 4> const face) { // 15
#ifndef PERFORMANCE
   (void) face; //Avoid compiler warning
#endif
   return 1.0;
}

/**
 * Function lookup table for aperture calculation.
 */
std::array<std::function<double (std::array<double, 4> const)>, 16> ApertureFunctionLookup = {ApertureMMMM, AperturePMMM, ApertureMPMM, AperturePPMM, ApertureMMPM, AperturePMPM, ApertureMPPM, AperturePPPM, ApertureMMMP, AperturePMMP, ApertureMPMP, AperturePPMP, ApertureMMPP, AperturePMPP, ApertureMPPP, AperturePPPP};


/**
 * @brief Calculates the aperture of a single cell patch from the level-set values at the corners using triangulation.
 * @param face Level-set values at the corners of the cell patch.
 * @return Cell-face aperture.
 */
double CellFaceAperture(std::array<double, 4> const face){

   std::bitset<4> positive_corners;
   for(unsigned int b = 0; b < 4; ++b) {
      positive_corners.set(b, face[b] > 0.0);
   }
   return ApertureFunctionLookup[positive_corners.to_ulong()](face);
}

/**
 * @brief Calculates the cell face apertures according to \cite Lauer2012 .
 * @param levelset The level-set field.
 * @param i Index i of the cell of interest.
 * @param j Index j of the cell of interest.
 * @param k Index k of the cell of interest.
 * @param material_sign Indicates for which phase the cell-face apertures are calculated ( > 0 -> positive phase, < 0 -> negative phase). Default: positive material.
 * @return Cell-face apertures.
 */
std::array<double, 6> ComputeCellFaceApertureCellBasedImplementation(double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
   unsigned int const i, unsigned int const j, unsigned int const k,
   std::int8_t const material_sign){

   double cell_levelset_cube[cell_box_size_x][cell_box_size_y][cell_box_size_z];
   for(unsigned int r = 0; r < cell_box_size_x; ++r) {
      for(unsigned int s = 0; s < cell_box_size_y; ++s) {
         for(unsigned int t = 0; t < cell_box_size_z; ++t) {
            cell_levelset_cube[r][s][t] = 0.0;
         }
      }
   }

   GetLevelsetAtCellCorners(cell_levelset_cube,levelset,i,j,k);

   unsigned int r = 0;
   unsigned int s = 0;
   unsigned int t = 0;

   unsigned int const r_offset = 1;
   unsigned int const s_offset = CC::DIM() != Dimension::One   ? 1 : 0;
   unsigned int const t_offset = CC::DIM() == Dimension::Three ? 1 : 0;

   std::array<double, 4> faces = {0.0, 0.0, 0.0, 0.0};

   std::array<double, 6> cell_apertures = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

   //face (i-1/2,j,k)
   r = 0;
   s = 0;
   t = 0;
   faces[0] = cell_levelset_cube[r][s]         [t]         ;
   faces[1] = cell_levelset_cube[r][s+s_offset][t]         ;
   faces[2] = cell_levelset_cube[r][s+s_offset][t+t_offset];
   faces[3] = cell_levelset_cube[r][s]         [t+t_offset];
   cell_apertures[0] = CellFaceAperture(faces);

   //face (i+1/2,j,k)
   r = 0 + r_offset;
   s = 0;
   t = 0;
   faces[0] = cell_levelset_cube[r][s]         [t]         ;
   faces[1] = cell_levelset_cube[r][s+s_offset][t]         ;
   faces[2] = cell_levelset_cube[r][s+s_offset][t+t_offset];
   faces[3] = cell_levelset_cube[r][s]         [t+t_offset];
   cell_apertures[1] = CellFaceAperture(faces);

   //only needed in 2D/3D
   if constexpr(CC::DIM() != Dimension::One) {

      //face (i,j-1/2,k)
      r = 0;
      s = 0;
      t = 0;
      faces[0] = cell_levelset_cube[r]         [s][t]         ;
      faces[1] = cell_levelset_cube[r+r_offset][s][t]         ;
      faces[2] = cell_levelset_cube[r+r_offset][s][t+t_offset];
      faces[3] = cell_levelset_cube[r]         [s][t+t_offset];
      cell_apertures[2] = CellFaceAperture(faces);

      //face (i,j+1/2,k)
      r = 0;
      s = 0 + s_offset;
      t = 0;
      faces[0] = cell_levelset_cube[r]         [s][t]         ;
      faces[1] = cell_levelset_cube[r+r_offset][s][t]         ;
      faces[2] = cell_levelset_cube[r+r_offset][s][t+t_offset];
      faces[3] = cell_levelset_cube[r]         [s][t+t_offset];
      cell_apertures[3] = CellFaceAperture(faces);

   } else {
      //pad with 0.0 to length 6
      cell_apertures[2] = 0.0;
      cell_apertures[3] = 0.0;
   }

   //only needed in 3D
   if constexpr(CC::DIM() == Dimension::Three) {

      //face (i,j,k-1/2)
      r = 0;
      s = 0;
      t = 0;
      faces[0] = cell_levelset_cube[r]         [s]         [t];
      faces[1] = cell_levelset_cube[r+r_offset][s]         [t];
      faces[2] = cell_levelset_cube[r+r_offset][s+s_offset][t];
      faces[3] = cell_levelset_cube[r]         [s+s_offset][t];
      cell_apertures[4] = CellFaceAperture(faces);

      //face (i,j,k+1/2)
      r = 0;
      s = 0;
      t = 0 + t_offset;
      faces[0] = cell_levelset_cube[r]         [s]         [t];
      faces[1] = cell_levelset_cube[r+r_offset][s]         [t];
      faces[2] = cell_levelset_cube[r+r_offset][s+s_offset][t];
      faces[3] = cell_levelset_cube[r]         [s+s_offset][t];
      cell_apertures[5] = CellFaceAperture(faces);

   } else {
      //pad with 0.0 to length 6
      cell_apertures[4] = 0.0;
      cell_apertures[5] = 0.0;
   }


   for(double& cell_apertures_component : cell_apertures) {
      cell_apertures_component = std::min(cell_apertures_component,1.0);
      cell_apertures_component = std::max(cell_apertures_component,0.0);
   }

   //consider right phase
   if(material_sign < 0) {
      for(double& cell_apertures_component : cell_apertures) {
         cell_apertures_component = 1.0 - cell_apertures_component;
      }
   }

   return cell_apertures;
}

/**
 * @brief Calculates the cell face apertures according to \cite Lauer2012 .
 * @param levelset The level-set field.
 * @param i Index i of the cell of interest.
 * @param j Index j of the cell of interest.
 * @param k Index k of the cell of interest.
 * @param material_sign Indicates for which phase the cell-face apertures are calculated ( > 0 -> positive phase, < 0 -> negative phase). Default: positive material.
 * @return Cell-face apertures.
 */
std::array<double, 6> ComputeCellFaceApertureSubCellBasedImplementation(double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
   unsigned int const i, unsigned int const j, unsigned int const k,
   std::int8_t const material_sign){

   double levelset_cube[subcell_box_size_x][subcell_box_size_y][subcell_box_size_z];
   for(unsigned int r = 0; r < subcell_box_size_x; ++r) {
      for(unsigned int s = 0; s < subcell_box_size_y; ++s) {
         for(unsigned int t = 0; t < subcell_box_size_z; ++t) {
            levelset_cube[r][s][t] = 0.0;
         }
      }
   }

   //refine the cell by one level and compute level-set values on the corner of these subcells by linear interpolation
   GetLevelsetAtSubcellCorners(levelset_cube, levelset, i, j, k);

   unsigned int r = 0;
   unsigned int s = 0;
   unsigned int t = 0;
   unsigned int const r_max = 2;
   unsigned int const s_max = CC::DIM() != Dimension::One   ? 2 : 1;
   unsigned int const t_max = CC::DIM() == Dimension::Three ? 2 : 1;
   unsigned int const r_offset = 1;
   unsigned int const s_offset = CC::DIM() != Dimension::One   ? 1 : 0;
   unsigned int const t_offset = CC::DIM() == Dimension::Three ? 1 : 0;

   double const multiplier = CC::DIM() == Dimension::Three ? 0.25 : (CC::DIM() != Dimension::One ? 0.5 : 1.0);

   std::vector<double> subcell_apertures;
   std::array<double, 6> cell_apertures = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   std::array<double, 4> faces = {0.0, 0.0, 0.0, 0.0};

   //compute cell-face apertures. They are computed as average from the 4/2/1 subcell-face apertures in 3D/2D/1D.
   //face (i-1/2,j,k)
   r = 0;
   for(s = 0; s < s_max; ++s) {
      for(t = 0; t < t_max; ++t) {
         faces[0] = levelset_cube[r][s]         [t]         ;
         faces[1] = levelset_cube[r][s+s_offset][t]         ;
         faces[2] = levelset_cube[r][s+s_offset][t+t_offset];
         faces[3] = levelset_cube[r][s]         [t+t_offset];
         subcell_apertures.emplace_back(CellFaceAperture(faces));
      }
   }
   //apertures are computed on the subcells and then averaged. Reorder by size to decrease floating point errors.
   std::sort( subcell_apertures.begin(), subcell_apertures.end() );
   cell_apertures[0] = multiplier * (std::accumulate(subcell_apertures.begin(), subcell_apertures.end(), 0.0));
   subcell_apertures.clear();

   //face (i+1/2,j,k)
   r = r_max;
   for(s = 0; s < s_max; ++s) {
      for(t = 0; t < t_max; ++t) {
         faces[0] = levelset_cube[r][s]         [t]         ;
         faces[1] = levelset_cube[r][s+s_offset][t]         ;
         faces[2] = levelset_cube[r][s+s_offset][t+t_offset];
         faces[3] = levelset_cube[r][s]         [t+t_offset];
         subcell_apertures.emplace_back(CellFaceAperture(faces));
      }
   }
   std::sort( subcell_apertures.begin(), subcell_apertures.end() );
   cell_apertures[1] = multiplier * (std::accumulate(subcell_apertures.begin(), subcell_apertures.end(), 0.0));
   subcell_apertures.clear();

   //only needed in 2D/3D
   if constexpr(CC::DIM() != Dimension::One) {
      //face (i,j-1/2,k)
      s = 0;
      for(r = 0; r < r_max; ++r) {
         for(t = 0; t < t_max; ++t) {
            faces[0] = levelset_cube[r]         [s][t]         ;
            faces[1] = levelset_cube[r+r_offset][s][t]         ;
            faces[2] = levelset_cube[r+r_offset][s][t+t_offset];
            faces[3] = levelset_cube[r]         [s][t+t_offset];
            subcell_apertures.emplace_back(CellFaceAperture(faces));
         }
      }
      std::sort( subcell_apertures.begin(), subcell_apertures.end() );
      cell_apertures[2] = multiplier * (std::accumulate(subcell_apertures.begin(), subcell_apertures.end(), 0.0));
      subcell_apertures.clear();

      //face (i,j+1/2,k)
      s = s_max;
      for(r = 0; r < r_max; ++r) {
         for(t = 0; t < t_max; ++t) {
            faces[0] = levelset_cube[r]         [s][t]         ;
            faces[1] = levelset_cube[r+r_offset][s][t]         ;
            faces[2] = levelset_cube[r+r_offset][s][t+t_offset];
            faces[3] = levelset_cube[r]         [s][t+t_offset];
            subcell_apertures.emplace_back(CellFaceAperture(faces));
         }
      }
      std::sort( subcell_apertures.begin(), subcell_apertures.end() );
      cell_apertures[3] = multiplier * (std::accumulate(subcell_apertures.begin(), subcell_apertures.end(), 0.0));
      subcell_apertures.clear();
   } else {
      //pad with 0.0 to length 6
      cell_apertures[2] = 0.0;
      cell_apertures[3] = 0.0;
   }

   //only needed in 3D
   if constexpr(CC::DIM() == Dimension::Three) {
      //face (i,j,k-1/2)
      t = 0;
      for(r = 0; r < r_max; ++r) {
         for(s = 0; s < s_max; ++s) {
            faces[0] = levelset_cube[r]         [s]         [t];
            faces[1] = levelset_cube[r+r_offset][s]         [t];
            faces[2] = levelset_cube[r+r_offset][s+s_offset][t];
            faces[3] = levelset_cube[r]         [s+s_offset][t];
            subcell_apertures.emplace_back(CellFaceAperture(faces));
         }
      }
      std::sort( subcell_apertures.begin(), subcell_apertures.end() );
      cell_apertures[4] =  multiplier * (std::accumulate(subcell_apertures.begin(), subcell_apertures.end(), 0.0));
      subcell_apertures.clear();

      //face (i,j,k+1/2)
      t = t_max;
      for(r = 0; r < r_max; ++r) {
         for(s = 0; s < s_max; ++s) {
            faces[0] = levelset_cube[r]         [s]         [t];
            faces[1] = levelset_cube[r+r_offset][s]         [t];
            faces[2] = levelset_cube[r+r_offset][s+s_offset][t];
            faces[3] = levelset_cube[r]         [s+s_offset][t];
            subcell_apertures.emplace_back(CellFaceAperture(faces));
         }
      }
      std::sort( subcell_apertures.begin(), subcell_apertures.end() );
      cell_apertures[5] = multiplier * (std::accumulate(subcell_apertures.begin(), subcell_apertures.end(), 0.0));
      subcell_apertures.clear();

   } else {
      //pad with 0.0 to length 6
      cell_apertures[4] = 0.0;
      cell_apertures[5] = 0.0;
   }

   //consider right phase
   if(material_sign < 0) {
      for(double& cell_apertures_component : cell_apertures) {
         cell_apertures_component = 1.0 - cell_apertures_component;
      }
   }

   return cell_apertures;
}

/**
 * @brief Calculates the volume fraction of a single cell based on \cite Lauer2012 .
 * @param levelset The level-set field.
 * @param i Index i of the cell of interest.
 * @param j Index j of the cell of interest.
 * @param k Index k of the cell of interest.
 * @param material_sign Indicates for which phase the volume fraction is calculated ( > 0 -> positive phase, < 0 -> negative phase). Default: positive material.
 * @return The calculated volume fraction.
 */
double ComputeVolumeFractionCellBasedImplementation(double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
   unsigned int const i, unsigned int const j, unsigned int const k,
   std::int8_t const material_sign){

   double cell_levelset_cube[cell_box_size_x][cell_box_size_y][cell_box_size_z];
   for(unsigned int r = 0; r < cell_box_size_x; ++r) {
      for(unsigned int s = 0; s < cell_box_size_y; ++s) {
         for(unsigned int t = 0; t < cell_box_size_z; ++t) {
            cell_levelset_cube[r][s][t] = 0.0;
         }
      }
   }

   GetLevelsetAtCellCorners(cell_levelset_cube,levelset,i,j,k);

   unsigned int r = 0;
   unsigned int s = 0;
   unsigned int t = 0;

   constexpr unsigned int r_offset = 1;
   constexpr unsigned int s_offset = CC::DIM() != Dimension::One   ? 1 : 0;
   constexpr unsigned int t_offset = CC::DIM() == Dimension::Three ? 1 : 0;

   constexpr double one_third = 1.0/3.0;

   double A11 = 0.0;
   double A12 = 0.0;
   double A21 = 0.0;
   double A22 = 0.0;
   double A31 = 0.0;
   double A32 = 0.0;

   double levelset_center = 0.0;
   double volume_fraction =  0.0;
   double temp_i = 0.0;
   double temp_j = 0.0;
   double temp_k = 0.0;
   double delta_gamma = 0.0;

   std::array<double, 4> faces = {0.0, 0.0, 0.0, 0.0};

   //only for 3D the entire algorithm is required. For 2D/1D, only one patch needs to be computed.
   if constexpr(CC::DIM() == Dimension::Three) {
      std::array<double, 6> const cell_face_apertures = ComputeCellFaceApertureCellBasedImplementation(levelset,i,j,k,material_sign);
      A11 = cell_face_apertures[0];
      A12 = cell_face_apertures[1];
      A21 = cell_face_apertures[2];
      A22 = cell_face_apertures[3];
      A31 = cell_face_apertures[4];
      A32 = cell_face_apertures[5];
   } else {
      //face (i,j,k+1/2)
      //for 2D/1D, the volume fraction is equal to this single face aperture
      r = 0;
      s = 0;
      t = 0 + t_offset;
      faces[0] = cell_levelset_cube[r]          [s]          [t];
      faces[1] = cell_levelset_cube[r+r_offset][s]          [t];
      faces[2] = cell_levelset_cube[r+r_offset][s+s_offset][t];
      faces[3] = cell_levelset_cube[r]          [s+s_offset][t];
      A32 = CellFaceAperture(faces);
   }


   //the full algorithm is only required for 3D simulations
   if constexpr(CC::DIM() == Dimension::Three) {

      std::vector<double> corner_values;
      for(unsigned int r = 0; r < cell_box_size_x; ++r) {
         for(unsigned int s = 0; s < cell_box_size_y; ++s) {
            for(unsigned int t = 0; t < cell_box_size_z; ++t) {
               corner_values.push_back(cell_levelset_cube[r][s][t]);
            }
         }
      }
      std::sort( corner_values.begin(), corner_values.end() );
      levelset_center = 0.125 * (std::accumulate(corner_values.begin(), corner_values.end(), 0.0));

      temp_i = (A12-A11);
      temp_j = (A22-A21);
      temp_k = (A32-A31);

      delta_gamma = std::sqrt(temp_i*temp_i+temp_j*temp_j+temp_k*temp_k);

      volume_fraction = one_third * (0.5*A11 + 0.5*A12 + 0.5*A21 + 0.5*A22 + 0.5*A31 + 0.5*A32 + delta_gamma*levelset_center);
   } else {
      //for 2D/1D, the volume fraction is equal to this single face aperture. Be aware that this makes 1D/2D
      //simulations not fully comparable to 3D simulations
      volume_fraction = A32;
   }

   // this algorithm does not guarantee volume fractions between 0 and 1 and requires limiter
   volume_fraction = std::min(volume_fraction,1.0);
   volume_fraction = std::max(volume_fraction,0.0);

   return material_sign > 0 ? volume_fraction : 1.0 - volume_fraction;
}

/**
 * @brief Calculates the volume fraction of a single cell based on \cite Lauer2012 .
 * @param levelset The level-set field.
 * @param i Index i of the cell of interest.
 * @param j Index j of the cell of interest.
 * @param k Index k of the cell of interest.
 * @param material_sign Indicates for which phase the volume fraction is calculated ( > 0 -> positive phase, < 0 -> negative phase). Default: positive material.
 * @return The calculated volume fraction.
 */
double ComputeVolumeFractionSubCellBasedImplementation(double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
   unsigned int const i, unsigned int const j, unsigned int const k,
   std::int8_t const material_sign){

   double levelset_cube[subcell_box_size_x][subcell_box_size_y][subcell_box_size_z];
   for(unsigned int r = 0; r < subcell_box_size_x; ++r) {
      for(unsigned int s = 0; s < subcell_box_size_y; ++s) {
         for(unsigned int t = 0; t < subcell_box_size_z; ++t) {
            levelset_cube[r][s][t] = 0.0;
         }
      }
   }

   GetLevelsetAtSubcellCorners(levelset_cube,levelset,i,j,k);

   unsigned int r = 0;
   unsigned int s = 0;
   unsigned int t = 0;

   double A11 = 0.0;
   double A12 = 0.0;
   double A21 = 0.0;
   double A22 = 0.0;
   double A31 = 0.0;
   double A32 = 0.0;

   double levelset_center = 0.0;
   double subvolume_fraction =  0.0;
   double temp_i = 0.0;
   double temp_j = 0.0;
   double temp_k = 0.0;
   double delta_gamma = 0.0;

   unsigned int const r_max = 2;
   unsigned int const s_max = CC::DIM() != Dimension::One   ? 2 : 1;
   unsigned int const t_max = CC::DIM() == Dimension::Three ? 2 : 1;
   unsigned int const r_offset = 1;
   unsigned int const s_offset = CC::DIM() != Dimension::One   ? 1 : 0;
   unsigned int const t_offset = CC::DIM() == Dimension::Three ? 1 : 0;

   std::array<double, 4> faces = {0.0, 0.0, 0.0, 0.0};
   std::vector<double> subvolume_fraction_vector;
   std::vector<double> corner_values;

   constexpr double one_third = 1.0/3.0;

   //the cell volume is computed by averaging the volume fractions of the 8/4/2 subcells in 3D/2D/1D
   for(unsigned int l = 0; l < r_max; ++l) {
      for(unsigned int m = 0; m < s_max; ++m) {
         for(unsigned int n = 0; n < t_max; ++n) {

            //only for 3D the entire algorithm is required. For 2D/1D, only one patch needs to be computed.
            if constexpr(CC::DIM() == Dimension::Three) {
               //face (i-1/2,j,k)
               r = l;
               s = m;
               t = n;
               faces[0] = levelset_cube[r][s]         [t]         ;
               faces[1] = levelset_cube[r][s+s_offset][t]         ;
               faces[2] = levelset_cube[r][s+s_offset][t+t_offset];
               faces[3] = levelset_cube[r][s]         [t+t_offset];
               A11 = CellFaceAperture(faces);

               //face (i+1/2,j,k)
               r = l + r_offset;
               s = m;
               t = n;
               faces[0] = levelset_cube[r][s]         [t]         ;
               faces[1] = levelset_cube[r][s+s_offset][t]         ;
               faces[2] = levelset_cube[r][s+s_offset][t+t_offset];
               faces[3] = levelset_cube[r][s]         [t+t_offset];
               A12 = CellFaceAperture(faces);

               //face (i,j-1/2,k)
               r = l;
               s = m;
               t = n;
               faces[0] = levelset_cube[r]         [s][t]         ;
               faces[1] = levelset_cube[r+r_offset][s][t]         ;
               faces[2] = levelset_cube[r+r_offset][s][t+t_offset];
               faces[3] = levelset_cube[r]         [s][t+t_offset];
               A21 = CellFaceAperture(faces);

               //face (i,j+1/2,k)
               r = l;
               s = m + s_offset;
               t = n;
               faces[0] = levelset_cube[r]         [s][t]         ;
               faces[1] = levelset_cube[r+r_offset][s][t]         ;
               faces[2] = levelset_cube[r+r_offset][s][t+t_offset];
               faces[3] = levelset_cube[r]         [s][t+t_offset];
               A22 = CellFaceAperture(faces);

               //face (i,j,k-1/2)
               r = l;
               s = m;
               t = n;
               faces[0] = levelset_cube[r]         [s]         [t];
               faces[1] = levelset_cube[r+r_offset][s]         [t];
               faces[2] = levelset_cube[r+r_offset][s+s_offset][t];
               faces[3] = levelset_cube[r]         [s+s_offset][t];
               A31 = CellFaceAperture(faces);
            }

            //face (i,j,k+1/2)
            //for 2D/1D, the volume fraction is equal to this single face aperture
            r = l;
            s = m;
            t = n + t_offset;
            faces[0] = levelset_cube[r]         [s]         [t];
            faces[1] = levelset_cube[r+r_offset][s]         [t];
            faces[2] = levelset_cube[r+r_offset][s+s_offset][t];
            faces[3] = levelset_cube[r]         [s+s_offset][t];
            A32 = CellFaceAperture(faces);

            //the full algorithm is only required for 3D simulations
            if constexpr(CC::DIM() == Dimension::Three) {
               for(unsigned int o=0; o<r_max; ++o) {
                  for(unsigned int p=0; p<s_max; ++p) {
                     for(unsigned int q=0; q<t_max; ++q) {
                        corner_values.emplace_back(levelset_cube[l+o][m+p][n+q]);
                     }
                  }
               }
               std::sort( corner_values.begin(), corner_values.end() );
               levelset_center = 0.125*(std::accumulate(corner_values.begin(), corner_values.end(), 0.0));
               corner_values.clear();

               temp_i = ((A12-A11)*0.25);
               temp_j = ((A22-A21)*0.25);
               temp_k = ((A32-A31)*0.25);

               delta_gamma = std::sqrt(temp_i*temp_i+temp_j*temp_j+temp_k*temp_k);

               subvolume_fraction = one_third * (0.5*A11 + 0.5*A12 + 0.5*A21 + 0.5*A22 + 0.5*A31 + 0.5*A32 +
                  ((delta_gamma*levelset_center) * 8));
            } else {
               //for 2D/1D, the volume fraction is equal to this single face aperture. Be aware that this makes 1D/2D
               //simulations not fully comparable to 3D simulations
               subvolume_fraction = A32;
            }

            // this algorithm does not guarantee volume fractions between 0 and 1 and requires limiter
            subvolume_fraction = std::min(subvolume_fraction,1.0);
            subvolume_fraction = std::max(subvolume_fraction,0.0);

            subvolume_fraction_vector.emplace_back(subvolume_fraction);
         }
      }
   }

   //normalization by the number of subvolume which are used
   double const multiplier = CC::DIM() == Dimension::Three ? 0.125 : (CC::DIM() != Dimension::One ? 0.25 : 0.5);

   std::sort( subvolume_fraction_vector.begin(), subvolume_fraction_vector.end() );
   double volume_fraction = multiplier*(std::accumulate(subvolume_fraction_vector.begin(), subvolume_fraction_vector.end(), 0.0));

   return material_sign > 0 ? volume_fraction : 1.0 - volume_fraction;
}

}

/**
 * @brief Calculates the cell face apertures according to \cite Lauer2012 .
 * @param levelset The level-set field.
 * @param i Index i of the cell of interest.
 * @param j Index j of the cell of interest.
 * @param k Index k of the cell of interest.
 * @param material_sign Indicates for which phase the cell-face apertures are calculated ( > 0 -> positive phase, < 0 -> negative phase). Default: positive material.
 * @return Cell-face apertures.
 */
std::array<double, 6> GeometryCalculatorMarchingCubes::ComputeCellFaceApertureImplementation(double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
   unsigned int const i, unsigned int const j, unsigned int const k,
   std::int8_t const material_sign) const {
   if(cell_based_geometry_calculations_) {
      return ComputeCellFaceApertureCellBasedImplementation(levelset, i, j, k, material_sign);
   } else {
      return ComputeCellFaceApertureSubCellBasedImplementation(levelset, i, j, k, material_sign);
   }
}

/**
 * @brief Calculates the volume fraction of a single cell based on \cite Lauer2012 .
 * @param levelset The level-set field.
 * @param i Index i of the cell of interest.
 * @param j Index j of the cell of interest.
 * @param k Index k of the cell of interest.
 * @param material_sign Indicates for which phase the volume fraction is calculated ( > 0 -> positive phase, < 0 -> negative phase). Default: positive material.
 * @return The calculated volume fraction.
 */
double GeometryCalculatorMarchingCubes::ComputeVolumeFractionImplementation(double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
   unsigned int const i, unsigned int const j, unsigned int const k,
   std::int8_t const material_sign) const {
   if(cell_based_geometry_calculations_) {
      return ComputeVolumeFractionCellBasedImplementation(levelset, i, j, k, material_sign);
   } else {
      return ComputeVolumeFractionSubCellBasedImplementation(levelset, i, j, k, material_sign);
   }
}
