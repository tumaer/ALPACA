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
#ifndef SPATIAL_RECONSTRUCTION_STENCIL_H
#define SPATIAL_RECONSTRUCTION_STENCIL_H

#include <vector>
#include <limits>

/**
 * @brief The SpatialReconstructionStencil class provides an interface for computational Stencils needed throughout the simulation.
 */
template<typename DerivedSpatialReconstructionStencil>
class SpatialReconstructionStencil {

   friend DerivedSpatialReconstructionStencil;

   explicit SpatialReconstructionStencil() = default;

protected:
   static constexpr double epsilon_ = std::numeric_limits<double>::epsilon();

public:
   virtual ~SpatialReconstructionStencil() = default;
   SpatialReconstructionStencil( SpatialReconstructionStencil const& ) = delete;
   SpatialReconstructionStencil& operator=( SpatialReconstructionStencil const& ) = delete;
   SpatialReconstructionStencil( SpatialReconstructionStencil&& ) = delete;
   SpatialReconstructionStencil& operator=( SpatialReconstructionStencil&& ) = delete;

   /**
    * @brief Applies the SpatialReconstructionStencil to the provided Array
    * @param evaluation_properties[0] Gives an offset to be used to weight the data in upwind direction stronger.
    * @param evaluation_properties[1] Indicates from which direction the solution should be calculated. Positive (1) or negative (-1) direction is possible.
    * @param cell_size .
    * @return Value at the position of interest.
    * @note Hotpath function.
    */
   double Apply(const std::vector<double>& array, const int evaluation_properties[0], const int evaluation_properties[1], const double cell_size) const {
      return static_cast<DerivedSpatialReconstructionStencil const&>(*this).ApplyImplementation(array, evaluation_properties[0], evaluation_properties[1], cell_size);
   }
   /**
    * @brief Gives the number of cells needed for a single stencil evaluation.
    * @return Size of the stencil, i.e. number of data cells the stencil works on.
    */
   static constexpr unsigned int StencilSize() {return DerivedSpatialReconstructionStencil::stencil_size_;}
   /**
    * @brief Gives the size of the stencil in down stream direction.
    * @return Size of the stencil arm reaching down stream, i.e. number of data cells that lay down stream the stencil works on.
    */
   static constexpr unsigned int DownstreamStencilSize() {return DerivedSpatialReconstructionStencil::downstream_stencil_size_;}
};

#endif // SPATIAL_RECONSTRUCTION_STENCIL_H
