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
#ifndef MATERIAL_SIGN_CAPSULE_H
#define MATERIAL_SIGN_CAPSULE_H

#include "materials/material_definitions.h"

/**
 * @brief Work-around class to reduce cyclic dependencies and still query the material sign where needed.
 *        It is the user's responsibility to call functions only after initialization (currently done in the constructor of the material manager
 *        class that holds all information for the two materials).
 */
class MaterialSignCapsule {

private:
   static MaterialName positive_material_;
   static MaterialName negative_material_;

public:
   MaterialSignCapsule() = delete;
   ~MaterialSignCapsule() = default;
   MaterialSignCapsule( MaterialSignCapsule const& ) = delete;
   MaterialSignCapsule& operator=( MaterialSignCapsule const& ) = delete;
   MaterialSignCapsule( MaterialSignCapsule&& ) = delete;
   MaterialSignCapsule& operator=( MaterialSignCapsule&& ) = delete;

   MaterialSignCapsule( MaterialName const positive_material, MaterialName const negative_material ) {
      positive_material_ = positive_material;
      negative_material_ = negative_material;
   }

   /**
    * @brief Static function to get the positive material material identifier in a single-level-set simulation.
    *        $Always MATERIAL ONE in inputfile. This function can be uninitialized if called too early! Must not be called as long as no object is available.$
    * @return Positive material material identifier .
    */
   static inline MaterialName PositiveMaterial() { return positive_material_; }

   /**
    * @brief Static function to get the negative material material identifier in a single-level-set simulation.
    *        $Always MATERIAL TWO in inputfile. This function can be uninitialized if called too early! Must not be called as long as no object is available.$
    * @return Negative material material identifier.
    */
   static inline MaterialName NegativeMaterial() { return negative_material_; }

   /**
    * @brief Gives the sign of the given material used in the signed levelset and signed interface tag description.
    * @param material Material of interest.
    * @return return Sign of the material.
    */
   static inline std::int8_t SignOfMaterial( MaterialName const material ) { return material == PositiveMaterial() ? 1 : -1; }
};

#endif // MATERIAL_SIGN_CAPSULE_H
