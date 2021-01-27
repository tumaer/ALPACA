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
#ifndef MATERIAL_MANAGER_H
#define MATERIAL_MANAGER_H

#include <memory>
#include <vector>

#include "materials/material_type_definitions.h"
#include "materials/material_definitions.h"
#include "materials/material.h"
#include "materials/material_pairing.h"
#include "unit_handler.h"

/**
 * @brief The MaterialManager class provides access to all materials and material pairings present in the current simulation and forwards the appropriate
 *        object to the caller. The MaterialManager does not change any data. It provides the functionality to map material names and indices to the correct
 *        material or pairing class.
 */
class MaterialManager {
   // Vector with all materials
   std::vector<std::tuple<MaterialType, Material>> const materials_;
   // Vector with all material pairings
   std::vector<MaterialPairing> const material_pairings_;

   // offset vector to provide proper mapping from material names to its material pairings
   std::vector<int> const pairing_offset_;

   // local function to map a pairing of two materials to  corresponding vector index of the material pairings
   unsigned int MapPairingToIndex( MaterialName const first_material, MaterialName const second_material ) const;
   // factory function to generate the pairing offset vector (allows constness of it)
   std::vector<int> GenerateMaterialPairingOffset( std::size_t const number_of_materials ) const;

public:
   MaterialManager() = delete;
   explicit MaterialManager( std::vector<std::tuple<MaterialType, Material>> materials,
                             std::vector<MaterialPairing> material_pairings );
   ~MaterialManager()                        = default;
   MaterialManager( MaterialManager const& ) = delete;
   MaterialManager& operator=( MaterialManager const& ) = delete;
   MaterialManager( MaterialManager&& )                 = delete;
   MaterialManager& operator=( MaterialManager&& ) = delete;

   // Get the number of all materials contained in the simulation
   std::size_t GetNumberOfMaterials() const;

   // Get all materials used in the simulation
   std::vector<MaterialName> GetMaterialNames() const;

   // provides a material for a given index
   Material const& GetMaterial( std::size_t const index ) const;

   // provides a material for a given identifier
   Material const& GetMaterial( MaterialName const material ) const;

   // provides the material type for a given identifier
   MaterialType GetMaterialType( MaterialName const material ) const;

   // provides the pairing of two material identifier
   MaterialPairing const& GetMaterialPairing( MaterialName const first_material, MaterialName const second_material ) const;

   // checks whether the given material is a solid boundary
   bool IsSolidBoundary( MaterialName const material ) const;
};

#endif// MATERIAL_MANAGER_H
