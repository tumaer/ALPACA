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
#ifndef HALO_MANAGER_H
#define HALO_MANAGER_H

#include "communication/internal_halo_manager.h"
#include "boundary_condition/external_halo_manager.h"
#include "communication/communication_manager.h"

class HaloManager {

private:
   Tree & tree_;
   ExternalHaloManager const& external_halo_manager_;
   InternalHaloManager & internal_halo_manager_;
   CommunicationManager const& communication_manager_;
   unsigned int const maximum_level_;

public:
   HaloManager() = delete;
   explicit HaloManager( Tree& tree, ExternalHaloManager const& external_halo_manager, InternalHaloManager & internal_halo_manager,
                         CommunicationManager const& communication_manager, unsigned int const maximum_level );
   ~HaloManager() = default;
   HaloManager( HaloManager const& ) = delete;
   HaloManager& operator=( HaloManager const& ) = delete;
   HaloManager( HaloManager&& ) = delete;
   HaloManager& operator=( HaloManager&& ) = delete;

   void MaterialHaloUpdate( std::vector<unsigned int> const& levels_ascending, MaterialFieldType const field_type, bool const cut_jumps = false ) const;
   void MaterialHaloUpdateOnLevel( unsigned int const level, MaterialFieldType const field_type, bool const cut_jumps = false ) const;
   void MaterialHaloUpdateOnLmax( MaterialFieldType const field_type, bool const cut_jumps = true ) const;
   void MaterialHaloUpdateOnLmaxMultis( MaterialFieldType const field_type ) const;

   void MaterialInternalHaloUpdateOnLevel( unsigned int const level, MaterialFieldType const field_type, bool const cut_jumps = false ) const;
   void MaterialExternalHaloUpdateOnLevel( unsigned int const level, MaterialFieldType const field_type ) const;

   void InterfaceTagHaloUpdateOnLevelList( std::vector<unsigned int> const& updated_levels ) const;
   void InterfaceTagHaloUpdateOnLmax() const;

   void InterfaceHaloUpdateOnLevelList( std::vector<unsigned int> const updated_levels, InterfaceBlockBufferType const halo_type ) const;
   void InterfaceHaloUpdateOnLmax( InterfaceBlockBufferType const type ) const;
};

#endif //HALO_MANAGER_H