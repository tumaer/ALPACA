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
* 5. ApprovalTests.cpp  : See LICENSE_APPROVAL_TESTS.txt for more information            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, February 10th, 2021                                                            *
*                                                                                        *
*****************************************************************************************/
#include "interface_tag_functions.h"

#include "enums/interface_tag_definition.h"
#include "utilities/mathematical_functions.h"
#include "levelset/geometry/geometry_calculator.h"
#include "user_specifications/two_phase_constants.h"

namespace InterfaceTagFunctions {

   /**
 * @brief Initializes all tags by default as cut-cells. This is necessary for the function InterfaceTagManager::SetInternalCutCellTagsFromLevelset.
 * @param interface_tags Indirect return parameter for the interface tags.
 */
   void InitializeInternalInterfaceTags( std::int8_t ( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] ) {

      for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
               interface_tags[i][j][k] = ITTI( IT::OldCutCell );
            }//k
         }   //j
      }      //i
   }

   /**
 * @brief Sets the internal cut-cell tags according to the internal cells of the given levelset buffer.
 * @param levelset Reference to the levelset values.
 * @param interface_tags Indirect return parameter for the derived interface tags.
 */
   void SetInternalCutCellTagsFromLevelset( double const ( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()], std::int8_t ( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] ) {

      //find cut-cells in internal cells - due to the CFL condition, it is enough to test only those cells which were cut cells or neighbors before
      //of cut cells
      for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
               //cell was close to the interface in previous step
               if( std::abs( interface_tags[i][j][k] ) < ITTI( IT::ExtensionBand ) || std::abs( interface_tags[i][j][k] ) == ITTI( IT::ScaleSeparatedCell ) ) {
                  if( IsCutCell<GeometryCalculationSettings::CutCellCriteria>( levelset, i, j, k ) ) {
                     //this cut-cell was a cut-cell in the previous step
                     if( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::NewCutCell ) )
                        interface_tags[i][j][k] = ITTI( IT::OldCutCell );
                     //this cut-cell was not a cut-cell before, but is now
                     else
                        interface_tags[i][j][k] = Signum( levelset[i][j][k] ) * ITTI( IT::NewCutCell );
                  } else {
                     //this is not a cut-cell -> reset to bulk phase. Later, this cell interface tag is correctly set
                     //in function InterfaceTagManager::SetTotalInterfaceTagsFromCutCells.
                     interface_tags[i][j][k] = Signum( levelset[i][j][k] ) * ITTI( IT::BulkPhase );
                  }
               } else {
                  //this cannot be a cut-cell, as it was too far away from the interface -> reset it to bulk phase. Later, this
                  //cell interface tag is correctly set in function InterfaceTagManager::SetTotalInterfaceTagsFromCutCells.
                  interface_tags[i][j][k] = Signum( levelset[i][j][k] ) * ITTI( IT::BulkPhase );
               }
            }//k
         }   //j
      }      //i
   }

   /**
 * @brief Updates the interface tags with respect to all cut cells within the block.
 * @param interface_tags Reference to the interface tag buffer that has to be updated and already holds the cut cell tags.
 */
   void SetTotalInterfaceTagsFromCutCells( std::int8_t ( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] ) {

      //these type-casts are necessary to allow for checking for cut-cells in the outermost halo cell
      constexpr int total_cells_x                   = static_cast<int>( CC::TCX() );
      constexpr int total_cells_y                   = static_cast<int>( CC::TCY() );
      constexpr int total_cells_z                   = static_cast<int>( CC::TCZ() );
      constexpr int cut_cell_neighbour_band_width_x = 1;
      constexpr int cut_cell_neighbour_band_width_y = CC::DIM() != Dimension::One ? 1 : 0;
      constexpr int cut_cell_neighbour_band_width_z = CC::DIM() == Dimension::Three ? 1 : 0;
      constexpr int extension_band_width_x          = static_cast<int>( CC::EBW() );
      constexpr int extension_band_width_y          = CC::DIM() != Dimension::One ? static_cast<int>( CC::EBW() ) : 0;
      constexpr int extension_band_width_z          = CC::DIM() == Dimension::Three ? static_cast<int>( CC::EBW() ) : 0;
      constexpr int reinitialization_band_width_x   = static_cast<int>( CC::RBW() );
      constexpr int reinitialization_band_width_y   = CC::DIM() != Dimension::One ? static_cast<int>( CC::RBW() ) : 0;
      constexpr int reinitialization_band_width_z   = CC::DIM() == Dimension::Three ? static_cast<int>( CC::RBW() ) : 0;

      //set interface tags from cut-cells in full block
      for( int i = 0; i < total_cells_x; ++i ) {
         for( int j = 0; j < total_cells_y; ++j ) {
            for( int k = 0; k < total_cells_z; ++k ) {

               //find all cells that are now cut cells
               if( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::NewCutCell ) ) {

                  //tag cells which are neighbors to cut cells
                  for( int r = std::max( i - cut_cell_neighbour_band_width_x, 0 ); r <= std::min( i + cut_cell_neighbour_band_width_x, total_cells_x - 1 ); ++r ) {
                     for( int s = std::max( j - cut_cell_neighbour_band_width_y, 0 ); s <= std::min( j + cut_cell_neighbour_band_width_y, total_cells_y - 1 ); ++s ) {
                        for( int t = std::max( k - cut_cell_neighbour_band_width_z, 0 ); t <= std::min( k + cut_cell_neighbour_band_width_z, total_cells_z - 1 ); ++t ) {
                           if( std::abs( interface_tags[r][s][t] ) > ITTI( IT::NewCutCell ) )
                              interface_tags[r][s][t] = Signum( interface_tags[r][s][t] ) * ITTI( IT::CutCellNeighbor );
                        }//t
                     }   //s
                  }      //r

                  //tag the extension band
                  for( int r = std::max( i - extension_band_width_x, 0 ); r <= std::min( i + extension_band_width_x, total_cells_x - 1 ); ++r ) {
                     for( int s = std::max( j - extension_band_width_y, 0 ); s <= std::min( j + extension_band_width_y, total_cells_y - 1 ); ++s ) {
                        for( int t = std::max( k - extension_band_width_z, 0 ); t <= std::min( k + extension_band_width_z, total_cells_z - 1 ); ++t ) {
                           if( std::abs( interface_tags[r][s][t] ) > ITTI( IT::CutCellNeighbor ) )
                              interface_tags[r][s][t] = Signum( interface_tags[r][s][t] ) * ITTI( IT::ExtensionBand );
                        }//t
                     }   //s
                  }      //r

                  // tag the reinitialization band
                  for( int r = std::max( i - reinitialization_band_width_x, 0 ); r <= std::min( i + reinitialization_band_width_x, total_cells_x - 1 ); ++r ) {
                     for( int s = std::max( j - reinitialization_band_width_y, 0 ); s <= std::min( j + reinitialization_band_width_y, total_cells_y - 1 ); ++s ) {
                        for( int t = std::max( k - reinitialization_band_width_z, 0 ); t <= std::min( k + reinitialization_band_width_z, total_cells_z - 1 ); ++t ) {
                           if( std::abs( interface_tags[r][s][t] ) > ITTI( IT::ExtensionBand ) )
                              interface_tags[r][s][t] = Signum( interface_tags[r][s][t] ) * ITTI( IT::ReinitializationBand );
                        }//t
                     }   //s
                  }      //r
               }
            }//k
         }   //j
      }      //i
   }

   /**
 * @brief Gives whether or not the given interface tags are uniform (hence, a bulk phase).
 * @param interface_tags Reference to the interface tags to be checked.
 * @return True if all interface tags are uniform, false otherwise.
 */
   bool TotalInterfaceTagsAreUniform( std::int8_t const ( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] ) {
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               if( std::abs( interface_tags[i][j][k] ) != ITTI( IT::BulkPhase ) ) {
                  return false;
               }
            }//k
         }   //j
      }      //i
      return true;
   }

}// namespace InterfaceTagFunctions
