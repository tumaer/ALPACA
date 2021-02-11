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
#include "id_information.h"

#include <bitset>
#include "user_specifications/compile_time_constants.h"

namespace {
   /**
 * Gives the initial id for the most bottom-south-west node on all levels (number must equal maximum number of allowed levels CC::AMNL())
 */
   // NH: DO NOT TOUCH UNLESS YOU HAVE MY WRITTEN PERMISSION!
   constexpr std::array<nid_t, CC::AMNL()> headbits =
         { { 0x1400000, 0xA000000, 0x50000000, 0x280000000, 0x1400000000, 0xA000000000,
             0x50000000000, 0x280000000000, 0x1400000000000, 0xA000000000000, 0x50000000000000, 0x280000000000000,
             0x1400000000000000, 0xA000000000000000 } };
}// namespace

/**
 * @brief Gives the inital id seed. I. e. the id of the most bottom-south-west node on level zero.
 * @return seed.
 */
nid_t IdSeed() {
   return headbits[0];
}

/**
 * @brief Transforms the id into a normal Morton Order id (Z-Curve).
 * @param id Node id that is transformed.
 * @param level Level on which the node lies.
 * @return Morton Id.
 */
nid_t CutHeadBit( nid_t const id, unsigned int const level ) {
   return id - headbits[level];
}

/**
 * @brief Transforms a normal Morton Order id (Z-Curve) into an ALPACA id.
 * @param id Node id that is transformed.
 * @param level Level on which the node lies.
 * @return ALPACA Id.
 */
nid_t AddHeadBit( nid_t const id, unsigned int const level ) {
   return id + headbits[level];
}

/**
 * @brief Calculates the id of the eastern neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @return Id of eastern neighbor.
 */
nid_t EastNeighborOfNodeWithId( nid_t const id ) {

   unsigned int const level = LevelOfNode( id );

   std::bitset<1> indicator = 0;
   std::bitset<64> working_id( CutHeadBit( id, level ) );

   for( std::uint8_t i = 0; i < 60; i += 3 ) {
      working_id[i] = ( ( ~working_id[i] & ~indicator[0] ) | ( working_id[i] & indicator[0] ) );
      indicator[0]  = ( working_id[i] | indicator[0] );
   }

   return AddHeadBit( working_id.to_ullong(), level );
}

/**
 * @brief Calculates the id of the western neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @return Id of western neighbor.
 */
nid_t WestNeighborOfNodeWithId( nid_t const id ) {

   unsigned int const level = LevelOfNode( id );

   std::bitset<1> indicator = 0;
   std::bitset<64> working_id( CutHeadBit( id, level ) );

   for( std::uint8_t i = 0; i < 60; i += 3 ) {
      working_id[i] = ( ( ~working_id[i] & ~indicator[0] ) | ( working_id[i] & indicator[0] ) );
      indicator[0]  = ( ~working_id[i] | indicator[0] );
   }

   return AddHeadBit( working_id.to_ullong(), level );
}

/**
 * @brief Calculates the id of the northern neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @return Id of northern neighbor.
 */
nid_t NorthNeighborOfNodeWithId( nid_t const id ) {

   unsigned int const level = LevelOfNode( id );

   std::bitset<1> indicator = 0;
   std::bitset<64> working_id( CutHeadBit( id, level ) );

   for( std::uint8_t i = 1; i < 60; i += 3 ) {
      working_id[i] = ( ( ~working_id[i] & ~indicator[0] ) | ( working_id[i] & indicator[0] ) );
      indicator[0]  = ( working_id[i] | indicator[0] );
   }

   return AddHeadBit( working_id.to_ullong(), level );
}

/**
 * @brief Calculates the id of the southern neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @return Id of southern neighbor.
 */
nid_t SouthNeighborOfNodeWithId( nid_t const id ) {

   unsigned int const level = LevelOfNode( id );

   std::bitset<1> indicator = 0;
   std::bitset<64> working_id( CutHeadBit( id, level ) );

   for( std::uint8_t i = 1; i < 60; i += 3 ) {
      working_id[i] = ( ( ~working_id[i] & ~indicator[0] ) | ( working_id[i] & indicator[0] ) );
      indicator[0]  = ( ~working_id[i] | indicator[0] );
   }

   return AddHeadBit( working_id.to_ullong(), level );
}

/**
 * @brief Calculates the id of the top neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @return Id of top neighbor.
 */
nid_t TopNeighborOfNodeWithId( nid_t const id ) {

   unsigned int const level = LevelOfNode( id );

   std::bitset<1> indicator = 0;
   std::bitset<64> working_id( CutHeadBit( id, level ) );

   for( std::uint8_t i = 2; i < 60; i += 3 ) {
      working_id[i] = ( ( ~working_id[i] & ~indicator[0] ) | ( working_id[i] & indicator[0] ) );
      indicator[0]  = ( working_id[i] | indicator[0] );
   }

   return AddHeadBit( working_id.to_ullong(), level );
}

/**
 * @brief Calculates the id of the bottom neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @return Id of bottom neighbor.
 */
nid_t BottomNeighborOfNodeWithId( nid_t const id ) {

   unsigned int const level = LevelOfNode( id );

   std::bitset<1> indicator = 0;
   std::bitset<64> working_id( CutHeadBit( id, level ) );

   for( std::uint8_t i = 2; i < 60; i += 3 ) {
      working_id[i] = ( ( ~working_id[i] & ~indicator[0] ) | ( working_id[i] & indicator[0] ) );
      indicator[0]  = ( ~working_id[i] | indicator[0] );
   }

   return AddHeadBit( working_id.to_ullong(), level );
}

/**
 * @brief Gives the id of a neighbor at the provided direction.
 * @param id The id of the node whose neighbor is to be found.
 * @param location Direction in which the neighbor is located.
 * @return Id of the neighbor.
 */
nid_t GetNeighborId( nid_t const id, BoundaryLocation const location ) {

   switch( location ) {
      // Natural (planes)
      case BoundaryLocation::East:
         return EastNeighborOfNodeWithId( id );
      case BoundaryLocation::West:
         return WestNeighborOfNodeWithId( id );
      case BoundaryLocation::North:
         return NorthNeighborOfNodeWithId( id );
      case BoundaryLocation::South:
         return SouthNeighborOfNodeWithId( id );
      case BoundaryLocation::Top:
         return TopNeighborOfNodeWithId( id );
      case BoundaryLocation::Bottom:
         return BottomNeighborOfNodeWithId( id );

      // Sticks
      case BoundaryLocation::BottomNorth:
         return BottomNeighborOfNodeWithId( NorthNeighborOfNodeWithId( id ) );
      case BoundaryLocation::BottomSouth:
         return BottomNeighborOfNodeWithId( SouthNeighborOfNodeWithId( id ) );
      case BoundaryLocation::TopNorth:
         return TopNeighborOfNodeWithId( NorthNeighborOfNodeWithId( id ) );
      case BoundaryLocation::TopSouth:
         return TopNeighborOfNodeWithId( SouthNeighborOfNodeWithId( id ) );

      case BoundaryLocation::BottomEast:
         return BottomNeighborOfNodeWithId( EastNeighborOfNodeWithId( id ) );
      case BoundaryLocation::BottomWest:
         return BottomNeighborOfNodeWithId( WestNeighborOfNodeWithId( id ) );
      case BoundaryLocation::TopEast:
         return TopNeighborOfNodeWithId( EastNeighborOfNodeWithId( id ) );
      case BoundaryLocation::TopWest:
         return TopNeighborOfNodeWithId( WestNeighborOfNodeWithId( id ) );

      case BoundaryLocation::NorthEast:
         return NorthNeighborOfNodeWithId( EastNeighborOfNodeWithId( id ) );
      case BoundaryLocation::NorthWest:
         return NorthNeighborOfNodeWithId( WestNeighborOfNodeWithId( id ) );
      case BoundaryLocation::SouthEast:
         return SouthNeighborOfNodeWithId( EastNeighborOfNodeWithId( id ) );
      case BoundaryLocation::SouthWest:
         return SouthNeighborOfNodeWithId( WestNeighborOfNodeWithId( id ) );
      // Cubes
      case BoundaryLocation::EastNorthTop:
         return EastNeighborOfNodeWithId( NorthNeighborOfNodeWithId( TopNeighborOfNodeWithId( id ) ) );
      case BoundaryLocation::EastNorthBottom:
         return EastNeighborOfNodeWithId( NorthNeighborOfNodeWithId( BottomNeighborOfNodeWithId( id ) ) );
      case BoundaryLocation::EastSouthTop:
         return EastNeighborOfNodeWithId( SouthNeighborOfNodeWithId( TopNeighborOfNodeWithId( id ) ) );
      case BoundaryLocation::EastSouthBottom:
         return EastNeighborOfNodeWithId( SouthNeighborOfNodeWithId( BottomNeighborOfNodeWithId( id ) ) );

      case BoundaryLocation::WestNorthTop:
         return WestNeighborOfNodeWithId( NorthNeighborOfNodeWithId( TopNeighborOfNodeWithId( id ) ) );
      case BoundaryLocation::WestNorthBottom:
         return WestNeighborOfNodeWithId( NorthNeighborOfNodeWithId( BottomNeighborOfNodeWithId( id ) ) );
      case BoundaryLocation::WestSouthTop:
         return WestNeighborOfNodeWithId( SouthNeighborOfNodeWithId( TopNeighborOfNodeWithId( id ) ) );
#ifndef PERFORMANCE
      case BoundaryLocation::WestSouthBottom:
         return WestNeighborOfNodeWithId( SouthNeighborOfNodeWithId( BottomNeighborOfNodeWithId( id ) ) );
      default:
         throw std::invalid_argument( "Boundary Location not Found in Node::GetNeighborId() - Impossible Error" );
#else
      default: /* BoundaryLocation::WestSouthBottom */
         return WestNeighborOfNodeWithId( SouthNeighborOfNodeWithId( BottomNeighborOfNodeWithId( id ) ) );
#endif
   }
}

/**
 * @brief Determines whether a node face is also an external or domain edge.
 * @param location The direction of the edge under consideration.
 * @param id The id of the node under investigation.
 * @param number_of_nodes_on_level_zero Number of nodes on level zero in the three Cartesian directions.
 * @return True if the edge is a domain edge, false otherwise, i.e. internal edge.
 * @note Does not check for dimensionality! I. e. callers responsibility to only call on existing locations (e. g. NOT Top in 1D).
 */
bool IsNaturalExternalBoundary( BoundaryLocation const location, nid_t const id, std::array<unsigned int, 3> const number_of_nodes_on_level_zero ) {

   std::bitset<64> input( CutHeadBit( id, LevelOfNode( id ) ) );

   switch( location ) {
      case BoundaryLocation::West: {
         std::bitset<64> west_mask( 0x249249249249249 );
         std::bitset<64> result;
         result = input & west_mask;
         if( !result.to_ullong() ) {
            return true;
         }
      } break;
      case BoundaryLocation::South: {
         std::bitset<64> south_mask( 0x492492492492492 );
         std::bitset<64> result;
         result = input & south_mask;
         if( !result.to_ullong() ) {
            return true;
         }
      } break;
      case BoundaryLocation::Bottom: {
         std::bitset<64> bottom_mask( 0x924924924924924 );
         std::bitset<64> result;
         result = input & bottom_mask;
         if( !result.to_ullong() ) {
            return true;
         }
      } break;
      case BoundaryLocation::East: {
         std::bitset<1> indicator( 1 );
         std::uint32_t tows_exponent = 1;
         tows_exponent <<= LevelOfNode( id );
         std::bitset<20> east_mask( ( number_of_nodes_on_level_zero[0] * tows_exponent ) - 1 );
         for( unsigned int i = 0; i < 20; ++i ) {
            indicator[0] = ( ~( input[i * 3] ^ east_mask[i] ) & indicator[0] );
         }
         if( indicator[0] ) {
            return true;
         }
      } break;
      case BoundaryLocation::North: {
         std::bitset<1> indicator( 1 );
         std::uint32_t tows_exponent = 1;
         tows_exponent <<= LevelOfNode( id );
         ;
         std::bitset<20> north_mask( ( number_of_nodes_on_level_zero[1] * tows_exponent ) - 1 );
         for( unsigned int i = 0; i < 20; ++i ) {
            indicator[0] = ( ~( input[( i * 3 ) + 1] ^ north_mask[i] ) & indicator[0] );
         }
         if( indicator[0] ) {
            return true;
         }
      } break;
#ifndef PERFORMANCE
      case BoundaryLocation::Top: {
         std::bitset<1> indicator( 1 );
         std::uint32_t tows_exponent = 1;
         tows_exponent <<= LevelOfNode( id );
         ;
         std::bitset<20> top_mask( ( number_of_nodes_on_level_zero[2] * tows_exponent ) - 1 );
         for( unsigned int i = 0; i < 20; ++i ) {
            indicator[0] = ( ~( input[( i * 3 ) + 2] ^ top_mask[i] ) & indicator[0] );
         }
         if( indicator[0] ) {
            return true;
         }
      } break;
      default: {
         throw std::invalid_argument( "Boundary Type in IsExternal does not exist" );
      }
#else
      default: /* BoundaryLocation::Top */ {
         std::bitset<1> indicator( 1 );
         std::uint32_t tows_exponent = 1;
         tows_exponent <<= LevelOfNode( id );
         ;
         std::bitset<20> top_mask( ( number_of_nodes_on_level_zero[2] * tows_exponent ) - 1 );
         for( unsigned int i = 0; i < 20; ++i ) {
            indicator[0] = ( ~( input[( i * 3 ) + 2] ^ top_mask[i] ) & indicator[0] );
         }
         if( indicator[0] ) {
            return true;
         }
      }
#endif
   }

   return false;
}

/**
 * @brief Determines whether a location (including edges and corners) of a block is at the edge of the computational domain.
 * @param location The  direction of the edge under consideration.
 * @param id The id of the node under investigation.
 * @param number_of_nodes_on_level_zero Number of nodes on level zero in the three Cartesian directions.
 * @return True if the edge is a domain edge, false otherwise, i.e. internal edge.
 * @note Does not check for dimensionality! I. e. callers responsibility to only call on existing locations (e. g. NOT Top in 1D).
 */
bool IsExternalBoundary( BoundaryLocation const location, nid_t const id, std::array<unsigned int, 3> const number_of_nodes_on_level_zero ) {
   //natural | NH Such comparison are okay by (enforced) definiton of BoundaryLocation
   if( LTI( location ) <= LTI( BoundaryLocation::Bottom ) ) {
      return IsNaturalExternalBoundary( location, id, number_of_nodes_on_level_zero );
   }

   switch( location ) {
      // Sticks
      case BoundaryLocation::BottomNorth:
         return ( IsNaturalExternalBoundary( BoundaryLocation::Bottom, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::North, id, number_of_nodes_on_level_zero ) );
      case BoundaryLocation::BottomSouth:
         return ( IsNaturalExternalBoundary( BoundaryLocation::Bottom, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::South, id, number_of_nodes_on_level_zero ) );
      case BoundaryLocation::TopNorth:
         return ( IsNaturalExternalBoundary( BoundaryLocation::Top, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::North, id, number_of_nodes_on_level_zero ) );
      case BoundaryLocation::TopSouth:
         return ( IsNaturalExternalBoundary( BoundaryLocation::Top, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::South, id, number_of_nodes_on_level_zero ) );

      case BoundaryLocation::BottomEast:
         return ( IsNaturalExternalBoundary( BoundaryLocation::Bottom, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::East, id, number_of_nodes_on_level_zero ) );
      case BoundaryLocation::BottomWest:
         return ( IsNaturalExternalBoundary( BoundaryLocation::Bottom, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::West, id, number_of_nodes_on_level_zero ) );
      case BoundaryLocation::TopEast:
         return ( IsNaturalExternalBoundary( BoundaryLocation::Top, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::East, id, number_of_nodes_on_level_zero ) );
      case BoundaryLocation::TopWest:
         return ( IsNaturalExternalBoundary( BoundaryLocation::Top, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::West, id, number_of_nodes_on_level_zero ) );

      case BoundaryLocation::NorthEast:
         return ( IsNaturalExternalBoundary( BoundaryLocation::North, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::East, id, number_of_nodes_on_level_zero ) );
      case BoundaryLocation::NorthWest:
         return ( IsNaturalExternalBoundary( BoundaryLocation::North, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::West, id, number_of_nodes_on_level_zero ) );
      case BoundaryLocation::SouthEast:
         return ( IsNaturalExternalBoundary( BoundaryLocation::South, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::East, id, number_of_nodes_on_level_zero ) );
      case BoundaryLocation::SouthWest:
         return ( IsNaturalExternalBoundary( BoundaryLocation::South, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::West, id, number_of_nodes_on_level_zero ) );

      // Cubes
      case BoundaryLocation::EastNorthTop:
         return ( IsNaturalExternalBoundary( BoundaryLocation::East, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::North, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::Top, id, number_of_nodes_on_level_zero ) );
      case BoundaryLocation::EastNorthBottom:
         return ( IsNaturalExternalBoundary( BoundaryLocation::East, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::North, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::Bottom, id, number_of_nodes_on_level_zero ) );
      case BoundaryLocation::EastSouthTop:
         return ( IsNaturalExternalBoundary( BoundaryLocation::East, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::South, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::Top, id, number_of_nodes_on_level_zero ) );
      case BoundaryLocation::EastSouthBottom:
         return ( IsNaturalExternalBoundary( BoundaryLocation::East, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::South, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::Bottom, id, number_of_nodes_on_level_zero ) );

      case BoundaryLocation::WestNorthTop:
         return ( IsNaturalExternalBoundary( BoundaryLocation::West, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::North, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::Top, id, number_of_nodes_on_level_zero ) );
      case BoundaryLocation::WestNorthBottom:
         return ( IsNaturalExternalBoundary( BoundaryLocation::West, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::North, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::Bottom, id, number_of_nodes_on_level_zero ) );
      case BoundaryLocation::WestSouthTop:
         return ( IsNaturalExternalBoundary( BoundaryLocation::West, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::South, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::Top, id, number_of_nodes_on_level_zero ) );
#ifndef PERFORMANCE
      case BoundaryLocation::WestSouthBottom:
         return ( IsNaturalExternalBoundary( BoundaryLocation::West, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::South, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::Bottom, id, number_of_nodes_on_level_zero ) );
      default:
         throw std::invalid_argument( "Boundary Location not Found in GetNeighborId() - Impossible Error" );
#else
      default: /* BoundaryLocation::WestSouthBottom */
         return ( IsNaturalExternalBoundary( BoundaryLocation::West, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::South, id, number_of_nodes_on_level_zero ) ||
                  IsNaturalExternalBoundary( BoundaryLocation::Bottom, id, number_of_nodes_on_level_zero ) );
#endif
   }
}

/**
 * @brief Gives the equivalaent morton index of the id.
 * @note Morton indices are only unique within one level.
 */
nid_t RawMortonIndexOfId( nid_t const id ) {
   return CutHeadBit( id, LevelOfNode( id ) );
}
