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

#include <catch.hpp>
#include "topology/space_filling_curve_index.h"
#include <cstdint>
#include <iterator>
#include <vector>
#include "topology/node_id_type.h"
#include "topology/id_information.h"
#include "utilities/container_operations.h"

namespace {
   static nid_t const origin_id = IdSeed();

   /**
    * @brief Gives the node ids that make up a 2x2x2 cube on level zero.
    */
   std::vector<nid_t> TwoByTwoByTwoCube() {
      return { origin_id, EastNeighborOfNodeWithId( origin_id ), NorthNeighborOfNodeWithId( origin_id ),
               EastNeighborOfNodeWithId( NorthNeighborOfNodeWithId( origin_id ) ), TopNeighborOfNodeWithId( origin_id ),
               TopNeighborOfNodeWithId( EastNeighborOfNodeWithId( origin_id ) ), TopNeighborOfNodeWithId( NorthNeighborOfNodeWithId( origin_id ) ),
               TopNeighborOfNodeWithId( NorthNeighborOfNodeWithId( EastNeighborOfNodeWithId( origin_id ) ) ) };
   }

   /**
    * @brief Gives the node ids that make up a 4x2x1 channel on level zero.
    */
   std::vector<nid_t> FourByTwoByOneChannel() {
      return { origin_id, EastNeighborOfNodeWithId( origin_id ), EastNeighborOfNodeWithId( EastNeighborOfNodeWithId( origin_id ) ), EastNeighborOfNodeWithId( EastNeighborOfNodeWithId( EastNeighborOfNodeWithId( origin_id ) ) ),
               NorthNeighborOfNodeWithId( origin_id ), NorthNeighborOfNodeWithId( EastNeighborOfNodeWithId( origin_id ) ), NorthNeighborOfNodeWithId( EastNeighborOfNodeWithId( EastNeighborOfNodeWithId( origin_id ) ) ), NorthNeighborOfNodeWithId( EastNeighborOfNodeWithId( EastNeighborOfNodeWithId( EastNeighborOfNodeWithId( origin_id ) ) ) ) };
   }

   /*
    * @brief Gives the node ids that make up a 1x1x512 channel on level two.
    */
   std::vector<nid_t> LevelTwoChannelInTopBottomOrientation() {
      nid_t const start_on_level_two = IdsOfChildren( IdsOfChildren( origin_id ).front() ).front();
      std::vector<nid_t> nodes( 512, 0 );
      nid_t temp = start_on_level_two;
      for( auto& n : nodes ) {
         n    = temp;
         temp = TopNeighborOfNodeWithId( temp );
      }
      return nodes;
   }
}// namespace

SCENARIO( "Hilbert index is correctly computed", "[1rank]" ) {
   GIVEN( "Ids forming a 2x2x2 cube" ) {
      auto cube = TwoByTwoByTwoCube();
      WHEN( "We compute the index of each id in the cube" ) {
         std::vector<sfcidx_t> indices;
         indices.reserve( cube.size() );
         std::transform( std::cbegin( cube ), std::cend( cube ), std::back_inserter( indices ), []( auto const id ) { return HilbertIndex( id ); } );
         THEN( "The indices match the expected hilbert traversal thorugh the cube" ) {
            std::vector<sfcidx_t> expected_hilbert_order = { 0, 3, 7, 4, 1, 2, 6, 5 };
            REQUIRE( indices == expected_hilbert_order );
         }
         THEN( "The indices are unique" ) {
            REQUIRE( ContainerOperations::HoldsUniqueElements( indices ) );
         }
      }
   }
   GIVEN( "Ids forming 4x2x1 " ) {
      auto channel = FourByTwoByOneChannel();
      WHEN( "We compute the index of each id in the channel" ) {
         std::vector<sfcidx_t> indices;
         indices.reserve( channel.size() );
         std::transform( std::cbegin( channel ), std::cend( channel ), std::back_inserter( indices ), []( auto const id ) { return HilbertIndex( id ); } );
         THEN( "The indices match the expected hilbert traversal thorugh the cube" ) {
            std::vector<sfcidx_t> expected_hilbert_order = { 0, 3, 60, 63, 7, 4, 59, 56 };
            REQUIRE( indices == expected_hilbert_order );
         }
         THEN( "The indices are unique" ) {
            REQUIRE( ContainerOperations::HoldsUniqueElements( indices ) );
         }
      }
   }
   GIVEN( "A channel in Z direction on level two" ) {
      auto ltwo_channel_z = LevelTwoChannelInTopBottomOrientation();
      WHEN( "We compute the index of each id in the channel" ) {
         std::vector<sfcidx_t> indices;
         indices.reserve( ltwo_channel_z.size() );
         std::transform( std::cbegin( ltwo_channel_z ), std::cend( ltwo_channel_z ), std::back_inserter( indices ), []( auto const id ) { return HilbertIndex( id ); } );
         THEN( "The indices match the expected hilbert traversal thorugh the cube" ) {
            std::vector<sfcidx_t> expected_hilbert_order = { 0, 1, 26, 25, 486, 485, 510, 511, 512, 519,
                                                             520, 523, 724, 723, 716, 715, 13620, 13621, 13578, 13579,
                                                             13556, 13557, 13514, 13515, 13108, 13109, 13066, 13067, 13044, 13045,
                                                             13002, 13003, 249140, 249141, 249098, 249099, 249076, 249077, 249034, 249035,
                                                             248628, 248629, 248586, 248587, 248564, 248565, 248522, 248523, 261428, 261427,
                                                             261420, 261419, 261620, 261623, 261624, 261631, 261632, 261633, 261658, 261657,
                                                             262118, 262117, 262142, 262143, 262144, 262151, 262152, 262155, 262356, 262355,
                                                             262348, 262347, 266036, 266035, 266028, 266027, 266228, 266231, 266232, 266239,
                                                             266240, 266243, 266300, 266303, 266304, 266305, 266330, 266329, 267942, 267937,
                                                             267934, 267929, 267878, 267873, 267870, 267865, 371110, 371111, 371112, 371115,
                                                             370772, 370775, 370776, 370777, 370598, 370599, 370600, 370603, 370260, 370263,
                                                             370264, 370265, 367014, 367015, 367016, 367019, 366676, 366679, 366680, 366681,
                                                             366502, 366503, 366504, 366507, 366164, 366167, 366168, 366169, 6973862, 6973861,
                                                             6973886, 6973887, 6973888, 6973891, 6973948, 6973951, 6973952, 6973959, 6973960, 6973963,
                                                             6974164, 6974163, 6974156, 6974155, 6952244, 6952243, 6952236, 6952235, 6952436, 6952439,
                                                             6952440, 6952447, 6952448, 6952451, 6952508, 6952511, 6952512, 6952513, 6952538, 6952537,
                                                             6941094, 6941093, 6941118, 6941119, 6941120, 6941123, 6941180, 6941183, 6941184, 6941191,
                                                             6941192, 6941195, 6941396, 6941395, 6941388, 6941387, 6919476, 6919475, 6919468, 6919467,
                                                             6919668, 6919671, 6919672, 6919679, 6919680, 6919683, 6919740, 6919743, 6919744, 6919745,
                                                             6919770, 6919769, 6711718, 6711717, 6711742, 6711743, 6711744, 6711747, 6711804, 6711807,
                                                             6711808, 6711815, 6711816, 6711819, 6712020, 6712019, 6712012, 6712011, 6690100, 6690099,
                                                             6690092, 6690091, 6690292, 6690295, 6690296, 6690303, 6690304, 6690307, 6690364, 6690367,
                                                             6690368, 6690369, 6690394, 6690393, 6678950, 6678949, 6678974, 6678975, 6678976, 6678979,
                                                             6679036, 6679039, 6679040, 6679047, 6679048, 6679051, 6679252, 6679251, 6679244, 6679243,
                                                             6657332, 6657331, 6657324, 6657323, 6657524, 6657527, 6657528, 6657535, 6657536, 6657539,
                                                             6657596, 6657599, 6657600, 6657601, 6657626, 6657625, 127560102, 127560101, 127560126, 127560127,
                                                             127560128, 127560131, 127560188, 127560191, 127560192, 127560199, 127560200, 127560203, 127560404, 127560403,
                                                             127560396, 127560395, 127538484, 127538483, 127538476, 127538475, 127538676, 127538679, 127538680, 127538687,
                                                             127538688, 127538691, 127538748, 127538751, 127538752, 127538753, 127538778, 127538777, 127527334, 127527333,
                                                             127527358, 127527359, 127527360, 127527363, 127527420, 127527423, 127527424, 127527431, 127527432, 127527435,
                                                             127527636, 127527635, 127527628, 127527627, 127505716, 127505715, 127505708, 127505707, 127505908, 127505911,
                                                             127505912, 127505919, 127505920, 127505923, 127505980, 127505983, 127505984, 127505985, 127506010, 127506009,
                                                             127297958, 127297957, 127297982, 127297983, 127297984, 127297987, 127298044, 127298047, 127298048, 127298055,
                                                             127298056, 127298059, 127298260, 127298259, 127298252, 127298251, 127276340, 127276339, 127276332, 127276331,
                                                             127276532, 127276535, 127276536, 127276543, 127276544, 127276547, 127276604, 127276607, 127276608, 127276609,
                                                             127276634, 127276633, 127265190, 127265189, 127265214, 127265215, 127265216, 127265219, 127265276, 127265279,
                                                             127265280, 127265287, 127265288, 127265291, 127265492, 127265491, 127265484, 127265483, 127243572, 127243571,
                                                             127243564, 127243563, 127243764, 127243767, 127243768, 127243775, 127243776, 127243779, 127243836, 127243839,
                                                             127243840, 127243841, 127243866, 127243865, 133851558, 133851559, 133851560, 133851563, 133851220, 133851223,
                                                             133851224, 133851225, 133851046, 133851047, 133851048, 133851051, 133850708, 133850711, 133850712, 133850713,
                                                             133847462, 133847463, 133847464, 133847467, 133847124, 133847127, 133847128, 133847129, 133846950, 133846951,
                                                             133846952, 133846955, 133846612, 133846615, 133846616, 133846617, 133949862, 133949857, 133949854, 133949849,
                                                             133949798, 133949793, 133949790, 133949785, 133951398, 133951397, 133951422, 133951423, 133951424, 133951427,
                                                             133951484, 133951487, 133951488, 133951495, 133951496, 133951499, 133951700, 133951699, 133951692, 133951691,
                                                             133955380, 133955379, 133955372, 133955371, 133955572, 133955575, 133955576, 133955583, 133955584, 133955585,
                                                             133955610, 133955609, 133956070, 133956069, 133956094, 133956095, 133956096, 133956103, 133956104, 133956107,
                                                             133956308, 133956307, 133956300, 133956299, 133969204, 133969205, 133969162, 133969163, 133969140, 133969141,
                                                             133969098, 133969099, 133968692, 133968693, 133968650, 133968651, 133968628, 133968629, 133968586, 133968587,
                                                             134204724, 134204725, 134204682, 134204683, 134204660, 134204661, 134204618, 134204619, 134204212, 134204213,
                                                             134204170, 134204171, 134204148, 134204149, 134204106, 134204107, 134217012, 134217011, 134217004, 134217003,
                                                             134217204, 134217207, 134217208, 134217215, 134217216, 134217217, 134217242, 134217241, 134217702, 134217701,
                                                             134217726, 134217727 };

            REQUIRE( indices == expected_hilbert_order );
         }
         THEN( "The indices are unique" ) {
            REQUIRE( ContainerOperations::HoldsUniqueElements( indices ) );
         }
      }
   }
}

SCENARIO( "Lebesgue index is correctly computed", "[1rank]" ) {
   GIVEN( "Ids forming a 2x2x2 cube" ) {
      auto cube = TwoByTwoByTwoCube();
      WHEN( "We compute the index of each id in the cube" ) {
         std::vector<sfcidx_t> indices;
         indices.reserve( cube.size() );
         std::transform( std::cbegin( cube ), std::cend( cube ), std::back_inserter( indices ), []( auto const id ) { return LebesgueIndex( id ); } );
         THEN( "The indices match the expected hilbert traversal thorugh the cube" ) {
            std::vector<sfcidx_t> expected_lebesgue_order = { 0, 1, 2, 3, 4, 5, 6, 7 };
            REQUIRE( indices == expected_lebesgue_order );
         }
         THEN( "The indices are unique" ) {
            REQUIRE( ContainerOperations::HoldsUniqueElements( indices ) );
         }
      }
   }
   GIVEN( "Ids forming 4x2x1 " ) {
      auto channel = FourByTwoByOneChannel();
      WHEN( "We compute the index of each id in the cube" ) {
         std::vector<sfcidx_t> indices;
         indices.reserve( channel.size() );
         std::transform( std::cbegin( channel ), std::cend( channel ), std::back_inserter( indices ), []( auto const id ) { return LebesgueIndex( id ); } );
         THEN( "The indices match the expected hilbert traversal thorugh the cube" ) {
            std::vector<sfcidx_t> expected_lebesgue_order = { 0, 1, 8, 9, 2, 3, 10, 11 };
            REQUIRE( indices == expected_lebesgue_order );
         }
         THEN( "The indices are unique" ) {
            REQUIRE( ContainerOperations::HoldsUniqueElements( indices ) );
         }
      }
   }
}