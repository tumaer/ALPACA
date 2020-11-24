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
#include "topology/topology_manager.h"
#include "topology/tree.h"
#include "communication/mpi_utilities.h"
#include "topology/id_information.h"
#include "communication/communication_manager.h"
#include "communication/internal_halo_manager.h"
#include "multiresolution/averager.h"

namespace {
   // Some useful constants
   constexpr nid_t root_id = 0x1400000;
   auto const child_ids    = IdsOfChildren( root_id );
}// namespace

namespace {
   /**
    * @brief Fills all cells in a three dimensional array with the provided value.
    * @param buffer The array to be filled.
    * @param value The value to be set in every cell.
    * @tparam T value type.
    */
   template<typename T>
   void FillThreeDimensionalBufferWithValue( T ( &buffer )[CC::TCX()][CC::TCY()][CC::TCZ()], T value ) {
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               buffer[i][j][k] = value;
            }
         }
      }
   }

   /**
    * @brief Fills all cells in a three dimensional array with some arbitrary (not random!) values.
    * @param buffer The array to be filled.
    * @tparam T value type.
    */
   template<typename T>
   void FillThreeDimensionalBufferWithArbitratyValues( T ( &buffer )[CC::TCX()][CC::TCY()][CC::TCZ()] ) {
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               buffer[i][j][k] = -1.0 * double( i ) * double( i ) + -2.0 * double( j ) * double( j ) + -3.0 * double( k ) * double( k );
            }
         }
      }
   }

   /**
    * @brief Adds the given material into the first eight nodes on level one of the provided topology.
    * @param topology The topology whose nodes receive an additional fluid.
    * @param material The material top be added.
    * @note Undefined if topology does not have the first eight nodes on level one.
    */
   void AddMaterialToAllNodesInOneParentEightChildrenTopologyAndUpdate( TopologyManager& topology, MaterialName material ) {
      topology.AddMaterialToNode( root_id, material );
      for( auto const c : child_ids ) {
         topology.AddMaterialToNode( c, material );
      }
      topology.UpdateTopology();
   }

}// namespace

namespace {
   /**
    * @brief Inserts a given value into all equation right-hand-side buffers of the given node.
    * @param node Single-phase node whose buffers are to be filled.
    * @param value Value to poulate the buffers with.
    */
   void FillNodesRightHandSideEquationBuffersWithValue( Node& node, double const value ) {
      for( auto const eq : MF::ASOE() ) {
         FillThreeDimensionalBufferWithValue( node.GetSinglePhase().GetRightHandSideBuffer( eq ), value );
      }
   }

   /**
    * @brief Inserts arbitrary (not random!) values into all equation right-hand-side buffers of the given node.
    * @param node Single-phase node whose buffers are to be filled.
    */
   void FillNodesRightHandSideEquationBuffersWithArbitraryValues( Node& node ) {
      for( auto const eq : MF::ASOE() ) {
         FillThreeDimensionalBufferWithArbitratyValues( node.GetSinglePhase().GetRightHandSideBuffer( eq ) );
      }
   }

}// namespace

namespace {
   /**
    * @brief Fills the interface tag buffer of the given node with the provided value.
    * @param node The node whose interface are to be filled.
    * @param The interface tag value to put into every cell of the interface tag buffer.
    */
   void SetNodesInterfaceTags( Node& node, InterfaceTag const value ) {
      FillThreeDimensionalBufferWithValue( node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>(), ITTI( value ) );
   }

   /*
    * @brief Sets the internal corner-cell interface tags on the given node to the provided value.
    * @param node The node whose interface are to be filled.
    * @param The interface tag value to put into the corner cells.
    */
   void SetCornerTagsToValue( Node& node, InterfaceTag const value ) {
      auto& interface_tags                               = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
      interface_tags[CC::FICX()][CC::FICY()][CC::FICZ()] = ITTI( value );
      interface_tags[CC::LICX()][CC::FICY()][CC::FICZ()] = ITTI( value );
      interface_tags[CC::FICX()][CC::LICY()][CC::FICZ()] = ITTI( value );
      interface_tags[CC::LICX()][CC::LICY()][CC::FICZ()] = ITTI( value );
      interface_tags[CC::FICX()][CC::FICY()][CC::LICZ()] = ITTI( value );
      interface_tags[CC::LICX()][CC::FICY()][CC::LICZ()] = ITTI( value );
      interface_tags[CC::FICX()][CC::LICY()][CC::LICZ()] = ITTI( value );
      interface_tags[CC::LICX()][CC::LICY()][CC::LICZ()] = ITTI( value );
   }

   /**
    * @brief Verifies if the corner tags from all children have been averaged to the correct value and into the correct location in the provided node.
    * @param parent The node which has received the averaged tags.
    */
   void VerifyCornerTagsAveragedIntoParent( Node const& parent ) {
      auto parent_tags = parent.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
      REQUIRE( parent_tags[CC::FICX()][CC::FICY()][CC::FICZ()] == ITTI( IT::OldCutCell ) );
      REQUIRE( parent_tags[CC::FICX() + CC::ICX() / 2 - 1][CC::FICY()][CC::FICZ()] == ITTI( IT::OldCutCell ) );
      REQUIRE( parent_tags[CC::FICX()][CC::FICY() + CC::ICX() / 2 - 1][CC::FICZ()] == ITTI( IT::OldCutCell ) );
      REQUIRE( parent_tags[CC::FICX() + CC::ICX() / 2 - 1][CC::FICY() + CC::ICX() / 2 - 1][CC::FICZ()] == ITTI( IT::OldCutCell ) );
      REQUIRE( parent_tags[CC::FICX()][CC::FICY()][CC::FICZ() + CC::ICX() / 2 - 1] == ITTI( IT::OldCutCell ) );
      REQUIRE( parent_tags[CC::FICX() + CC::ICX() / 2 - 1][CC::FICY()][CC::FICZ() + CC::ICX() / 2 - 1] == ITTI( IT::OldCutCell ) );
      REQUIRE( parent_tags[CC::FICX()][CC::FICY() + CC::ICX() / 2 - 1][CC::FICZ() + CC::ICX() / 2 - 1] == ITTI( IT::OldCutCell ) );
      REQUIRE( parent_tags[CC::FICX() + CC::ICX() / 2 - 1][CC::FICY() + CC::ICX() / 2 - 1][CC::FICZ() + CC::ICX() / 2 - 1] == ITTI( IT::OldCutCell ) );

      REQUIRE( parent_tags[CC::FICX() + 1][CC::FICY() + 1][CC::FICZ() + 1] != ITTI( IT::OldCutCell ) );

      REQUIRE( parent_tags[CC::FICX() + CC::ICX() / 2][CC::FICY() + CC::ICX() / 2][CC::FICZ() + CC::ICX() / 2] == ITTI( IT::OldCutCell ) );
      REQUIRE( parent_tags[CC::LICX()][CC::FICY() + CC::ICX() / 2][CC::FICZ() + CC::ICX() / 2] == ITTI( IT::OldCutCell ) );
      REQUIRE( parent_tags[CC::FICX() + CC::ICX() / 2][CC::LICY()][CC::FICZ() + CC::ICX() / 2] == ITTI( IT::OldCutCell ) );
      REQUIRE( parent_tags[CC::LICX()][CC::LICY()][CC::FICZ() + CC::ICX() / 2] == ITTI( IT::OldCutCell ) );
      REQUIRE( parent_tags[CC::FICX() + CC::ICX() / 2][CC::FICY() + CC::ICX() / 2][CC::LICZ()] == ITTI( IT::OldCutCell ) );
      REQUIRE( parent_tags[CC::LICX()][CC::FICY() + CC::ICX() / 2][CC::LICZ()] == ITTI( IT::OldCutCell ) );
      REQUIRE( parent_tags[CC::FICX() + CC::ICX() / 2][CC::LICY()][CC::LICZ()] == ITTI( IT::OldCutCell ) );
      REQUIRE( parent_tags[CC::LICX()][CC::LICY()][CC::LICZ()] == ITTI( IT::OldCutCell ) );

      REQUIRE( parent_tags[CC::FICX() + CC::ICX() / 2 + 1][CC::FICY() + CC::ICX() / 2 + 1][CC::FICZ() + CC::ICX() / 2 + 1] != ITTI( IT::OldCutCell ) );
   }

   /*
    * @brief Verifies that all interface tags in the provided node are set to the bulk phase tag.
    * @param node The node whose interface tag buffer is to be examined.
    */
   void VerifyHomogenousBulkTags( Node const& node ) {
      auto const& tags = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
      for( unsigned int i = CC::LICX(); i < CC::FICX(); ++i ) {
         for( unsigned int j = CC::LICY(); j < CC::FICY(); ++j ) {
            for( unsigned int k = CC::LICZ(); k < CC::FICZ(); ++k ) {
               REQUIRE( tags[i][j][k] == ITTI( IT::BulkPhase ) );
            }
         }
      }
   }
}// namespace

namespace {

   /**
    * @brief Gives the desired result of an averaging of the eight cells around the one with the given array indices.
    * @param buffer The array whose values are averaged.
    * @param x,y,z The indices in x-,y- and z-direction.
    * @tparam buffer value type.
    */
   template<typename T>
   T AverageAroundGivenIndices( T const ( &buffer )[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const x, unsigned int const y, unsigned int const z ) {
      return buffer[x][y][z] + buffer[x + 1][y][z] + buffer[x][y + 1][z] + buffer[x + 1][y + 1][z] + buffer[x][y][z + 1] + buffer[x + 1][y][z + 1] + buffer[x][y + 1][z + 1] + buffer[x + 1][y + 1][z + 1];
   }

   /**
    * @brief Verifies if the parent buffer contains the expected averaged values from its eight child buffers.
    * @param parent_values The buffer holding the averaged values.
    * @param child0_values,child1_values,child2_values,child3_values,child4_values,child5_values,child6_values,child7_values The buffers that were
    * averaged. The name indicates the respective child position among its siblings. In Morton indexing stating at zero for the x- then y- then z-axis.
    * @tparam buffer value type.
    */
   template<class T>
   void VerifyAveragedMaterialValues( T const ( &parent_values )[CC::TCX()][CC::TCY()][CC::TCZ()], T const ( &child0_values )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                      T const ( &child1_values )[CC::TCX()][CC::TCY()][CC::TCZ()], T const ( &child2_values )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                      T const ( &child3_values )[CC::TCX()][CC::TCY()][CC::TCZ()], T const ( &child4_values )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                      T const ( &child5_values )[CC::TCX()][CC::TCY()][CC::TCZ()], T const ( &child6_values )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                      T const ( &child7_values )[CC::TCX()][CC::TCY()][CC::TCZ()] ) {

      unsigned int const x_child_start = CC::FICX();
      unsigned int const y_child_start = CC::FICY();
      unsigned int const z_child_start = CC::FICZ();

      unsigned int x_child = x_child_start;
      for( unsigned x = CC::FICX(); x < CC::FICX() + CC::ICX() / 2; x++ ) {
         unsigned int y_child = y_child_start;
         for( unsigned y = CC::FICY(); y < CC::FICY() + CC::ICY() / 2; y++ ) {
            unsigned int z_child = z_child_start;
            for( unsigned z = CC::FICZ(); z < CC::FICZ() + CC::ICZ() / 2; z++ ) {
               REQUIRE( 8 * parent_values[x][y][z] == Approx( AverageAroundGivenIndices( child0_values, x_child, y_child, z_child ) ) );
               REQUIRE( 8 * parent_values[x + CC::ICX() / 2][y][z] == Approx( AverageAroundGivenIndices( child1_values, x_child, y_child, z_child ) ) );
               REQUIRE( 8 * parent_values[x][y + CC::ICY() / 2][z] == Approx( AverageAroundGivenIndices( child2_values, x_child, y_child, z_child ) ) );
               REQUIRE( 8 * parent_values[x + CC::ICX() / 2][y + CC::ICY() / 2][z] == Approx( AverageAroundGivenIndices( child3_values, x_child, y_child, z_child ) ) );
               REQUIRE( 8 * parent_values[x][y][z + CC::ICZ() / 2] == Approx( AverageAroundGivenIndices( child4_values, x_child, y_child, z_child ) ) );
               REQUIRE( 8 * parent_values[x + CC::ICX() / 2][y][z + CC::ICZ() / 2] == Approx( AverageAroundGivenIndices( child5_values, x_child, y_child, z_child ) ) );
               REQUIRE( 8 * parent_values[x][y + CC::ICY() / 2][z + CC::ICZ() / 2] == Approx( AverageAroundGivenIndices( child6_values, x_child, y_child, z_child ) ) );
               REQUIRE( 8 * parent_values[x + CC::ICX() / 2][y + CC::ICY() / 2][z + CC::ICZ() / 2] == Approx( AverageAroundGivenIndices( child7_values, x_child, y_child, z_child ) ) );
               z_child += 2;
            }
            y_child += 2;
         }
         x_child += 2;
      }
   }
}// namespace

SCENARIO( "Averaging material data yields correct results on a simple topology", "[1rank],[2rank]" ) {
   GIVEN( "A topology and a corresponding tree with one parent on level zero and its eight children holding just one material" ) {
      constexpr unsigned int maximum_level = 1;
      constexpr MaterialName material      = MaterialName::MaterialOne;
      int const my_rank                    = MpiUtilities::MyRankId();
      TopologyManager topology             = TopologyManager( { 1, 1, 1 }, maximum_level, 0 );
      topology.UpdateTopology();

      Tree tree = Tree( topology, maximum_level, 1.0 );
      if( topology.NodeIsOnRank( root_id, my_rank ) ) {
         tree.CreateNode( root_id, { material } );
      }

      topology.RefineNodeWithId( root_id );
      topology.UpdateTopology();
      AddMaterialToAllNodesInOneParentEightChildrenTopologyAndUpdate( topology, material );

      for( auto const id : topology.LocalLeafIds() ) {
         tree.CreateNode( id, { material } );
      }

      WHEN( "Material values are set differently between children and parent and averaged into parent" ) {
         if( topology.NodeIsOnRank( root_id, my_rank ) ) { FillNodesRightHandSideEquationBuffersWithValue( tree.GetNodeWithId( root_id ), 42.0 ); }
         if( topology.NodeIsOnRank( child_ids[1], my_rank ) ) { FillNodesRightHandSideEquationBuffersWithValue( tree.GetNodeWithId( child_ids[1] ), 1.0 ); }
         if( topology.NodeIsOnRank( child_ids[2], my_rank ) ) { FillNodesRightHandSideEquationBuffersWithValue( tree.GetNodeWithId( child_ids[2] ), 2.0 ); }
         if( topology.NodeIsOnRank( child_ids[7], my_rank ) ) { FillNodesRightHandSideEquationBuffersWithValue( tree.GetNodeWithId( child_ids[7] ), -3.0 ); }
         if( topology.NodeIsOnRank( child_ids[5], my_rank ) ) { FillNodesRightHandSideEquationBuffersWithArbitraryValues( tree.GetNodeWithId( child_ids[5] ) ); }

         CommunicationManager communicator( topology, maximum_level );
         Averager averager( topology, communicator, tree );
         averager.AverageMaterial( { maximum_level } );

         THEN( "Material averages in the parent are according to child values" ) {
            double empty_buffer[CC::TCX()][CC::TCY()][CC::TCZ()];
            double one_buffer[CC::TCX()][CC::TCY()][CC::TCZ()];
            double two_buffer[CC::TCX()][CC::TCY()][CC::TCZ()];
            double negative_three_buffer[CC::TCX()][CC::TCY()][CC::TCZ()];
            double arbitrary_buffer[CC::TCX()][CC::TCY()][CC::TCZ()];
            FillThreeDimensionalBufferWithValue( empty_buffer, 0.0 );
            FillThreeDimensionalBufferWithValue( one_buffer, 1.0 );
            FillThreeDimensionalBufferWithValue( two_buffer, 2.0 );
            FillThreeDimensionalBufferWithValue( negative_three_buffer, -3.0 );
            FillThreeDimensionalBufferWithArbitratyValues( arbitrary_buffer );
            if( topology.NodeIsOnRank( root_id, my_rank ) ) {
               for( auto const eq : MF::ASOE() ) {
                  VerifyAveragedMaterialValues( tree.GetNodeWithId( root_id ).GetSinglePhase().GetRightHandSideBuffer( eq ), empty_buffer, one_buffer,
                                                two_buffer, empty_buffer, empty_buffer, arbitrary_buffer, empty_buffer, negative_three_buffer );
               }
            }
         }
      }

      WHEN( "Corner interfaces tags in children are cut cells, the rest bulk and an average is created" ) {
         if( topology.NodeIsOnRank( root_id, my_rank ) ) { SetNodesInterfaceTags( tree.GetNodeWithId( root_id ), IT::BulkPhase ); }
         if( topology.NodeIsOnRank( child_ids[0], my_rank ) ) { SetNodesInterfaceTags( tree.GetNodeWithId( child_ids[0] ), IT::BulkPhase ); }
         if( topology.NodeIsOnRank( child_ids[7], my_rank ) ) { SetNodesInterfaceTags( tree.GetNodeWithId( child_ids[7] ), IT::BulkPhase ); }
         if( topology.NodeIsOnRank( child_ids[0], my_rank ) ) { SetCornerTagsToValue( tree.GetNodeWithId( child_ids[0] ), IT::NewCutCell ); }
         if( topology.NodeIsOnRank( child_ids[7], my_rank ) ) { SetCornerTagsToValue( tree.GetNodeWithId( child_ids[7] ), IT::OldCutCell ); }

         CommunicationManager communicator( topology, maximum_level );
         Averager averager( topology, communicator, tree );

         THEN( "The parent interface tags remain all bulk given a single-phase topology when averaged" ) {
            averager.AverageInterfaceTags( { maximum_level } );
            if( topology.NodeIsOnRank( root_id, my_rank ) ) { VerifyHomogenousBulkTags( tree.GetNodeWithId( root_id ) ); }
         }

         THEN( "The cut cells are propagated down into the parent on a (faked) multi-phase topology when averaged" ) {
            AddMaterialToAllNodesInOneParentEightChildrenTopologyAndUpdate( topology, material );
            averager.AverageInterfaceTags( { maximum_level } );
            if( topology.NodeIsOnRank( root_id, my_rank ) ) { VerifyCornerTagsAveragedIntoParent( tree.GetNodeWithId( root_id ) ); }
         }
      }
   }
}