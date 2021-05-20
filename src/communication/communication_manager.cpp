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
#include "communication_manager.h"

#include <bitset>

#include "topology/id_information.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "enums/interface_tag_definition.h"
#include "boundary_condition/material_boundary_condition.h"
#include "mpi_utilities.h"
#include "communication/communication_statistics.h"

/**
 * @brief Default constructor.
 * @param topology Instance that provides node data on a global level.
 * @param maximum_level Maximum present level for the simulation.
 */
CommunicationManager::CommunicationManager( TopologyManager& topology, unsigned int const maximum_level ) :                      // Start initializer list
                                                                                                            CommunicationTypes(),// For allocation of the MPI Datatypes
                                                                                                            topology_( topology ),
                                                                                                            maximum_level_( maximum_level ),
                                                                                                            my_rank_id_( MpiUtilities::MyRankId() ),
                                                                                                            mpi_tag_ub_( MpiUtilities::MpiTagUb() ),
                                                                                                            partner_tag_map_( MpiUtilities::NumberOfRanks(), 0 ),
                                                                                                            internal_boundaries_( maximum_level_ + 1 ),
                                                                                                            internal_boundaries_mpi_( maximum_level_ + 1 ),
                                                                                                            internal_multi_boundaries_(),
                                                                                                            internal_multi_boundaries_mpi_(),
                                                                                                            internal_boundaries_jump_( maximum_level_ + 1 ),
                                                                                                            internal_boundaries_jump_mpi_( maximum_level_ + 1 ),
                                                                                                            external_boundaries_( maximum_level_ + 1 ),
                                                                                                            external_multi_boundaries_(),
                                                                                                            boundaries_valid_( maximum_level_ + 1, false ) {
   // Initialize cache for Halo Update
   for( unsigned int level = 0; level <= maximum_level_; level++ ) {
      jump_send_count_.emplace_back( std::array<unsigned int, 3>( { 0, 0, 0 } ) );
   }
}

/**
 * @brief Counts the necessary amount of planes, sticks and cubes for jump sends.
 * @param loc Location necessary to determine type.
 * @param jump_send_counters [0] plane, [1] stick and [2] cube count.
 */
void IncrementJumpSendCounter( BoundaryLocation const loc, std::array<unsigned int, 3>& jump_send_counters ) {
   if( LTI( loc ) <= LTI( BoundaryLocation::Bottom ) ) {
      jump_send_counters[ETTI( ExchangeType::Plane )]++;//Plane
   } else if( LTI( loc ) <= LTI( BoundaryLocation::SouthWest ) ) {
      jump_send_counters[ETTI( ExchangeType::Stick )]++;//Stick
   } else {
      jump_send_counters[ETTI( ExchangeType::Cube )]++;//Cube
   }
}

/**
 * @brief Computes and caches the relations between ranks (who is whose communication partner).
 * @param level The level on which the relations shall be obtained and cached.
 */
void CommunicationManager::GenerateNeighborRelationForHaloUpdate( unsigned int const level ) {
   if( boundaries_valid_[level] )
      return;

   //clear existing cache
   internal_boundaries_[level].clear();
   internal_boundaries_mpi_[level].clear();
   internal_boundaries_jump_[level].clear();
   internal_boundaries_jump_mpi_[level].clear();
   external_boundaries_[level].clear();
   if( level == maximum_level_ ) {
      internal_multi_boundaries_.clear();
      internal_multi_boundaries_mpi_.clear();
      external_multi_boundaries_.clear();
   }

   // Declare temporary variables and reserve maximum possible space for vectors
   jump_send_count_[level] = { 0, 0, 0 };
   std::vector<std::tuple<nid_t, BoundaryLocation>> tmp_neighbor_location_vector;
   std::vector<std::tuple<nid_t, BoundaryLocation>> tmp_external_id_location_vector;
   tmp_neighbor_location_vector.reserve( 26 );
   tmp_external_id_location_vector.reserve( 6 );
   // Get List of Halos on this Level
   // All (globally existing) nodes are viewed by all ranks in order to have the same order on all ranks, so that the tagging system works properly.
   for( auto const global_id : topology_.IdsOnLevel( level ) ) {
      NeighborsOfNode( global_id, tmp_neighbor_location_vector, tmp_external_id_location_vector );
      if( topology_.NodeIsOnRank( global_id, my_rank_id_ ) ) {
         external_boundaries_[level].insert( external_boundaries_[level].end(), tmp_external_id_location_vector.begin(), tmp_external_id_location_vector.end() );
         if( level == maximum_level_ && topology_.IsNodeMultiPhase( global_id ) ) {
            external_multi_boundaries_.insert( external_multi_boundaries_.end(), tmp_external_id_location_vector.begin(), tmp_external_id_location_vector.end() );
         }
      }
      for( auto const& neighbor_id_location_element : tmp_neighbor_location_vector ) {
         //sort boundary type into vector
         if( !topology_.NodeExists( std::get<0>( neighbor_id_location_element ) ) ) {
            //jump boundary
            auto const parent_id = ParentIdOfNode( global_id );
            if( topology_.NodeIsOnRank( global_id, my_rank_id_ ) ) {
               // The node belongs to me -> I must Update it
               if( topology_.NodeIsOnRank( parent_id, my_rank_id_ ) ) {
                  // I have the child and the parent
                  internal_boundaries_jump_[level].push_back( std::make_tuple( global_id, std::get<1>( neighbor_id_location_element ), InternalBoundaryType::JumpBoundaryLocal ) );
               } else {
                  // I have the child, but not the parent
                  internal_boundaries_jump_mpi_[level].push_back(
                        std::make_tuple( global_id, std::get<1>( neighbor_id_location_element ), InternalBoundaryType::JumpBoundaryMpiRecv ) );
               }
            } else if( topology_.NodeIsOnRank( parent_id, my_rank_id_ ) ) {
               // I have the parent, but not the child -> I Send the needed Info
               internal_boundaries_jump_mpi_[level].push_back(
                     std::make_tuple( global_id, std::get<1>( neighbor_id_location_element ), InternalBoundaryType::JumpBoundaryMpiSend ) );
               IncrementJumpSendCounter( std::get<1>( neighbor_id_location_element ), jump_send_count_[level] );
            }
         } else {
            // No Jump
            bool const my_node     = topology_.NodeIsOnRank( global_id, my_rank_id_ );
            bool const my_neighbor = topology_.NodeIsOnRank( std::get<0>( neighbor_id_location_element ), my_rank_id_ );
            if( my_node ) {
               if( my_neighbor ) {
                  internal_boundaries_[level].push_back( std::make_tuple( global_id, std::get<1>( neighbor_id_location_element ), InternalBoundaryType::NoJumpBoundaryLocal ) );
               } else {
                  internal_boundaries_mpi_[level].push_back(
                        std::make_tuple( global_id, std::get<1>( neighbor_id_location_element ), InternalBoundaryType::NoJumpBoundaryMpiRecv ) );
               }
            } else {
               if( my_neighbor ) {
                  internal_boundaries_mpi_[level].push_back(
                        std::make_tuple( std::get<0>( neighbor_id_location_element ), OppositeDirection( std::get<1>( neighbor_id_location_element ) ),
                                         InternalBoundaryType::NoJumpBoundaryMpiSend ) );
               }
            }
            if( level == maximum_level_ && topology_.IsNodeMultiPhase( global_id ) && topology_.IsNodeMultiPhase( std::get<0>( neighbor_id_location_element ) ) ) {
               if( my_node ) {
                  if( my_neighbor ) {
                     internal_multi_boundaries_.push_back( std::make_tuple( global_id, std::get<1>( neighbor_id_location_element ), InternalBoundaryType::NoJumpBoundaryLocal ) );
                  } else {
                     internal_multi_boundaries_mpi_.push_back(
                           std::make_tuple( global_id, std::get<1>( neighbor_id_location_element ), InternalBoundaryType::NoJumpBoundaryMpiRecv ) );
                  }
               } else {
                  if( my_neighbor ) {
                     internal_multi_boundaries_mpi_.push_back(
                           std::make_tuple( std::get<0>( neighbor_id_location_element ), OppositeDirection( std::get<1>( neighbor_id_location_element ) ),
                                            InternalBoundaryType::NoJumpBoundaryMpiSend ) );
                  }
               }
            }
         }
      }
      tmp_neighbor_location_vector.clear();
      tmp_external_id_location_vector.clear();
   }
   boundaries_valid_[level] = true;
}

namespace {

   /**
 * Bitset used to convert natural boundaries into it's 9 Halo Boundary sides.
 * e.g. std::bitset<26> (CC::ANBL()[LTI(EAST)] gives all halo positions contained in east natural boundary side, e.g. east-south-bottom.
 * The bitset fills the tables with the configuration below with 0 or 1.
 * wsb, wst, wnb, wnt, esb, est, enb, ent, sw, se, nw, ne, tw, te, bw, be, ts, tn, bs, bn, b, t, s, n, w, e.
 * returns an unsigned int that has all bits set for a specific natural boundary location.
 */
   constexpr std::array<nid_t, 6> bitsets_for_natural_boundary_locations = {
         0x3d5401, //east
         0x3c2a802,//west
         0xccc144, //north
         0x3330288,//south
         0x1543310,//top
         0x2a80ce0,//bottom
   };
}// namespace

/**
 * @brief Finds the neighbors of a specific node and sorts them into the array. Domain Boundaries are inserted only as natural (e,w,n,s,t,b), internals are inserted per direction e.g diagonal wnb.
 * @param global_id global node's id.
 * @param nodes_internal_boundaries output array for internal boundaries, all HaloBoundarySides are included that are not part of externals BC.
 * @param external_boundaries all external boundaries are included as natural Boundary Side e.g. east, west...
 */
void CommunicationManager::NeighborsOfNode( nid_t const global_id, std::vector<std::tuple<nid_t, BoundaryLocation>>& nodes_internal_boundaries, std::vector<std::tuple<nid_t, BoundaryLocation>>& external_boundaries ) {
   std::bitset<26> sides( 0 );
   for( BoundaryLocation loc : CC::ANBS() ) {
      if( topology_.IsExternalTopologyBoundary( loc, global_id ) ) {
         sides |= std::bitset<26>( bitsets_for_natural_boundary_locations[LTI( loc )] );
         external_boundaries.push_back( std::make_tuple( global_id, loc ) );
      }
   }
   for( BoundaryLocation loc : CC::HBS() ) {
      if( !sides.test( LTI( loc ) ) ) {
         nid_t neighbor_id = topology_.GetTopologyNeighborId( global_id, loc );
         nodes_internal_boundaries.push_back( std::make_tuple( neighbor_id, loc ) );
      }
   }
}

/**
 * @brief Tells the communication manager that there was a change in the topology and it needs to generate the material boundaries from scratch.
 */
void CommunicationManager::InvalidateCache() {
   for( unsigned i = 0; i < boundaries_valid_.size(); i++ ) {
      boundaries_valid_[i] = false;
   }
}

/**
 * @brief Wrapper for MPI_Send or MPI_Isend, use like MPI_Send.
 * @param buffer initial address of send buffer (choice).
 * @param count number of elements in send buffer (integer).
 * @param datatype datatype of each send buffer element (handle).
 * @param destination_rank rank of destination (integer).
 * @param requests vector of communication request (handle), new handle will be added at the end of the vector.
 * @return Error value see MPI_Isend for details.
 * @note Must not be called in single-core mode (will hang).
 */
int CommunicationManager::Send( void const* buffer, int const count, MPI_Datatype const datatype, int const destination_rank, std::vector<MPI_Request>& requests ) {
   int const tag = TagForRank( destination_rank );
   requests.push_back( MPI_Request() );
   int const status = MPI_Isend( buffer, count, datatype, destination_rank, tag, MPI_COMM_WORLD, &requests.back() );
#ifndef PERFORMANCE
   if( status != MPI_SUCCESS ) {
      throw std::logic_error( "Send error" );
   }
#endif
   return status;
}

/**
 * @brief Wrapper for MPI_Recv or MPI_Irecv, use like MPI_Recv.
 * @param buffer initial address of send buffer (choice).
 * @param count number of elements in send buffer (integer).
 * @param datatype datatype of each send buffer element (handle).
 * @param source_rank MPI rank of source.
 * @param requests vector of communication request (handle), new handle will be added at the end of the vector.
 * @return Error value see MPI_Irecv for details.
 * @note Must not be called in single-core mode (will hang).
 */
int CommunicationManager::Recv( void* buffer, int const count, MPI_Datatype const datatype, int const source_rank, std::vector<MPI_Request>& requests ) {
   int const tag = TagForRank( source_rank );
   requests.push_back( MPI_Request() );
   int const status = MPI_Irecv( buffer, count, datatype, source_rank, tag, MPI_COMM_WORLD, &requests.back() );
#ifndef PERFORMANCE
   if( status != MPI_SUCCESS ) {
      throw std::logic_error( "Receive error" );
   }
#endif
   return status;
}

/**
 * @brief Returns a tag, that can be used for one communication with another rank. TagForRank must be called in the same order on all participating ranks. If replacing the communication with MPI_Send/MPI_Recv is valid, it's OK.
 * @param rank of the communication partner.
 * @return tag for communication.
 */
int CommunicationManager::TagForRank( unsigned int const rank ) {
   partner_tag_map_[rank] = ( partner_tag_map_[rank]++ ) % mpi_tag_ub_;
   return partner_tag_map_[rank];
}

/**
 * @brief Resets the tags used for communication between ranks to reduce tag overflows over the mpi tag upper bound.
 */
void CommunicationManager::ResetTagsForPartner() {
   std::fill( std::begin( partner_tag_map_ ), std::end( partner_tag_map_ ), 0 );
}

/**
 * @brief Gives a reference to the list for a given level holding all internal jump boundary relations that require mpi communication.
 * @param level Level for which the list should be returned.
 * @return List with relations for all internal mpi jump boundaries.
 */
std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const& CommunicationManager::InternalBoundariesJumpMpi( unsigned int const level ) const {
   return internal_boundaries_jump_mpi_[level];
}

/**
 * @brief Gives a reference to the list for a given level holding all internal jump boundary relations that do not require mpi communication.
 * @param level Level for which the list should be returned.
 * @return List with relations for all internal no-mpi jump boundaries.
 */
std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const& CommunicationManager::InternalBoundariesJump( unsigned int const level ) const {
   return internal_boundaries_jump_[level];
}

/**
 * @brief Gives a reference to the list for a given level holding all internal non-jump boundary relations that require mpi communication.
 * @param level Level for which the list should be returned.
 * @return List with relations for all internal mpi non-jump boundaries.
 */
std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const& CommunicationManager::InternalBoundariesMpi( unsigned int const level ) const {
   return internal_boundaries_mpi_[level];
}

/**
 * @brief Gives a reference to the list for a given level holding all internal non-jump boundary relations that do not require mpi communication.
 * @param level Level for which the list should be returned.
 * @return List with relations for all internal no-mpi non-jump boundaries.
 */
std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const& CommunicationManager::InternalBoundaries( unsigned int const level ) const {
   return internal_boundaries_[level];
}

/**
 * @brief Gives a reference to the list of multi-material nodes holding all internal boundary relations that require mpi communication.
 * @return List with relations for all internal mpi multi-material boundaries.
 */
std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const& CommunicationManager::InternalMultiBoundariesMpi() const {
   return internal_multi_boundaries_mpi_;
}

/**
 * @brief Gives a reference to the list of multi-material nodes holding all internal boundary relations that do not require mpi communication.
 * @return List with relations for all internal no-mpi multi-material boundaries.
 */
std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const& CommunicationManager::InternalMultiBoundaries() const {
   return internal_multi_boundaries_;
}

/**
 * @brief Gives a reference to the list for a given level holding all external boundary relations.
 * @param level Level for which the list should be returned.
 * @return List with relations for all external boundary nodes.
 */
std::vector<std::tuple<nid_t, BoundaryLocation>> const& CommunicationManager::ExternalBoundaries( unsigned int const level ) const {
   return external_boundaries_[level];
}

/**
 * @brief Gives a reference to the list of multi-material nodes on the maximum level holding all external boundary relations.
 * @return List with relations for all external boundary nodes.
 */
std::vector<std::tuple<nid_t, BoundaryLocation>> const& CommunicationManager::ExternalMultiBoundaries() const {
   return external_multi_boundaries_;
}

/**
 * @brief Gives a reference to the list for a given level holding all internal jump boundary relations that require mpi communication.
 * @param level Level for which the list should be returned.
 * @return List with internal mpi jump boundaries.
 */
bool CommunicationManager::AreBoundariesValid( unsigned const level ) const {
   return boundaries_valid_[level];
}

/**
 * @brief Gives a reference to the list for a given level holding all internal jump boundaries that require mpi communication.
 * @param level Level for which the list should be returned.
 * @param type Type of the exchange buffer that is used.
 * @return List with internal mpi jump boundaries.
 */
unsigned CommunicationManager::JumpSendCount( unsigned const level, ExchangeType const type ) {
   return jump_send_count_[level][ETTI( type )];
}
