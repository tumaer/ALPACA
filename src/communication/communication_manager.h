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
#ifndef COMMUNICATION_MANAGER_H
#define COMMUNICATION_MANAGER_H

#include <vector>
#include <numeric>
#include "topology/tree.h"
#include "topology/topology_manager.h"
#include "communication/communication_types.h"
#include "communication/exchange_types.h"
#include "internal_boundary_types.h"

/**
 * @brief The CommunicationManager class provides the functionality for communicating data between nodes and ranks. Furthermore, it holds the neighbor relations
 *        between nodes and external relations for external node boundaries.
 */
class CommunicationManager : public CommunicationTypes {
   TopologyManager& topology_;
   unsigned int const maximum_level_;
   int const my_rank_id_;
   int const mpi_tag_ub_;
   std::vector<unsigned int> partner_tag_map_;

   //Cache for Halo Update Pattern
   std::vector<std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>>> internal_boundaries_;
   std::vector<std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>>> internal_boundaries_mpi_;
   std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> internal_multi_boundaries_;
   std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> internal_multi_boundaries_mpi_;
   std::vector<std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>>> internal_boundaries_jump_;
   std::vector<std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>>> internal_boundaries_jump_mpi_;
   std::vector<std::vector<std::tuple<std::uint64_t, BoundaryLocation>>> external_boundaries_;
   std::vector<std::tuple<std::uint64_t, BoundaryLocation>> external_multi_boundaries_;

   // Vector holding flags dor each level that the lists have been created successfully
   std::vector<bool> boundaries_valid_;

   // Three values for three dimensions, even if only one dimension is simulated
   std::vector<std::array<unsigned int, 3>> jump_send_count_; // [0]: Plane, [1]: Stick, [2]: cube

   // Function that gives all neighbor-location and external-location relations for a given global node
   void NeighborsOfNode( std::uint64_t const global_id, std::vector<std::tuple<std::uint64_t, BoundaryLocation>>& nodes_internal_boundaries, std::vector<std::tuple<std::uint64_t, BoundaryLocation>>& external_boundaries );

public:
   CommunicationManager() = delete;
   explicit CommunicationManager( TopologyManager& topology, unsigned int const maximum_level );
   ~CommunicationManager() = default;
   CommunicationManager( CommunicationManager const& ) = delete;
   CommunicationManager& operator=( CommunicationManager const& ) = delete;
   CommunicationManager( CommunicationManager&& ) = delete;
   CommunicationManager& operator=( CommunicationManager&& ) = delete;

   // Function to fill the lists holding the relation to neighbor nodes and external boundaries
   void GenerateNeighborRelationForHaloUpdate( unsigned int const level );

   // return functions for the relation lists
   std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> const& InternalBoundariesJumpMpi( unsigned level ) const;
   std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> const& InternalBoundariesJump( unsigned level ) const;
   std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> const& InternalBoundariesMpi( unsigned level ) const;
   std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> const& InternalBoundaries( unsigned level ) const;
   std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> const& InternalMultiBoundariesMpi() const;
   std::vector<std::tuple<std::uint64_t, BoundaryLocation, InternalBoundaryType>> const& InternalMultiBoundaries() const;
   std::vector<std::tuple<std::uint64_t, BoundaryLocation>> const& ExternalBoundaries( unsigned level ) const;
   std::vector<std::tuple<std::uint64_t, BoundaryLocation>> const& ExternalMultiBoundaries() const;

   // Functions to get the status of the list creations and to empty the flags to regenerate the lists
   bool AreBoundariesValid( unsigned level ) const;
   void InvalidateCache();

   // Returns the counter for jump boundaries for the different exchange types
   unsigned int JumpSendCount( unsigned int const level, ExchangeType const type );

   // Send and receive function to buffer data between nodes and ranks
   int Send( void const* buffer, int const count, MPI_Datatype const datatype, int const destination_rank, std::vector<MPI_Request>& requests );
   int Recv( void *buf, int count, MPI_Datatype datatype, int source, std::vector<MPI_Request>& requests );

   // Helping functions to provide current rank and partner tags (MyRankId as member variable to avoid multiple calls of Mpi library)
   int TagForRank( unsigned int const partner );
   void ResetTagsForPartner();
   inline int MyRankId() { return my_rank_id_; }
};

#endif /* COMMUNICATION_MANAGER_H */
