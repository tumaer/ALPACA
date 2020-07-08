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
#ifndef MPI_UTILITIES_H
#define MPI_UTILITIES_H

#include <mpi.h>
#include <numeric>
#include <vector>

namespace MpiUtilities {

   /**
    * @brief Reduces a bool across MPI ranks.
    * @param input local bool.
    * @param operation The redcution operation. Default: or.
    * @return Global results.
    */
   inline bool GloballyReducedBool( bool const input, MPI_Op const operation = MPI_LOR ) {
      bool result = input;
      MPI_Allreduce( MPI_IN_PLACE, &result, 1, MPI_CXX_BOOL, operation, MPI_COMM_WORLD );
      return result;
   }

   /**
    * @brief Gives the rank id from MPI directly as int. Avoids handle creation, e.g. for const members in initializer list.
    * @return Rank id.
    */
   inline int MyRankId() {
      int rank_id = -1;
      MPI_Comm_rank( MPI_COMM_WORLD, &rank_id );
      return rank_id;
   }

   /**
    * @brief Gives the number of ranks in the MPI communicator "MPI_COMM_WORLD". Avoids handle creation, e.g. for const members in initializer list.
    * @return Communicator Size which is the number of ranks.
    */
   inline int NumberOfRanks() {
      int communicator_size = -1;
      MPI_Comm_size( MPI_COMM_WORLD, &communicator_size );
      return communicator_size;
   }

   /**
   * @brief Gives the MPI_TAG_UB. Avoids handle creation, e.g. for const members in initializer list.
   * @return MPI_TAG_UB
   */
   inline int MpiTagUb() {
      int* tag_ub;
      int flag;
      MPI_Comm_get_attr( MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &flag );
      return *tag_ub;
   }

   /**
    * @brief Wrapper function to collect data from a local (= on one MPI rank) vector into a large global (= data of all ranks) one via a gatherv operation.
    * @param local_data The data present at this rank.
    * @param type The MPI datatype to be used in the gather calls.
    * @param number_of_ranks The number of ranks in the communicator.
    * @param global_data Vector holding the collected data from all ranks (indirect return parameter).
    * @tparam Data type.
    * @note Does not perform sanity checks. If the template type and the MPI datatype do not match the results will be corrupted.
    *       Uses MPI_COMM_WORLD as communicator.
    *       Overrides the provided global_data array.
    */
   template<class T>
   void LocalToGlobalData( std::vector<T> const& local_data, MPI_Datatype const type, int const number_of_ranks, std::vector<T>& global_data ) {
      int length = local_data.size();// Must be int due to MPI standard.
      std::vector<int> all_lengths( number_of_ranks );
      MPI_Allgather( &length, 1, MPI_INT, all_lengths.data(), 1, MPI_INT, MPI_COMM_WORLD );

      std::vector<int> offsets( number_of_ranks );
      int insert_key = 0;
      for( int i = 0; i < number_of_ranks; ++i ) {
         offsets[i] = insert_key;
         insert_key += all_lengths[i];
      }

      global_data.resize( std::accumulate( all_lengths.begin(), all_lengths.end(), 0 ) );
      MPI_Allgatherv( local_data.data(), length, type, global_data.data(), all_lengths.data(), offsets.data(), type, MPI_COMM_WORLD );
   }
}// namespace MpiUtilities

#endif// MPI_UTILITIES_H
