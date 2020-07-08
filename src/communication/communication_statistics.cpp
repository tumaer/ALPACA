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
#include <mpi.h>
#include "communication_statistics.h"
#include "mpi_utilities.h"

long CommunicationStatistics::no_jump_halos_recv_ = 0;
long CommunicationStatistics::no_jump_halos_send_ = 0;
long CommunicationStatistics::jump_halos_recv_    = 0;
long CommunicationStatistics::jump_halos_send_    = 0;
long CommunicationStatistics::balance_send_       = 0;
long CommunicationStatistics::balance_recv_       = 0;
long CommunicationStatistics::average_level_send_ = 0;
long CommunicationStatistics::average_level_recv_ = 0;

/**
 * @brief Sums the MPI statistic over all MPI ranks and gives a string respresentation of the result.
 * @return .
 */
std::string SummedCommunicationStatisticsString() {

   //Logging Stats
   std::string statistics;
   statistics.append( " Ranks: " + std::to_string( MpiUtilities::NumberOfRanks() ) + " | " );

   long global_statistic;
   MPI_Allreduce( &CommunicationStatistics::balance_send_, &global_statistic, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
   statistics.append( " Balance Send: " + std::to_string( global_statistic ) + " | " );

   MPI_Allreduce( &CommunicationStatistics::balance_recv_, &global_statistic, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
   statistics.append( " Balance Recv: " + std::to_string( global_statistic ) + " | " );

   MPI_Allreduce( &CommunicationStatistics::no_jump_halos_send_, &global_statistic, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
   statistics.append( " No Jump Halos Send: " + std::to_string( global_statistic ) + " | " );

   MPI_Allreduce( &CommunicationStatistics::no_jump_halos_recv_, &global_statistic, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
   statistics.append( " No Jump Halos Recv: " + std::to_string( global_statistic ) + " | " );

   MPI_Allreduce( &CommunicationStatistics::jump_halos_send_, &global_statistic, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
   statistics.append( " Jump Halos Send: " + std::to_string( global_statistic ) + " | " );

   MPI_Allreduce( &CommunicationStatistics::jump_halos_recv_, &global_statistic, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
   statistics.append( " Jump Halos Recv: " + std::to_string( global_statistic ) + " | " );

   MPI_Allreduce( &CommunicationStatistics::average_level_send_, &global_statistic, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
   statistics.append( " Proj.lvl-send: " + std::to_string( global_statistic ) + " | " );

   MPI_Allreduce( &CommunicationStatistics::average_level_recv_, &global_statistic, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
   statistics.append( " Proj.lvl-recv: " + std::to_string( global_statistic ) + " | " );
   return statistics;
}