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
#ifndef FINEST_LEVEL_MESH_GENERATOR_H
#define FINEST_LEVEL_MESH_GENERATOR_H

#include <hdf5.h>
#include "input_output/output_writer/mesh_generator.h"

/**
 * @brief The StandardFinestLevelMeshGenerator generates a mesh for the output (currently xdmf + hdf5) using a virtual mapping of the mesh on the finest level. 
 * 
 * It provides the functionality of a non-ambiguous mesh (vertex IDs and coordinates) using the currently finest level as a 
 * reference. All blocks are virtually mapped onto the current maximum level. Therefore, the final hdf5 file is 
 * overloaded with vertex coordinates of the finest level, where a lot of are not used in the reading, e.g. using ParaView.
 * Nevertheless, the final mesh represents the current multi-resolution situation including jumps between different blocks. Only 
 * leaf nodes are written.  
 * 
 * Example for specification of vertex IDs:
 * 
 * Actual mesh situation:                         Virtual mapping on finest level 
 * 
 *  4____5________6_______________7                8____9____10___11___12___13___14___15
 *  |____|________|_______________|                |____|____|____|____|____|____|____|
 *  0    1        2               3                0    1    2    3    4    5    6    7  
 * 
 * Final Output: 
 * 
 *  8____9________11______________15               
 *  |____|________|_______________|               
 *  0    1        3               7               
 * 
 * Coordinates written to hdf5: all from 0 to 15, even only 8 are used, to provide direct mapping from vertex IDs to coordinates (required for reading)
 * 
 */
class StandardFinestLevelMeshGenerator : public MeshGenerator {

   // Variable specification from the base class
   using MeshGenerator::topology_;
   using MeshGenerator::tree_;
   using MeshGenerator::dimensionalized_node_size_on_level_zero_;

   // Self defined variables 
   std::array<unsigned int, 3> const number_of_nodes_on_level_zero_;

   // local functions required for the computation
   std::array<unsigned long long int, 3> GlobalNumberOfVerticesPerDimension() const;
   unsigned long long int TotalNumberOfVerticesOnFinestLevel() const;
   std::array<unsigned long long int, 2> GetLocalNumberOfVerticesAndStartIndex() const;

   // virtual functions required from the base class to compute data to hdf5 file
   void DoComputeVertexIDs( std::vector<unsigned long long int> & vertex_ids ) const override;
   void DoComputeVertexCoordinates( std::vector<double> & vertex_coordinates ) const override;   

   // virtual dimension functions required from base class 
   std::vector<std::reference_wrapper<Node const>> DoGetLocalNodes() const override;
   hsize_t DoGetGlobalNumberOfCells() const override;
   hsize_t DoGetLocalNumberOfCells() const override;
   hsize_t DoGetLocalCellsStartIndex() const override;
   std::vector<hsize_t> DoGetGlobalDimensionsOfVertexCoordinates() const override;
   std::vector<hsize_t> DoGetLocalDimensionsOfVertexCoordinates() const override;
   hsize_t DoGetLocalVertexCoordinatesStartIndex() const override;

public:
   StandardFinestLevelMeshGenerator() = delete;
   explicit StandardFinestLevelMeshGenerator( TopologyManager const& topology_manager, 
                                              Tree const& flower, 
                                              double const dimensionalized_node_size_on_level_zero, 
                                              std::array<unsigned int, 3> const number_of_nodes_on_level_zero );
   virtual ~StandardFinestLevelMeshGenerator() = default;
   StandardFinestLevelMeshGenerator( StandardFinestLevelMeshGenerator const& ) = delete;
   StandardFinestLevelMeshGenerator& operator=( StandardFinestLevelMeshGenerator const& ) = delete;
   StandardFinestLevelMeshGenerator( StandardFinestLevelMeshGenerator&& ) = delete;
   StandardFinestLevelMeshGenerator& operator=( StandardFinestLevelMeshGenerator&& ) = delete;
};

#endif // FINEST_LEVEL_MESH_GENERATOR_H