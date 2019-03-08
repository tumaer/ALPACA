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
#ifndef MESH_QUANTITY_H
#define MESH_QUANTITY_H

#include <hdf5.h>
#include "topology/node.h"
#include "topology/tree.h"
#include "topology/topology_manager.h"
#include "unit_handler.h"
#include "materials/material_definitions.h"

/**
 * @brief The MeshGenerator class handles the output to the filesystem in Xdmf+HDF5 file format for ParaView. MeshGenerator must not change any data.
 *
 * The mesh generator class serves as a base class to generate the mesh (vertex IDs and coordinates) of the current simulation. In general three different
 * mesh generations are distinguished (standard, interface and debug). In the standard mesh, the full MultiResolution mesh is written. In the interface mesh
 * only the interface nodes are written and in debug mode all nodes with all cells (internal + halos) are provided.
 */
class MeshGenerator {

private:
   // Naming of the vertex IDs and coordinates in the final file
   std::string const vertex_ids_name_ = "cell_vertex_IDs";
   std::string const vertex_coordinates_name_ = "cell_vertex_coordinates";

protected:
   // topology manager containing the global information of nodes
   TopologyManager const& topology_;
   // tree containing all local information of the nodes
   Tree const& tree_;
   // block size on level zero (already dimensionalized)
   double const dimensionalized_node_size_on_level_zero_;

   /**
    * @brief Appends the vertex IDs to the given vector
    * @param vertex_ids Vector for the vertex IDs (indirect return)
    * @note Pure virtual function that is required from derived class
    */
   virtual void DoComputeVertexIDs( std::vector<unsigned long long int> & vertex_ids ) const = 0;
   /**
    * @brief Appends the vertex coordinates to the given vector
    * @param vertex_coordinates Vector for the vertex coordinates (indirect return)
    * @note Pure virtual function that is required from derived class
    */
   virtual void DoComputeVertexCoordinates( std::vector<double> & vertex_coordinates ) const = 0;

   /**
    * @brief See public function for reference
    */
   virtual std::vector<std::reference_wrapper<Node const>> DoGetLocalNodes() const = 0;

   /**
    * @brief See public function for reference
    */
   virtual hsize_t DoGetGlobalNumberOfCells() const = 0;
   /**
    * @brief See public function for reference
    */
   virtual hsize_t DoGetLocalNumberOfCells() const = 0;
   /**
    * @brief See public function for reference
    */
   virtual hsize_t DoGetLocalCellsStartIndex() const = 0;

   /**
    * @brief See public function for reference
    */
   virtual std::vector<hsize_t> DoGetGlobalDimensionsOfVertexCoordinates() const = 0;
   /**
    * @brief See public function for reference
    */
   virtual std::vector<hsize_t> DoGetLocalDimensionsOfVertexCoordinates() const = 0;
   /**
    * @brief See public function for reference
    */
   virtual hsize_t DoGetLocalVertexCoordinatesStartIndex() const = 0;

   // Constructor can only be called from derived classes
   explicit MeshGenerator( TopologyManager const& topology_manager, Tree const& flower, double const dimensionalized_node_size_on_level_zero );

public:
   MeshGenerator() = delete;
   virtual ~MeshGenerator() = default;
   MeshGenerator( MeshGenerator const& ) = delete;
   MeshGenerator& operator=( MeshGenerator const& ) = delete;
   MeshGenerator( MeshGenerator&& ) = delete;
   MeshGenerator& operator=( MeshGenerator&& ) = delete;

   /**
    * @brief Returns the nodes to be filled for the given mesh.
    *
    * A mesh generator creates the output mesh and stores vertex IDs and coordinates in an appropriate order.
    * Therefore, all subsequent classes that write cell data needs to store the data in the correct order for the
    * correct alignment of cell data and cell vertices. The nodes contain all required information for the output (material fields and interface fields).
    *
    * @return vector with references to all nodes considered for the computed mesh. The nodes are in the order in which
    *         they should be filled by output quantities.
    */
   std::vector<std::reference_wrapper<Node const>> GetLocalNodes() const {
      return DoGetLocalNodes();
   }

   /**
    * @brief Gives the global number of all cells (on all ranks) for the given mesh generator
    * @return Global number of cells (hsize_t: hdf5 specific unsigned long long int)
    */
   hsize_t GetGlobalNumberOfCells() const {
      return DoGetGlobalNumberOfCells();
   }

   /**
    * @brief Gives the local number of all cells (on current rank) for the given mesh generator
    * @return Local number of cells (hsize_t: hdf5 specific unsigned long long int)
    */
   hsize_t GetLocalNumberOfCells() const {
      return DoGetLocalNumberOfCells();
   }

   /**
    * @brief Gives the local start index of cells (on current rank) for the given mesh generator. The Start index
    *        is used to access the correct position in a hdf5 hyperslab.
    * @return Local start index of cells (hsize_t: hdf5 specific unsigned long long int)
    */
   hsize_t GetLocalCellsStartIndex() const {
      return DoGetLocalCellsStartIndex();
   }

   /**
    * @brief Gives the global dimensions (on all ranks) of the vertex coordinates for the given mesh generator
    * @return 2-dimensional vector with the specific dimensions of vertex coordinates (hsize_t: hdf5 specific unsigned long long int)
    */
   std::vector<hsize_t> GetGlobalDimensionsOfVertexCoordinates() const {
      return DoGetGlobalDimensionsOfVertexCoordinates();
   }

   /**
    * @brief Gives the local dimensions (on current rank) of the vertex coordinates for the given mesh generator
    * @return 2-dimensional vector with the specific dimensions of vertex coordinates (hsize_t: hdf5 specific unsigned long long int)
    */
   std::vector<hsize_t> GetLocalDimensionsOfVertexCoordinates() const {
      return DoGetLocalDimensionsOfVertexCoordinates();
   }

   /**
    * @brief Gives the start index (on current rank) of the vertex coordinates
    *        for the given mesh generator on the current rank
    * @return Local start index of vertex coordinates (hsize_t: hdf5 specific unsigned long long int)
    */
   hsize_t GetLocalVertexCoordinatesStartIndex() const {
      return DoGetLocalVertexCoordinatesStartIndex();
   }

   /**
    * @brief Gives the global dimensions of the vertex IDs for the given mesh generator
    * @return 2-dimensional vector with the specific dimensions of vertex IDs (hsize_t: hdf5 specific unsigned long long int)
    */
   std::vector<hsize_t> GetGlobalDimensionsOfVertexIDs() const {
      return { DoGetGlobalNumberOfCells(), hsize_t( 8 ) };
   };

   /**
    * @brief Gives the local dimensions of the vertex IDs for the given mesh generator on the current rank
    * @return 2-dimensional vector with the specific dimensions of vertex IDs (hsize_t: hdf5 specific unsigned long long int)
    */
   std::vector<hsize_t> GetLocalDimensionsOfVertexIDs() const {
      return { DoGetLocalNumberOfCells(), hsize_t( 8 ) };
   };

   /**
    * @brief Gives the start index for each local dimension of the vertex IDs
    *        for the given mesh generator on the current rank
    * @return Local start index of vertex IDs (hsize_t: hdf5 specific unsigned long long int)
    */
   hsize_t GetLocalVertexIDsStartIndex() const {
      return DoGetLocalCellsStartIndex();
   };

   /**
    * @brief Return the name of the vertex IDs to be used in the files
    * @return .
    */
   std::string GetVertexIDsName() const {
      return vertex_ids_name_;
   }

   /**
    * @brief Return the name of the vertex coordinates to be used in the files
    * @return .
    */
   std::string GetVertexCoordinatesName() const {
      return vertex_coordinates_name_;
   }

   // Functions to append vertex IDs and coordinatesS
   void ComputeVertexIDs( std::vector<unsigned long long int> & vertex_ids ) const;
   void ComputeVertexCoordinates( std::vector<double> & vertex_coordinates ) const;

   // Creates the appropriate strings for the topology (vertex ids) and geometry (vertex coordinates)
   std::string GetXdmfTopologyString( std::string const& filename, std::string const& group_name ) const;
   std::string GetXdmfGeometryString( std::string const& filename, std::string const& group_name ) const;
};

#endif // MESH_QUANTITY_H