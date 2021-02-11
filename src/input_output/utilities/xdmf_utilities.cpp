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
#include "input_output/utilities/xdmf_utilities.h"
#include "utilities/string_operations.h"

namespace XdmfUtilities {

   namespace {
      /**
       * @brief Retruns a string for the dimensions where all components are split by spaces.
       * @param dimensions The Dimensions to be formatted.
       * @return the string.
       */
      template<typename T>
      std::string DimensionsToString( std::array<T, 2> const& dimensions ) {
         // If dimensions is empty return also an empty string
         if( dimensions.empty() ) {
            return "";
         } else if( dimensions.front() == 1 && dimensions.back() == 1 ) {// scalar quantity also return empty string
            return "";
         } else if( dimensions.back() == 1 && dimensions.front() <= 3 ) {// vectorial (up to three components) return only the first
            return " " + std::to_string( dimensions.front() );
         } else {// multi-dimensional return all
            std::string dimensions_string;
            for( auto const& element : dimensions ) {
               dimensions_string += " " + std::to_string( element );
            }
            return dimensions_string;
         }
      }
   }// namespace

   /**
    * @brief Generates a properly formated string for the DataItem node for the temporal information in the Xdmf file.
    * @param output_time Time to be added to the data item.
    * @return string for the time data item.
    */
   std::string TimeDataItem( double const output_time ) {
      return StringOperations::Indent( 6 ) + "<Time TimeType=\"Single\" Value=\"" + StringOperations::ToScientificNotationString( output_time ) + "\" />\n";
   }

   /**
    * @brief Generates a properly formated string for DataItem nodes in the Xdmf file.
    * @param hdf5_filename The name of the hdf5 file (without path).
    * @param item_name The name of the dataitem where the data can be found in the hdf5 file.
    * @param number_of_cells The total number of cells in the data item.
    * @param dimensions Dimensions the data item consists of ({1,1} : scalar, {3,1} : vector, {n,m} : matrix/tensor).
    * @return string for data item item.
    */
   std::string DataItemString( std::string const& hdf5_filename, std::string const& item_name, unsigned long long int const number_of_cells, std::array<unsigned int, 2> const& dimensions ) {

      return "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"" + std::to_string( number_of_cells ) + DimensionsToString( dimensions ) + "\"> " + hdf5_filename + ":/" + item_name + " </DataItem>\n";
   }

   /**
    * @brief Returns the attribute string used for the description of the topology in the Xdmf file.
    *        The topology describes the vertex IDs forming a given cell. In this framework only cubic hexahedrons are used.
    *        One cell is build by 8 vertices (IDs). The actual information of the vertex IDs (where they are placed in hdf5 file) is contained in
    *        the data_item inside of the topology string.
    *
    * @param data_item data_item to be placed in the topology string (where the vertex ID information can be found in hdf5 file).
    * @param number_of_cells total number of cells of the complete domain.
    * @return Attribute string for the topology.
    */
   std::string TopologyString( std::string const& data_item, unsigned long long int const number_of_cells ) {
      return StringOperations::Indent( 6 ) + "<Topology TopologyType=\"Hexahedron\" NumberOfElements=\"" + std::to_string( number_of_cells ) + "\">\n" + StringOperations::Indent( 8 ) + data_item + StringOperations::Indent( 6 ) + "</Topology>\n";
   }

   /**
    * @brief Returns the attribute string used for the description of the geometry in the Xdmf file.
    *        The geometry describes the vertex coordinates for each vertex. Therefore, for each vertex ID (ID=0 corresponds to first
    *        position in the coordinates data item) three coordinates must be specified. This information is placed in the data item
    *        inside the geometry attribute string (where it can be found in the hdf5 file).
    *
    * @param data_item data_item to be placed in the geometry string (where the vertex coordinates can be found in hdf5 file).
    * @param number_of_vertices total number of vertices of the complete domain.
    * @return Attribute string for the geometry.
    */
   std::string GeometryString( std::string const& data_item, unsigned long long int const number_of_vertices ) {
      return StringOperations::Indent( 6 ) + "<Geometry name=\"geometry\" GeometryType=\"XYZ\" NumberOfElements=\"" + std::to_string( number_of_vertices ) + "\">\n" + StringOperations::Indent( 8 ) + data_item + StringOperations::Indent( 6 ) + "</Geometry>\n";
   }

   /**
    * @brief Generates a properly formated string for an Attribute node of a scalar quantity in the Xdmf file.
    * @param attribute_name The name of the attribute.
    * @param data_item The data item information placed into the attribute (where the data vector can be found in the hdf5 file).
    * @return the complete attribute string.
    */
   std::string ScalarAttributeString( std::string const& attribute_name, std::string const& data_item ) {
      return StringOperations::Indent( 6 ) + "<Attribute Name=\"" + attribute_name + "\" AttributeType=\"Scalar\" Center=\"Cell\">\n" + StringOperations::Indent( 8 ) + data_item + StringOperations::Indent( 6 ) + "</Attribute>\n";
   }

   /**
    * @brief Generates a properly formated string for an Attribute node of a vectorial quantity in the Xdmf file ( 3 x 1 maximum).
    * @param attribute_name The name of the attribute.
    * @param data_item The data item information placed into the attribute (where the data vector can be found in the hdf5 file).
    * @return the complete attribute string.
    */
   std::string VectorAttributeString( std::string const& attribute_name, std::string const& data_item ) {
      return StringOperations::Indent( 6 ) + "<Attribute Name=\"" + attribute_name + "\" AttributeType=\"Vector\" Center=\"Cell\">\n" + StringOperations::Indent( 8 ) + data_item + StringOperations::Indent( 6 ) + "</Attribute>\n";
   }

   /**
    * @brief Generates a properly formated string for an Attribute node of a matrix quantity in the Xdmf file ( n x m ).
    * @param attribute_name The name of the attribute.
    * @param data_item The data item information placed into the attribute (where the data vector can be found in the hdf5 file).
    * @return the complete attribute string.
    */
   std::string MatrixAttributeString( std::string const& attribute_name, std::string const& data_item ) {
      return StringOperations::Indent( 6 ) + "<Attribute Name=\"" + attribute_name + "\" AttributeType=\"Matrix\" Center=\"Cell\">\n" + StringOperations::Indent( 8 ) + data_item + StringOperations::Indent( 6 ) + "</Attribute>\n";
   }

   /**
    * @brief Generates a properly formated string for an Attribute node of a tensor quantity in the Xdmf file ( nxn with maximum 3x3 dimensions).
    * @param attribute_name The name of the attribute.
    * @param data_item The data item information placed into the attribute (where the data vector can be found in the hdf5 file).
    * @return the complete attribute string.
    */
   std::string TensorAttributeString( std::string const& attribute_name, std::string const& data_item ) {
      return StringOperations::Indent( 6 ) + "<Attribute Name=\"" + attribute_name + "\" AttributeType=\"Tensor\" Center=\"Cell\">\n" + StringOperations::Indent( 8 ) + data_item + StringOperations::Indent( 6 ) + "</Attribute>\n";
   }

   /**
    * @brief Gives a properly formatted string for the complete spatial data information.
    * @param spatial_data_name Name of the grid to be used.
    * @param spatial_data_information Complete data information (topology and data item).
    * @return Complete formatted spatial data string.
    */
   std::string SpatialDataInformation( std::string const& spatial_data_name, std::string const& spatial_data_information ) {
      return StringOperations::Indent( 4 ) + "<Grid Name=\"" + spatial_data_name + "\" GridType=\"Uniform\">\n" + spatial_data_information + StringOperations::Indent( 4 ) + "</Grid>\n";
   }

   /**
   * @brief Gives a properly formatted string for the XDMF file header.
   * @param data_name Name of the complete data stored for this file.
   * @return string with header information.
   */
   std::string HeaderInformation( std::string const& data_name ) {
      return StringOperations::Indent( 0 ) + "<?xml version=\"1.0\" ?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n" + StringOperations::Indent( 0 ) + "<Xdmf Version=\"3.0\">\n" + StringOperations::Indent( 1 ) + "<Domain>\n" + StringOperations::Indent( 2 ) + "<Grid Name=\"" + data_name + "\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";
   }

   /**
   * @brief Give a properly formatted string for the XDMF file footer (closes all tags opened by the header).
   * @return The footer string.
   */
   std::string FooterInformation() {
      return StringOperations::Indent( 2 ) + "</Grid>\n" + StringOperations::Indent( 1 ) + "</Domain>\n" + StringOperations::Indent( 0 ) + "</Xdmf>";
   }
}// namespace XdmfUtilities
