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

#include <array>

#include "input_output/utilities/xdmf_utilities.h"

SCENARIO( "TimeDataItems can be correctly populated", "[1rank]" ) {
   GIVEN( "Output time t = 0.5 " ) {
      double const output_time = 0.5;
      WHEN( "Created with output time t = 0.5" ) {
         REQUIRE( XdmfUtilities::TimeDataItem( output_time ) == "      <Time TimeType=\"Single\" Value=\"5.0000000000000000e-01\" />\n" );
      }
   }
}

SCENARIO( "DataItems can be correctly populated", "[1rank]" ) {
   GIVEN( "An hdf5 filename (test.h5) + 42 global cells + a dataset named ds" ) {
      std::string const filename         = "test.h5";
      unsigned int const number_of_cells = 42;
      std::string const dataset          = "ds";
      WHEN( "Created with dimensions of scalar" ) {
         std::array<unsigned int, 2> const dimensions = { 1, 1 };
         REQUIRE( XdmfUtilities::DataItemString( filename, dataset, number_of_cells, dimensions ) == "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"42\"> test.h5:/ds </DataItem>\n" );
      }
      WHEN( "Created with dimensions of vector (first component < 3) " ) {
         std::array<unsigned int, 2> const dimensions = { 2, 1 };
         REQUIRE( XdmfUtilities::DataItemString( filename, dataset, number_of_cells, dimensions ) == "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"42 2\"> test.h5:/ds </DataItem>\n" );
      }
      WHEN( "Created with dimensions of matrix (first component = 1, second component > 1)" ) {
         std::array<unsigned int, 2> const dimensions = { 1, 2 };
         REQUIRE( XdmfUtilities::DataItemString( filename, dataset, number_of_cells, dimensions ) == "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"42 1 2\"> test.h5:/ds </DataItem>\n" );
      }
      WHEN( "Created with dimensions of matrix (first component > 1, second component > 1)" ) {
         std::array<unsigned int, 2> const dimensions = { 2, 2 };
         REQUIRE( XdmfUtilities::DataItemString( filename, dataset, number_of_cells, dimensions ) == "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"42 2 2\"> test.h5:/ds </DataItem>\n" );
      }
   }
}

SCENARIO( "Attributes can be properly created", "[1rank]" ) {
   GIVEN( "A scalar attribute name (scalar)" ) {
      std::string const name = "scalar";
      WHEN( "Created with empty data item string" ) {
         REQUIRE( XdmfUtilities::ScalarAttributeString( name, "" ) == "      <Attribute Name=\"scalar\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
                                                                      "        "
                                                                      "      </Attribute>\n" );
      }
      WHEN( "Created with a data item string (data_item)" ) {
         std::string const data_item = "data_item";
         REQUIRE( XdmfUtilities::ScalarAttributeString( name, data_item ) == "      <Attribute Name=\"scalar\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
                                                                             "        data_item"
                                                                             "      </Attribute>\n" );
      }
   }
   GIVEN( "A vector attribute name (vector)" ) {
      std::string const name = "vector";
      WHEN( "Created with empty data item string" ) {
         REQUIRE( XdmfUtilities::VectorAttributeString( name, "" ) == "      <Attribute Name=\"vector\" AttributeType=\"Vector\" Center=\"Cell\">\n"
                                                                      "        "
                                                                      "      </Attribute>\n" );
      }
      WHEN( "Created with a data item string (data_item)" ) {
         std::string const data_item = "data_item";
         REQUIRE( XdmfUtilities::VectorAttributeString( name, data_item ) == "      <Attribute Name=\"vector\" AttributeType=\"Vector\" Center=\"Cell\">\n"
                                                                             "        data_item"
                                                                             "      </Attribute>\n" );
      }
   }
   GIVEN( "A matrix attribute name (matrix)" ) {
      std::string const name = "matrix";
      WHEN( "Created with empty data item string" ) {
         REQUIRE( XdmfUtilities::MatrixAttributeString( name, "" ) == "      <Attribute Name=\"matrix\" AttributeType=\"Matrix\" Center=\"Cell\">\n"
                                                                      "        "
                                                                      "      </Attribute>\n" );
      }
      WHEN( "Created with a data item string (data_item)" ) {
         std::string const data_item = "data_item";
         REQUIRE( XdmfUtilities::MatrixAttributeString( name, data_item ) == "      <Attribute Name=\"matrix\" AttributeType=\"Matrix\" Center=\"Cell\">\n"
                                                                             "        data_item"
                                                                             "      </Attribute>\n" );
      }
   }
   GIVEN( "A tensor attribute name (tensor)" ) {
      std::string const name = "tensor";
      WHEN( "Created with empty data item string" ) {
         REQUIRE( XdmfUtilities::TensorAttributeString( name, "" ) == "      <Attribute Name=\"tensor\" AttributeType=\"Tensor\" Center=\"Cell\">\n"
                                                                      "        "
                                                                      "      </Attribute>\n" );
      }
      WHEN( "Created with a data item string (data_item)" ) {
         std::string const data_item = "data_item";
         REQUIRE( XdmfUtilities::TensorAttributeString( name, data_item ) == "      <Attribute Name=\"tensor\" AttributeType=\"Tensor\" Center=\"Cell\">\n"
                                                                             "        data_item"
                                                                             "      </Attribute>\n" );
      }
   }
}

SCENARIO( "Spatial data information string can be properly created", "[1rank]" ) {
   GIVEN( "A spatial data grid name (spatial_data)" ) {
      std::string const name = "spatial_data";
      WHEN( "Created with empty data item string" ) {
         REQUIRE( XdmfUtilities::SpatialDataInformation( name, "" ) == "    <Grid Name=\"spatial_data\" GridType=\"Uniform\">\n"
                                                                       ""
                                                                       "    </Grid>\n" );
      }
      WHEN( "Created with a data item string (spatial_data_information)" ) {
         std::string const spatial_data_information = "spatial_data_information";
         REQUIRE( XdmfUtilities::SpatialDataInformation( name, spatial_data_information ) == "    <Grid Name=\"spatial_data\" GridType=\"Uniform\">\n"
                                                                                             "spatial_data_information"
                                                                                             "    </Grid>\n" );
      }
   }
}

SCENARIO( "Topology string can be properly created", "[1rank]" ) {
   GIVEN( "Number of global cells 42" ) {
      unsigned int const cells = 42;
      WHEN( "Created with empty data item string" ) {
         REQUIRE( XdmfUtilities::TopologyString( "", cells ) == "      <Topology TopologyType=\"Hexahedron\" NumberOfElements=\"42\">\n"
                                                                "        "
                                                                "      </Topology>\n" );
      }
      WHEN( "Created with a data item string (topology_data)" ) {
         std::string const topology_data = "topology_data";
         REQUIRE( XdmfUtilities::TopologyString( topology_data, cells ) == "      <Topology TopologyType=\"Hexahedron\" NumberOfElements=\"42\">\n"
                                                                           "        topology_data"
                                                                           "      </Topology>\n" );
      }
   }
}

SCENARIO( "Geometry string can be properly created", "[1rank]" ) {
   GIVEN( "Number of global vertices 42" ) {
      unsigned int const cells = 42;
      WHEN( "Created with empty data item string" ) {
         REQUIRE( XdmfUtilities::GeometryString( "", cells ) == "      <Geometry name=\"geometry\" GeometryType=\"XYZ\" NumberOfElements=\"42\">\n"
                                                                "        "
                                                                "      </Geometry>\n" );
      }
      WHEN( "Created with a data item string (geometry_data)" ) {
         std::string const geometry_data = "geometry_data";
         REQUIRE( XdmfUtilities::GeometryString( geometry_data, cells ) == "      <Geometry name=\"geometry\" GeometryType=\"XYZ\" NumberOfElements=\"42\">\n"
                                                                           "        geometry_data"
                                                                           "      </Geometry>\n" );
      }
   }
}

SCENARIO( "Header string can be properly created", "[1rank]" ) {
   GIVEN( "A data grid name (data_name)" ) {
      std::string const name = "data_name";
      WHEN( "Created with given value" ) {
         REQUIRE( XdmfUtilities::HeaderInformation( name ) == "<?xml version=\"1.0\" ?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
                                                              "<Xdmf Version=\"3.0\">\n"
                                                              " <Domain>\n"
                                                              "  <Grid Name=\"data_name\" GridType=\"Collection\" CollectionType=\"Temporal\">\n" );
      }
   }
}

SCENARIO( "Footer string can be properly created", "[1rank]" ) {
   GIVEN( "Default values (nothing" ) {
      WHEN( "Created with default values" ) {
         REQUIRE( XdmfUtilities::FooterInformation() == "  </Grid>\n"
                                                        " </Domain>\n"
                                                        "</Xdmf>" );
      }
   }
}