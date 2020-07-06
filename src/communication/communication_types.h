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
#ifndef COMMUNICATION_TYPES_H_
#define COMMUNICATION_TYPES_H_

#include <mpi.h>
#include <array>
#include "user_specifications/compile_time_constants.h"
#include "enums/dimension_definition.h"
#include "enums/datatype_for_mpi_definition.h"

#include "block_definitions/field_material_definitions.h"

/**
 * @brief The CommunicationTypes class collects all information about the different types used for the MPI-communication. It specifies
 *        the sizes, start indices and other information to allow specific buffer communication in one single chunk.
 */
class CommunicationTypes{

   static constexpr std::array<int,3 > block_size_ = { CC::TCX(), CC::TCY(), CC::TCZ() };
   // 26 Elements for 3 Dimensions, numbering like boundary_specification.h -> BoundaryLocation
   static constexpr std::array<std::array<int, 3>,26> start_indices_halo_send_ = { {
                                                                                     { { CC::ICX(),  CC::FICY(), CC::FICZ() } }, //e  -planes
                                                                                     { { CC::FICX(), CC::FICY(), CC::FICZ() } }, //w
                                                                                     { { CC::FICX(), CC::ICY(),  CC::FICZ() } }, //n
                                                                                     { { CC::FICX(), CC::FICY(), CC::FICZ() } }, //s
                                                                                     { { CC::FICX(), CC::FICY(), CC::ICZ()  } }, //t
                                                                                     { { CC::FICX(), CC::FICY(), CC::FICZ() } }, //b
                                                                                     { { CC::FICX(), CC::ICY(),  CC::FICZ() } }, //bn -sticks
                                                                                     { { CC::FICX(), CC::FICY(), CC::FICZ() } }, //bs
                                                                                     { { CC::FICX(), CC::ICY(),  CC::ICZ()  } }, //tn
                                                                                     { { CC::FICX(), CC::FICY(), CC::ICZ()  } }, //ts
                                                                                     { { CC::ICX(),  CC::FICY(), CC::FICZ() } }, //be
                                                                                     { { CC::FICX(), CC::FICY(), CC::FICZ() } }, //bw
                                                                                     { { CC::ICX(),  CC::FICY(), CC::ICZ()  } }, //te
                                                                                     { { CC::FICX(), CC::FICY(), CC::ICZ()  } }, //tw
                                                                                     { { CC::ICX(),  CC::ICY(),  CC::FICZ() } }, //ne
                                                                                     { { CC::FICX(), CC::ICY(),  CC::FICZ() } }, //nw
                                                                                     { { CC::ICX(),  CC::FICY(), CC::FICZ() } }, //se
                                                                                     { { CC::FICX(), CC::FICY(), CC::FICZ() } }, //sw
                                                                                     { { CC::ICX(),  CC::ICY(),  CC::ICZ()  } }, //ent -cubes
                                                                                     { { CC::ICX(),  CC::ICY(),  CC::FICZ() } }, //enb
                                                                                     { { CC::ICX(),  CC::FICY(), CC::ICZ()  } }, //est
                                                                                     { { CC::ICX(),  CC::FICY(), CC::FICZ() } }, //esb
                                                                                     { { CC::FICX(), CC::ICY(),  CC::ICZ()  } }, //wnt
                                                                                     { { CC::FICX(), CC::ICY(),  CC::FICZ() } }, //wnb
                                                                                     { { CC::FICX(), CC::FICY(), CC::ICZ()  } }, //wst
                                                                                     { { CC::FICX(), CC::FICY(), CC::FICZ() } }  //wsb
                                                                                   } };
   // 26 Elements for 3 Dimensions, numbering like boundary_specification.h -> BoundaryLocation
   static constexpr std::array<std::array<int, 3>,26> start_indices_halo_recv_ = { {
                                                                                     { { CC::FHHX(), CC::FICY(), CC::FICZ() } }, //e  -planes
                                                                                     { { 0,          CC::FICY(), CC::FICZ() } }, //w
                                                                                     { { CC::FICX(), CC::FHHY(), CC::FICZ() } }, //n
                                                                                     { { CC::FICX(), 0,          CC::FICZ() } }, //s
                                                                                     { { CC::FICX(), CC::FICY(), CC::FHHZ() } }, //t
                                                                                     { { CC::FICX(), CC::FICY(), 0          } }, //b
                                                                                     { { CC::FICX(), CC::FHHY(), 0          } }, //bn -sticks
                                                                                     { { CC::FICX(), 0,          0          } }, //bs
                                                                                     { { CC::FICX(), CC::FHHY(), CC::FHHZ() } }, //tn
                                                                                     { { CC::FICX(), 0,          CC::FHHZ() } }, //ts
                                                                                     { { CC::FHHX(), CC::FICY(), 0          } }, //be
                                                                                     { { 0,          CC::FICY(), 0          } }, //bw
                                                                                     { { CC::FHHX(), CC::FICY(), CC::FHHZ() } }, //te
                                                                                     { { 0,          CC::FICY(), CC::FHHZ() } }, //tw
                                                                                     { { CC::FHHX(), CC::FHHY(), CC::FICZ() } }, //ne
                                                                                     { { 0,          CC::FHHY(), CC::FICZ() } }, //nw
                                                                                     { { CC::FHHX(), 0,          CC::FICZ() } }, //se
                                                                                     { { 0,          0,          CC::FICZ() } }, //sw
                                                                                     { { CC::FHHX(), CC::FHHY(), CC::FHHZ() } }, //ent -cubes
                                                                                     { { CC::FHHX(), CC::FHHY(), 0          } }, //enb
                                                                                     { { CC::FHHX(), 0,          CC::FHHZ() } }, //est
                                                                                     { { CC::FHHX(), 0,          0          } }, //esb
                                                                                     { { 0,          CC::FHHY(), CC::FHHZ() } }, //wnt
                                                                                     { { 0,          CC::FHHY(), 0          } }, //wnb
                                                                                     { { 0,          0,          CC::FHHZ() } }, //wst
                                                                                     { { 0,          0,          0          } }  //wsb
                                                                                   } };
   static constexpr std::array<std::array<int, 3>,26> halo_size_ = { {
                                                                       { { CC::HSSX(), CC::ICY(),  CC::ICZ()  } }, //e  -planes
                                                                       { { CC::HSSX(), CC::ICY(),  CC::ICZ()  } }, //w
                                                                       { { CC::ICX(),  CC::HSSY(), CC::ICZ()  } }, //n
                                                                       { { CC::ICX(),  CC::HSSY(), CC::ICZ()  } }, //s
                                                                       { { CC::ICX(),  CC::ICY(),  CC::HSSZ() } }, //t
                                                                       { { CC::ICX(),  CC::ICY(),  CC::HSSZ() } }, //b
                                                                       { { CC::ICX(),  CC::HSSY(), CC::HSSZ() } }, //bn -sticks
                                                                       { { CC::ICX(),  CC::HSSY(), CC::HSSZ() } }, //bs
                                                                       { { CC::ICX(),  CC::HSSY(), CC::HSSZ() } }, //tn
                                                                       { { CC::ICX(),  CC::HSSY(), CC::HSSZ() } }, //ts
                                                                       { { CC::HSSX(), CC::ICY(),  CC::HSSZ() } }, //be
                                                                       { { CC::HSSX(), CC::ICY(),  CC::HSSZ() } }, //bw
                                                                       { { CC::HSSX(), CC::ICY(),  CC::HSSZ() } }, //te
                                                                       { { CC::HSSX(), CC::ICY(),  CC::HSSZ() } }, //tw
                                                                       { { CC::HSSX(), CC::HSSY(), CC::ICZ()  } }, //ne
                                                                       { { CC::HSSX(), CC::HSSY(), CC::ICZ()  } }, //nw
                                                                       { { CC::HSSX(), CC::HSSY(), CC::ICZ()  } }, //se
                                                                       { { CC::HSSX(), CC::HSSY(), CC::ICZ()  } }, //sw
                                                                       { { CC::HSSX(), CC::HSSY(), CC::HSSZ() } }, //ent -cubes
                                                                       { { CC::HSSX(), CC::HSSY(), CC::HSSZ() } }, //enb
                                                                       { { CC::HSSX(), CC::HSSY(), CC::HSSZ() } }, //est
                                                                       { { CC::HSSX(), CC::HSSY(), CC::HSSZ() } }, //esb
                                                                       { { CC::HSSX(), CC::HSSY(), CC::HSSZ() } }, //wnt
                                                                       { { CC::HSSX(), CC::HSSY(), CC::HSSZ() } }, //wnb
                                                                       { { CC::HSSX(), CC::HSSY(), CC::HSSZ() } }, //wst
                                                                       { { CC::HSSX(), CC::HSSY(), CC::HSSZ() } }, //wsb
                                                                     } };
   // Averaging
   static constexpr std::array<int,3> child_size_ = { { CC::PSOCICX(), CC::PSOCICY(), CC::PSOCICZ() } };
   static constexpr std::array<std::array<int,3>,8> start_index_child_ = { {
                                                                             { CC::FICX(),      CC::FICY(),      CC::FICZ()      },
                                                                             { CC::PIOHCFICX(), CC::FICY(),      CC::FICZ()      },
                                                                             { CC::FICX(),      CC::PIOHCFICY(), CC::FICZ()      },
                                                                             { CC::PIOHCFICX(), CC::PIOHCFICY(), CC::FICZ()      },
                                                                             { CC::FICX(),      CC::FICY(),      CC::PIOHCFICZ() },
                                                                             { CC::PIOHCFICX(), CC::FICY(),      CC::PIOHCFICZ() },
                                                                             { CC::FICX(),      CC::PIOHCFICY(), CC::PIOHCFICZ() },
                                                                             { CC::PIOHCFICX(), CC::PIOHCFICY(), CC::PIOHCFICZ() }
                                                                           } };
   // Datatypes for Boundaries ( doubles )
   // 26 Elements for 3 Dimensions, numbering like boundary_specification.h ->BoundaryLocation
   // 1D: East, West
   // 2D: East, West, North, South
   // 3D: East, West, North, South, Top, Bottom
   std::array<std::array<MPI_Datatype, 26>, 2> send_types_;
   std::array<std::array<MPI_Datatype, 26>, 2> recv_types_;
   // Datatypes for ProjectLevel Recv into block
   std::array<std::array<MPI_Datatype, CC::NOC()>, 2> averaging_send_;
   MPI_Datatype jump_plane_ew_, jump_plane_ns_, jump_plane_tb_;
   MPI_Datatype jump_stick_x_, jump_stick_y_, jump_stick_z_;
   MPI_Datatype jump_cube_;
   // Datatypes for Load Balancing, sending whole Struct
   MPI_Datatype single_conservatives_;
   MPI_Datatype single_boundary_jump_;

   /**
    * @brief CreateDataTypes Creates and commits MPI datatypes for in order to exchange boundary subarrays and such in one go as one subarray.
    */
   void CreateDataTypes();

   /**
    * @brief FreeTypes Frees the memory of the used MPI datatypes.
    */
   void FreeTypes();

public:
   //NH TODO-19 Rule-of-five omitted here.
   explicit CommunicationTypes();
   virtual ~CommunicationTypes();

   MPI_Datatype SendDatatype(BoundaryLocation const location, DatatypeForMpi const datatype) const;
   MPI_Datatype RecvDatatype(BoundaryLocation const location, DatatypeForMpi const datatype) const;
   MPI_Datatype JumpPlaneSendDatatype(BoundaryLocation const location) const;
   MPI_Datatype ConservativesDatatype() const;
   MPI_Datatype JumpSurfaceDatatype() const;
   MPI_Datatype AveragingSendDatatype(unsigned int const child_position, DatatypeForMpi const datatype) const;

   static constexpr std::array<int, 3> GetStartIndicesHaloSend(BoundaryLocation const location) {
      return start_indices_halo_send_[LTI(location)];
   }
   static constexpr std::array<int, 3> GetStartIndicesHaloRecv(BoundaryLocation const location) {
      return start_indices_halo_recv_[LTI(location)];
   }
   static constexpr std::array<int, 3> GetHaloSize(BoundaryLocation const location) {
      return halo_size_[LTI(location)];
   }

};

/**
 * @brief Gives the total number values need to be transferred per conservative variable (+buffer copies),
 *        which need to be transferred if the content of a whole block is to be sent from one rank to another.
 * @return The size of a full block (internal + halos).
 */
constexpr int FullBlockSendingSize() {return CC::TCX()*CC::TCY()*CC::TCZ();}

/**
 * @brief Gives the total number values need to be transferred per conservative variable (+buffer copies),
 *        if the content of all internal cells of a whole block is to be sent from one rank to another.
 * @return The size of the internal block.
 */
constexpr int DomainBlockSendingSize() {return CC::ICX()*CC::ICY()*CC::ICZ();}

/**
 * @brief Gives the total number of values to be transferred if the jump flux buffers of all conservative variables are to be sent form
 *        one rank to another at once.
 * @return The size of jump flux buffer.
 */
constexpr int JumpBufferSendingSize() {return MF::ANOE()*CC::ICY()*CC::ICZ();}

#endif /* COMMUNICATION_TYPES_H_ */
