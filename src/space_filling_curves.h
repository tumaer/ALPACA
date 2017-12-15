/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*  \\\\                                                                                  *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA                                                                                 *
* Copyright (c) 2017 Nikolaus A. Adams and contributors (see AUTHORS list)               *
* All rights reserved.                                                                   *
*                                                                                        *
* Chair of Aerodynamics and Fluid Mechanics                                              *
* Technical University of Munich                                                         *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
* This project has received funding from the European Reseach Council (ERC)              *
* under the European Union's Horizon 2020 research and innovation programme              *
* (grant agreement No 667483).                                                           *
*                                                                                        *
* ERC Advanced Grant No 667483, Prof. Dr. Nikolaus A. Adams:                             *
* "NANOSHOCK - Manufacturing Shock Interactions for Innovative Nanoscale Processes"      *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Redistribution and use in source and binary forms, with or without                     *
* modification, are permitted provided that the following conditions are met:            *
*                                                                                        *
* 1. Redistributions of source code must retain the above copyright notice,              *
*    this list of conditions and the following disclaimer.                               *
*                                                                                        *
* 2. Redistributions in binary form must reproduce the above copyright notice            *
*    this list of conditions and the following disclaimer in the documentation           *
*    and/or other materials provided with the distribution.                              *
*                                                                                        *
* 3. Neither the name of the copyright holder nor the names of its                       *
*    contributors may be used to endorse or promote products derived from this           *
*    software without specific prior written permission.                                 *
*                                                                                        *
* 4. Any redistribution of substantial fractions of the code as a                        *
*    different project should preserve the word ALPACA in the name                       *
*    of the code                                                                         *
*                                                                                        *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"            *
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE              *
* IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE            *
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE              *
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR                    *
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF                   *
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS               *
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN                *
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)                *
* ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE            *
* POSSIBILITY OF SUCH DAMAGE.                                                            *
*                                                                                        *
* Please note, several third-party tools are used within the ALPACA code under           *
* their own license agreement.                                                           *
*                                                                                        *
* 1. xdmf_writer        : Licensed by Technische Universitaet Muenchen                   *
*                         See 'COPYING_XDMF_WRITER' for more information.                *
*                                                                                        *
* 2. tiny_xml           : This software is provided 'as-is', without any express or      *
*                         implied warranty. In no event will the authors be held         *
*                         liable for any damages arising from the use of this software.  *
*                         See COPYING_TINY_XMLfor more information.                      *
*                                                                                        *
* 3. expression_toolkit : Free use of The C++ Mathematical Expression Toolkit Library is *
*                         permitted under the guidelines and in accordance with the most *
*                         current version of the Common Public License.                  *
*                         http://www.opensource.org/licenses/cpl1.0.php                  *
*                         See COPYING_EXPRESSION_TOOLKITfor more information.            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* AUTHORS                                                                                *
*                                                                                        *
*   Prof. Dr. Nikolaus A. Adams                                                          *
*                                                                                        *
*   Dr. Stefan Adami                                                                     *
*   Vladimir Bogdanov                                                                    *
*   Nico Fleischmann                                                                     *
*   Nils Hoppe                                                                           *
*   Naeimeh Hosseini                                                                     *
*   Jakob Kaiser                                                                         *
*   Aleksandr Lunkov                                                                     *
*   Thomas Paula                                                                         *
*   Josef Winter                                                                         *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
*   nanoshock@aer.mw.tum.de                                                              *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, December 15th 2017                                                             *
*                                                                                        *
*****************************************************************************************/

#ifndef SPACE_FILLING_CURVES_H
#define SPACE_FILLING_CURVES_H

#include <array>

/* Nomenclature according to Bader, M.
 * Space Filling curves: An Introduction with Applications in Scientifc Computing
 * ISBN 978-3-642-31045-4 (2013), P. 115
 * _x implies "x bar"
 */
enum class HilbertPosition : unsigned short {xyz, x_y_z, yzx, y_z_x, zxy, z_x_y, _xy_z, _x_yz, _yz_x, _y_zx, _zx_y,  _z_xy};

/**
 * @brief The SpaceFillingCurves class provides space filling curve Information %Currently only for Hilbert Curve% for the Load Balancing.
 * $PLEASE NOTE: SPACEFILLING CURVE LOAD BALANCING IS ONLY USEFUL IN 3D SIMULATIONS. HENCE THIS CLASS ONLY PROVIDES 3D VARIANTS$
 */
class SpaceFillingCurves {

  // Due to nomenclature in Litertaure no trailing underscores ...
  static constexpr std::array<HilbertPosition,8> Rxyz   ={{ HilbertPosition::yzx,    HilbertPosition::zxy,    HilbertPosition::zxy,    HilbertPosition::_x_yz,
                                                            HilbertPosition::_x_yz,  HilbertPosition::_zx_y,  HilbertPosition::_zx_y,  HilbertPosition::y_z_x}};
  static constexpr std::array<HilbertPosition,8> Rx_y_z ={{ HilbertPosition::y_z_x,  HilbertPosition::z_x_y,  HilbertPosition::z_x_y,  HilbertPosition::_xy_z,
                                                            HilbertPosition::_xy_z,  HilbertPosition::_z_xy,  HilbertPosition::_z_xy,  HilbertPosition::yzx  }};
  static constexpr std::array<HilbertPosition,8> Ryzx   ={{ HilbertPosition::zxy,    HilbertPosition::xyz,    HilbertPosition::xyz,    HilbertPosition::_yz_x,
                                                            HilbertPosition::_yz_x,  HilbertPosition::x_y_z,  HilbertPosition::x_y_z,  HilbertPosition::_z_xy}};
  static constexpr std::array<HilbertPosition,8> Ry_z_x ={{ HilbertPosition::z_x_y,  HilbertPosition::x_y_z,  HilbertPosition::x_y_z,  HilbertPosition::_y_zx,
                                                            HilbertPosition::_y_zx,  HilbertPosition::xyz,    HilbertPosition::xyz,    HilbertPosition::_zx_y}};

  static constexpr std::array<HilbertPosition,8> Rzxy   ={{ HilbertPosition::xyz,    HilbertPosition::yzx,    HilbertPosition::yzx,    HilbertPosition::z_x_y,
                                                            HilbertPosition::z_x_y,  HilbertPosition::_y_zx,  HilbertPosition::_y_zx,  HilbertPosition::_xy_z}};
  static constexpr std::array<HilbertPosition,8> Rz_x_y ={{ HilbertPosition::x_y_z,  HilbertPosition::y_z_x,  HilbertPosition::y_z_x,  HilbertPosition::zxy,
                                                            HilbertPosition::zxy,    HilbertPosition::_yz_x,  HilbertPosition::_yz_x,  HilbertPosition::_x_yz}};
  static constexpr std::array<HilbertPosition,8> R_xy_z ={{ HilbertPosition::_yz_x,  HilbertPosition::_zx_y,  HilbertPosition::_zx_y,  HilbertPosition::x_y_z,
                                                            HilbertPosition::x_y_z,  HilbertPosition::zxy,    HilbertPosition::zxy,    HilbertPosition::_y_zx}};
  static constexpr std::array<HilbertPosition,8> R_x_yz ={{ HilbertPosition::_y_zx,  HilbertPosition::_z_xy,  HilbertPosition::_z_xy,  HilbertPosition::xyz,
                                                            HilbertPosition::xyz,    HilbertPosition::z_x_y,  HilbertPosition::z_x_y,  HilbertPosition::_yz_x}};

  static constexpr std::array<HilbertPosition,8> R_yz_x ={{ HilbertPosition::_zx_y,  HilbertPosition::_xy_z,  HilbertPosition::_xy_z,  HilbertPosition::yzx,
                                                            HilbertPosition::yzx,    HilbertPosition::_x_yz,  HilbertPosition::_x_yz,  HilbertPosition::z_x_y}};
  static constexpr std::array<HilbertPosition,8> R_y_zx ={{ HilbertPosition::_z_xy,  HilbertPosition::_x_yz,  HilbertPosition::_x_yz,  HilbertPosition::y_z_x,
                                                            HilbertPosition::y_z_x,  HilbertPosition::_xy_z,  HilbertPosition::_xy_z,  HilbertPosition::zxy  }};
  static constexpr std::array<HilbertPosition,8> R_zx_y ={{ HilbertPosition::_xy_z,  HilbertPosition::_yz_x,  HilbertPosition::_yz_x,  HilbertPosition::_z_xy,
                                                            HilbertPosition::_z_xy,  HilbertPosition::y_z_x,  HilbertPosition::y_z_x,  HilbertPosition::xyz  }};
  static constexpr std::array<HilbertPosition,8> R_z_xy ={{ HilbertPosition::_x_yz,  HilbertPosition::_y_zx,  HilbertPosition::_y_zx,  HilbertPosition::_zx_y,
                                                            HilbertPosition::_zx_y,  HilbertPosition::yzx,    HilbertPosition::yzx,    HilbertPosition::x_y_z}};

  static constexpr std::array<unsigned short,8> Oxyz   ={{0,2,3,1, 5,7,6,4}};
  static constexpr std::array<unsigned short,8> Ox_y_z ={{6,4,5,7, 3,1,0,2}};
  static constexpr std::array<unsigned short,8> Oyzx   ={{0,1,5,4, 6,7,3,2}};
  static constexpr std::array<unsigned short,8> Oy_z_x ={{6,7,3,2, 0,1,5,4}};

  static constexpr std::array<unsigned short,8> Ozxy   ={{0,4,6,2, 3,7,5,1}};
  static constexpr std::array<unsigned short,8> Oz_x_y ={{6,2,0,4, 5,1,3,7}};
  static constexpr std::array<unsigned short,8> O_xy_z ={{5,7,6,4, 0,2,3,1}};
  static constexpr std::array<unsigned short,8> O_x_yz ={{3,1,0,2, 6,4,5,7}};

  static constexpr std::array<unsigned short,8> O_yz_x ={{5,4,0,1, 3,2,6,7}};
  static constexpr std::array<unsigned short,8> O_y_zx ={{3,2,6,7, 5,4,0,1}};
  static constexpr std::array<unsigned short,8> O_zx_y ={{5,1,3,7, 6,2,0,4}};
  static constexpr std::array<unsigned short,8> O_z_xy ={{3,7,5,1, 0,4,6,2}};

  static constexpr std::array<unsigned short,8> O_z_curve = {{0,1,2,3, 4,5,6,7}};

  public:
  /**
     * @brief Gives the Hilbert Order for the given position.
     * @param position The current position of the cube in the Hilbert Curve.
     * @return The traversal order.
     */
    static std::array<unsigned short,8> GetHilbertOrder(const HilbertPosition position) {
        switch(position){
            case HilbertPosition::xyz:   return Oxyz;
            case HilbertPosition::x_y_z: return Ox_y_z;
            case HilbertPosition::yzx:   return Oyzx;
            case HilbertPosition::y_z_x: return Oy_z_x;

            case HilbertPosition::zxy:   return Ozxy;
            case HilbertPosition::z_x_y: return Oz_x_y;
            case HilbertPosition::_xy_z: return O_xy_z;
            case HilbertPosition::_x_yz: return O_x_yz;

            case HilbertPosition::_yz_x: return O_yz_x;
            case HilbertPosition::_y_zx: return O_y_zx;
            case HilbertPosition::_zx_y: return O_zx_y;
            case HilbertPosition::_z_xy: return O_z_xy;
        }
    return O_z_curve;
    }

    /**
     * @brief Gives the Hilbert Replacement for the given Position.
     * @param position The current position of the cube in the Hilbert Curve.
     * @return The replacement order.
     */
    static std::array<HilbertPosition,8> GetHilbertReplacement(const HilbertPosition position) {
        switch(position){
            case HilbertPosition::xyz:   return Rxyz;
            case HilbertPosition::x_y_z: return Rx_y_z;
            case HilbertPosition::yzx:   return Ryzx;
            case HilbertPosition::y_z_x: return Ry_z_x;

            case HilbertPosition::zxy:   return Rzxy;
            case HilbertPosition::z_x_y: return Rz_x_y;
            case HilbertPosition::_xy_z: return R_xy_z;
            case HilbertPosition::_x_yz: return R_x_yz;

            case HilbertPosition::_yz_x: return R_yz_x;
            case HilbertPosition::_y_zx: return R_y_zx;
            case HilbertPosition::_zx_y: return R_zx_y;
            case HilbertPosition::_z_xy: return R_z_xy;
        }
        return Rxyz;
    }
};

#endif // SPACE_FILLING_CURVES_H
