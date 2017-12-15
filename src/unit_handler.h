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

#ifndef UNIT_HANDLER
#define UNIT_HANDLER

#include <limits>
#include "inputfile_parser.h"
#include "log_writer.h"

/**
 * @brief The UnitHandler class takes care of the conversions between unitless and non-unitless representations of quantities.
 *        Within the Simulation Kernel only non-dimensional quantities are used. In the user input and output,
 *        however, values are given in unit repesentation.
 */
class UnitHandler{

  const double length_reference_;
  const double velocity_reference_;
  const double density_reference_;
  const double mu_reference_;
  const double gravity_reference_;
  const double time_reference_;

  public:
    UnitHandler() = delete;
    UnitHandler(const InputFileParser& parser, LogWriter &logger);

    double ComputeReynoldsNumber() const;
    double ComputeFroudeNumber() const;

    /**
     * @brief Translates a unit based density value into a unitless one, using the provided reference density.
     * @param rho_dimensional Unit based density.
     * @return Unitless density value.
     */
    inline double NonDimensionalizeDensity(const double density_dimensional) const {return density_dimensional/density_reference_;}

    /**
     * @brief Translates a unit based length/size/position value into a unitless one, using the provided reference length.
     * @param length_dimensional Unit based length.
     * @return Unitless length value.
     */
    inline double NonDimensionalizeLength(const double length_dimensional) const {return length_dimensional/length_reference_;}

    /**
     * @brief Translates a unit based velocity value into a unitless one, using the provided reference velocity.
     * @param velocity_dimensional Unit based velocity.
     * @return Unitless velocity value.
     */
    inline double NonDimensionalizeVelocity(double velocity_dimensional) const {return velocity_dimensional/velocity_reference_;}

    /**
     * @brief Translates a unit based viscosity value into a unitless one, using the provided reference mu.
     * @param mu_dimensional Unit based viscosity.
     * @return Unitless viscosity value.
     */
    inline double NonDimensionalizeViscosity(double mu_dimensional) const {return mu_dimensional/mu_reference_;}

    /**
     * @brief Translates a unit based time value into a unitless one, using the provided reference time.
     * @param time_dimensional Unit based time.
     * @return Unitless time value.
     */
    inline double NonDimensionalizeTime(const double time_dimensional) const {return time_dimensional/time_reference_;}

    /**
     * @brief Translates a unit based pressure value into a unitless one, using the provided reference parameters: velocity and density.
     * @param pressure_dimensional Unit pressure density.
     * @return Unitless pressure value.
     */
    inline double NonDimensionalizePressure(const double pressure_dimensional) const {
        return pressure_dimensional/(velocity_reference_*velocity_reference_*density_reference_);
    }

    /**
     * @brief Translates a unit based gravity value into a unitless one, using the provided reference gravity.
     * @param gravity Unit based gravity.
     * @return Unitless gravity value.
     */
    inline std::array<double,3> NonDimensionalizeGravity(const std::array<double,3> gravity) const {
        return {gravity[0]/gravity_reference_,gravity[1]/gravity_reference_,gravity[2]/gravity_reference_};
    }

    /**
     * @brief Translates a unitless density value into a unit based one, unsing the provided reference density.
     * @param density_unitless Unitless density.
     * @return Unit based density value.
     */
    inline double DimensionalizeDensity(const double density_unitless) const {return density_unitless*density_reference_;}

    /**
     * @brief Translates a unitless momentum value into a unit based one, unsing the provided reference density and velocity.
     * @param momentum_unitless Unitless momentum.
     * @return Unit based momentum value.
     */
    inline double DimensionalizeMomentum(const double momentum_unitless) const {return momentum_unitless*density_reference_*velocity_reference_;}

    /**
     * @brief Translates a unitless length/size/position value into a unit based one, unsing the provided reference length.
     * @param length_unitless Unitless length/size/position.
     * @return Unit based length/size/position value.
     */
    inline double DimensionalizeLength(const double length_unitless) const {return length_unitless*length_reference_;}

    /**
     * @brief Translates a unitless time value into a unit based one, unsing the provided reference time.
     * @param time_unitless Unitless time.
     * @return Unit based time value.
     */
    inline double DimensionalizeTime(const double time_unitless) const {return time_unitless*time_reference_;}

    /**
     * @brief Translates a unitless energy value into a unit based one, unsing the provided reference energy.
     * @param energy_unitless Unitless energy.
     * @return Unit based energy value.
     */
    inline double DimensionalizeEnergy(const double energy_unitless) const {
        return energy_unitless*(velocity_reference_*velocity_reference_*density_reference_);
    }

};

#endif // UNIT_HANDLER
