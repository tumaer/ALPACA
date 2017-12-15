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

#ifndef INPUTFILE_PARSER_H
#define INPUTFILE_PARSER_H

#include <string>
#include <array>
#include <vector>
#include "3rdParty/tiny_xml/tinyxml2.h"
#include "boundary_condition/boundary_specifications.h"
#include "materials/material_names.h"
#include "output/output_types.h"
#include "compile_time_constants.h"

/**
 * @brief The InputFileParser class provides functions to read out a xml-inputfile. The input is parsed and converted into machine
 *        readable formats (strings, enums, doubles, etc.). The input is neither checked for plausibility nor converted, e.g. into conservative form.
 * @note Uses the TinyXML library by Lee Thomason (www.grinninglizard.com). See respective files for License and Copyright information.
 */
class InputFileParser{

  tinyxml2::XMLDocument               inputfile_;

  double ReadDouble(const tinyxml2::XMLElement *node, const char* tag, const double default_value) const;
  int ReadInt(const tinyxml2::XMLElement *node, const char* tag, const int default_value) const;
  std::string ReadString(const tinyxml2::XMLElement *node) const;

  BoundaryType SelectBoundaryCondition(const std::string boundary_string) const;
  std::string ReadInitialConditionString(const std::string name) const;

  MaterialName ReadMaterialName(const std::string fluid) const;
  MaterialName SelectMaterial(const std::string material_name) const;
  std::vector<double> ReadMaterialProperties(const std::string fluid) const;

  public:

    InputFileParser(const std::string & filename);

    double ReadBlockSize() const;
    int ReadNumberOfBlocksX() const;
    int ReadNumberOfBlocksY() const;
    int ReadNumberOfBlocksZ() const;
    std::array<BoundaryType, 6> ReadBoundaryConditions() const;
    std::array<double,CC::NoEq()> ReadFixedValueBoundaryCondition(const BoundaryLocation location) const;

    std::vector<std::string> ReadInitialConditionOfFluids() const;

    int ReadNumberOfFluids() const;
    std::vector<std::pair<MaterialName, std::vector<double>>> ReadParameterOfAllFluids() const;

    std::array<double, 3> ReadGravity() const;

    int ReadMaximumLevel() const;
    int ReadLevelOfReferenceMultiResolutionEpsilon() const;
    double ReadReferenceMultiResolutionEpsilon() const;

    double ReadStartTime() const;
    double ReadEndTime() const;
    double ReadCFLnumber() const;

    double ReadReferenceLength() const ;
    double ReadReferenceVelocity() const;
    double ReadReferenceDensity() const;
    double ReadReferenceMu() const;
    double ReadReferenceGravity() const;
    double ReadReferenceTime() const;

    OutputType ReadOutputFormat() const;
    std::string ReadOutputTimesType() const;
    double ReadOutputPeriod() const;
    int ReadNumberOfTimestamps() const;
    std::vector<double> ReadTimestamps() const;
    double ReadTimeNamingFactor() const;

};
#endif // INPUTFILE_PARSER_H
