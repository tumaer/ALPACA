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

#include "inputfile_parser.h"

#include <stdexcept>
#include <iostream>
#include <algorithm>

/**
 * @brief Standard constructor opening and registering the xml file for later parsings.
 * @param filename XML-File including its path to be opend and read-out.
 */
InputFileParser::InputFileParser(const std::string& filename){

  // Open Config file
  inputfile_.LoadFile( filename.c_str());
  if (inputfile_.FirstChildElement() == NULL){
    throw std::logic_error("Error parsing the XML inputfile file");
  }

}

/**
 * @brief Reads out a numeric value from an XML node, treats and converts it into a double value.
 * @param node The XML node holding the desired information.
 * @param tag Name of the variable to be read-out
 * @param default_value Emergency return value.
 * @return The read-out and converted value.
 */
double InputFileParser::ReadDouble(const tinyxml2::XMLElement *node, const char* tag, const double default_value) const {

  double value;
  if (node->QueryDoubleText( &value) != tinyxml2::XML_NO_ERROR){
    std::cerr << "Error while reading user defined input argument" << tag << "' set to default value: " << default_value << std::endl;
    return default_value;
  } else {
    return (double)value;
  }

}

/**
 * @brief Reads out a numeric value from an XML node, treats and converts it into an int value.
 * @param node The XML node holding the desired information.
 * @param tag Name of the variable to be read-out.
 * @param default_value Emergency return value.
 * @return The read-out and converted value.
 */
int InputFileParser::ReadInt(const tinyxml2::XMLElement *node, const char *tag, const int default_value) const {

  int value;
  if (node->QueryIntText(&value) != tinyxml2::XML_NO_ERROR){
      std::cerr << "Error while reading argument (int)" << std::endl;
      std::cerr << "'" << tag << "' set to default value: " << default_value << std::endl;
      return default_value;
  } else {
      return value;
  }

}

/**
 * @brief Reads out a string value from an XML node.
 * @param node The XML node holding the desired information.
 * @return Found String
 */
std::string InputFileParser::ReadString(const tinyxml2::XMLElement *node) const{

  std::string value = node->GetText();

  if (value.empty() || !value.compare("")){
    const std::string nodename = node->Name();
    throw std::logic_error("ERROR in inputfile file in Node " + nodename);
  }
  return value;
}

/**
 * @brief Assosciates a boundary type to the string describing the boundary condition.
 * @param boundary_string The boundary description as found in the XML file.
 * @return Machine readable boundary definition.
 */
BoundaryType InputFileParser::SelectBoundaryCondition(const std::string boundary_string) const {
  // NH 2016-10-24: Switch Statements do not work with String

  //make entire string upper case
  std::string boundary_string_upper_case = boundary_string;
  std::transform(boundary_string_upper_case.begin(), boundary_string_upper_case.end(),boundary_string_upper_case.begin(), ::toupper);

  if      (boundary_string_upper_case == "ZEROGRADIENT") {return BoundaryType::eZeroGradient;}
  else if (boundary_string_upper_case == "SYMMETRY") {return BoundaryType::eSymmetry;}
  else if (boundary_string_upper_case == "FIXEDVALUE") {return BoundaryType::eFixedValue;}
  else if (boundary_string_upper_case == "WALL") {return BoundaryType::eWall;}
  else { throw std::invalid_argument("Boundary condition in input file is undefined / not yet implmented");}
}

/**
 * @brief Reads out the initial condition of the simulation. Gives the initial condition as string (equation form).
 *        No conversion to machine-readbale inputs is performed.
 * @param name The name of the object (first/second/n-th fluid, phi, ..) whose initial condition is to be retrieved.
 * @return Initial Conditions in equation form.
 */
std::string InputFileParser::ReadInitialConditionString(const std::string name) const {
  const tinyxml2::XMLElement *node = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("initialConditions")->FirstChildElement(name.c_str());
  return ReadString(node);
}

/**
 * @brief Reads-out the types of the materials, i. e. Equation of State.
 * @param fluid The name of the considered fluid as in XML file.
 * @return The material identifer enum associated with the equation of state.
 */
MaterialName InputFileParser::ReadMaterialName(const std::string fluid) const {

  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("fluids")->FirstChildElement(fluid.c_str())->FirstChildElement("type");
  std::string  material_name = ReadString(node);
  return SelectMaterial(material_name);
}

/**
 * @brief Associates a material identifer to the string describing the material.
 * @param material_name The material description as found in the XML file.
 * @return Machine readable material.
 */
MaterialName InputFileParser::SelectMaterial(const std::string material_name) const {

  //make entire string upper case
  std::string material_name_upper_case = material_name;
  std::transform(material_name_upper_case.begin(), material_name_upper_case.end(),material_name_upper_case.begin(), ::toupper);

  if (material_name_upper_case == "STIFFENEDGAS")  {return MaterialName::eStiffenedGas;}
  else if (material_name_upper_case == "WATERLIKEFLUID") {return MaterialName::eWaterlikeFluid;}
  else { throw std::invalid_argument("Material name in input file is undefined / not yet implmented");}
}

/**
 * @brief Reads out the fluid properties, i.e. material parameters.
 *        ATTENTION: do not reorder the parameters from 0 to 5, as they need to be in this order for material initialization!
 * @param fluid The name of the considered fluid as in the xml file.
 * @return The obtained values in the order
 *         0: Gamma Default: 4.4,
 *         1: A Default: 0.0,
 *         2: B Default: 6.0x10^8,
 *         3: Rho_zero Default: 998,
 *         4: dynamic shear viscosity Default 0.0,
 *         5: dynamic bulk viscosity Default 0.0.
 */
std::vector<double> InputFileParser::ReadMaterialProperties(const std::string fluid) const {

  std::vector<double> material_properties;
  const tinyxml2::XMLElement* node1 = inputfile_.FirstChildElement()->FirstChildElement("fluids")->FirstChildElement(fluid.c_str())->FirstChildElement("gamma");
  material_properties.push_back(ReadDouble(node1, "gamma", 4.4));
  const tinyxml2::XMLElement* node2 = inputfile_.FirstChildElement()->FirstChildElement("fluids")->FirstChildElement(fluid.c_str())->FirstChildElement("A");
  material_properties.push_back(ReadDouble(node2, "A", 0.0));
  const tinyxml2::XMLElement* node3 = inputfile_.FirstChildElement()->FirstChildElement("fluids")->FirstChildElement(fluid.c_str())->FirstChildElement("B");
  material_properties.push_back(ReadDouble(node3, "B", 6.0e8));
  const tinyxml2::XMLElement* node4 = inputfile_.FirstChildElement()->FirstChildElement("fluids")->FirstChildElement(fluid.c_str())->FirstChildElement("rho0");
  material_properties.push_back(ReadDouble(node4, "rho0", 998));
  const tinyxml2::XMLElement* node5 = inputfile_.FirstChildElement()->FirstChildElement("fluids")->FirstChildElement(fluid.c_str())->FirstChildElement("viscosity")->FirstChildElement("dynamicShear");
  material_properties.push_back(ReadDouble(node5, "dynamicShear", 0.0));
  const tinyxml2::XMLElement* node6 = inputfile_.FirstChildElement()->FirstChildElement("fluids")->FirstChildElement(fluid.c_str())->FirstChildElement("viscosity")->FirstChildElement("dynamicBulk");
  material_properties.push_back(ReadDouble(node6, "dynamicBulk", 0.0));

  return material_properties;
}

/**
 * @brief Reads out the size of blocks on level zero.
 * @return Size of the blocks. Default: 1.
 */
double InputFileParser::ReadBlockSize() const{
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("blockSize");
  return ReadDouble(node, "Block Size", 1);
}

/**
 * @brief Reads out how many blocks on level zero will be present in x-direction.
 * @return Number of blocks in x-direction. Default: 1.
 */
int InputFileParser::ReadNumberOfBlocksX() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("blockRatio")->FirstChildElement("x");
  return ReadInt(node, "Block Ratio x", 1);
}

/**
 * @brief Reads out how many blocks on level zero will be present in y-direction.
 * @return Number of blocks in y-direction. Default: 1.
 */
int InputFileParser::ReadNumberOfBlocksY() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("blockRatio")->FirstChildElement("y");
  return ReadInt(node, "Block Ratio y", 1);
}

/**
 * @brief Reads out how many blocks on level zero will be present in z-direction.
 * @return Number of blocks in z-direction. Default: 1.
 */
int InputFileParser::ReadNumberOfBlocksZ() const{
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("blockRatio")->FirstChildElement("z");
  return ReadInt(node, "Block Ratio z", 1);
}

/**
 * @brief Reads out the type of the external boundary conditions and converts them into a machine-readable BoundaryType enum.
 * @return Boundary types at the respective edges of the domain. 0: West, 1: East, 2: South, 3: North, 4: Bottom, 5: Top.
 */
std::array<BoundaryType, 6> InputFileParser::ReadBoundaryConditions() const {
  std::array<BoundaryType, 6> boundaries;
  const tinyxml2::XMLElement* node0 = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("boundaryConditions")->FirstChildElement("fluid")->FirstChildElement("east");
  boundaries[0] = SelectBoundaryCondition(ReadString(node0));
  const tinyxml2::XMLElement* node1 = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("boundaryConditions")->FirstChildElement("fluid")->FirstChildElement("west");
  boundaries[1] =  SelectBoundaryCondition(ReadString(node1));
  const tinyxml2::XMLElement* node2 = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("boundaryConditions")->FirstChildElement("fluid")->FirstChildElement("north");
  boundaries[2] =  SelectBoundaryCondition(ReadString(node2));
  const tinyxml2::XMLElement* node3 = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("boundaryConditions")->FirstChildElement("fluid")->FirstChildElement("south");
  boundaries[3] = SelectBoundaryCondition(ReadString(node3));
  const tinyxml2::XMLElement* node4 = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("boundaryConditions")->FirstChildElement("fluid")->FirstChildElement("top");
  boundaries[4] =  SelectBoundaryCondition(ReadString(node4));
  const tinyxml2::XMLElement* node5 = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("boundaryConditions")->FirstChildElement("fluid")->FirstChildElement("bottom");
  boundaries[5] = SelectBoundaryCondition(ReadString(node5));

  return boundaries;
}

/**
 * @brief Reads out prime states specified in the fixed value boundary condition.
 * @return Boundary types at the respective edges of the domain. 0: West, 1: East, 2: South, 3: North, 4: Bottom, 5: Top.
 */
std::array<double,CC::NoEq()> InputFileParser::ReadFixedValueBoundaryCondition(const BoundaryLocation location) const {

    std::array<double,CC::NoEq()> fixed_values;
    for (unsigned int i = 0; i<CC::NoEq(); ++i) {fixed_values[i] = 0.0;}
    std::string location_name;

    //select name of BC that needs to be read in from inputfile
    if (location == BoundaryLocation::eEast){
        location_name = "valuesEast";
    } else if (location == BoundaryLocation::eWest){
        location_name = "valuesWest";
    } else if (location == BoundaryLocation::eNorth){
        location_name = "valuesNorth";
    } else if (location == BoundaryLocation::eSouth){
        location_name = "valuesSouth";
    } else if (location == BoundaryLocation::eTop){
        location_name = "valuesTop";
    } else if (location == BoundaryLocation::eBottom){
        location_name = "valuesBottom";
    }

    //read values from inputfile
    const tinyxml2::XMLElement* node0 = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("boundaryConditions")->FirstChildElement("fluid")->FirstChildElement(location_name.c_str())->FirstChildElement("density");
    fixed_values[0] = ReadDouble(node0, "fixedValueDensity", 0);
    const tinyxml2::XMLElement* node1 = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("boundaryConditions")->FirstChildElement("fluid")->FirstChildElement(location_name.c_str())->FirstChildElement("pressure");
    fixed_values[1] = ReadDouble(node1, "fixedValuePressure", 0);
    const tinyxml2::XMLElement* node2 = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("boundaryConditions")->FirstChildElement("fluid")->FirstChildElement(location_name.c_str())->FirstChildElement("velocityX");
    fixed_values[2] = ReadDouble(node2, "fixedValueVelocityX", 0);

    //read y-velocity only when required, otherwise leave 0 from initialization
    if (CC::DIM() != Dimension::One ) {
        const tinyxml2::XMLElement* node3 = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("boundaryConditions")->FirstChildElement("fluid")->FirstChildElement(location_name.c_str())->FirstChildElement("velocityY");
        fixed_values[3] = ReadDouble(node3, "fixedValueVelocityY", 0);
    }
    //read z-velocity only when required, otherwise leave 0 from initialization
    if (CC::DIM() == Dimension::Three ) {
        const tinyxml2::XMLElement* node4 = inputfile_.FirstChildElement()->FirstChildElement("domain")->FirstChildElement("boundaryConditions")->FirstChildElement("fluid")->FirstChildElement(location_name.c_str())->FirstChildElement("velocityZ");
        fixed_values[4] = ReadDouble(node4, "fixedValueVelocityZ", 0);
    }

  return fixed_values;
}

/**
 * @brief Reads out the initial condition of all fluids specifed in the inputfile.
 * @return Equation form of the initial condition of all fluids.
 */
std::vector<std::string>  InputFileParser::ReadInitialConditionOfFluids() const {

  std::vector<std::string> fluid_initialisation;
  std::string fluid;
  std::string base_name("fluid");

  for(int i = 1; i <= ReadNumberOfFluids(); ++i){
    fluid = base_name + std::to_string(i);
    fluid_initialisation.push_back(ReadInitialConditionString(fluid));
  }

  return fluid_initialisation;
}

/**
 * @brief Reads out the number of fluids present in the simulation.
 * @return Number of fluids. Default: 1.
 */
int InputFileParser::ReadNumberOfFluids() const {
    const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("fluids")->FirstChildElement("numberOfFluids");
    return ReadInt(node, "numberOfFluids", 1);
}

/**
 * @brief Reads out the material parameters for all fluids defined in the inputfile.
 * @return Vector of pairs defining a fluid by its material type, i.e. equation of state, and its properties, e.g. gamma, viscostiy, etc.
 */
std::vector<std::pair<MaterialName, std::vector<double>>>  InputFileParser::ReadParameterOfAllFluids() const {

  std::string           fluid;
  std::string base_name("fluid");
  std::vector<double> properties;
  std::vector<std::pair<MaterialName, std::vector<double>>>  fluids;

  for(int i = 1; i <= ReadNumberOfFluids(); ++i){
    fluid = base_name + std::to_string(i);
    properties = ReadMaterialProperties(fluid);

    //genrate pair out of data above
    std::pair<MaterialName, std::vector<double>> fluid_data_pair;
    fluid_data_pair = std::make_pair(ReadMaterialName(fluid),properties);

    // add to vector
    fluids.push_back(fluid_data_pair);
    }

  return fluids;
}

/**
 * @brief Reads out the gravitational pull as three dimensional array, one for the pull in each cartesian direction.
 * @return Gravity in x-, y- and z-direction.
 */
std::array<double, 3> InputFileParser::ReadGravity() const {
  std::array<double, 3> gravity;
  const tinyxml2::XMLElement* node0 = inputfile_.FirstChildElement()->FirstChildElement("sourceTerms")->FirstChildElement("gravity")->FirstChildElement("x");
  gravity[0] = ReadDouble(node0, "gravity x", 0);
  const tinyxml2::XMLElement* node1 = inputfile_.FirstChildElement()->FirstChildElement("sourceTerms")->FirstChildElement("gravity")->FirstChildElement("y");
  gravity[1] = ReadDouble(node1, "gravity y", 0);
  const tinyxml2::XMLElement* node2 = inputfile_.FirstChildElement()->FirstChildElement("sourceTerms")->FirstChildElement("gravity")->FirstChildElement("z");
  gravity[2] = ReadDouble(node2, "gravity z", 0);
  return gravity;
}

/**
 * @brief Reads out the maximum level (of refinement of the simulation).
 * @return Finest level. Default: 0.
 */
int InputFileParser::ReadMaximumLevel() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("multiResolution")->FirstChildElement("maximumLevel");
  return ReadInt(node, "maximum level", 0);
}

/**
 * @brief Reads-out epsilon_ref, which defines the threshold for coarsening of levels.
 *        Epsilon_ref is thereby not the finally used threshold directly, but is modified according to Roussel et al. (2003)- "A conservative
 *        fully adaptive multiresolution algorithm for parabolic PDEs".
 * @return Epsilon_ref Reference epsilon. Default: 0.05.
 */
double InputFileParser::ReadReferenceMultiResolutionEpsilon() const {
  const tinyxml2::XMLElement *node = inputfile_.FirstChildElement()->FirstChildElement("multiResolution")->FirstChildElement("refinementCriterion")->FirstChildElement("epsilonReference");
  return ReadDouble(node, "epsilonReference", 0.01);
}

/**
 * @brief Reads out the level on which the reference epsilon (see ReadReferenceMultiResolutionEpsilon()) is enforced.
 * @return Level on which reference epsilon is to be enforced. Default: 1.
 */
int InputFileParser::ReadLevelOfReferenceMultiResolutionEpsilon() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("multiResolution")->FirstChildElement("refinementCriterion")->FirstChildElement("levelOfEpsilonReference");
  return ReadInt(node, "levelOfEpsilonReference", 1);
}

/**
 * @brief Reads-out the time at which the simulation is supposed to start.
 * @return Start time. Default: 0.0.
 */
double InputFileParser::ReadStartTime() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("timeControl")->FirstChildElement("startTime");
  return ReadDouble(node, "startTime", 0.0);
}

/**
 * @brief Reads out the time at which the simulation is supposed to end. I. e. EndTime - StartTime = SimulatedTime.
 * @return End time. Default: 0.0.
 */
double InputFileParser::ReadEndTime() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("timeControl")->FirstChildElement("endTime");
  return ReadDouble(node, "endTime", 0.0);
}

/**
 * @brief Read-out the Courant–Friedrichs–Lewy number.
 * @return CFL number CFL number. Default 0.6.
 */
double InputFileParser::ReadCFLnumber() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("CFLNumber");
  return ReadDouble(node, "CFLNumber", 0.6);
}

/**
 * @brief Reads-out the reference length to (non-)dimensionalize the (inputs) outputs.
 * @return Reference length. Default: 1.0.
 */
double InputFileParser::ReadReferenceLength() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("referenceParameter")->FirstChildElement("lRef");
  return ReadDouble(node, "lRef", 1.0);
}

/**
 * @brief ReadReferenceVelocity Reads-out the reference velocity to (non-)dimensionalize the (inputs) outputs.
 * @return Reference veloctiy. Default: 1.0.
 */
double InputFileParser::ReadReferenceVelocity() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("referenceParameter")->FirstChildElement("uRef");
  return ReadDouble(node, "U_ref", 1.0);
}

/**
 * @brief Reads-out the reference density to (non-)dimensionalize the (inputs) outputs.
 * @return Reference density. Default: 1.0.
 */
double InputFileParser::ReadReferenceDensity() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("referenceParameter")->FirstChildElement("rhoRef");
  return ReadDouble(node, "Rho_ref", 1.0);
}

/**
 * @brief Reads-out the reference viscosity to (non-)dimensionalize the (inputs) outputs.
 * @return Reference Mu. Default: 1.0
 */
double InputFileParser::ReadReferenceMu() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("referenceParameter")->FirstChildElement("muRef");
  return ReadDouble(node, "Mu_ref", 1.0);
}

/**
 * @brief Reads-out the reference gravity to (non-)dimensionalize the (inputs) outputs.
 * @return Reference gravity. Default: 1.0
 */
double InputFileParser::ReadReferenceGravity() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("referenceParameter")->FirstChildElement("gRef");
  return ReadDouble(node, "g_ref", 1.0);
}

/**
 * @brief Reads-out the reference time to (non-)dimensionalize the (inputs) outputs.
 * @return Reference time. Default: 1.0
 */
double InputFileParser::ReadReferenceTime() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("calculation")->FirstChildElement("referenceParameter")->FirstChildElement("tRef");
  return ReadDouble(node, "t_ref", 1.0);
}

/**
 * @brief Reads out the format in which the output should be written.
 * @return The output format in machine-readable format. Default: ASCII.
 */
OutputType InputFileParser::ReadOutputFormat() const {
    const tinyxml2::XMLElement *node = inputfile_.FirstChildElement()->FirstChildElement("output")->FirstChildElement("outputFileType");
    std::string type = ReadString(node);
    // NH 2017-02-08: Switch Statements do not work with String
    if (type == "XDMF") {return OutputType::eXdmf;}
    else {return OutputType::eAscii;}
}

/**
 * @brief Reads out the way output times are determined. Inputs should be 'Interval' or 'Timestamps'.
 * @return Output times type.
 */
std::string InputFileParser::ReadOutputTimesType() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("output")->FirstChildElement("outputTimesType");
  return ReadString(node);
}

/**
 * @brief Reads out the output period. The period defines the time interval in between two outputs.
 * @return Intervall after which the next output is to be written. Default: 1.0.
 */
double InputFileParser::ReadOutputPeriod() const {
  const tinyxml2::XMLElement *node = inputfile_.FirstChildElement()->FirstChildElement("output")->FirstChildElement("outputPeriod");
  return ReadDouble(node, "outputPeriod", 1.0);
}

/**
 * @brief Reads out the number of timestamps.
 * @return Number of timestamps specifed.
 */
int InputFileParser::ReadNumberOfTimestamps() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("output")->FirstChildElement("numberOfTimestamps");
  return ReadInt(node, "numberOfTimestamps", 0);
}

/**
 * @brief Reads out the timestamps to be for outputs.
 * @return Timestamps.
 */
std::vector<double> InputFileParser::ReadTimestamps() const {

  std::string timestamp_name;
  std::string base_name("ts");

  std::vector<double> timestamps;

  for(int i = 1; i <= ReadNumberOfTimestamps(); ++i) {
    timestamp_name = base_name + std::to_string(i);
    const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("output")->FirstChildElement("timestamps")->FirstChildElement(timestamp_name.c_str());
    timestamps.emplace_back(ReadDouble(node,timestamp_name.c_str(),0.0));
  }

  return timestamps;
}

/**
 * @brief Reads out the detail to which the timestep is printed, i. e. number of decimals shown.
 * @return Number with the corresponding amount of decimals. Default: 1.
 */
double InputFileParser::ReadTimeNamingFactor() const {
  const tinyxml2::XMLElement* node = inputfile_.FirstChildElement()->FirstChildElement("output")->FirstChildElement("timeNamingFactor");
  return ReadDouble(node, "timeNamingFactor", 1.0);
}
