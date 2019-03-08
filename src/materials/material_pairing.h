#ifndef MATERIAL_PAIRING_H
#define MATERIAL_PAIRING_H

#include "unit_handler.h"

/**
 * @brief The MaterialPairing class defines an interface for all pairing properties used in the simulation (e.g. surface tension coefficient). 
 *        Furthermore, it stores parameter models for properties that allow such computations. The MaterialPairing class does not manipulate any data.
 */
class MaterialPairing {

protected:

   // material pairing properties (fixed values). 
   double const surface_tension_coefficient_;

public:
   explicit MaterialPairing( double const surface_tension_coefficient, 
                             UnitHandler const& unit_handler );
   explicit MaterialPairing();
   virtual ~MaterialPairing() = default;
   MaterialPairing( MaterialPairing const& ) = delete;
   MaterialPairing& operator=( MaterialPairing const& ) = delete;
   // default move constructor (in case pointers to a parameter model are used the move constructor must be self written (see materials))
   MaterialPairing( MaterialPairing&& ) = default;
   MaterialPairing& operator=( MaterialPairing&& ) = delete;

   // return functions of member variables
   double GetSurfaceTensionCoefficient() const;
};

#endif //MATERIAL_PAIRING_H