#ifndef STATE_RECONSTRUCTION_SETTINGS_H
#define STATE_RECONSTRUCTION_SETTINGS_H

enum class StateReconstructionType { Conservative,
                                     Primitive,
                                     RoeCharacteristic };

constexpr StateReconstructionType state_reconstruction_type = StateReconstructionType::RoeCharacteristic;

/**
    * @brief Provides a string representation of the state reconstruction type.
    * @param reconstruction State reconstruction set to be stringified.
    * @return String of the given state reconstruction type.
    */
inline std::string SetToString( StateReconstructionType const reconstruction ) {
   switch( reconstruction ) {
      case StateReconstructionType::Conservative:
         return "Conservative";
      case StateReconstructionType::Primitive:
         return "Primitive";
      case StateReconstructionType::RoeCharacteristic:
         return "RoeCharacteristic";
      default:
         return "ERROR: This reconstruction type is not (yet) defined!";
   }
}

#endif// STATE_RECONSTRUCTION_SETTINGS_H
