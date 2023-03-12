#include "utilities/random_number_generator.h"

/**
 * @brief This function is used to get the Logger. If no Logger exists yet it is
 * created, otherwise the existing logger is passed back. "Singleton
 * Constructor"
 * @param save_all_ranks Decider whether or not all ranks are to be logged.
 * (Only relevant at first call!)
 * @return The logger instance.
 */
RandomNumberGenerator &RandomNumberGenerator::Instance(int const rank_id) {
  static RandomNumberGenerator instance_(rank_id);
  return instance_;
}