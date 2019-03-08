#ifndef RANDOM_NUMBER_GENERATOR_H
#define RANDOM_NUMBER_GENERATOR_H

#include <random>

class RandomNumberGenerator {


    std::mt19937_64 gen_; //Standard mersenne_twister_engine
    std::uniform_real_distribution<> dis_;

    //Singleton has only privat Constructor
    explicit RandomNumberGenerator(int const rank_id) :
    gen_(rank_id),
    dis_(-1.0, +1.0){
    }


public:
    //Singelton "Constructor":
    static RandomNumberGenerator& Instance(int const rank_id = 0);

    //Singeltons may never call these methods.
    RandomNumberGenerator(const RandomNumberGenerator&) = delete;
    void operator=(const RandomNumberGenerator&) = delete;

    double GiveRandomNumber() {
       return dis_(gen_);
    }
};


#endif //RANDOM_NUMBER_GENERATOR_H
