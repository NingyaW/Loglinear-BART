#ifndef GIG_H
#define GIG_H

#include <random>

using std::mt19937_64;

class Random {
private:
  mt19937_64 generator;

public:
  Random(unsigned int seed) : generator(seed) {}
  Random(mt19937_64 &generator) : generator(generator) {}
  double gig(double lambda, double chi, double psi);
};

#endif
