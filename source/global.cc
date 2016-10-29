#ifndef RandomB
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
gsl_rng * gslRand;
#endif

bool VERBOSITY_HIGH = false;
int WARNINGS = 0;

template<typename T> void safe_delete(T& obj) {
    typename T::iterator it =obj.begin();
    for( ; it != obj.end() ; ++it ) {
        delete *it;
    }
    obj.clear();
}

