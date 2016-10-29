///////////////////////////////////////////////////////////////////////
//  This file is part of ToFeT.
//  
//  ToFeT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//  
//  ToFeT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//  
//  You should have received a copy of the GNU Lesser General Public License
//  along with ToFeT.  If not, see <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////
#include "RandomB.h"
long SEED = -1;                 // Same sequence
//long SEED=-time(NULL);        // Unique sequence
long *p_SEED = &SEED;           // Declaring and initialising pointer to seed


/*============================================================================================
JK: keep the functions etc, kill the class wrap

  RandomB class
  ------------

  This version of the random number class generates various types of random number using
  routines obtained from "Numerical Recipes in C"

  Constructor:	RandomB()

  Destructor:	~RandomB()

  Methods:		double Rndm(long *idum)
					Generates a uniform random number using a SEED
				double Uniform()
					Returns a random number between 0 and 1 (not inclusive) 
				double RandLog()
					Returns -logarithm of a uniform random number
				double RandConstLog(double constant)
					Returns -log(1 -(uniform()/constant))
				int RandPos(int n)
					Returns a random integer between 0 and n-1
				int Dice()
					Returns a random integer between 1 and 6

  Accessor methods:

  Created:		22/10/02 by R Chandler

  Modified:		25/02/03 by R Chandler
				Variable and method names changed so that methods are capitalised and
				variables start with lower case

============================================================================================*/



// Random number generation using SEED
// -----------------------------------
double Rndm(long *idum)
/*-------------------------------------------------------------------------------------------------
"Minimal" random number generator of Park and Miller with Bays-Durham shuffle and added safeguards. 

Inputs:	long *idum		'SEED' for generator
 
Returns:	A uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values).

Notes:		Call with idum a negative integer to initialise; thereafter do not alter idum between 
			successive deviates in a sequence.  RNMX should approximate the largest floating value 
			that is less than 1.

 From:		Based around "Numerical Recipes in C - The Art of Scientific Computing, Second Edition".
			W H Press, S A Tekkolsky, W T Vetterling, B P Flannery.

Created:	07/10/02	R Chandler

Modified:
---------------------------------------------------------------------------------------------------*/
{
// Constant local variables
	const long IA = 16807;
	const long IM = 2147483647;
	const double AM = (1.0 / IM);
	const long IQ = 127773;
	const long IR = 2836;
	const int NTAB = 32;
	const double NDIV = (1 + (IM - 1) / NTAB);
	const double EPS = 1.2e-7;
	const double RNMX = 1.0-EPS;
	
// Other local variables
	int j;
	long k;
	static long iy = 0;
	static long iv[NTAB];
	double temp;

// This section executes if *idum is a negative integer.  It serves to initialise the random 
	//number generator
	if (*idum <= 0 || !iy) 
	{
		if (-(*idum) < 1) *idum = 1;		// Prevents idum = 0
		else *idum = -(*idum);
		for (j = NTAB + 7; j >= 0; j--)		// Load shuffle table after 8 warm-ups)
		{
			k = (*idum) / IQ;
			*idum = IA * (*idum - k * IQ) - IR * k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}

// Execution starts here when not initialising SEED
	k = (*idum) / IQ;
	*idum = IA * (*idum - k * IQ) - IR * k; // Compute idum=(IA*idum) % IM without overflows by
												// Schrage's method
	if (*idum < 0) *idum += IM;				// Will be in the range 0..NTAB-1
	j = int(iy / NDIV);						
	iy = iv[j];								// Output previously stored value and refill shuffle table
	iv[j] = *idum;
	if ((temp = AM*iy) > RNMX) return RNMX;	// Because users don't expect endpoint values
	else return temp;
}


// Uniform random number
// ---------------------
// Returns uniform random number between 0 and 1
double Uniform()
{
	double x = 0.0;
	x = Rndm(p_SEED);
	return x;
}

double UniformPos()
{
    double x=0.0;
    while (x==0) {
	x = Uniform();
    }
    return x;
}


// - logarithm
// -----------
// Returns -logarithm of a uniform random number 
double RandLog()
{
	double rnd = -log(Uniform());
	if (rnd == 0) rnd = 1;
	return rnd;
} // End randlog method


// -log(const * uniform)
// ---------------------
double RandConstLog(double constant)
{
	double rnd = -log(1 - (Uniform() / constant));
	if (rnd == 0) rnd = 1;
	return rnd;
} // End randConstTimesLog method


// Random positive integer
// -----------------------
// Returns a random positive integer between 0 and n-1
int RandPos(int n)
{
	int number = int(n*Uniform());
	if (number < 0) return 0;			// number should never be < 0 but check just incase
	else if (number >= n) return n-1;   // Ensure number is less than n
	else return number;
} // End randpos method


// Random positive double
// ======================
// Returns a random positive double between 0 and maxSize
double RandDouble(double maxSize)
{
	double number = maxSize*Uniform();
	if (number < 0) return 0;
	else if (number > maxSize) return maxSize;
	else return number;
} // End of RandDouble method


// Simulates dice roll
// -------------------
// Returns a random integer between 1 and 6
int Dice()
{
	int number = RandPos(6) + 1;
	return number;
}
