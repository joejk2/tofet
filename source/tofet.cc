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
//long SEED=-time(NULL);        // Unique
#include "RandomB.h"
#include "global.h"
#include "graph.h"
#include "hoppers.h"
#include "kmc.h"

int main(int argc, char * argv[]) {
    // DETERMINE INPUT FILES
    if(argc < 4) {				
		cout << "*** ERROR ***: Expect at least three input files: .sim, .xyz, .edge\n";
		exit(-1);
    }
    char sim[128], xyz[128], edge[128], occ[128];
    for (int i=1; i<argc; i++) {
		if (strstr(argv[i],".sim"))  strcpy(sim,argv[i]);	
		if (strstr(argv[i],".xyz"))  strcpy(xyz,argv[i]);	
		if (strstr(argv[i],".edge")) strcpy(edge,argv[i]);	
		if (strstr(argv[i],".occ")) strcpy(occ,argv[i]);	
    }
    cout << "Taking input from " << sim << ", " << xyz << ", " << edge << " ..." << endl;
    cout << "Read simulation parameters ..." << endl;
    PrintAll(sim);
    // Check for incompatibilities in sim file
    if ( Read(sim,"hopperInteractions","0")=="1" && Read(sim,"mode","tof")=="tof" ) {
        cout << "!!! ERROR !!!: hopperInteractions aren't currently implemented in the 'tof' mode\n";
        exit(-1);
    }
    // Determine verbosity of output
    if (Read(sim, "verbosity", "low") == "high") VERBOSITY_HIGH = true;

    // SETUP GSL RANDOM NUMBER GENERATOR (IF NEEDED)
    #ifndef RandomB
    cout << "Setting up GSL random number generator...\n";
    const gsl_rng_type * T;				
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gslRand = gsl_rng_alloc (T);
    printf ("\tGenerator type: %s\n", gsl_rng_name ( gslRand ));
    printf ("\tSeed = %lu\n", gsl_rng_default_seed);
    #else
    cout << "Setting up RandomB random number generator...\n"; 
    #endif

    // INITIALISE GRAPH
    if ( VERBOSITY_HIGH ) cout << "Initialising Graph...\n";		
    graph Graph(sim, xyz, edge);  

    // INITIALISE HOPPERS
    int totalHoppers=0;
	if ( VERBOSITY_HIGH ) cout << "Initialising Hoppers...\n";
	hoppers Hoppers(&Graph,sim);
    if ( Read(sim,"mode","tof")!="fet" ) { 
        totalHoppers = atoi( Read(sim,"hoppers").c_str() ); 
    }
    else {
        if (strstr(occ,".occ")) {
            totalHoppers = Hoppers.GenerateOnPreviouslyOccupied(occ, 0.0);
            cout << "Generated " << totalHoppers << " charges, as read from " << occ << endl;
        } 
        if (totalHoppers == 0) {
            totalHoppers = Hoppers.SetSourceDrainOccupation(0.0);
            cout << "No '.occ' file successfully read, so generated " << totalHoppers 
                 << " hoppers at random\n";
        }
    }

    // INITIALISE KMC
	if ( VERBOSITY_HIGH ) cout << "Initialising KMC...\n";
    kmc KMC(sim, &Hoppers, totalHoppers, &Graph);
    cout << "All systems go!  Beginning KMC...\n"
         << ".................................\n"
         << ".................................\n";

    // RUN FET SIMULATIONS
    if ( Read(sim,"mode","tof")=="fet" ) { 
		if ( VERBOSITY_HIGH ) cout << "Using algorithm KMC::FRM_FET()\n";
		KMC.FRM_FET(); 
	}

    // RUN TOF SIMULATIONS (DEFAULT)
	else {
        if ( VERBOSITY_HIGH ) cout << "Using algorithm KMC::FRM()\n";
        KMC.FRM();	
    }

    // OUTPUT
    cout << "Simulation finished with " << WARNINGS << " warnings\n" 
         << ".................................\n"
         << ".................................\n";
    if (Read(sim,"mode","tof")=="fet") {
        cout << "> TIME = " << KMC.GetTime() << endl
             << "> NUMBER OF HOPPERS LEFT = " << Hoppers.GetActive() << endl
             << "> TOTAL NUMBER OF HOPPERS COLLECTED AT DRAIN = " << Hoppers.GetCollectorCurrent() << endl
             << "> TOTAL NUMBER OF HOPPERS INJECTED AT SOURCE = " << Hoppers.GetCollectorCurrent() << endl
             << "> CURRENT (A) = " << Hoppers.GetFETCurrent() << endl
             << "> OCCUPIED MOLECULES AT END OF SIMULATION" << endl
             << "\tmolecule_ID\n";
        Hoppers.PrintOccupiedVertices();
    }
    // for regenerate and tof modes...
    else if (Read(sim, "mode", "tof")=="regenerate" || Read(sim, "mode", "tof")=="tof") { 
        cout << "> PHOTOCURRENT TRANSIENT\n"
			 << "\ttime (s)\t-\tcurrent (A)\n";
		KMC.PrintCurrent();
	    cout << "> TOTAL RUNS = " << KMC.GetnRuns() << endl
		     << "> MOBILITY FROM COLLECTION TIMES (cm^2/V.s)= " 
	         << Hoppers.GetSumReciprocalCollTimes() / (double (Hoppers.GetTotalCollectionEvents())) * 1e-16 
                * (Graph.GetDepth() / -Graph.GetFieldZ()) << endl
             << "> MOBILITY FROM TOTAL DISPLACEMENT AND TOTAL TIME (cm^2/V.s)= " 
             << KMC.GetSumDz() * 1e-16 / (KMC.GetTotalTimeOverAllRuns() * totalHoppers * -Graph.GetFieldZ()) 
             << endl 
             << "> AVERAGE NUMBER OF HOPPERS COLLECTED PER RUN = "
             << double(Hoppers.GetTotalCollectionEvents()) / KMC.GetnRuns() << endl 
             << "> PROBABILITY OF HOPPER BEING COLLECTED DURING RUN = "
	         << double(Hoppers.GetTotalCollectionEvents()) / (KMC.GetnRuns() * totalHoppers) << endl;
	}

    // for tofetOccupation simulations...
    #ifdef printTotalOccupation
    if (Read(sim,"printOccupation","0") == "1") {
        Graph.NormaliseOccupationTimes( KMC.GetTime(), totalHoppers );
        Graph.PrintTotalOccupationTimes();
    }
	if (Read(sim,"printEnergies","0") == "1") { 
        Graph.PrintEnergies();
    }
    #endif 

    #ifndef RandomB
    gsl_rng_free(gslRand);
    #endif

    return 0;
}
