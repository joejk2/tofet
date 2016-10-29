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

/**********************************************************************
 * 'kmc' is the central kinetic Monte Carlo class.
 * At the moment it contains two 'First Reaction Method' algorithms 
 * to cater for time-of-flight and variants, and FETs. 
 *********************************************************************/

#ifndef _KMC_H
#define	_KMC_H
#include "hoppers.h"
#include "graph.h"

using namespace std;

class kmc{
    private:
        graph * _graph;
        hoppers * _Hoppers;
        double _time;   // the simulation time 
        double _totalTimeOverAllRuns;   // the total time, summed over all runs 
        int _geometricBin;  // time bins for photocurrent transients (tof)
        int _oldGeometricBin;
        double _maxTime; 		
        int _nRuns;  // the number of runs in each simulation
        double _maxRuns;  // set to double so that inf can be represented
        double _tol;  // tolerance of results of simulation
        double _lowerTol, _upperTol;
        double _dt;  // width of first time bin
        double _logDt;				
        double _alpha;  // subsequent log time bins are dt * [(alpha ^ n) - (alpha ^ (n-1))] wide
        double _logAlpha;
        double _sum_dz;  // the total distance moved along z, summed over all hoppers
        int _nLogTimeBins;  // number of geometric time bins
        vector <double> _current;  // photocurrent
        int _nHoppers;  // initial number of hoppers.  NOTE: this is not updated as hoppers are collected
        string _mode;  // mode of simulation (FET, tof, regenerate...)
        bool _hopperInteractions;  // Coulombic interactions?
    
        void UpdatePhotocurrent(const double & ); 
        double (hoppers::*moveFastest)();  // pointer to appropriate MoveFastest_* function 
    //end of private:

    public:
        kmc(){}
        kmc(char * sim, hoppers * Hoppers, int totalHoppers, graph * Graph ){
            _totalTimeOverAllRuns = 0.0;
            _sum_dz = 0.0;
            _graph = Graph;
            _Hoppers = Hoppers;
            _maxTime=atof(Read(sim,"maxTime").c_str());
            _mode = Read(sim, "mode", "tof");
            if (Read(sim, "hopperInteractions", "0") == "1") {
                _hopperInteractions = true;
            }
            else {
                _hopperInteractions = false;
            }
            if (_mode == "tof" || _mode == "regenerate") {
                _dt=atof(Read(sim,"deltaTime").c_str());
                if (_dt > _maxTime) {
                    cout << "*** ERROR *** : deltaTime > maxTime!\n";
                    exit(-1);
                }
                _alpha=atof(Read(sim,"alpha").c_str());
                _nHoppers=totalHoppers;
                _tol=atof(Read(sim,"tol","0.0").c_str());
                _upperTol=1.0+_tol;
                _lowerTol=1.0-_tol;
                _maxRuns=atof(Read(sim,"maxRuns","inf").c_str()); 
                _logAlpha = log(_alpha);
                _logDt = log(_dt);
                _nLogTimeBins = int ( (log(_maxTime)  - _logDt ) / _logAlpha );
                _current.resize(_nLogTimeBins);
            }
            if (_mode=="tof") { 
                if (VERBOSITY_HIGH) {
                    cout << "Setting to MoveFastest_C" << endl;
                }
                moveFastest=&hoppers::MoveFastest_C;
            }
            if (_mode=="regenerate") {
                if (_hopperInteractions) { 
                    if (VERBOSITY_HIGH) {
                        cout << "Setting moveFastest to 'MoveFastest_RCI'" << endl;
                    }
                    moveFastest=&hoppers::MoveFastest_RCI;
                }
                else {
                    if (VERBOSITY_HIGH) {
                        cout << "Setting moveFastest to 'MoveFastest_R'" << endl;
                    }
                    moveFastest=&hoppers::MoveFastest_R;
                }
            }
            if (_mode=="fet") {
                if ( VERBOSITY_HIGH ) {
                    cout << "Setting moveFasest to 'MoveFastest_F'" << endl;
                }
                moveFastest=&hoppers::MoveFastest_F; 
                if (Read(sim,"converged","0")=="1") {
                    _Hoppers->SetActiveHoppersConverged();
                    cout << "Assuming the charge density is already converged\n";
                }
            }
        }
        ~kmc(){
            _current.clear();
        }

        /***************************************************
         * DO'S
         * TODO: These functions could probably be merged...
         **************************************************/
        void FRM();  // simple First Reaction Method	
        void FRM_FET();  // FRM with all necessary add-ons for FET simulations
        
        /***************************************************
         * GET'S
         **************************************************/
        const double & GetSumDz() const {return _sum_dz;}
        const double & GetTotalTimeOverAllRuns() const {return _totalTimeOverAllRuns;}
        const double & GetDt() const	{return _dt;}
        const double & GetTime() const 	{return _time;}
        const double & GetAlpha() const	{return _alpha;}
        vector <double> & GetTimeBins()	{return _current;}
        const int & GetnRuns() const	{return _nRuns;}
        void PrintCurrent(string dest="") {
            double t1, t2, dt, depth;
            // If necessary, redirect 'cout' to 'fout'
            streambuf* cout_sbuf = std::cout.rdbuf();
            ofstream   fout;
            if(dest=="file") {
                fout.open("occVert.out");
                cout.rdbuf(fout.rdbuf());
            }
            for ( int i = 0; i < _geometricBin - 1; i++){ 
                    if (i == 0 ) {
                        t1 = 0;
                    }
                    else {
                        t1 = _dt * pow(_alpha, i-1);
                    }
                    t2 = _dt * pow(_alpha, i);
                    dt =  t2 - t1;
                    depth = (*_graph).GetDepth();
                    if (_current.at(i) != 0) {
                        cout << '\t' << t2 << "\t\t" << e * _current.at(i) / (dt*depth) <<endl;
                    }
            }
            if (dest=="file") {
                fout.close();
                // Restore the original stream buffer  
                cout.rdbuf(cout_sbuf);
            }
        }
    // end of public:
};
#endif	/* _KMC_H */

