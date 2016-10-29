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
 * 'hoppers' is the containing class for all hoppers, and all their 
 * associated functions.  It is responsible for deciding which hopper 
 * should be moved first, and executing that hop.
 * Most comments are given in 'hoppers.cc'.
 *********************************************************************/

#ifndef _HOPPERS_H
#define	_HOPPERS_H
#include "graph.h"
#include "vertex.h"
#include "hopper.h"
#include "global.h"
#include "vec.h"

using namespace std;

class hoppers{
    private:
        int _nHoppers;  // number of active hoppers
        list <hopper *> _hoppers; 
        vector <double > _reciprocalCollectionTimes; 
        double _totalReciprocalCollectionTimes;
        list <hopper *>::iterator _fastest;  // hopper with most imminent hop time
        double _fastestTime;  // time of most imminent hop
        map < vertex *, hopper * >  _mapVertexToHopper;
        graph * _graph;
        int _printOccupation;  // track occupation of vertices?	
        bool _track;  // track the movement of charges?
        int _hopperInteractions;  // Coulombic interactions?
        // These are just used in FET simulations
        vector <vertex *> _generators;  
        vector <vertex *> _collectors; 
        int _collectorCurrent;  // drain current
        int _generatorCurrent;  // source current
        int _totalCurrent;  
        deque <double> _currentStore;  // store current (not geometric time bins!)
        int _moves;  // number of MC moves
        int _movesCycle;  // check convergence, update current every '_movesCycle'
        int _cyclesForConvergence;  //  current must be stable over this many cycles for convergence
        double _maxTime;
        bool _activeHoppersConverged;  
        double _activeHoppersConvergedTime;
        double _tol, _lowerTol, _upperTol;  // tolerance of FET simulation result	

        void SetMap(){
             _mapVertexToHopper.clear();
             list <hopper *>::iterator it_hop;
             for (it_hop = _hoppers.begin(); it_hop != _hoppers.end(); ++it_hop){
                 vertex * from = (*it_hop)->GetFrom();
                _mapVertexToHopper[from] = *it_hop;
             }
         } 
    // end of private:
    
    public:
        bool _run;  // run FET simulations whilst(_run).  UGLY!	
        hoppers(){}
        hoppers(graph * Graph, char * sim){
            _activeHoppersConverged=false;
            _printOccupation=atoi(Read(sim, "printOccupation","0.0").c_str());
            _hopperInteractions =atoi(Read(sim, "hopperInteractions", "0.0").c_str());
            _graph = Graph;
            _nHoppers=0;
            _generatorCurrent=0;
            _collectorCurrent=0;
            _totalReciprocalCollectionTimes=0.0;
            if (Read(sim, "track", "0") == "1") {
                _track = true;
            }
            else {
                _track = false;
            }
            if (Read(sim, "mode","tof")=="fet") {
                _generators= _graph->GetGenerators();
                _collectors= _graph->GetCollectors();
                if (VERBOSITY_HIGH) {
                    cout << "Initialised " << _generators.size() 
                         << " generators and " << _collectors.size() << " collectors\n";
                }
                _currentStore.assign(15, 0);
                _maxTime=atof(Read(sim,"maxTime").c_str());
                _activeHoppersConvergedTime=0.0;
                _activeHoppersConverged=false;
                _run=true;
                _tol=atof(Read(sim,"tol","0.0").c_str());
                _upperTol=1.0+_tol;
                _lowerTol=1.0-_tol;
                _moves=0;
                _movesCycle = (int) atof(Read(sim,"movesCycle","2e4").c_str());
                _cyclesForConvergence = atoi(Read(sim,"cyclesForConvergence","15").c_str());
                if (Read(sim,"converged","0")=="1") {
                    _activeHoppersConverged=true;
                    _activeHoppersConvergedTime=0.0;
                }
            }	
            _fastestTime=0.0;
        }
        ~hoppers(){
            softClear();
        }
        void softClear() {
            map < vertex *, hopper *>::iterator it_map = _mapVertexToHopper.begin();
            for (; it_map != _mapVertexToHopper.end(); ++it_map){
                delete (it_map -> second );            
            }
            _hoppers.clear();
            _nHoppers=0;
            _mapVertexToHopper.clear();
            if (_hoppers.size() != 0 || _mapVertexToHopper.size() != 0 ) {
                cerr << "*** ERROR *** : softClear failed in hoppers.h\n"; 
                exit(-1);
            }
        }
            
        /***********************************
        * DO'S
        ************************************/
        void Generate(vertex * , const double &);
        int GenerateOnPreviouslyOccupied(char *, const double &); 
        int GenerateOccProb(const double & time);
        void GenerateAll(const int & nHoppers, const double & time);	
        void GenerateRandom_F(const int & nHoppers, const double & time);
        int SetSourceDrainOccupation(const double & time);
        void Remove(list <hopper *>::iterator, const double & );	 	
        double Move(list <hopper *>::iterator , vertex *, double & );
        double MoveFastest_C() ;						
        double MoveFastest_RCI() ;					
        double MoveFastest_R() ;				
        double MoveFastest_F() ;			
        void FindFastest();				
        void SetWaitTimes(double time);
        void SetActiveHoppersConverged() {_activeHoppersConverged=true;}
        void FETConvergence();
        void activeHoppersConvergence();
        void SetHops_C(const double &);	
        void AddCoulomb(vertex *, int sign=1 );	
        void UpdateCoulomb_all(vertex *,int);
        void UpdateCoulomb_single(vertex *, vertex *, int);
        void DeleteCoulomb(vertex *);		
        double const GetSingleCoulombEnergy(vertex *, vertex *);
        double GetAllCoulombEnergies(vertex *, vertex *);

        /***********************************
         * GET'S
        ************************************/
        double GetFETCurrent()  {return _currentStore.back();}
        hopper * GetHopper(vertex * v);
        list <hopper *>::iterator GetHopperIterator(vertex * v);
        int GetHopperNumber(vertex * v);					
        const int GetActive() const  {return _nHoppers;}
        const double & GetFastestTime () const	{return (*_fastest)->GetWaitTime();}
        vec GetFastestPos  () const  {return (*_fastest)->GetFrom()->GetPos();}
        double  GetFastestZ  () const  {return (*_fastest)->GetFrom()->GetZ();}
        const double  & GetFastestDz  () const  {return (*_fastest)->GetDz();}
        double GetSumReciprocalCollTimes()  {return _totalReciprocalCollectionTimes;}
        unsigned int GetTotalCollectionEvents()  {return _reciprocalCollectionTimes.size();}
        void PrintOccupiedVertices(string dest="");
        int GetCollectorCurrent()  {return _collectorCurrent;}
        int GetGeneratorCurrent()  {return _generatorCurrent;}
    // end of public:
};
#endif	/* _HOPPER_H */

