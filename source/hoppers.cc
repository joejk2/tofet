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
#include "hoppers.h"
    
/*******************
 * OUTPUT FUNCTIONS
 *******************/
// Called from the FET mode...
void hoppers::PrintOccupiedVertices(string dest) {
    // If necessary, redirect 'cout' to 'fout'
    streambuf* cout_sbuf = std::cout.rdbuf();
    ofstream   fout;
    if(dest=="file") {
        fout.open("occVert.out");
        cout.rdbuf(fout.rdbuf()); 
    }
    map <vertex *, hopper *> ::iterator it_vert = _mapVertexToHopper.begin();
    for (; it_vert!=_mapVertexToHopper.end(); ++it_vert) {                      
        cout << "\t" <<  _graph->GetVertex(it_vert->first) << endl;
    }
    if (dest=="file") {
        fout.close();
        // Restore the original stream buffer  
        cout.rdbuf(cout_sbuf);
    }
}

/***************************************************************************
 * COULOMBIC INTERACTIONS 
 * When a hopper is created or deleted, then need to add or remove
 *     it's contribution to the Coulombic forces felt by all other hoppers.
 *   Every time a hopper is moved, it is removed, and generated.
 *   All Coulomb energies are either updated by:
 *     1) AddCoulomb, or
 *     2) DeleteCoulomb
 *
 *  TODO:  There are multiple, and substantial, improvements that can 
 *         be made to improve the speed and memory usage of these functions.
 *         However, I haven't got round to these because this section has 
 *         never yet been the bottleneck in my simulations.
 ***************************************************************************/
// Given a 'newlyOccupied' vertex, update all the necessary DC's
void hoppers::AddCoulomb(vertex * newlyOccupied, int sign) {
    map <vertex *, hopper *> ::iterator it_vert = _mapVertexToHopper.begin();
    for (; it_vert!=_mapVertexToHopper.end(); ++it_vert) {		 		
        // For the hopper that has just been added, need to calculate 
        //   Coulombic interactions with *all* other hoppers:
    	if ( it_vert->first == newlyOccupied )	{  
            UpdateCoulomb_all(newlyOccupied, sign); 
        }
        // For other hoppers, only need to update the Coulombic 
        //   energy with the contribution from the most recently 
        //   added hopper
        else { 
            UpdateCoulomb_single(it_vert->first, newlyOccupied, sign); 
        }
    }
}
// Given a 'newlyUnoccupied' vertex, update all the necessary DC's
void hoppers::DeleteCoulomb(vertex * newlyUnoccupied) {
    AddCoulomb(newlyUnoccupied,-1);
    // If you're worried that something untoward is happening, uncomment
    //   this to check that DC for each hop is always reset to zero for 'newlyOccupied'
    /*for (int i=0; i < int (newlyUnoccupied->GetNumberNeighbours()); i++) {
        if (newlyUnoccupied->GetDC(i) > 1e-10) {
            cout << "Oops!  Expect DCs to be zero for unoccupied vertices\n"; 
            cout << "DC = " << newlyUnoccupied->GetDC(i) << " at "; newlyUnoccupied->PrintPos(); cout << endl;
            exit(-1);
        }
    }*/
}
// For a hopper on 'newlyOccupied', calculate the Coulombic interactions with 
//   *all* other hoppers
void hoppers::UpdateCoulomb_all(vertex * newlyOccupied, int sign) {
    if (sign==-1) {  // deleting a hopper...
        newlyOccupied -> ClearDCs(); 
        newlyOccupied -> SetEC(0.0, _fastestTime);
    }
    else {  // adding a hopper...
        double deltaCurrentCoulomb, deltaNeighbourCoulomb;
        deltaCurrentCoulomb = GetAllCoulombEnergies(newlyOccupied, newlyOccupied); 
        #ifdef printTotalOccupation
        newlyOccupied->SetEC(deltaCurrentCoulomb,_fastestTime);
        #endif
        // Update energetics for all reactions from 'newlyOccupied'
        vector <vertex *> neighbours = newlyOccupied -> GetNeighbours();
        vector <vertex *>::iterator neighbour = neighbours.begin();
        for (int i=0; neighbour!=neighbours.end(); neighbour++,i++) {  			
            deltaNeighbourCoulomb = GetAllCoulombEnergies(newlyOccupied, *neighbour);
            newlyOccupied -> IncrementDCs(i, (deltaNeighbourCoulomb - deltaCurrentCoulomb));
        }
    }
}
// Given a new hopper on 'newlyOccupied', update the Coulombic interactions of
//   all other hoppers.
void hoppers::UpdateCoulomb_single(vertex * interacting, vertex * newlyOccupied, int sign){
    double deltaCurrentCoulomb, deltaNeighbourCoulomb; 
    deltaCurrentCoulomb = GetSingleCoulombEnergy(interacting, newlyOccupied);
    #ifdef printTotalOccupation
	interacting->IncrementEC(sign*deltaCurrentCoulomb, _fastestTime);
    #endif
    // Update energetics for all reactions from 'interacting'
    vector <vertex *> neighbours = interacting -> GetNeighbours();
    vector <vertex *>::iterator neighbour = neighbours.begin();
    for (int i=0; neighbour!=neighbours.end(); neighbour++,i++) {  			
        deltaNeighbourCoulomb = GetSingleCoulombEnergy(*neighbour, newlyOccupied);
        interacting -> IncrementDCs(i, sign*(deltaNeighbourCoulomb - deltaCurrentCoulomb) );
    }
}
// Get the Coulomb energy between two vertices
double const hoppers::GetSingleCoulombEnergy(vertex * v1, vertex * v2) {
    // If you've got a very small grid, you might prefer to use a look-up table.
    //   However, I found that for ~7000 molecules that the 10% speed-up wasn't 
    //   worth the massive increase in memory costs
    // If you use this, you'll also have to uncomment the function 
    //   'MakeCoulombEnergyGrid()' in 'graph.h'
    //return _graph->GetCoulomb( v1, v2 );			

    // On the fly calculations could have a cutoff built in, but I've not dealt with 
    //   large enough morphologies yet for this to be worthwhile.
    if (v1 != v2 ) { // don't think we need this check....?
        return _graph->_coulombPrefactor/_graph->GetDistance(v1, v2);  
    }	
    else return 0.0;

    /** TODO:  Use trees to do this *really* fast... **/
}
// Get the Coulomb energy between 'interacting' and every other occupied vertex except 'ignore'
double hoppers::GetAllCoulombEnergies(vertex * ignore, vertex * interacting) {
    double coulomb=0;
    map <vertex *, hopper *> ::iterator occupied = _mapVertexToHopper.begin();
    for (; occupied!=_mapVertexToHopper.end(); ++occupied) { 						
        if ( occupied->first != ignore ) {  // ignore interactions with self...
            coulomb += GetSingleCoulombEnergy(interacting, occupied->first);	
        }
    }
    return coulomb;
}
// Once all the Coulomb energies have been updated, need to 
//   recalculate rates and reset all hops
void hoppers::SetHops_C(const double & fastestTime) {
    // C++ doesn't necessarily iterate through maps in the same order.
    //   Rather, it depends on the size of the objects in the map.
    //   If you're worried about this, use the following code instead...
    /*list <hopper *>::iterator it_hop = _hoppers.begin();
    for (; it_hop!=_hoppers.end(); ++it_hop) {
        ( (*it_hop)->GetFrom() ) -> UpdateRates_C(_graph->_reorg, _graph->_kT);
	    (*it_hop) -> SetHop(fastestTime);
    }*/

    // If you're not worried about the order of iteration, use this.
    //   Note: If you change the size of your hoppers / vertex object
    //         you may change the order of iteration and therefore your 
    //         results...
    map <vertex *, hopper *> ::iterator it_vert = _mapVertexToHopper.begin();
    for (; it_vert!=_mapVertexToHopper.end(); ++it_vert) { 	
        (it_vert->first)  -> UpdateRates_C(_graph->_reorg, _graph->_kT);	
        (it_vert->second) -> SetHop(fastestTime);
    }
}

/********************************
 * GENERATION / REMOVAL FUNCTIONS
 ********************************/ 
// Generate on a given vertex at a given time
void hoppers::Generate(vertex * V, const double & time){
    hopper * newhopper;
    newhopper = new hopper(V,time);
    _hoppers.push_back(newhopper);
    _mapVertexToHopper[V]=newhopper;
    _nHoppers++;
    if (_hopperInteractions) {
        AddCoulomb(V);
    }
    else {
        newhopper->SetHop(time);
    }
}
// Generate on previously occupied vertices
int hoppers::GenerateOnPreviouslyOccupied(char * filename, const double & time) {
    vector <vertex *> generateOnMe;
    generateOnMe = _graph->GetPreviouslyOccupied(filename);
    vector <vertex *>::iterator it = generateOnMe.begin();
    for (; it!=generateOnMe.end(); ++it) {
        Generate( (*it), time );
    }
    return generateOnMe.size();
}
// Generate randomly on unoccupied vertices that are electrodes
void hoppers::GenerateAll(const int & nHoppers, const double & time) {
    vertex * emptyGenerator;
    for (int i=0; i<nHoppers; i++) {
        emptyGenerator = _graph -> GetEmptyGenerator();
        Generate(emptyGenerator, time);
    }
}
// Remove hopper 'H' at 'time'
void hoppers::Remove(list <hopper *>::iterator H, const double & time){
    vertex * from=(*H)->GetFrom();
    if (_hopperInteractions) {
        DeleteCoulomb(from);
    }
    _mapVertexToHopper.erase(from);
    (*H)->SetWaitTime(time); 	
    delete *H;
    _hoppers.erase(H);		
    _nHoppers--;
}
// Set the Fermi-level of the source and drain in FETs. 
//   Called at every MC step
//   TODO: This could be a lot nicer....
int hoppers::SetSourceDrainOccupation(const double & time) {
    double energy;
    list <hopper *>::iterator it_hop;
    vector <vertex *>::iterator it;
    // Set occupation of the generators (the source)
    for (it=_generators.begin(); it!=_generators.end(); ++it) {
        energy = (*it)->GetE() + GetAllCoulombEnergies(*it, *it); 

        #ifdef RandomB
        if ( exp( (_graph->_sourceFermiEnergy-energy)/_graph->_kT ) > Uniform() ) { // vertex should be occupied...
        #else
        if ( exp( (_graph->_sourceFermiEnergy-energy)/_graph->_kT ) > gsl_rng_uniform(gslRand) ) { 
        #endif
            if ( !(*it)->IsOccupied() ) { 
                Generate(*it, time); 
                _generatorCurrent++;
            }
        }
        else {  // vertex shouldn't be occupied...
            if ( (*it)->IsOccupied() ) {	
                Remove(GetHopperIterator(*it),time);  // BODGE: This ain't fast!
                _generatorCurrent--;
            }
        }
    }
    // Do the same thing for the collectors (the drain)
    for (it=_collectors.begin(); it!=_collectors.end(); ++it) {
        energy = (*it)->GetE() + GetAllCoulombEnergies(*it, *it);
        #ifdef RandomB
        if ( exp( (_graph->_drainFermiEnergy - energy)/_graph->_kT ) > Uniform() ) {
        #else
        if ( exp( (_graph->_drainFermiEnergy - energy)/_graph->_kT ) > gsl_rng_uniform(gslRand) ) {
        #endif
            if ( !(*it)->IsOccupied() ) {
                Generate(*it, time); 
                _collectorCurrent--;
            }
        }
        else {
            if ( (*it)->IsOccupied() ) {
                Remove(GetHopperIterator(*it),time); 
                _collectorCurrent++;
            }
        }
    }
    return _nHoppers;
}

/*************************************************
 * DETERMINE FET CONVERGENCE 
 * Do two things to test FET convergence:
 *   1) Test that the number of hoppers is steady
 *   2) Test that the current is steady
 *************************************************/
void hoppers::activeHoppersConvergence() {
    if (VERBOSITY_HIGH) {
        cout << "Hoppers not converged: time = " << GetFastestTime() << "; hoppers = " << _nHoppers 
             << "; hoppers collected/injected = " << _collectorCurrent 
             << "/" << _generatorCurrent << endl;
    }
    if (_generatorCurrent>10) {
        if ( _generatorCurrent >= 0.7*_collectorCurrent && _generatorCurrent <= 1.3*_collectorCurrent ) {
            _activeHoppersConvergedTime=GetFastestTime();
            _activeHoppersConverged=true;
            cout << "Charge density converged at time = " << _activeHoppersConvergedTime 
                 << "\n\tNumber of charges = " << _nHoppers
                 << "\n\tCharges collected / injected since last reset = " << _collectorCurrent << " / " << _generatorCurrent << endl;
        }
        else {
            _generatorCurrent=0;
            _collectorCurrent=0;
            if (VERBOSITY_HIGH) {
                cout << "Hoppers not converged: resetting collection/injection events. Time = " 
                     << GetFastestTime() << endl; 
            }
        }
    }
}
void hoppers::FETConvergence() {
    // Update the current
    _currentStore.pop_front();
    _currentStore.push_back(e * (double(_collectorCurrent + _generatorCurrent))
                              / (2.0 * (GetFastestTime() - _activeHoppersConvergedTime))) ;
    if (VERBOSITY_HIGH) {
        cout << "Hoppers converged: time (s) = " << GetFastestTime() << "; hoppers = " << _nHoppers
             << "; hoppers collected/injected = " << _collectorCurrent
             << "/" << _generatorCurrent
             << "; current (A) = " <<  _currentStore.back() << endl;
    }
    // See if the current has converged
    if (_currentStore.back() != 0.0) {  // make sure there is a current.... 
        _run = false;
        for (int i = 0; i < _cyclesForConvergence; i++) {
            if (_currentStore[0]/_currentStore[i] < _lowerTol || _currentStore[0]/_currentStore[i] > _upperTol) {
                _run = true;  // hasn't converged yet...
                break;
            }
        }
    }
    if (!_run) {
         cout << "Current converged at time = " << GetFastestTime()
              << "\n\tNumber of charges = " << _nHoppers
              << "\n\tCurrent (A) = " <<  _currentStore.back() << endl;
    }
}

/**********************************************************************
 * MOVE FUNCTIONS
 * The actual move is executed by 'Move', but this is wrapped by a 
 * 'MoveFastest' function which also takes care of everything else
 *  that has to be done before / after a 'Move'
 *
 * TODO: Some of the 'MoveFastest' functions should certainly be merged
 ***********************************************************************/
// Actually move the charge.
double hoppers::Move( list <hopper *>::iterator H, vertex * to, double &fastestTime){
    vertex * from = (*H)->GetFrom() ;
    #ifdef printTotalOccupation
    if (_track) {
        cout << GetHopperNumber(from) 
             << '\t' << fastestTime
             << '\t' << from->GetX() 
             << '\t' << from->GetY()
             << '\t' << from->GetZ() << endl;
    }
    #endif
    double dz = GetFastestDz();
    if ( !to->IsOccupied() ) {     	
        if (_hopperInteractions) {
            DeleteCoulomb(from);
        }
        from->SetUnoccupied(fastestTime);  // Note: do this after DeleteCoulomb
        _mapVertexToHopper.erase(from); 
        to->SetOccupied(fastestTime);  // Note: do this before AddCoulomb
        _mapVertexToHopper[to] = *H;
        if(_hopperInteractions)	{
            (*H) -> Move(to);
            AddCoulomb(to);  // If 'to' is generator, shouldn't be here!
        }
        else {
            (*H) -> SetHop(to,fastestTime);
        }
        return dz;
    }
    // If 'to' is occupied, don't move but just give a new waitTime
    //   (taking into account the disabled reaction).
    else {				
        (*H) -> SetHopOccNeigh(from, fastestTime);	
        return 0.0;
    }						
}
// The MoveFastest function for simple ToF
double hoppers::MoveFastest_C() {
    vertex * to = (*_fastest)->GetTo();
    double fastestTime=GetFastestTime();
    double dz;
    if ( to->IsCollector() ) {
        _reciprocalCollectionTimes.push_back(1.0/fastestTime);
        _totalReciprocalCollectionTimes+=(1.0/fastestTime);
        dz = GetFastestDz();
        Remove(_fastest,fastestTime);
    }
    else dz = Move(_fastest, to, fastestTime);
    FindFastest();        
    return dz;
}
// The MoveFastest function for simple Regenerate mode
double hoppers::MoveFastest_R() {
    vertex * to = (*_fastest)->GetTo();
    double fastestTime = GetFastestTime();
    double dz;
    if ( to->IsCollector() ) {			
        double transitTime = fastestTime - (*_fastest)->GetGenerationTime();
        _reciprocalCollectionTimes.push_back(1.0 / transitTime);
        _totalReciprocalCollectionTimes += 1.0 / transitTime;
        dz = GetFastestDz();
        Remove(_fastest,fastestTime);
        GenerateAll(1,fastestTime);
    }
    else dz = Move(_fastest, to, fastestTime);
    FindFastest();
    return dz;
}
// The 'MoveFastest' function for simple Regnerate mode, with 
//   Coulombic interactions.  !!! WARNING !!! Poorly tested.
//   TODO: Merge with above function
double hoppers::MoveFastest_RCI() {
    vertex * to = (*_fastest)->GetTo();	
    double fastestTime=GetFastestTime();
    double dz;
    if ( to->IsCollector() ) {
        double transitTime = fastestTime - (*_fastest)->GetGenerationTime();
        _reciprocalCollectionTimes.push_back(1.0 / transitTime);
        _totalReciprocalCollectionTimes += 1.0 / transitTime;
        dz = GetFastestDz();
        Remove(_fastest,fastestTime);
        GenerateAll(1,fastestTime);
    }
    else dz = Move(_fastest, to, fastestTime);
    SetHops_C(fastestTime);
    FindFastest();        
    return dz;
}
// The 'MoveFastest' function for the FET mode
double hoppers::MoveFastest_F() {
    vertex * to = (*_fastest)->GetTo();	
    double fastestTime=GetFastestTime();
    double dz;
    dz = Move(_fastest, to, fastestTime);
    _totalCurrent = _generatorCurrent + _collectorCurrent;
    SetSourceDrainOccupation(fastestTime); 
    SetHops_C(fastestTime);
    // Check there are still hoppers 
    if (_nHoppers == 0) { 
        cout << "!!! WARNING !!! : Ran out of hoppers\n";
        WARNINGS++;
        _run = false;
        return .0;
    }
    FindFastest();        
    ++_moves;
    // Check if maxTime has been exceeded
    if (GetFastestTime() > _maxTime) {
        cout << "!!! WARNING !!! : maxTime exceeded!  Time = " << GetFastestTime() << endl;
        WARNINGS++;
        _run=false;
        return .0;
    }
    // Check for convergence
    if ( _moves==_movesCycle ) {
	    _moves=0;
	    if (_activeHoppersConverged) {
            FETConvergence(); 
        }
	    else { 
            activeHoppersConvergence();
        }
    }
    return dz;
}

/******************
 * TIME FUNCTIONS
 ******************/
// Find the fastest hopper
void hoppers::FindFastest() {
    if (_hoppers.empty()) {
        _fastestTime = 1e50;
    }
    else{
        list <hopper *>::iterator it_hop = _hoppers.begin();
        _fastest = it_hop;
        for ( ; it_hop != _hoppers.end(); ++it_hop){
            if ( (*it_hop)->GetWaitTime() < (*_fastest)->GetWaitTime()) _fastest = it_hop;
        }
        _fastestTime=(*_fastest)->GetWaitTime();
    }
}
// Set all hoppers' waitTimes to 'time'.  
//   Used at end of simulations, needed for occupation times.
void hoppers::SetWaitTimes(double time) {
    list <hopper *>::iterator it_hop = _hoppers.begin();
    for (; it_hop!=_hoppers.end(); ++it_hop) {
        (*it_hop)->SetWaitTime(time);
    }
}

/***************************************
 * GET HOPPER FUNCTIONS
 *
 * TODO: Should merge these functions... 
 ***************************************/ 
// Find the hopper pointer, given the vertex.  Return Hopper.
hopper * hoppers::GetHopper(vertex * v) {
    map < vertex *, hopper * > :: iterator it = _mapVertexToHopper.find(v);
    if (it != _mapVertexToHopper.end()){ 
        return (it -> second);
    }
    else {
        cout << "***ERROR*** : Thought vertex " << _graph->GetVertex(v) << " was occupied but can't find hopper\n";
        cout << "              " << v->IsOccupied() << '\t' << v->IsCollector() << '\t' << v->IsGenerator() << endl;
        cout << "Active hoppers = " << GetActive() << endl; 
        cout << "Occupied vertices = "; _graph->PrintOccupied();
        exit(-1);
    }
}
// Find the hopper iterator to pointer, given the vertex.  Return iterator.
list <hopper *>::iterator hoppers::GetHopperIterator(vertex *v) {
    map < vertex *, hopper * > :: iterator it_vert = _mapVertexToHopper.find(v);
    list <hopper *>::iterator it_hop = _hoppers.begin();
    for (int i=0; it_hop!=_hoppers.end(); i++, it_hop++) {
        if ( *it_hop == it_vert->second ){ return it_hop; }
    }
    cout << "***ERROR*** : Thought vertex " << _graph->GetVertex(v) << " was occupied but can't find hopper\n";
    cout << "              " << v->IsOccupied() << '\t' << v->IsCollector() << '\t' << v->IsGenerator() << endl;
    cout << "Active hoppers = " << GetActive() << endl;
    cout << "Occupied vertices = "; _graph->PrintOccupied();
    exit(-1);
}
// Find the hopper index, given the vertex. Return index.
int hoppers::GetHopperNumber(vertex *v) {
    map < vertex *, hopper * > :: iterator it_vert = _mapVertexToHopper.find(v);
    list <hopper *>::iterator it_hop = _hoppers.begin();
    for (int i=0; it_hop!=_hoppers.end(); i++, it_hop++) {
        if ( *it_hop == it_vert->second ){ return i; }
    }
    cout << "***ERROR*** : Thought vertex " << _graph->GetVertex(v) << " was occupied but can't find hopper\n";
    cout << "              " << v->IsOccupied() << '\t' << v->IsCollector() << '\t' << v->IsGenerator() << endl;
    cout << "Active hoppers = " << GetActive() << endl;
    cout << "Occupied vertices = "; _graph->PrintOccupied();
    exit(-1);
}
