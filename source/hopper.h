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

/*********************************************************************
 * 'hopper' describes a single hopper (charge), where it is, where it 
 * will go and when. 'hopper' is also able to make the actual move of 
 * a hopper.
 ********************************************************************/

#ifndef _HOPPER_H
#define	_HOPPER_H
#include "global.h"
#include "vertex.h"

class hopper{
    private:
        vertex * _from;  // where the hopper currently is (and will hop from)
        vertex * _to;  // where the hopper will hop to
        double _waitTime;  // when the hopper will hop
        double _dZ;  // how far along the 'z' axis the hopper will hop
        double _timeGenerated;  // the time at which the hopper was generated
     
    public:
        hopper() {
            _waitTime=0.0;
            _dZ=0.0;
            _timeGenerated = 0.0;
        }
        hopper(vertex * V, const double & time) {
            _from = V;
            _from->SetOccupied(time);
            _timeGenerated = time;
        }
        ~hopper() {
            _from->SetUnoccupied(_waitTime);
            _waitTime=2e10;
            _timeGenerated = 2e10;
        }
    
    /**********
     * DO'S
     *********/
    void Move(vertex *V) {
        _from=V;
    }
    void SetHop(vertex * V, const double &time) {
        Move(V);
        SetHop(time);
    }
    // SetHop when a neighbour is occupied
    void SetHopOccNeigh(vertex * V, const double &time) {
        Move(V);
        double totalRate=_from->CalcTotalRateToUnoccupied();  // recalculate total rate
        #ifdef RandomB
        _waitTime = time - log( UniformPos() ) / totalRate;
        #else
        _waitTime = time - log( gsl_rng_uniform_pos(gslRand) ) / totalRate;
        #endif
        if (totalRate>0) {
            _to = _from->ChooseToUnoccupied(totalRate);
            _dZ = (_to -> GetZ()) - (_from ->GetZ());
        }
        else {
            _to = _from;
            _dZ = 0.0;
        }
    }
    void SetHop(const double &time){
        #ifdef RandomB
        _waitTime = time - log( UniformPos() ) / _from->GetTotalRate();	
        #else
        _waitTime = time - log( gsl_rng_uniform_pos(gslRand) ) / _from->GetTotalRate();	
        #endif
        if (_from->GetTotalRate()>0) {
            _to = _from->ChooseTo();
            _dZ = (_to -> GetZ()) - (_from ->GetZ());
        }
        else {
            _to = _from;
            _dZ = 0.0;
        }
    }
    void SetWaitTime(double time) {
        _waitTime=time;
    }

    /**********
     * GET'S
     *********/
    const double  & GetDz () const {
        return _dZ;
    }
    const double & GetWaitTime () const	{
        return _waitTime;
    }
    const double & GetGenerationTime () const {
        return _timeGenerated;  
    }
    vertex * GetTo () const	{
        return _to;
    }
    vertex * GetFrom () const {return 
        _from;
    }
};
#endif	/* _HOPPER_H */
