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
 * 'vertex' is the object that describes a single molecule, including 
 * position, list of neighbours, and rates to those neighbours.
 ********************************************************************/
#ifndef _VERTEX_H
#define	_VERTEX_H
#include "vec.h"
#include "global.h"
#include <algorithm>

using namespace std;

class vertex{

    private:
        int _ID;  // ID of vertex in graph::_vertices
        double _posZ;  // position along the 'z' axis
        vector <vertex *> _neighbours;
        vector <double> _Js;
        vector <double> _DEs;  // deltaE between vertices
        double _E;  // site energy, as read in from ***.xyz
        vector <double> _rates;
        double _totalRate; 
        bool _occupied;
        string _type;		// generator (g), collector (c), other (-)
        // The following are only used when 'hopperInteractions'
        //   are enabled.
        vector <double> _DCs;  // difference in Coulomb energies 
                               //   between vertices 
        vector <double> _ratesPrefactor; 
        double _EC;  // Coulomb energy of a charge on this vertex 
        double _EC_time;  // a running sum of _EC * time (used for calculating potentials)
        double _oldTime;  // last time the occupation status of this vertex changed 
        // The following are used when 'printTotalOccupation' is enabled
        double _totalOccupationTime;  // total time this vertex is occupied by a hopper
        double _timeOfOccupation;  // when a hopper last moved onto the vertex
        unsigned int _timesOccupied;
    
    public:
        vec _pos;  // public because needed so often for Coulombic calculations
        bool _electrode;  // likewise, a bit messy, but speeds up Coulombic calculations.
        vertex(){
            _occupied = false;
            _totalOccupationTime=0.0;
            _ID=-1;
            _E=0.0;
            #ifdef printTotalOccupation
            _EC=0.0;
            _EC_time=0.0;
            _oldTime=0.0;
            _timesOccupied=0;
            #endif
            _electrode=false;
        }
        ~vertex(){
            _neighbours.clear();
            _rates.clear();
            _Js.clear();
            _DEs.clear();
            _DCs.clear();
        }
        
        /*******************************
         * SETUP
         ******************************/
        void AddNeighbour(vertex *, const double &); 
        void AddNeighbour(vertex *, const double &, const double &); 
        void SetPos(const vec & pos);
        void SetType (string);
        void SetID(int i) 	{_ID = i;}
        void SetE(double E);
        void SetDEs();
        void SetField_E(const double &_fieldZ);
        void SetField_DE(const double &_fieldZ);
        
        /*******************************
         * SET RATES
         ******************************/
        void SetRates_DE(const double &, const double &, const double &, const double &);
        void SetRatesPrefactor_C(const double &, const double &);
        void UpdateRates_C(const double &, const double &);

        /*******************************
         * MISCELLANEOUS
         *******************************/
        void SetOccupied(const double & time);
        void SetUnoccupied(double time);
        void IncrementEC(const double newEC, const double time);
        void SetEC(const double newEC, const double time);
        double CalcTotalRateToUnoccupied(); 
        vertex * ChooseTo() const;
        vertex * ChooseToUnoccupied(double ) const;
        void IncrementDCs(int, double);
        void IncrementTotalOccupationTime(const double & time) {_totalOccupationTime+=time;}
        void NormaliseTotalOccupationTime(const double maxTime, int totalHoppers);
        void ClearDCs();

        /*******************
         * PRINTS and GETS 
         ******************/
        void PrintEdges();
        void PrintPos() {cout << "(" << _pos.getX() << ", " << _pos.getY() << ", " << GetZ() <<  "), " << _type ;} 
        const int &GetID() const {return _ID;}
        const double &GetTotalRate() const {return _totalRate;} 
        const double &GetRate(int i) const {return _rates.at(i);}
        const double &GetJ(const int & i) {return _Js.at(i);}
        const double &GetDE(const int & i) {return _DEs.at(i);}
        const double &GetDC(const int & i) {return _DCs.at(i);}
        const double &GetE() {return _E;}
        const double &GetEC_time() {return _EC_time;}
        const double &GetX() const {return _pos.getX();}
        const double &GetY() const {return _pos.getY();}
        const double &GetZ() const {return _posZ;}
        const vec &GetPos() const {return _pos;}
        const string &GetType() const {return _type;}
        unsigned int GetNumberNeighbours() {return _neighbours.size();}
        vector <vertex *> & GetNeighbours()	{return _neighbours;}
        const bool & IsOccupied() const	{return _occupied;}
        vector < double > & GetRates() {return _rates;}
        const double & GetTotalOccupationTime() {return _totalOccupationTime;}
        unsigned int GetTimesOccupied()	{return _timesOccupied;}
        const bool IsCollector() const {if (_type=="c") return true; else return false;}
        const bool IsGenerator() const {if (_type=="g") return true; else return false;}
};
#endif	/* _VERTEX_H */
