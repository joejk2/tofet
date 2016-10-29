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
 * 'graph' is the class that defines the graph through which hoppers 
 * (charges) move.  The graph is constituted of vertices (molecules) 
 * and edges (connections between molecules).
 *********************************************************************/
#ifndef _GRAPH_H
#define	_GRAPH_H
#include "vertex.h"
#include "global.h"
#include "IO.h"

class graph{
    private:
        vector <vertex *> _vertices;
        double _Vg;  // V.Ang^-1 (sorry!)
        double _fieldZ;  // V.Ang^-1
        double _temp;  // K
        vector <vector <double> > _CoulombGrid;
        bool _hopperInteractions; 
        double _tmpX, _tmpY, _tmpZ;
    // end of private:
    
    public:
        double _reorg; // eV
        double _kT;    // eV
        double _sourceFermiEnergy;
        double _drainFermiEnergy;
        double _coulombPrefactor;
        graph(){}
        graph(char * sim, char *xyz, char *edge){
            if (Read(sim, "mode", "tof")!="fet") {
                _fieldZ = atof(Read(sim, "fieldZ").c_str()); 
            }
            else {
                _fieldZ = 1e50;
                _Vg = atof(Read(sim, "Vg").c_str()); 
                double Vds = atof(Read(sim, "Vds").c_str());
                _sourceFermiEnergy = _Vg;
                _drainFermiEnergy  = _Vg + Vds;
                cout << "Source Fermi energy = " << _sourceFermiEnergy
                     << ", drain Fermi energy = " << _drainFermiEnergy << endl;
            }
            _reorg = atof(Read(sim, "reorg").c_str()); 
            _temp = atof(Read(sim, "temp").c_str());
            _kT = _temp*k_eVK;
            if (Read(sim, "hopperInteractions", "0") == "1") {
                _hopperInteractions = true;
            }
            else {
                _hopperInteractions = false;
            }
            // Read site energies:
            if (Read(sim, "siteEnergies", "0") == "1"\
             || _hopperInteractions\
             || Read(sim, "mode", "tof") == "fet") {
                if (VERBOSITY_HIGH) {
                    cout << "Reading E's from ***.xyz\n";
                }
                ReadVertices(xyz, _vertices, true);  // read E's 
                ReadEdges(edge, _vertices, false);  // don't read deltaE's
                if (Read(sim, "mode", "tof") == "fet") {
                    SetDEs();
                    SetRatesPrefactor_C();  // doesn't set the field!
                    // The Coulomb prefactor in eV.Ang/e^2:
                    _coulombPrefactor=14.3996442/atof(Read(sim, "dielectric").c_str());
                    // The following should be uncommented if you want to use a look-up 
                    //   table for the Coulombic interactions.  See also 'GetSingleCoulomb'
                    //   in hoppers.cc
                    // MakeCoulombEnergyGrid(); 
                }
                if (Read(sim, "mode", "tof") == "regenerate") {
                    SetField_E();
                    SetDEs();
                    if (_hopperInteractions) {
                        SetRatesPrefactor_C();
                        _coulombPrefactor=14.3996442/atof(Read(sim, "dielectric").c_str());
                        // See above notes about MakeCoulombEnergyGrid()
                    }
                    else {
                        SetRates_DE();
                    }
                }
                if (Read(sim, "mode", "tof") == "tof") {
                    SetField_E();
                    SetDEs();
                    SetRates_DE();
                }
                if (Read(sim,"printVertices","0")=="1") {
                    PrintVertices("E");
                }
                if (Read(sim,"printEdges","0")=="1") {
                    PrintEdges("E", _hopperInteractions);
                } 
            }
            // Read delta E's
            else {
                if (VERBOSITY_HIGH) {
                    cout << "Reading delta E's from ***.edge\n";
                }
                ReadVertices(xyz,_vertices,false);  // don't read E
                ReadEdges(edge,_vertices,true);	  // do read deltaE
                SetField_DE();
                SetRates_DE();
                if (Read(sim,"printVertices","0")=="1") PrintVertices("DE");
                if (Read(sim,"printEdges","0")=="1") {
                    PrintEdges("DE", _hopperInteractions);
                }
            }
        }
        ~graph(){
            vector < vertex * >::iterator it = _vertices.begin();
            for (int i=0; it!=_vertices.end(); i++,it++) {
                delete (*it);
            }
            _vertices.clear();
        }

    /*****************************
     * SETS and ADDS
     ****************************/
    void ClearDCs(); // reset _DCs to 0
    //void AddEdge(const unsigned int &, const unsigned int &, const double &);
    void SetDEs();  // set deltaE's from E's
    void SetField_E();  // modify E's to reflect an applied field...
    void SetField_DE();  // ... likewise for deltaE's
    void SetRatesPrefactor_C();  // set prefactors (no energies) 
    void SetRates_DE();  // set rates from deltaE's 
    void NormaliseOccupationTimes(const double, int);  
    void MakeCoulombEnergyGrid();  
    double const &GetCoulomb(vertex *, vertex *);  // ... from a grid

    /*****************************
     * PRINTS AND READS 
     ****************************/
    void ReadEdges(char *, vector <vertex *> &, bool);
    void ReadVertices(char *, vector <vertex *> &, bool);
    void PrintEdges(string, bool);
    void PrintVertices(string);
    void PrintEnergies();  // print average sum of static + coulomb energies
    void PrintOccupied();  // print all occupied molecules	
    void PrintTotalOccupationTimes();

    /*****************************
     * GETS 
     ****************************/
    vector <vertex *> GetPreviouslyOccupied(char *);  // read occupied vertices from file
    vertex * GetEmptyGenerator();  // returns random empty generator
    double GetDepth();  // get the depth of the graph in the z direction
    double GetDistance(vertex *, vertex *);  // get the distance between two vertices
    int GetVertex( vertex *);			// return index of vertex in vertices
    int CountTotalElectrodes();
    const double & GetFieldZ() 	const {return _fieldZ;}
    vector <vertex *> GetCollectors();
    vector <vertex *> GetGenerators(); 
};
#endif	/* _GRAPH_H */
