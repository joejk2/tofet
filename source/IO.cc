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

/***********************
 * Lesser General input/output 
 ***********************/

#include "IO.h"
void open(char *filename, ifstream &in) {
    in.open(filename);
    if (!in) {
        cout << "***ERROR***: Unable to open " << filename << endl;
        exit(-1);
    }
}
//
void PrintAll(char *filename) {
    ifstream in;
    open(filename,in);
    string line;
    while (in) {
        getline(in,line);
        cout << "  : " << line << endl;
    }
}
//
string Read(char *filename, string name, string value) {
    ifstream in;
    open(filename, in);
    string word;
    while (in) {
        in >> word; 
        if (word==name) {
            in >> value;
            return value;
        }
    }
    if (value.empty()) {
        cout << "***ERROR***: Didn't find " << name << " in " << filename << endl;
        exit(-1);
    }
    in.close();
    return value;
}
