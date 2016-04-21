#include "std.h"
#include "data.h"
#include "io.h"
#include "functions.h"
#include "fileio.h"

void fileio::read_map (string filename, unordered_map<string,double> *m, int idx, int vidx){ 
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
        vector<string> toks;
        functions::tokenize(line.c_str(),toks," \t");
        string val = line;
        if (vidx >= 0 && vidx < toks.size() )
            val = toks[vidx];
        if (idx < toks.size()){
            (*m)[toks[idx]] = atof(val.c_str());
        } else {
            cerr << "Invalid column number in file " << filename << endl;
            exit(1);
        }
	}	
	inp.close ();
}

void fileio::read_map (string filename, unordered_map<string,int> *m, int idx, int vidx){ 
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
        vector<string> toks;
        functions::tokenize(line.c_str(),toks," \t");
        string val = line;
        if (vidx >= 0 && vidx < toks.size() )
            val = toks[vidx];
        if (idx < toks.size()){
            (*m)[toks[idx]] = atoi(val.c_str());
        } else {
            cerr << "Invalid column number in file " << filename << endl;
            exit(1);
        }
	}	
	inp.close ();
}

void fileio::read_map (string filename, unordered_map<string,string> *m, int idx, int vidx){ 
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
		if (idx>=0){
			vector<string> toks;
			functions::tokenize(line.c_str(),toks," \t");
			string val = line;
			if (vidx >= 0 && vidx < toks.size() )
				val = toks[vidx];
			if (idx < toks.size()){
				(*m)[toks[idx]] = val;
			} else {
				cerr << "Invalid column number in file " << filename << endl;
				exit(1);
			}
		} else {
			(*m)[line] = line;
		}
	}	
	inp.close ();
}


void fileio::read_vector (string filename, vector<string> *m, int idx){ 
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
		if (idx>=0){
			vector<string> toks;
			functions::tokenize(line.c_str(),toks," \t");
			if (idx < toks.size()){
                m->push_back(toks[idx]);
			} else {
				cerr << "Invalid column number in file " << filename << endl;
				exit(1);
			}
		} else {
            m->push_back(line);
		}
	}	
	inp.close ();
}




