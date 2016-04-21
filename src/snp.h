#ifndef SNP_H

#define SNP_H

#include "std.h"
#include "printable.h"

class snp : public printable {
	public :
		string id;
		string chr;
		double genpos;
		double physpos;	
		char var;
		char ref;

		vector<int> gtype;
		int vcount, rcount;
		double freq;
		bool ignore;
        bool transition;
        bool transversion;

		snp (string, string, double, double genpos=0, char var='X', char ref='X', bool ignore=false);
		snp () {}
		string to_string () const;
		string stats () const;
		string getphyspos () const;
        string getgeneticpos () const ;
		string getgtype () const;
        bool istransversion () ;
        bool istransition () ;
        void setattributes () ;
        static string tophyspos (double) ;
};

std::size_t hash_value (const snp& s) ;
bool operator < (const snp& s, const snp& t) ; 
bool operator == (const snp &s, const snp& t) ; 
#endif

