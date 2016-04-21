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
		bool ignore;
		bool genetic;
		string var;
		string ref;
		double d1;
		double p1, p2;
		// c2 = total counts between this and the previous snp
		// c1 = total counts between this and the previous map snp|this snp is a map snp
		int c1, c2;

		// Dirichlet
		double u1;
		
		// z2 = genetic distance between this and previous snp
		// z1 = genetic distance between this and previous map snp|this snp is a map snp
		double z1, z2;

		int truec;

		double truez;
		double scale;


		int geno;
		double af;

		bool insideinterval;
		bool start;
		bool end;

		vector<int> gtype;
		int vcount, rcount;
		double freq;

        bool transition;
        bool transversion;

		snp (string, string, double, double genpos=0, string var="", string ref="", bool ignore=false, bool genetic=true);
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

