#ifndef INTERVAL_H

#define INTERVAL_H

#include "std.h"
#include "printable.h"

class interval:public printable  {
		
	public:
		interval () {}
		interval (int id, string chr, double start, double end, string info = "");
        void set (int id,string chr, double start, double end, string info) ;
		string to_string() const;

		int id ;
		string chr;
		double start;
		double end;
        string info;

};

std::size_t hash_value (const interval& s) ;
bool operator < (const interval& s, const interval& t) ; 
bool operator == (const interval &s, const interval& t) ; 
bool operator && (const interval &s, const interval &t) ;
interval operator + (const interval &s, const interval &t) ; 
interval operator * (const interval &s, const interval &t) ; 

#endif
