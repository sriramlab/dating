#include "printable.h"
#include "io.h"

int io::log = 0 ; 
int io::debug = 0 ;
ofstream io::logfs;
/*
inline
void io::print (const string &s, int level )  {
	if ( log == 0 || logfs == NULL ) { 
		if (level <= debug)
			cout << s ;
	} else {
		if (level <= io::debug)
			logfs << s;
	}
}

inline
void io::println (const string &s, int level) {
	print (s+"\n",level);
}

inline
void io::print (const string &s, int level, string file) { 
	if (level <= debug) { 
		ofstream ofs (file.c_str());
		ofs << s;
		ofs.close();
	}
}

inline
void io::print (printable &p, int level) { 
	print (p.to_string(), level);
}

inline
void io::println (const printable &p, int level) {
	io::print (p.to_string()+"\n", level);
}

inline
void io::print (printable &p, int level, string file) { 
	print (p.to_string(), level, file);
}*/
