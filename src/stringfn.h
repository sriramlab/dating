#ifndef STRINGFN_H
#define STRINGFN_H

#include "std.h"
#include <sstream>

class stringfn{
	public:
//	static void tokenize(const string &str, vector<string> &tokens, const string& delimiters = " ");
};

template <class T>
inline std::string tostring (const T& t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
}


template <class T>
inline std::string tostring (const vector<T>& t, string delim = " ")
{
	std::stringstream ss;
    for (int i = 0;  i  <  t.size(); i++)
    	ss << t[i] << delim;
	return ss.str();
}
/*
void stringfn::tokenize(const string& str,
		                      vector<string>& tokens,
				      const string& delimiters )
{
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}
*/


/*
template <class T>
std::string operator+ (const T& t, const char* rhs) {
	std::stringstream ss;
	ss << t;
	ss << rhs;
	return ss.str();
		
}

template <class T>
std::string operator+ (const char* &lhs, const T &t) {
	std::stringstream ss;
	ss << lhs;
	ss << t;
	return ss.str();
		
}
*/


//std::string operator+ (const string &lhs, const double &t) ;

//std::string operator+ (const string &lhs, const int &t) ;
#endif
