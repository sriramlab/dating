#include "interval.h"
#include "stringfn.h"



interval::interval (int id,string chr, double start, double end, string info) { 
	this->id = id;
	this -> chr = chr;
	this -> start = start;
	this -> end = end;
    this->info = info;
}

void interval::set (int id,string chr, double start, double end, string info) { 
	this->id = id;
	this -> chr = chr;
	this -> start = start;
	this -> end = end;
    this->info = info;
}
string interval::to_string () const{ 
	std::stringstream ss;
	ss.precision(0);
	ss.setf(ios::fixed);
	ss << "(" << tostring(id) << "," <<chr << ",";
	ss << noshowpoint <<  start << "," << end <<")";
	return ss.str();
}



bool operator < (const interval& s, const interval& t) { 
	if ( s.chr.compare(t.chr) == 0 ) {
		if (s.start == t.start && s.end == t.end)
			return (s.id < t.id);

		if (s.start == t.start)
			return (s.end < t.end);
		return (s.start < t.start);
	} else {
		return s.chr.compare(t.chr)<0;
	}
}


bool operator == (const interval &s, const interval& t) { 
	return ((s.id==t.id) && s.chr.compare(t.chr)==0 && (s.start==t.start) && (s.end==t.end));
}


bool operator && (const interval &s, const interval &t) {
	return (s.chr.compare(t.chr)==0 && 
            ((s.start<=t.start && s.end>=t.start) || 
             (s.start<=t.end && s.end>=t.end) ||
             (t.start<=s.start && t.end>=s.start)||
             (t.start<=s.end && t.end>=s.end))
            );
}

interval operator + (const interval &s, const interval &t) { 
    if (s && t) {
        double x = s.start<t.start?s.start:t.start;
        double y  = s.end<t.end?t.end:s.end;
        interval t2 = interval (0,s.chr,x,y);
        return (t2);
    } else {
        cerr << "Result not an interval " << endl;
        exit(1);
    }
}

interval operator * (const interval &s, const interval &t) { 
    if (s&&t) { 
        double x = s.start>t.start?s.start:t.start;
        double y  = s.end>t.end?t.end:s.end;
        interval t2 = interval (0,s.chr,x,y);
        return (t2);
    } else { 
        cerr << "Result not an interval " << endl;
        exit(1);
    }
}

std::size_t hash_value (const interval& s) {
	std::size_t seed = 0;
	boost::hash_combine (seed, s.id);
	boost::hash_combine (seed, s.chr);
	boost::hash_combine (seed, s.start);
	boost::hash_combine (seed, s.end);
	return seed;
}



