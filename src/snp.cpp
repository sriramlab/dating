#include "std.h"
#include "stringfn.h"
#include "snp.h"


snp::snp (string id, string chr,  double physpos, double genpos, char var, char ref, bool ignore) { 
	this->id = id;
	this->chr = chr;
	this->genpos = genpos;
	this->physpos = physpos;
	this->ignore =  ignore;
	this->var = var;
	this->ref = ref;
    setattributes ();
}


string snp::to_string () const{ 
	ostringstream oss;
	oss << id << ","<<chr<<","<<genpos<<","<<getphyspos()<<",I="<<ignore;
	return oss.str();
}

string snp::stats() const  {
    return "";
}

string snp::getgtype () const {
	string s = "";
	for (int i  = 0 ; i < gtype.size(); i++) {
		s += tostring(gtype[i]);
	}
	return s;
}

string snp::getphyspos () const { 
	std::stringstream ss;
	ss.precision(0);
	ss.setf(ios::fixed);
	ss << noshowpoint <<  physpos;
	return ss.str();

}

string snp::tophyspos (double physpos) { 
	std::stringstream ss;
	ss.precision(0);
	ss.setf(ios::fixed);
	ss << noshowpoint <<  physpos;
	return ss.str();

}

string snp::getgeneticpos () const {
	std::stringstream ss;
    ss.setf(ios::fixed,ios::floatfield); 
	ss.precision(7);
    ss << genpos;
    return ss.str();
}

void snp::setattributes () {
    char v = var;
    char r = ref;
    if (v>r) {
        char t = v;
        v =  r;
        r = t;
    }

    transition = false;
    transversion = false;
    if ((v=='A' && r=='G')||(v=='C'&&r=='T'))
        transition = true;
    if ((v=='A' && r=='T')||(v=='A' && r=='C')
            || (v=='C' && r=='G') || (v=='G'&&r=='T'))
        transversion = true;
}

bool snp::istransition () {
    return transition;
}

bool snp::istransversion () {
    return transversion;
}

bool operator < (const snp& s, const snp& t) { 
	if ( s.chr.compare(t.chr) == 0 ) {
		return (s.physpos <= t.physpos);
	} else {
		return s.chr.compare(t.chr)<0;
	}
}


bool operator == (const snp &s, const snp& t) { 
	return (s.chr.compare(t.chr)==0 && (s.physpos==t.physpos));
}


std::size_t hash_value (const snp& s) {
	std::size_t seed = 0;
	boost::hash_combine (seed, s.chr);
	boost::hash_combine (seed, s.physpos);
	return seed;
}


