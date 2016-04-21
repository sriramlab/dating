#include "std.h"
#include "ind.h"

ind::ind (string id, char gender, string group) { 
	this->id = id;
	this->gender = gender;
	this->group = group;
}

string ind::to_string ()  const{
	ostringstream oss;
	oss << id << ","<<gender<<","<<group;
	return oss.str();
}
