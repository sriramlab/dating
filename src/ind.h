#ifndef IND_H

#define IND_H

#include "std.h"
#include "snp.h"
#include "printable.h"
#include "data.h"
class ind : public printable {
	public:
		string id;
		char gender;
		string group;

		ind (string, char , string);
		string to_string () const;
};
#endif

