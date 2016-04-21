#ifndef CLDATA_H

#define CLDATA_H

#include "std.h"
#include "data.h"

class cldata:public data {

	public:
	cldata (int argc, char *argv[], const char *optString = NULL, const option *longOpts = NULL) ;

};
#endif
