#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "manip.h"



int main(int argc, char *argv[])
{
	unordered_map<string, string > *pmapsnps = NULL;
    pmapsnps = new unordered_map <string, string> ();
    string s (argv[1]);
    fileio::read_map (s,  pmapsnps, 0 , 1);

}

