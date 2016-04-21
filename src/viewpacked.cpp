#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "manip.h"


int main (int argc, char* argv[]) {
    string s = argv[1];
    packedgtype::print_packedes (s);
}
