#include "std.h"


class covariance {
    int numeg, numblock, numh2;
    vector<string> eglist;
};

void covariance::read_f3blocks (string filename ) {
    int blocksize  = 40;    
    char buf[blocksize];
    memset (buf, 0,blocksize); 
    ifstream ifs (filename, ios::in|ios::binary);
    FILE* fptr =fopen (filename.c_str(), "rb");
    fread (buf, sizeof(char), blocksize, fptr);
    sscanf (buf, "numeg: %d nblock: %d nh2: %d", &numeg, &numblock, &numh2) ;
    
    for (int k = 0 ; k < numeg; k++) {
        memset (buf, 0,blocksize); 
        fread (buf, sizeof(char), blocksize, fptr);
        string s (buf);
        cout << s << endl;
    }

    for (int k = 0 ; k < numblock; k++) { 

    }

    fclose (fptr);
    cout << numeg << "\t" << numblock << "\t" << numh2 << endl;

}

int main (int argc, char *argv[]) { 
    covariance c();
    read_f3blocks(string(argv[1]));

}
