#include "std.h"
#include "intervalmap.h"


int main (int argc, char *argv[]) {
    if (argc>2) { 
        string s1=  string(argv[1]);
        string s2 = string(argv[2]);
        intervalmap i1 = intervalmap (s1, NULL, NULL);
        intervalmap i2 = intervalmap (s2, NULL, NULL);
//        cout << i1.to_string () << endl;
//      cout << i2.to_string () << endl;
        i1.add (i2);
        i1.merge ();
        
       	for (int i = 0 ; i < i1.chrs.size(); i++) {
		unordered_map<string, vector<interval> >::const_iterator it  = i1.intervals.find(i1.chrs[i]);
		const vector<interval>& tmp = it->second;
		for (int j = 0; j < tmp.size();j++) {
            cout << tmp[j].to_string() << endl;
        }
	    }
    }
}
