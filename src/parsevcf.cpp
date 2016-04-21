#include "io.h"
#include "stringfn.h"
#include "mathfn.h"
#include "ind.h"
#include "gzstream.h"
#include "std.h"
#include "functions.h"
#include "fileio.h"

// Very basic program
// Parses a VCF file
// Outputs SNPs that pass some filters
class parsevcf  {
    public:
		parsevcf (int argc, char *argv[]) ;
        void read_vcf (string ) ;

		// Parameters passed to this program
		data *d;

		// Parameters specific to the data
		data *p;
        int nind;

        string filterpolicy;

};

parsevcf::parsevcf (int argc, char *argv[]) {
	static const char *optString = "vh";
    filterpolicy = "160812";

	static const struct option longOpts[] = {
		{ "parameter", required_argument, NULL, 'p' },
		{ "verbose", no_argument, NULL, 'v' },
		{ "help", no_argument, NULL, 'h' },
		{ NULL, no_argument, NULL, 0 }
	};
	d = new data (argc, argv, optString, longOpts);
	string pfile;
	d->get_string ("parameter",pfile,true);
	p = new data (pfile);
	p->print_parameters();
	
    string vcffile;
	p->get_string ("vcffile",vcffile,true);
    p->get_int ("debug", io::debug, 0 );
    p->get_string ("filterpolicy", filterpolicy, false);
    nind = -1;

//    cout << "debug =  " << io::debug << endl;

    read_vcf (vcffile);

}



void parsevcf::read_vcf (string vcffile) { 
    if (!functions::fileExists(vcffile)) {
        cerr  << "ERROR: cannot read " << vcffile << endl;
        exit(1);
    }
    igzstream ig (vcffile.c_str());

    // A hashtable that stores key,value pairs for a SNP
    // Use this to apply filters on attributes
    //
    unordered_map <string, string> attributemap;

    int linenum = 0 ;
    string line;
    bool y = std::getline (ig,line).good();
    vector<int> gtype;


    while (y) { 
        linenum ++;
        char c = line[0];
        bool flag = false;

        // Strip out header
        if (c=='#') { 
            // Get the header line
            // Use it to get the number of chromosomes
            if (line.find ("CHROM")!=string::npos) {
                vector<string> toks;
                functions::tokenize (line.c_str(), toks, " \t");
                nind = toks.size() - 9 ;
            }
            y = std::getline (ig,line).good();
            flag = true;
        }

        if (io::debug >=2 )  {
            cout << line << endl;
        }

        if (flag)
            continue;

        // Parse the fields separated by whitespace
        vector<string> toks;
        functions::tokenize (line.c_str(), toks, " \t");
        string chr = toks[0];
        double physpos =  atof(toks[1].c_str());
        double genpos = physpos*1.3/1e8;
        //
        // Use the id if given, else use chr:pos
        string id = toks[2];
        if (id.compare (".")==0) { 
            id = chr+":"+toks[1];
        }

        string ref = toks[3];
        string var = toks[4];
        string pass = toks[6];

        // Parse these fields to filter fields and get genotypes
        string info = toks[7];
        string format = toks[8];
        string geno = toks[9];

        string anc = "";
        string der = "";

        // Filter out indels
        if (ref.length() > 1 || var.length() > 1) {
            if (io::debug >= 1)
                cout << "Dropping position " << id << ": indel" <<endl;
            y = std::getline (ig,line).good();
            continue;
        }


        vector<string> tok1;
        vector<string> tok2;
        vector<string> tok3;

        functions::tokenize (info.c_str(), tok1, ";");

        // attributemap is a hash table 
        // (key=value) pairs are extracted from the info field
        // key,value pairs are also extracted from the format and geno fields
        attributemap.clear ();
        for (int i =  0 ; i < tok1.size(); i++){ 
            tok2.resize(0);
            functions::tokenize (tok1[i], tok2, "=");
            if ( tok2.size() > 1)
                attributemap[tok2[0]] = tok2[1];
        }
        tok3.resize(0);tok2.resize(0);
        functions::tokenize (format.c_str(), tok3, ":");
        functions::tokenize (geno.c_str(), tok2, ":");
        for (int i  = 0 ; i < tok3.size(); i++) { 
            attributemap[tok3[i]] = tok2[i];
        }

        // Apply all the filters here
        // To drop a SNP, set flag = false
        // To retrieve a values for a key, use attributemap[key]
        flag = true;
        if (filterpolicy.compare ("160812")==0) { 
            if (attributemap.find("Map20")!=attributemap.end()){
                if (atof(attributemap["Map20"].c_str()) !=1) {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": below Map20 cutoff"<<endl;
                    flag = false;
                }
            } else  {
                if (io::debug >= 1)
                    cout << "Dropping position " << id << ": no Map20 field"<<endl;
                flag = false;
            }

            if (attributemap.find("GQ")!=attributemap.end()){
                if (atof(attributemap["GQ"].c_str())>=30);
                else {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": below GQ cutoff"<<endl;
                    flag = false;
                }
            } else {
                if (io::debug >= 1)
                    cout << "Dropping position " << id << ": no GQ field"<<endl;
                flag = false;
            }

            if (attributemap.find ("DP")!=attributemap.end()){
                double rd = atof(attributemap["DP"].c_str());
                if (rd >= 31 && rd <= 79);
                else {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": outside DP range"<<endl;
                    flag = false;
                }
            } else  {
                if (io::debug >= 1)
                    cout << "Dropping position " << id << ": no DP field"<<endl;
                flag = false;
            }
        } else if (filterpolicy.compare ("")==0) {
            if (attributemap.find("Map20")!=attributemap.end()){
                if (atof(attributemap["Map20"].c_str()) !=1) {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": below Map20 cutoff"<<endl;
                    flag = false;
                }
            } else  {
                if (io::debug >= 1)
                    cout << "Dropping position " << id << ": no Map20 field"<<endl;
                flag = false;
            }

            if (attributemap.find("GQ")!=attributemap.end()){
                if (atof(attributemap["GQ"].c_str())>=30);
                else {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": below GQ cutoff"<<endl;
                    flag = false;
                }
            } else {
                if (io::debug >= 1)
                    cout << "Dropping position " << id << ": no GQ field"<<endl;
                flag = false;
            }

            if (attributemap.find ("DP")!=attributemap.end()){
                double rd = atof(attributemap["DP"].c_str());
                if (rd >= 31 && rd <= 79);
                else {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": outside DP range"<<endl;
                    flag = false;
                }
            } else  {
                if (io::debug >= 1)
                    cout << "Dropping position " << id << ": no DP field"<<endl;
                flag = false;
            }
        }

        /*
        if (attributemap.find("CAnc")!=attributemap.end()) {
            int hcount =  0; int pcount = 0 ;
            if (attributemap.find("TS")!=attributemap.end()){
                string val = attributemap["TS"];
                for (int i = 0  ; i < val.length();i++) {
                    hcount += (val[i]=='H');
                    pcount += (val[i]=='P');
                }
            }
            if ( hcount ==1 && pcount==1)
                anc = attributemap["CAnc"];
            else 
                flag = false;
        } else {
            flag = false;
        }*/


        if (ref.compare(".")==0) 
            ref="X";
        if (var.compare(".")==0)
            var="X";

        //
        // Fill in the genotypes
        //
        // Get the index of the genotype in the format field (the GT key)
        tok1.resize (0);
        functions::tokenize (format.c_str(), tok1,":");
        int gind = -1;

        for (int i = 0 ; i < tok1.size(); i++){
            if (tok1[i].compare ("GT")==0) { 
                gind = i;
                break;
            }
        }
        if (gind == -1)
            flag = false;
        if (nind < 0 ) 
            nind = toks.size () - 9;


        // Fill up the vector of genotypes
        gtype.resize(nind,9);

        if (!flag) {
            y = std::getline (ig,line).good();
            continue;
        }

        for (int i = 9, k = 0; i < toks.size(); i++, k++){ 
            tok1.resize(0);
            functions::tokenize(toks[i], tok1,",");
            int a1 =  tok1[gind][0];
            int a2 =  tok1[gind][2];
            gtype[k] = (a1-48) + (a2-48);
        }


        // Use this to polarize genotypes e.g.
        if ( ref.compare(anc)==0 ) {
        } else  if ( var.compare("X") == 0 || var.compare (anc) == 0 ) { 
        } else
            flag = false;


        if (io::debug>=1)  
            cout << "Accepted " << id << endl;

        // Do any processing on the filtered SNPs
        //
        cout << "Accepted " << id <<  endl;


        // Pull the next SNP
        y = std::getline (ig,line).good();
    }

    ig.close ();

}

int main (int argc, char* argv[]) {
	parsevcf p = parsevcf(argc,argv);
}

