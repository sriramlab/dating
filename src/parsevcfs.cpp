#include "packedgtype.h"
#include "genotype.h"
#include "io.h"
#include "stringfn.h"
#include "mathfn.h"
#include "ind.h"
#include "gzstream.h"
#include "std.h"
#include "functions.h"
#include "fileio.h"

class parsevcf  {
    public:
		parsevcf (int argc, char *argv[]) ;
        void read_vcf (string nvcf, string dvcf, string kgfile) ;

		// Parameters passed to this program
		data *d;

		// Parameters specific to the data
		data *p;
        string outfile;


};

parsevcf::parsevcf (int argc, char *argv[]) {
	static const char *optString = "vh";

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
	
	string geno1, snp1, ind1;
    string inputformat = "eigenstrat";

    string nvcf;
    string dvcf;
    string kgfile;
	p->get_string ("nvcf",nvcf,true);
	p->get_string ("dvcf",dvcf,true);
    p->get_string ("kgfile", kgfile,true);
    p->get_int ("debug", io::debug, 0 );
    p->get_string ("outfile", outfile, true);

    cout << "debug =  " << io::debug << endl;

    read_vcf (nvcf, dvcf, kgfile);

}



void parsevcf::read_vcf (string nvcf, string dvcf, string kgfile) { 
    igzstream ig1 (nvcf.c_str());
    igzstream ig2 (dvcf.c_str());
    string line;
	vector<snp> tmp;
	int j = 0 ;
    unordered_map <string, string> formatmap;
    unordered_map <string, string> infomap;
    unordered_map <string, string> attributemap;

    unordered_map<string, string> *tmpmap = new unordered_map<string, string> ();
    fileio::read_map (kgfile,tmpmap,0,-1);
    int npops;
    vector<int> nchr;
    for (unordered_map<string,string>::iterator i = tmpmap->begin(); i != tmpmap->end(); i++) { 
        string val = i->second;
        vector<string> tmptoks;
        functions::tokenize ( val.c_str(), tmptoks, " \t");
        npops  =(tmptoks.size() - 3)/2; 
        nchr.resize(npops);
        for (int i  = 0  ; i < npops; i++) { 
            int a = atoi(tmptoks[1+2*i].c_str());
            int b = atoi(tmptoks[2+2*i].c_str());
            nchr[i] = a+b;
        }
    }
    int total  =1;
    for (int i  =  0; i  < nchr.size(); i++) {
        total *= (nchr[i]+1);
    }
    total *= 9;
    vector<int> afs ( total);

	double prevgpos ;
	double prevppos;
	string prevchr="";
    int size =  10000;
    int inc = 10000;
    
    unordered_map <string, string> * gsmap  = NULL;
    bool givengoodsnps  = false;

    /*
    if (functions::fileExists (goodsnpsname) ) { 
        gsmap =new unordered_map <string, string> ();
        fileio::read_map (goodsnpsname, gsmap, 0);
        givengoodsnps  = true;
        cout << "goodsnpname =  " << goodsnpsname  << endl;
    }*/
	
    int linenum = 0 ;
    bool setind = false;
    int numhap = 2;

    string line1;
    string line2;
    bool y1 = std::getline (ig1,line1).good();
    bool y2 = std::getline (ig2,line2).good();
    while (y1 && y2) { 
        linenum ++;
        char c1 = line1[0];
        char c2 = line2[0];
        bool flag = false;
        if (c1=='#') { 
            y1 = std::getline (ig1,line1).good();
            flag = true;
        }
        if (c2=='#'){
            y2 = std::getline (ig2,line2);
            flag = true;
        }


        if (flag)
            continue;

        vector<string> toks1;
        functions::tokenize (line1.c_str(), toks1, " \t");
        string id1 = toks1[2];
        string chr1 = toks1[0];
        double physpos1 =  atof(toks1[1].c_str());
        double genpos1 = physpos1*1.3/1e8;
        if (id1.compare (".")==0) { 
            id1 = chr1+":"+toks1[1];
        }

        vector<string> toks2;
        functions::tokenize (line2.c_str(), toks2, " \t");
        string id2 = toks2[2];
        string chr2= toks2[0];
        double physpos2 =  atof(toks2[1].c_str());
        double genpos2 = physpos2*1.3/1e8;
        if (id2.compare (".")==0) { 
            id2 = chr2+":"+toks2[1];
        }


        if ( chr1.compare(chr2)==0 ){

            if ( physpos1 == physpos2 ) {

                if (io::debug>=1 ) {
                    cout << "Processing " << id1  << endl;
                }
                string ref = toks1[3];
                string var = toks1[4];
                string pass = toks1[6];
                string info = toks1[7];
                string format = toks1[8];
                string geno = toks1[9];
                string anc = "";
                string der = "";
                int ngtype; 
                int dgtype;

                if (ref.length() > 1 || var.length() > 1) {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id1 << ": indel" <<endl;
                    y1 = std::getline (ig1,line1).good();
                    continue;
                }


                vector<string> tok1;
                vector<string> tok2;
                vector<string> tok3;
                functions::tokenize (info.c_str(), tok1, ";");


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

                flag = true;
                if (attributemap.find("Map20")!=attributemap.end()){
                    if (atof(attributemap["Map20"].c_str()) !=1) {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id1 << ": below Map20 cutoff"<<endl;
                        flag = false;
                    }
                } else  {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id1 << ": no Map20 field"<<endl;
                    flag = false;
                }
                if (attributemap.find("GQ")!=attributemap.end()){
                    if (atof(attributemap["GQ"].c_str())>=30);
                    else {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id1 << ": below GQ cutoff"<<endl;
                        flag = false;
                    }
                } else {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id1 << ": no GQ field"<<endl;
                     flag = false;
                }
                if (attributemap.find ("DP")!=attributemap.end()){
                    double rd = atof(attributemap["DP"].c_str());
                    if (rd >= 31 && rd <= 79);
                    else {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id1 << ": outside DP range"<<endl;
                        flag = false;
                    }
                } else  {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id1 << ": no DP field"<<endl;
                     flag = false;
                }

                if (attributemap.find("CAnc")!=attributemap.end()) {
                    int hcount =  0; int pcount = 0 ;
                    if (attributemap.find("TS")!=attributemap.end()){
                        string val = attributemap["TS"];
                        for (int i = 0  ; i < val.length();i++) {
                            hcount += (val[i]=='H');
                            pcount += (val[i]=='P');
                        }
                    }
                    anc = attributemap["CAnc"];
                    /*
                    if ( hcount ==1 && pcount==1)
                        anc = attributemap["CAnc"];
                    else 
                        flag = false;*/
                } else {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id1 << ": no CAnc field"<<endl;
                    flag = false;
                }


                if (ref.compare(".")==0) 
                    ref="X";
                if (var.compare(".")==0)
                    var="X";

                //
                // Fill in the genotypes
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

                if (!flag) {
                    y1 = std::getline (ig1,line1).good();
                    continue;
                }

                for (int i = 9, k = 0; i <= 9; i++, k++){ 
                    tok1.resize(0);
                    functions::tokenize(toks1[i], tok1,",");
                    int a1 =  tok1[gind][0];
                    int a2 =  tok1[gind][2];
                    ngtype = (a1-48) + (a2-48);
                }

                if (io::debug>=0)  {
                    cout << "Accepted n " << ngtype << endl;
                }

                if ( ref.compare(anc)==0 ) {
                } else  if ( var.compare("X") == 0 || var.compare (anc) == 0 ) { 
                    ngtype =  2 - ngtype;
                    der = ref; 
                } else
                    flag = false;

                // Denisova
                ref = toks2[3];
                var = toks2[4];
                pass = toks2[6];
                info = toks2[7];
                format = toks2[8];
                geno = toks2[9];

                if (ref.length() > 1 || var.length() > 1) {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id2 << ": indel" <<endl;
                    y2 = std::getline (ig2,line2);
                    continue;
                }

                tok1.clear();
                tok2.clear();
                tok3.clear();
                functions::tokenize (info.c_str(), tok1, ";");


                attributemap.clear ();
                for (int i =  0 ; i < tok1.size(); i++){ 
                    tok2.resize(0);
                    functions::tokenize (tok1[i], tok2, "=");
                     attributemap[tok2[0]] = tok2[1];
                }
                tok3.resize(0);tok2.resize(0);
                functions::tokenize (format.c_str(), tok3, ":");
                functions::tokenize (geno.c_str(), tok2, ":");
                for (int i  = 0 ; i < tok3.size(); i++) { 
                     attributemap[tok3[i]] = tok2[i];
                }

                flag = true;
                if (attributemap.find("Map20")!=attributemap.end()){
                    if (atof(attributemap["Map20"].c_str()) !=1) {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id2 << ": below Map20 cutoff"<<endl;
                        flag = false;
                    }
                } else  {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id2 << ": no Map20 field"<<endl;
                    flag = false;
                }
                if (attributemap.find("GQ")!=attributemap.end()){
                    if (atof(attributemap["GQ"].c_str())>=30);
                    else {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id2 << ": below GQ cutoff"<<endl;
                        flag = false;
                    }
                } else {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id2 << ": no GQ field"<<endl;
                     flag = false;
                }
                if (attributemap.find ("DP")!=attributemap.end()){
                    double rd = atof(attributemap["DP"].c_str());
                    if (rd >= 16 && rd <= 46);
                    else {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id2 << ": outside DP range"<<endl;
                        flag = false;
                    }
                } else  {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id2 << ": no DP field"<<endl;
                     flag = false;
                }

                //
                // Fill in the genotypes
                tok1.resize (0);
                functions::tokenize (format.c_str(), tok1,":");
                gind = -1;

                for (int i = 0 ; i < tok1.size(); i++){
                    if (tok1[i].compare ("GT")==0) { 
                        gind = i;
                        break;
                    }
                }
                if (gind == -1)
                    flag = false;

                if (!flag) {
                    y2 = std::getline (ig2,line1);
                    continue;
                }


                for (int i = 9, k = 0; i <= 9; i++, k++){ 
                    tok1.resize(0);
                    functions::tokenize(toks2[i], tok1,",");
                    int a1 =  tok1[gind][0];
                    int a2 =  tok1[gind][2];
                    dgtype = (a1-48) + (a2-48);
                }

                if (io::debug>=1)  {
                    cout << "Accepted d " << dgtype << endl;
                }

                bool flip = false;
                if ( ref.compare(anc)==0 ) {
                } else  if (var.compare("X") == 0 || var.compare (anc) == 0 ) { 
                    dgtype =  2 - dgtype;
                    der = ref; 
                    flip = true;
                } else
                    flag = false;


                if (tmpmap->find(id1) !=tmpmap->end()) { 
                    string val = (*tmpmap)[id1];
                    vector<string> tmptoks;
                    functions::tokenize ( val.c_str(), tmptoks, " \t");
                    int index = atoi (tmptoks[1].c_str()) ;
                    for (int i  = 1 ;  i < nchr.size(); i++) { 
                        index *= (nchr[i]+1);
                        index += atoi(tmptoks[2*i+1].c_str());
                    }

                    if (index == 3)  {
                        long tmppos = ((long)physpos1);
                        cout << "SNP with AC=3 " << chr1 << ":" << tmppos << endl;
                    }

                    index *= 3; index += dgtype;
                    index *= 3; index += ngtype;
                    afs[index]++;

                    if (io::debug >=0 )  {
                        cout << "Found SNP in 1kg" << endl;
                        cout << id1 << "\t";
                        for ( int i =  0 ; i < nchr.size(); i++) {
                            cout << nchr[i] << "\t";
                        }
                        cout  << dgtype << "\t" << ngtype <<"\t" << index <<"\t" << afs[index] << endl;
                    }
                } else { 
                    if (0) {
                    if (flip) {
                        int index = nchr[0] ;
                        for (int i  = 1 ;  i < nchr.size(); i++) { 
                            index *= (nchr[i]+1);
                            index += nchr[i];
                        }
                        index *= 3; index += dgtype;
                        index *= 3; index += ngtype;
                        afs[index]++;

                        if (io::debug >=0 )  {
                            cout << "Did not find SNP in 1kg :  flipped " << endl;
                            for ( int i =  0 ; i < nchr.size(); i++) {
                                cout << nchr[i] << "\t";
                            }
                            cout << id1 << "\t" << dgtype << "\t" << ngtype <<"\t" << index <<"\t" << afs[index] << endl;
                        }

                    } else { 
                        int index =  0 ;
                        index *= 3; index += dgtype;
                        index *= 3; index += ngtype;
                        afs[index]++;


                        if (io::debug >=0 )  {
                            cout << "Did not find SNP in 1kg  " << endl;
                            for ( int i =  0 ; i < nchr.size(); i++) {
                                cout << nchr[i] << "\t";
                            }
                            cout << id1 << "\t" << dgtype << "\t" << ngtype <<"\t" << index <<"\t" << afs[index] << endl;
                        }


                    }
                    }
                }
                if (io::debug>=0)  {
                    cout << "Accepted " << id1 << endl;
                }
                
                y1 = std::getline (ig1,line1).good();
                y2 = std::getline (ig2,line2).good();
            } else if (physpos1 < physpos2) {
                y1 = std::getline (ig1,line1).good();
            } else if (physpos2 < physpos1) { 
                y2 = std::getline (ig2,line2).good();
            }
        }

    }

    if (io::debug >= 1) { 
        cout << "Printing" << endl;
    }

    ig1.close ();
    ig2.close ();

    vector<int> counter (nchr.size());
    int counts = 0 ;
    ofstream os (outfile);
    while (1) {
        if (io::debug >= 1) {
            cout << "counts = "<< counts << endl;
        }
        counts ++;

        for (int j = 0 ; j < 3; j++){ 
            for (int k = 0 ; k < 3; k++) { 
                 int index = counter[0] ;
                 for (int i  = 1 ;  i < counter.size(); i++) { 
                     index *= (nchr[i]+1);
                     index += counter[i];
                 }
                 index *= 3; index += j;
                 index *= 3; index += k;

                 string s = "";
                 for (int i  = 0 ; i < counter.size(); i++)
                     s += tostring(counter[i]) + "\t";
                s += tostring(j) + "\t" + tostring (k);
                s += "\t" + tostring (afs[index]);
                os << s << endl;

            }
        }
        int i = nchr.size () - 1;
        for (; i>=0 && counter[i] == nchr[i] ; i--) ;
        if (i<0)
            break;
        counter[i]++;
        for (int j  = i+1; j < counter.size(); j++)
            counter[j] =  0;
            
    }
        os.close();
}

int main (int argc, char* argv[]) {
	parsevcf p = parsevcf(argc,argv);
}

