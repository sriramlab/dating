#include "packedgtype.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "ind.h"
#include "gzstream.h"
#include "std.h"
#include "functions.h"


packedgtype::packedgtype (int nind, int nsnp, vector<unsigned int> cumnumsnps,  bool isgenotype, bool hasmissing, data *d, vector<snp> *snps, string entryformat, int order) {
    
    if (entryformat.compare ("prob")==0 ) { 
        // Stores a probability ( 0<=x<=1, x=double)
        // with precision of 4 digits
        // needs 2 chars
        // 2 digits stored in 1 char
        wordsize = sizeof(char)  * 8;
        unitsize =  4;
        unitsperword = wordsize/unitsize;
        mask = 0;
        for (int i = 0 ; i < unitsize; i++)
            mask = mask |(0x1<<i);

        if (io::debug >= 2)
            cout << "Mask = " << hex << setfill('0') << setw(2) << int(mask) << dec << endl;

    } else { 
        if (!isgenotype && !hasmissing )  {
            wordsize = sizeof(char)  * 8;
            unitsize =  1;
            unitsperword = wordsize/unitsize;
            mask = 0;
            for (int i = 0 ; i < unitsize; i++)
                mask = mask |(0x1<<i);

            if (io::debug >= 2)
               cout << "Mask = " << hex << setfill('0') << setw(2) << int(mask) << dec << endl;
        } else if (isgenotype||hasmissing) {
            wordsize = sizeof(char)  * 8;
            unitsize = 2;
            unitsperword = wordsize/unitsize;
            mask = 0;
            for (int i = 0 ; i < unitsize; i++)
                mask = mask |(0x1<<i);

        } else {
            cerr << "Error in packedgenotype::packedgenotype" << endl;
            exit(1);
        }
    }

    this->nind = nind;
    this->nsnp = nsnp;
    this->isgenotype = isgenotype;
    this->hasmissing = hasmissing;
    this->entryformat =  entryformat;
    this->entrytype = format_to_type (entryformat);
    this->order = order;
    if (entryformat.compare ("prob")==0) { 
        ncol = 2*nind;
        nrow = nsnp;
    } else { 
        if (order == SNPMAJOR) {
            ncol = ceil(1.0*nind/unitsperword);
            nrow = nsnp;
        } else if (order == INDMAJOR){ 
            ncol = ceil(1.0*nsnp/unitsperword);
            nrow = nind;
        }
    }
    if (io::debug >= 1)
       cout <<  "# nrow = " << nrow << ", ncol = " << ncol << "\t unitsperword = " << unitsperword <<  endl;
    gtype =  new unsigned char[nrow*ncol];
    vcount = new int [nsnp];
    rcount = new int [nsnp];
    this->cumnumsnps = cumnumsnps;
    this->d =  d;
    this->snps = snps;
}

packedgtype::packedgtype ( int nsnp, vector<unsigned int> cumnumsnps, bool isgenotype, bool hasmissing, data *d, vector<snp> *snps, string entryformat, int order) {
    if (entryformat.compare ("prob")==0) { 
        wordsize = sizeof(char)  * 8;
        unitsize =  4;
        unitsperword = wordsize/unitsize;
        mask = 0;
        for (int i = 0 ; i < unitsize; i++)
            mask = mask |(0x1<<i);

        if (io::debug >= 2)
            cout << "Mask = " << hex << setfill('0') << setw(2) << int(mask) << dec << endl;


    } else { 
        if (!isgenotype && !hasmissing )  {
            wordsize = sizeof(char)  * 8;
            unitsize =  1;
            unitsperword = wordsize/unitsize;
            mask = 0;
            for (int i = 0 ; i < unitsize; i++)
                mask = mask |(0x1<<i);

            if (io::debug >= 2)
               cout << "Mask = " << hex << setfill('0') << setw(2) << int(mask) << dec << endl;
        } else if (isgenotype) {
            wordsize = sizeof(char)  * 8;
            unitsize =  2;
            unitsperword = wordsize/unitsize;
            mask = 0;
            for (int i = 0 ; i < unitsize; i++)
                mask = mask |(0x1<<i);

        } else {
            cerr << "Error in packedgenotype::packedgenotype" << endl;
            exit(1);
        }
    }

    this->nind = 0;
    this->nsnp = nsnp;
    this->entryformat =  entryformat;
    this->entrytype = format_to_type (entryformat);
    this->isgenotype = isgenotype;
    this->hasmissing = hasmissing;
    this->order = order;
    vcount = new int [nsnp];
    rcount = new int [nsnp];
    this->cumnumsnps = cumnumsnps;
    this->d =  d;
    this->snps = snps;
}

void packedgtype::read_eigenstrat (string filename, vector<bool> *indivfilter, vector<bool> *snpfilter) {
    read_es (filename, indivfilter, snpfilter);
}

void packedgtype::read_es (string filename, vector<bool> *indivfilter, vector<bool> *snpfilter) { 
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
    unsigned int rowindex = 0 ;
    int effectiveindex  = 0;

	bool setind = false;
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
        vector<string> toks;

		if (c=='#')
			continue;
		
		if (line.empty())
			continue;

        effectiveindex ++;
        if (snpfilter != NULL && (*snpfilter)[effectiveindex-1] == false)
            continue;


        if (entryformat.compare ("prob")==0) 
            functions::tokenize(line.c_str (), toks, " \t");

        if (nind <= 0) {
            if (entryformat.compare ("prob")==0)  {
                nind = toks.size();
                ncol = 2*nind;
                nrow = nsnp;
            } else {
                nind = line.length ();
                if (order == SNPMAJOR) {
                    ncol = ceil(1.0*nind/unitsperword);
                    nrow = nsnp;
                } else if (order == INDMAJOR){ 
                    ncol = ceil(1.0*nsnp/unitsperword);
                    nrow = nind;
                }
            }
            if (io::debug >= 2)
                cout <<  " nrow = " << nrow << ", ncol = " << ncol << "\t unitsperword = " << unitsperword <<  endl;
            gtype =  new unsigned char[nrow*ncol];
        }

        if (entryformat.compare ("prob")==0) {
            unsigned int i = 0 ;
            for (int i1  = 0; i1 < toks.size();i1++){
            //for (int i = 0 ; i < nind; i++) { 

                if (indivfilter!=NULL && (*indivfilter)[i1]==false)
                    continue;
                double d = atof (toks[i1].c_str());
                ostringstream os;
                os.setf(ios::fixed,ios::floatfield); 
                os.precision(4); os << d;
                string s = os.str();
                double e  = stod (s);

                string::size_type pos = s.find_first_of (".");
                unsigned int ind = 2*i;
                if ( pos == string::npos ) {
                    if (d>=1){
                        gtype[rowindex*ncol + ind] = 160;
                        gtype[rowindex*ncol + ind+1] = 0;
                    }
                } else if (e>=1) {
                    gtype[rowindex*ncol + ind] = 160;
                    gtype[rowindex*ncol + ind+1] = 0;
                } else {
                    gtype[rowindex*ncol + ind] = gtype[rowindex*ncol + ind+1] = 0;
                    for (unsigned int j = pos+1, k=0 ; j < s.length() && k < 4; j++, k++) {
                        unsigned int ind1 = rowindex*ncol + ind  + k/unitsperword;
                        unsigned char c = s.at(j); c-=48;
                        gtype[ind1] = (gtype[ind1] << unitsize) | c; 
                    }
                }
                i++;
//            cout << d << "\t" << s << "\t" << gtype[rowindex*ncol  + ind ]<<","<<gtype[rowindex*ncol + ind + 1] <<endl;
//            exit(1);
            }

            rowindex++;

        } else  {
            int vc = 0;
            int rc = 0;
            int shift =  0;
            unsigned int i =  0;
            for (int j = 0; j < line.length();j++){
                if (indivfilter!=NULL && (*indivfilter)[j]==false)
                    continue;
                if (i%unitsperword==0)
                    gtype[rowindex*ncol+i/unitsperword]=0;
                unsigned char c = line.at(j);
                unsigned char d = c;
                c-= 48;
                c=(c==9)?3:c;

                gtype[(rowindex)*ncol+i/unitsperword] = gtype[(rowindex)*ncol+i/unitsperword] |(c<< (unitsize*shift)) ;
                //            cout <<linenum << "\t" << hex << setfill('0') << setw(2)<<  int(gtype[(rowindex)*ncol+i/unitsperword]) << endl;

                if (c!=9) {
                    vc += c;
                    rc += isgenotype?(2-c):(1-c);
                }
                shift = (shift +1 ) % unitsperword; i ++;
            }
            vcount[rowindex] = vc;
            rcount[rowindex] = rc;
            rowindex++;
        }

    }
    inp.close ();
    /*
    if (entryformat.compare ("prob")==0) {
        cout << "genomatrix" << endl;
        for (int j = 0 ; j < nsnp && j < 100 ; j++) {
            for (int i = 0 ; i < nind; i++)  {
                int ind = j*ncol + 2 * i ;
                cout << gtype[ind] <<"," << gtype[ind+1] << "\t";
            }
            cout << endl;
        }
        exit(1);
    }*/
}

void packedgtype::set_freq ( int chr, vector<snp> &s) {
    int n = s.size();
    for (int i = 0 ; i < n  ; i++) {
        int vc = vcount [ cumnumsnps[chr]+i];
        int rc = rcount [ cumnumsnps[chr]+i];
        s[i].vcount = vc;
        s[i].rcount = rc;
        s[i].freq = (vc+rc>0)?((1.*vc)/(vc+rc)):0;
    }
    snps = &s;
}

ostream& packedgtype::write_packedes ( ostream &ofs) { 
    char format[20] = "packedes";
    int version = 1;
    cout << "# ncol = " << ncol <<"\tnrow = " << nrow << endl;

    fileio::binary_write ( ofs, format);
    fileio::binary_write ( ofs, version);
    fileio::binary_write ( ofs, entrytype);
    fileio::binary_write ( ofs, isgenotype);
    fileio::binary_write ( ofs, hasmissing);
    fileio::binary_write ( ofs, order);
    fileio::binary_write ( ofs, ncol);
    fileio::binary_write ( ofs, nrow);
    fileio::binary_write ( ofs, nind);
    fileio::binary_write ( ofs, nsnp);
    fileio::binary_write ( ofs, cumnumsnps);
    ofs.write (reinterpret_cast<const char*>(gtype), nrow*ncol*sizeof(unsigned char));
    return ofs;
    
}

void packedgtype::write_packedes(string filename) {
    ofstream ofs (filename, ios::out|ios::binary);
    write_packedes (ofs);
    ofs.close ();
}


void packedgtype::read_packedes(string filename, vector<bool> *indivfilter, vector<bool> *snpfilter) {
    int c = 0 ;
    for (int i = 0 ; indivfilter!=NULL && i < indivfilter->size(); i++)
        c+= ((*indivfilter)[i]==false);
    for (int i = 0 ; snpfilter!=NULL && i < snpfilter->size(); i++)
        c+= ((*snpfilter)[i]==false);

    if (c==0)
        read_packedes (filename);
    else {
        int version;
        char format[20] ;
        bool isgenotype;
        bool hasmissing; 
        int order;
        int ncol;
        int nrow;
        int nind;
        int nsnp;
        vector<unsigned int> cumnumsnps;
        etype entrytype;


        ifstream ifs (filename, ios::in|ios::binary);
        fileio::binary_read ( ifs, format);
        if (strcmp (format,"packedes") != 0 ) {
            cerr << "Unknown file format in " << filename << endl;
            exit (1);
        }

        fileio::binary_read ( ifs, version);
        if (version != 1) {
            cerr << "Incorrect version number in " << filename << endl;
            exit(1);
        }

        fileio::binary_read ( ifs, entrytype);
        if ( entrytype <0  || entrytype> 1 ) {
            cerr << "Incorrect entry type in " << filename << endl;
            exit(1);
        }
        if ( entrytype != this->entrytype ){ 
            cerr << "Metadata mismatch in " << filename << endl;
            exit(1);
        }
        string entryformat = get_format (entrytype);

        fileio::binary_read ( ifs, isgenotype);
        fileio::binary_read ( ifs, hasmissing);
        fileio::binary_read ( ifs, order);
        fileio::binary_read ( ifs, ncol);
        fileio::binary_read ( ifs, nrow);
        fileio::binary_read ( ifs, nind);
        fileio::binary_read ( ifs, nsnp);
        fileio::binary_read ( ifs, cumnumsnps);


        if (isgenotype!=this->isgenotype){
            cerr << "isgenotype entry does not match in file " << filename << endl;
            exit(1);
        }
        if (hasmissing!=this->hasmissing){
            cerr << "hasmissing entry does not match in file " << filename << endl;
            exit(1);
        }
        if (order!=this->order){
            cerr << "order entry does not match in file " << filename << endl;
            exit(1);
        }

        /*
        for (int i  = 0 ; i < cumnumsnps.size(); i++) {
            if (cumnumsnps[i] != this->cumnumsnps[i]){
                cerr << "snp vector does not match in file " << filename << endl;
                exit(1);
            }
        }*/

        set_metadata ();
        packedgtype *tmppg = new packedgtype (nind, 1, vector<unsigned int>(1,0), isgenotype, hasmissing, NULL, NULL, entryformat, order);
        vector<int> tmpg (nind, 0);
        vector<int> tmpnewg (this->nind,0);
        unsigned char *tmp = new unsigned char[ncol];
        for (int i = 0 , nr = 0; i < nrow; i++) {
            ifs.read ( reinterpret_cast<char *>(tmp), ncol*sizeof(unsigned char));
            if (snpfilter != NULL && (*snpfilter)[i] == false)
                continue;
            tmppg->gtype = tmp;
            tmppg->get_geno (0,tmpg);
            for (int j  =0, nc = 0 ; j < nind; j++) {
                if (indivfilter!=NULL && (*indivfilter)[j]==false)
                    continue;
                tmpnewg[nc] = tmpg[j];
                nc ++;
            }
            set_geno (nr, tmpnewg);
            nr ++;
        }
        delete[] tmp;
        ifs.close ();
    }
}

void packedgtype::read_packedes(string filename) {
    ifstream ifs (filename, ios::in|ios::binary);
    char format[20] ;
    int version;
    vector<unsigned int> cumnumsnps;

    fileio::binary_read ( ifs, format);
    if (strcmp(format,"packedes") != 0 ) {
        cerr << "Unknown file format in " << filename << endl;
        exit (1);
    }

    fileio::binary_read ( ifs, version);
    if (version != 1) {
        cerr << "Incorrect version number in " << filename << endl;
        exit(1);
    }

    int entrytype;
    fileio::binary_read ( ifs, entrytype);
    if ( entrytype <0 || entrytype > 1 ) {
        cerr << "Incorrect entry type in " << filename << endl;
    }

    fileio::binary_read ( ifs, isgenotype);
    fileio::binary_read ( ifs, hasmissing);
    fileio::binary_read ( ifs, order);
    fileio::binary_read ( ifs, ncol);
    fileio::binary_read ( ifs, nrow);
    fileio::binary_read ( ifs, nind);
    fileio::binary_read ( ifs, nsnp);
    fileio::binary_read ( ifs, cumnumsnps);
    ifs.read (reinterpret_cast<char*>(gtype), nrow*ncol*sizeof(unsigned char));
    if (io::debug>=1) {
        cout << isgenotype << "\t" << hasmissing << endl;
        cout << "gtype = " << gtype[0] << "," << gtype[1] << "," << gtype[2] << "," << gtype[3] << endl;
    }
    /*
    cout << "cumnumsnps = " << cumnumsnps.size() << endl;
    for (int i = 0 ; i  < cumnumsnps.size(); i++)
        cout << "chr= " << i << "\t"<< cumnumsnps[i] << endl;
        */
    cout << (*this)(0,0) << (*this)(0,1) << (*this)(0,2) << (*this)(0,3) << (*this)(0,4) <<endl;
    set_metadata ();
    ifs.close ();


}



void packedgtype::print_packedes(string filename) {
    ifstream ifs (filename, ios::in|ios::binary);
    char format[20] ;
    int version;
    bool isgenotype;
    bool hasmissing;
    int entrytype;
    int order;
    int nind; 
    int nsnp;
    int ncol;
    int nrow;
    vector<unsigned int> cumnumsnps;

    fileio::binary_read ( ifs, format);
    if (strcmp(format,"packedes") != 0 ) {
        cerr << "Unknown file format in " << filename << endl;
        exit (1);
    }

    fileio::binary_read ( ifs, version);
    if (version != 1) {
        cerr << "Incorrect version number in " << filename << endl;
        exit(1);
    }

    fileio::binary_read ( ifs, entrytype);
    if ( entrytype <0 || entrytype > 1 ) {
        cerr << "Incorrect entry type in " << filename << endl;
    }

    fileio::binary_read ( ifs, isgenotype);
    fileio::binary_read ( ifs, hasmissing);
    fileio::binary_read ( ifs, order);
    fileio::binary_read ( ifs, ncol);
    fileio::binary_read ( ifs, nrow);
    fileio::binary_read ( ifs, nind);
    fileio::binary_read ( ifs, nsnp);
    fileio::binary_read ( ifs, cumnumsnps);
    unsigned char *gtype = new unsigned char [nrow*ncol];
    ifs.read (reinterpret_cast<char*>(gtype), nrow*ncol*sizeof(unsigned char));
    if (io::debug>=1) {
        cout << isgenotype << "\t" << hasmissing << endl;
        cout << "gtype = " << gtype[0] << "," << gtype[1] << "," << gtype[2] << "," << gtype[3] << endl;
    }

    ifs.close ();
    vector <int> tmp (nind, 0);

    unsigned char mask;
    int wordsize;
    unsigned int unitsperword;
    int unitsize;
    switch (entrytype)  {
            case GENO:
            if (!isgenotype && !hasmissing)  {
                wordsize = sizeof(char)  * 8;
                unitsize =  1;
                unitsperword = wordsize/unitsize;
                mask = 0;
                for (int i = 0 ; i < unitsize; i++)
                    mask = mask |(0x1<<i);
            } else {
                wordsize = sizeof(char)  * 8;
                unitsize = 2;
                unitsperword = wordsize/unitsize;
                mask = 0;
                for (int i = 0 ; i < unitsize; i++)
                    mask = mask |(0x1<<i);

            }
            break;
        default:
            break;

    }

    for (int j = 0 ; j < nsnp; j ++ ){ 
        for (int i = 0 ;i < nind ; i++) {
            unsigned char c = gtype [ j *ncol + i/unitsperword];
            int y = (c>>unitsize*(i%unitsperword))&mask;
            y=(y==3)?9:y;
            cout << y;
        }
        cout << endl;
    }

}

void packedgtype::set_metadata () { 
    switch (entrytype)  {
        case PROB:
            wordsize = sizeof(char)  * 8;
            unitsize =  4;
            unitsperword = wordsize/unitsize;
            mask = 0;
            for (int i = 0 ; i < unitsize; i++)
                mask = mask |(0x1<<i);
            break;


        case GENO:
            if (!isgenotype && !hasmissing)  {
                wordsize = sizeof(char)  * 8;
                unitsize =  1;
                unitsperword = wordsize/unitsize;
                mask = 0;
                for (int i = 0 ; i < unitsize; i++)
                    mask = mask |(0x1<<i);
            } else {
                wordsize = sizeof(char)  * 8;
                unitsize = 2;
                unitsperword = wordsize/unitsize;
                mask = 0;
                for (int i = 0 ; i < unitsize; i++)
                    mask = mask |(0x1<<i);

            }
            break;


        default:
            break;

    }
}

packedgtype::etype packedgtype::format_to_type (string &entryformat) { 
    if (entryformat.compare ("geno")==0) { 
        return packedgtype::GENO;
    } else if (entryformat.compare ("geno") ==0) { 
        return packedgtype::PROB;
    }
}


string packedgtype::get_format (int e) { 
    string entryformat = ""; 
    switch (e) { 
        case GENO:  entryformat = "geno"; break;
        case PROB:  entryformat = "prob"; break;
        default:
            cerr << "Unknown format" << endl;
            exit(1);        
    }
    return entryformat;
}
