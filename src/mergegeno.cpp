mergegeno::mergegeno (int argc, char *argv[]) {
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
	p->get_string ("geno1",geno1,true);
	p->get_string ("snp1",snp1,true);
	p->get_string ("ind1",ind1,true);
    string inputformat1 = "eigenstrat";
    p->get_string ("inputformat1", inputformat, false);
    g1 = new genotype (snp1, ind1, geno1, inputformat1);


	string geno2, snp2, ind2;
	p->get_string ("geno2",geno2,true);
	p->get_string ("snp2",snp2,true);
	p->get_string ("ind2",ind2,true);
	g2 = new genotype (snp2, ind2, geno2);

}
