#include "cldata.h"

cldata::cldata (int argc, char *argv[], const char *optString, const option *longOpts) 
    : data(argc, argv, optString, longOpts) {
        string pfile;
        get_string ("parameter",pfile,true);
        data d (pfile, configMap);
        for (map<string,string>::iterator i = d.configMap.begin(); i!=d.configMap.end(); i++){
            configMap[i->first] = i->second;
        }
    }
