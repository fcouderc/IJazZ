#ifndef __eleIDMap__
#define __eleIDMap__

#include<iostream>
#include<map>

class eleIDMap {

 public:
  std::map<std::string,UInt_t> eleIDmap;

  eleIDMap() {

    eleIDmap["fiducial"]          =0x0001;
    eleIDmap["loose"]             =0x0002;
    eleIDmap["medium"]            =0x0004;
    eleIDmap["tight"]             =0x0008;

    eleIDmap["WP90PU"]            =0x0010;
    eleIDmap["WP80PU"]            =0x0020;
    eleIDmap["WP70PU"]            =0x0040;

    eleIDmap["loose25nsRun2"]     =0x0080;
    eleIDmap["medium25nsRun2"]    =0x0100;
    eleIDmap["tight25nsRun2"]     =0x0200;

    eleIDmap["loose50nsRun2"]     =0x0400;
    eleIDmap["medium50nsRun2"]    =0x0800;
    eleIDmap["tight50nsRun2"]     =0x1000;
  

    eleIDmap["eleID74X-veto"]    =  0x10000;
    eleIDmap["eleID74X-loose"]   =  0x20000;
    eleIDmap["eleID74X-medium"]  =  0x40000;
    eleIDmap["eleID74X-tight"]   =  0x80000;
    
    eleIDmap["eleID80X-veto" ]   =  0x4000000;
    eleIDmap["eleID80X-loose" ]  =  0x8000000;
    eleIDmap["eleID80X-medium"]  =  0x10000000;
    eleIDmap["eleID80X-tight" ]  =  0x20000000;

  }
  


};

#endif
