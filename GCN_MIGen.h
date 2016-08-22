/* 
 * File:   GCN_MIGen.h
 * Author: husain
 *
 * Created on April 12, 2012, 1:43 PM
 */

#ifndef GCN_MIGEN_H
#define	GCN_MIGEN_H

#include<vector>

using namespace std;

class GCN_MIGen {
public:
    GCN_MIGen();
    GCN_MIGen(const GCN_MIGen& orig);
    virtual ~GCN_MIGen();
    
    void ReadnWriteMIGen();
    vector<string> FindGenoFromString(string);
    void writeMIGen_chormosomeBychromosome();
    
    //---For Josefin Script------------
    void AddConnection_in_new_col();
private:

};

#endif	/* GCN_MIGEN_H */

