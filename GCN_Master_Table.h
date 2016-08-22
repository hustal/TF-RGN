/* 
 * File:   GCN_Master_Table.h
 * Author: husain
 *
 * Created on May 3, 2012, 5:02 PM
 */

#ifndef GCN_MASTER_TABLE_H
#define	GCN_MASTER_TABLE_H

using namespace std;

class GCN_Master_Table {
public:
    GCN_Master_Table();
    GCN_Master_Table(const GCN_Master_Table& orig);
    virtual ~GCN_Master_Table();
    
    void Retrieve_Patient_ID(string,string,int);
    
private:

};

#endif	/* GCN_MASTER_TABLE_H */

