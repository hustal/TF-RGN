/* 
 * File:   GCN_Background_Correction.h
 * Author: husain
 *
 * Created on October 7, 2011, 3:54 PM
 */

#ifndef GCN_BACKGROUND_CORRECTION_H
#define	GCN_BACKGROUND_CORRECTION_H

#include<vector>

using namespace std;

class GCN_Background_Correction {
public:
    GCN_Background_Correction();
    GCN_Background_Correction(const GCN_Background_Correction& orig);
    virtual ~GCN_Background_Correction();

    vector< vector<double> > CLR(vector<double>[], vector<string>);
    //--------row based----------------
    vector< vector<double> > Z_score_Regulator(vector<double>[]);
    //------column based---------------
    vector< vector<double> > Z_score_Target(vector<double>[]);

    
private:

};

#endif	/* GCN_BACKGROUND_CORRECTION_H */

