/* 
 * File:   GCN_Mutual_Information.h
 * Author: husain
 *
 * Created on October 3, 2011, 9:51 AM
 */

#ifndef GCN_MUTUAL_INFORMATION_H
#define	GCN_MUTUAL_INFORMATION_H

#include<vector>

using namespace std;

class GCN_Mutual_Information {
public:
    GCN_Mutual_Information();
    GCN_Mutual_Information(const GCN_Mutual_Information& orig);
    virtual ~GCN_Mutual_Information();

    int bin,spline_order;
    vector<int> knot_vector;
    vector< vector<double> > weighting_coefficient_X,weighting_coefficient_Y;

    void Set_bin(int);
    void Set_spline_order(int);
    void Set_knot_vector();
    double B_spline(int, int, double);

    double Entropy(vector<double>,char);
    double Entropy();
    double MI(vector<double>, vector<double>);
   

private:

};

#endif	/* GCN_MUTUAL_INFORMATION_H */

