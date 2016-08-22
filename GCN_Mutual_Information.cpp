/* 
 * File:   GCN_Mutual_Information.cpp
 * Author: husain
 * 
 * Created on October 3, 2011, 9:51 AM
 */

#include<iostream>
#include<algorithm>
#include<math.h>

#include "GCN_Mutual_Information.h"

using namespace std;

GCN_Mutual_Information::GCN_Mutual_Information() {
}

GCN_Mutual_Information::GCN_Mutual_Information(const GCN_Mutual_Information& orig) {
}

GCN_Mutual_Information::~GCN_Mutual_Information() {
}

void GCN_Mutual_Information::Set_bin(int p_bin) {
    this->bin = p_bin;
}

void GCN_Mutual_Information::Set_spline_order(int p_spline_order) {
    this->spline_order = p_spline_order;
}

void GCN_Mutual_Information::Set_knot_vector() {

    int M= this->bin;
    int k= this->spline_order;

    this->knot_vector.push_back(-1);//--To make 1-based vector----
    
    //------------AS new paper----------------------
    for(int i=1;i<=M+k;i++) {
        if(i<=k) {
            this->knot_vector.push_back(0);
        } else if((k<i) && (i<M+1)) {
            this->knot_vector.push_back(i-k);
        } else if(i >= M+1) {
            this->knot_vector.push_back(M-k+1);
        }
    }

    //*********To view the value of knot vector************
//    for(int i=1;i<=M+k;i++) {
//        cout << this->knot_vector.at(i) << '\n';
//    }
    //*****************************************************

}

double GCN_Mutual_Information::B_spline(int i, int k, double z) {

    double B,first,f_first,second,s_second;
    if(k==1) {
        if((this->knot_vector.at(i) <= z) && (z < this->knot_vector.at(i+1))) {
            B = 1;
        } else {
            B = 0;
        }
    } else {
        first = (knot_vector.at(i+k-1)-knot_vector.at(i));
        if(first == 0) {
            f_first = 0;
        } else {
            f_first = ((z-knot_vector.at(i))/first);
        }

        second = (knot_vector.at(i+k)-knot_vector.at(i+1));
        if(second == 0) {
            s_second = 0;
        } else {
            s_second = ((knot_vector.at(i+k)-z)/second);
        }

        B = f_first*B_spline(i,k-1,z) + s_second*B_spline(i+1,k-1,z);
    }

    return B;
}

double GCN_Mutual_Information::MI(vector<double> X, vector<double> Y) {

    weighting_coefficient_X.clear();
    weighting_coefficient_Y.clear();
    double mutual_information;
    mutual_information = Entropy(X,'x') + Entropy(Y,'y') - Entropy();
    //mutual_information = Entropy(X,'x');
    //cout << weighting_coefficient_Y.at(1).at(2);
    return mutual_information;
}

double GCN_Mutual_Information::Entropy() {

    int N = weighting_coefficient_X.at(0).size();
    double entropy_ij = 0;
    for(int i=0;i<bin;i++) {
        for(int j=0;j<bin;j++) {
            double sum=0;
            for(int k=0;k<N;k++) {
                sum = sum + (weighting_coefficient_X.at(i).at(k)*weighting_coefficient_Y.at(j).at(k));
            }
            //cout << sum << '\t';
            double prob_ij = sum/N;
            if(prob_ij==0) {
                entropy_ij=entropy_ij;
            } else {
                entropy_ij=entropy_ij + (prob_ij*log2(prob_ij));
            }
        }
        //cout << '\n';
    }
    entropy_ij = (-1)*entropy_ij;
    
    return entropy_ij;

}

double GCN_Mutual_Information::Entropy(vector<double> X,char char_tag) {

    vector<double>::const_iterator it_max,it_min;
    it_max = max_element(X.begin(),X.end());
    double Xmax = *it_max + 0.00001;
    it_min = min_element(X.begin(),X.end());
    double Xmin = *it_min;

    double entropy = 0;
    for(int i=1;i<=bin;i++) {
        vector<double> temp_weight_coef;
        temp_weight_coef.clear();
        int N = X.size();
        double sumX=0;
        for(int j=0;j<N;j++) {
            double z = ((X.at(j)-Xmin)*((bin-spline_order+1)/(Xmax-Xmin)));
            double weighting_coefficient = B_spline(i,spline_order,z);
            //cout << weighting_coefficient << '\t';
            temp_weight_coef.push_back(weighting_coefficient);
            sumX = sumX + weighting_coefficient;
        }
        if(char_tag == 'x') {
            weighting_coefficient_X.push_back(temp_weight_coef);
        } else if(char_tag == 'y') {
            weighting_coefficient_Y.push_back(temp_weight_coef);
        }
        double prob_i = sumX/N;
        if(prob_i==0) {
            entropy = entropy;
        } else {
            entropy = entropy + (prob_i*log2(prob_i));
        }
        //cout << '\n';
    }
    entropy = (-1)*entropy;
    
    return entropy;
}