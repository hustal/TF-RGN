/* 
 * File:   GCN_Background_Correction.cpp
 * Author: husain
 * 
 * Created on October 7, 2011, 3:54 PM
 */

#include<iostream>
#include<math.h>
#include<fstream>

#include "GCN_Background_Correction.h"
#include "GCN_Utility.h"
#include "GCN_Input_Output.h"

using namespace std;

GCN_Background_Correction::GCN_Background_Correction() {
};

GCN_Background_Correction::GCN_Background_Correction(const GCN_Background_Correction& orig) {
};

GCN_Background_Correction::~GCN_Background_Correction() {
};

vector< vector<double> > GCN_Background_Correction::CLR(vector<double> mi[], vector<string> listgene) {

    GCN_Utility GCN_Utility_;
    ofstream writeFile;
    writeFile.open(CLR_out_FileName);

    if(!writeFile) {
        cout << "Output File, CLR_out_FileName is Invalid !!!" << "\n";
        GCN_Utility_.Write_LogFile("Output File, CLR_out_FileName is Invalid !!!\n");
        //exit(0);
    }

    vector< vector<double> > mi_CLR;
    vector< vector<double> > z_score_reg = Z_score_Regulator(mi);
    vector< vector<double> > z_score_tar = Z_score_Target(mi);

    //--------for Normalization <for list of edge and their value>-----
    //vector<string> Edge;
    //vector<double> Edge_value;
    //--------------------------------------------------------

    cout << "\nCLR Finalizing ....\n";
    GCN_Utility_.Write_LogFile("\nCLR Finalizing ....\n");
    vector<double> temp_CLR;
    int N=z_score_reg.at(0).size();
    for(int i=0;i<N;i++) {
        temp_CLR.clear();
        for(int j=0;j<N;j++) {
            double dbl_clr_val = sqrt(pow(z_score_reg.at(i).at(j),2) + pow(z_score_tar.at(i).at(j),2));
            //-----For complete Matrix--------------
            //temp_CLR.push_back(dbl_clr_val);
            if(j>i) {
                writeFile << listgene.at(i) << '\t' << listgene.at(j) << '\t' << dbl_clr_val << '\n';
                //-----To store only unique value from Matrix----------
                //temp_CLR.push_back(dbl_clr_val);

                //----------Only for Normalization-------------------
                //string strEdge = listgene.at(i) + '\t' + listgene.at(j);
                //cout << strEdge << '\n';
                //Edge.push_back(strEdge);
                //Edge_value.push_back(dbl_clr_val);
            }
        }
        //mi_CLR.push_back(temp_CLR);
    }

    //--------------for normaliation value-----------------
    //GCN_Utility GCN_Utility_;
    //GCN_Utility_.Normalize_Data(Edge,Edge_value);

    writeFile.close();
    cout << "CLR Successfully Done...\n";
    GCN_Utility_.Write_LogFile("CLR Successfully Done...\n");
    return mi_CLR;
}

vector< vector<double> > GCN_Background_Correction::Z_score_Regulator(vector<double> mi_reg[]) {
    //-----------According to Row value--------------
    int N = mi_reg[0].size();
    vector< vector<double> > z_score_reg;
    z_score_reg.clear();
    vector<double> temp_score;
    double mean_val,std_dev;
    GCN_Utility GCN_Utility_;

    for(int i=0;i<N;i++) {
        temp_score.clear();
        mean_val = GCN_Utility_.A_Mean(mi_reg[i]);
        std_dev = GCN_Utility_.Standard_dev(mi_reg[i],mean_val);
        for(int j=0;j<N;j++) {
            double temp_val = (mi_reg[i].at(j) - mean_val)/std_dev;
            temp_score.push_back(max(0.0,temp_val));
        }
        z_score_reg.push_back(temp_score);
    }

    return z_score_reg;
}

vector< vector<double> > GCN_Background_Correction::Z_score_Target(vector<double> mi_target[]) {
    //-----------According to Column value--------------
    int N = mi_target[0].size();
    vector< vector<double> > z_score_tar;
    z_score_tar.clear();
    vector<double> temp_score;
    double mean_val,std_dev;
    GCN_Utility GCN_Utility_;

    vector<double> temp_column;
    for(int i=0;i<N;i++) {
        temp_column.clear();
        mean_val=0;
        for(int k=0;k<N;k++) {
            double temp_mi_target_val=mi_target[k].at(i);
            temp_column.push_back(temp_mi_target_val);
            mean_val = mean_val+temp_mi_target_val;
        }
        temp_score.clear();
        mean_val = mean_val/N;
        std_dev = GCN_Utility_.Standard_dev(temp_column,mean_val);
        for(int j=0;j<N;j++) {
            double temp_val = (temp_column.at(j) - mean_val)/std_dev;
            temp_score.push_back(max(0.0,temp_val));
        }
        
        if(i==0) {
            for(int l=0;l<N;l++) {
                z_score_tar.push_back(temp_column);
                z_score_tar.at(l).at(i)=temp_score.at(l);
            }
        } else {
            for(int l=0;l<N;l++) {
                z_score_tar.at(l).at(i)=temp_score.at(l);
            }
        }
    }

    return z_score_tar;
}

