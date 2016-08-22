/* 
 * File:   TF_Based_Cluster.h
 * Author: husain
 *
 * Created on October 17, 2011, 3:07 PM
 */

#ifndef GCN_TF_BASED_CLUSTER_H
#define	GCN_TF_BASED_CLUSTER_H

#include<vector>

using namespace std;

class GCN_TF_Based_Cluster {
public:
    GCN_TF_Based_Cluster();
    GCN_TF_Based_Cluster(const GCN_TF_Based_Cluster& orig);
    virtual ~GCN_TF_Based_Cluster();
   
    void Create_Cluster();
    
private:
    const char *path_Network_inputFile, *path_TF_inputFile, *path_TF_regulated_geneList, *path_TF_regulated_interactions, *path_TF_regulated_cluster;
    void Load_File_Path();

    vector<string> TF_Regulated_Gene(vector<string>);
    vector<string> TF_Interactions(vector<string>);
    vector<string> TF_Direct_Interactions_Only(vector<string>, vector<string>);
    //------Data Processing Inequality (DPI)------
    vector<string> DPI(vector<string>, vector<string>);
    vector< vector<double> > Fill_TF_interaction_Matrix(vector<string>, vector<string>, vector< vector<double> >);
    vector< vector<int> > Return_Discard_Edges_Index(vector<string>, vector< vector<double> >);
    vector< vector<string> > Index_to_Gene_of_Discard_Edges(vector<string>, vector< vector<int> >);
    vector<string> Return_TF_Reduced_Interactions(vector<string>, vector< vector<string> >);

};

#endif	/* TF_BASED_CLUSTER_H */

