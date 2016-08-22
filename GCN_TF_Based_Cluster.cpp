/* 
 * File:   TF_Based_Cluster.cpp
 * Author: husain
 * 
 * Created on October 17, 2011, 3:07 PM
 */

#include <stdlib.h>
#include <vector>
#include<iostream>
#include<fstream>

#include "GCN_TF_Based_Cluster.h"
#include "GCN_Utility.h"

using namespace std;

GCN_TF_Based_Cluster::GCN_TF_Based_Cluster() {

    Load_File_Path();
}

GCN_TF_Based_Cluster::GCN_TF_Based_Cluster(const GCN_TF_Based_Cluster& orig) {
}

GCN_TF_Based_Cluster::~GCN_TF_Based_Cluster() {
}

void GCN_TF_Based_Cluster::Load_File_Path() {

    path_TF_inputFile = "result/Top10_TF_50.txt";
    path_Network_inputFile = "result/Filtered_Unique_Normalized_Macro_Homologs50_gene_expression_0.50.txt";
    path_TF_regulated_geneList = "result/Top10_TF_regulated_gene_list_Filtered_Unique_Normalized_Macro_Homologs50_gene_expression_0.50.txt";
    path_TF_regulated_interactions = "result/Top10_TF_regulated_cluster_Filtered_Unique_Normalized_Macro_Homologs50_gene_expression_0.50.txt";
    path_TF_regulated_cluster = "result/.txt";
}

void GCN_TF_Based_Cluster::Create_Cluster() {

    vector<string> TF_list;
    vector<string> TF_regulated_gene_list;
    vector<string> TF_interaction_list;
    vector<string> TF_reduce_interaction_list;

    //----------------Read and Load TF----------------------------
    ifstream readFile;
    readFile.open(path_TF_inputFile);
    if(!readFile)
    {
        cout << "path_TF_inputFile is not found ...\n";
        exit(0);
    }
    string line;
    while(getline(readFile,line)) {
        TF_list.push_back(line);
    }
    readFile.close();
    //-------------------------------------------------------------

    TF_regulated_gene_list = TF_Regulated_Gene(TF_list);
    
    //-------------For TF_regulated Cluster up to First Neighbor------
//    TF_interaction_list = TF_Interactions(TF_regulated_gene_list);
//    TF_reduce_interaction_list = DPI(TF_regulated_gene_list,TF_interaction_list);
    
    //----------For Only TF_Direct_regulated Cluster--------------------
    TF_interaction_list = TF_Direct_Interactions_Only(TF_regulated_gene_list,TF_list);
    

    cout << "TF based Cluster Construction Successfully done.\n\n";
    cout << "Total number of TF                         : " << TF_list.size() << '\n';
    cout << "Total number of TF regulated gene          : " << TF_regulated_gene_list.size() << '\n';
    cout << "Total number of TF regulated interactions  : " << TF_interaction_list.size() << '\n';
    cout << "Total number of interactions after DPI     : " << TF_reduce_interaction_list.size() << "\n\n";
}

vector<string> GCN_TF_Based_Cluster::TF_Regulated_Gene(vector<string> TF_name) {

    vector<string> TF_regulated_gene_list;
    //-----Load TF as a TF regulated gene---------
    TF_regulated_gene_list=TF_name;

    ifstream readFile;
    readFile.open(path_Network_inputFile);
    if(!readFile)
    {
        cout << "File is not found ..." << '\n';
        exit(0);
    }

    ofstream writeFile;
    writeFile.open(path_TF_regulated_geneList);
    if(!writeFile) {
            cout << "output File is not found ..." << '\n';
            exit(0);
    }

    //------Write TF as a TF regulated Gene--------------
    for(int i=0;i<TF_regulated_gene_list.size();i++) {
        writeFile << TF_regulated_gene_list.at(i) << '\n';
    }

    //------------------Find, Load and write TF regulated Gene-------
    string line;
    int index, inner_index;
    vector<string> list_of_word;
    GCN_Utility GCN_Utility_;
    while(getline(readFile,line)) {
        list_of_word.clear();
        GCN_Utility_.Split(line,'\t',list_of_word);
        index = GCN_Utility_.Find_Gene(list_of_word.at(0),TF_name);
        if( index >=0) {
            inner_index = GCN_Utility_.Find_Gene(list_of_word.at(1),TF_regulated_gene_list);
            if(inner_index<0) {
                TF_regulated_gene_list.push_back(list_of_word.at(1));
                writeFile << list_of_word.at(1) << '\n';
            }
        } else {
            index = GCN_Utility_.Find_Gene(list_of_word.at(1),TF_name);
            if (index >=0) {
                inner_index = GCN_Utility_.Find_Gene(list_of_word.at(0),TF_regulated_gene_list);
                if(inner_index<0) {
                    TF_regulated_gene_list.push_back(list_of_word.at(0));
                    writeFile << list_of_word.at(0) << '\n';
                }
            }
        }
    }
    //---------------------------------------------------------------

    readFile.close();
    writeFile.close();
    
    return TF_regulated_gene_list;
}

vector<string> GCN_TF_Based_Cluster::TF_Interactions(vector<string> TF_regulated_gene) {

    vector<string> TF_interactions;
    ifstream readFile;
    readFile.open(path_Network_inputFile);
    if(!readFile)
    {
        cout << "File is not found ..." << '\n';
        exit(0);
    }

    ofstream writeFile;
    writeFile.open(path_TF_regulated_interactions);
    if(!writeFile) {
        cout << "Output File is Invalid !!!" << "\n";
        exit(0);
    }

    string line;
    vector<string> list_of_word;
    GCN_Utility GCN_Utility_;
    while(getline(readFile,line)) {
        list_of_word.clear();
        GCN_Utility_.Split(line,'\t',list_of_word);
        if((GCN_Utility_.Find_Gene(list_of_word.at(0),TF_regulated_gene) >= 0) && (GCN_Utility_.Find_Gene(list_of_word.at(1),TF_regulated_gene) >= 0)) {
            TF_interactions.push_back(line);
            writeFile << line << '\n';
        }
    }

    readFile.close();
    writeFile.close();
    
    return TF_interactions;
}

vector<string> GCN_TF_Based_Cluster::TF_Direct_Interactions_Only(vector<string> TF_regulated_gene, vector<string>TF_name) {
    
    vector<string> TF_interactions;
    ifstream readFile;
    readFile.open(path_Network_inputFile);
    if(!readFile)
    {
        cout << "File is not found ..." << '\n';
        exit(0);
    }

    ofstream writeFile;
    writeFile.open(path_TF_regulated_interactions);
    if(!writeFile) {
        cout << "Output File is Invalid !!!" << "\n";
        exit(0);
    }

    string line;
    vector<string> list_of_word;
    GCN_Utility GCN_Utility_;
    while(getline(readFile,line)) {
        list_of_word.clear();
        GCN_Utility_.Split(line,'\t',list_of_word);
        if(((GCN_Utility_.Find_Gene(list_of_word.at(0),TF_regulated_gene) >= 0) && (GCN_Utility_.Find_Gene(list_of_word.at(1),TF_name) >= 0)) || ((GCN_Utility_.Find_Gene(list_of_word.at(1),TF_regulated_gene) >= 0) && (GCN_Utility_.Find_Gene(list_of_word.at(0),TF_name) >= 0))) {
            TF_interactions.push_back(line);
            writeFile << line << '\n';
        }
    }

    readFile.close();
    writeFile.close();
    
    return TF_interactions;
    
}

vector<string> GCN_TF_Based_Cluster::DPI(vector<string> TF_regulated_gene, vector<string> TF_interaction) {
    
    vector< vector<double> > TF_interaction_matrix;
    vector<string> TF_reduced_interaction;

    //------Initialize TF_interaction_matrix with (-1)--------------
    int M = TF_regulated_gene.size();
    vector<double> temp_vector;
    for(int i=0;i<M;i++) {
        temp_vector.push_back(-1);
    }
    for(int i=0;i<M;i++) {
        TF_interaction_matrix.push_back(temp_vector);
    }

    //--------Fill MXM TF_interaction_matrix with respective edge value----
    TF_interaction_matrix = Fill_TF_interaction_Matrix(TF_regulated_gene,TF_interaction,TF_interaction_matrix);
    
    vector< vector<int> > discard_edge_pos;
    /*------Find the index of discard Edges from TF_interaction_matrix---
     * ----this is the part of DPI implementation-----------------------*/
    discard_edge_pos = Return_Discard_Edges_Index(TF_regulated_gene,TF_interaction_matrix);

    vector< vector<string> > discard_edge_gene;
    //------Discard Genes name from the index of TF_interaction_matrix---
    discard_edge_gene = Index_to_Gene_of_Discard_Edges(TF_regulated_gene,discard_edge_pos);

    //----------Find and write Reduced Gene Interaction list-----------------
    TF_reduced_interaction = Return_TF_Reduced_Interactions(TF_interaction,discard_edge_gene);

    return TF_reduced_interaction;

}

vector< vector<double> > GCN_TF_Based_Cluster::Fill_TF_interaction_Matrix(vector<string> TF_regulated_gene, vector<string> TF_interaction, vector<vector<double> > TF_interaction_matrix) {

    GCN_Utility GCN_Utility_;
    vector<string> list_of_word;
    string line;
    int N = TF_interaction.size();
    int x_index,y_index;
    for(int i=0;i<N;i++) {
        line = TF_interaction.at(i);
        list_of_word.clear();
        GCN_Utility_.Split(line,'\t',list_of_word);
        x_index = GCN_Utility_.Find_Gene(list_of_word.at(0),TF_regulated_gene);
        y_index = GCN_Utility_.Find_Gene(list_of_word.at(1),TF_regulated_gene);
        //------------convert string to double-------------
        char char_val[list_of_word.at(2).size()];
        for(int l=0;l<list_of_word.at(2).size();l++) {
            char_val[l] = list_of_word.at(2)[l];
        }
        double double_value = atof(char_val);

        TF_interaction_matrix.at(x_index).at(y_index) = double_value;
        TF_interaction_matrix.at(y_index).at(x_index) = double_value;
    }

    return TF_interaction_matrix;
}

vector< vector<int> > GCN_TF_Based_Cluster::Return_Discard_Edges_Index(vector<string> TF_regulated_gene, vector<vector<double> > TF_interaction_matrix) {

    int M = TF_regulated_gene.size();
    double x,y,z,minxy,min;
    int pos;
    vector< vector<int> > discard_edge_pos;
    discard_edge_pos.clear();
    vector<int> pair_pos;
    for(int i=0;i<M-2;i++) {
        for(int j=i+1;j<M-1;j++) {
            if(TF_interaction_matrix.at(i).at(j)>-1) {
                for(int k=j+1;k<M;k++) {
                    if(TF_interaction_matrix.at(i).at(k)>-1) {
                        if(TF_interaction_matrix.at(j).at(k)>-1) {
                            x=TF_interaction_matrix.at(i).at(j);
                            y=TF_interaction_matrix.at(i).at(k);
                            z=TF_interaction_matrix.at(j).at(k);
                            if(x<y) {
                                minxy = x;
                                pos = 0;
                            } else {
                                minxy = y;
                                pos = 1;
                            }

                            if(minxy < z) {
                                min = minxy;
                            } else {
                                min = z;
                                pos = 2;
                            }

                            pair_pos.clear();
                            if(pos==0) {
                                pair_pos.push_back(i);
                                pair_pos.push_back(j);
                            }
                            else if(pos==1) {
                                pair_pos.push_back(i);
                                pair_pos.push_back(k);
                            }
                            else if(pos==2) {
                                pair_pos.push_back(j);
                                pair_pos.push_back(k);
                            }

                            discard_edge_pos.push_back(pair_pos);
                        }

                    }

                }

            }
            //cout << TF_interaction_matrix.at(i).at(j) << '\t';
        }
        //cout << '\n';
    }

    return discard_edge_pos;
}

vector< vector<string> > GCN_TF_Based_Cluster::Index_to_Gene_of_Discard_Edges(vector<string> TF_regulated_gene, vector<vector<int> > discard_edge_pos) {

    vector< vector<string> > discard_edge_gene;
    discard_edge_gene.clear();
    vector<string> pair_gene;

    int discard_edge_num = discard_edge_pos.size();
    for(int i=0;i<discard_edge_num;i++) {
        pair_gene.clear();
        pair_gene.push_back(TF_regulated_gene.at(discard_edge_pos.at(i).at(0)));
        pair_gene.push_back(TF_regulated_gene.at(discard_edge_pos.at(i).at(1)));
        discard_edge_gene.push_back(pair_gene);
    }

    return discard_edge_gene;
}

vector<string> GCN_TF_Based_Cluster::Return_TF_Reduced_Interactions(vector<string> TF_interaction, vector<vector<string> > discard_edge_gene) {

    ofstream writeFile;
    writeFile.open(path_TF_regulated_cluster);
    if(!writeFile) {
        cout << "path_TF_regulated_cluster Output File is Invalid !!!\n";
        exit(0);
    }

    vector<string> TF_reduced_interaction;
    GCN_Utility GCN_Utility_;
    int N=TF_interaction.size();
    int discard_edge_num = discard_edge_gene.size();
    vector<string> list_of_word;
    string line;
    bool flag;
    for(int i=0;i<N;i++) {
        line = TF_interaction.at(i);
        list_of_word.clear();
        GCN_Utility_.Split(line,'\t',list_of_word);
        flag = false;
        for(int j=0;j<discard_edge_num;j++) {
            if(list_of_word.at(0) == discard_edge_gene.at(j).at(0)) {
                if(list_of_word.at(1) == discard_edge_gene.at(j).at(1)) {
                    flag = true;break;
                }
            } else if(list_of_word.at(0) == discard_edge_gene.at(j).at(1)) {
                if(list_of_word.at(1) == discard_edge_gene.at(j).at(0)) {
                    flag = true;break;
                }
            }
        }

        if(flag == false) {
            TF_reduced_interaction.push_back(TF_interaction.at(i));
            writeFile << TF_interaction.at(i) << '\n';
        }

    }

    writeFile.close();
    return TF_reduced_interaction;
}