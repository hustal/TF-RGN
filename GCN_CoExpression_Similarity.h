/* 
 * File:   GCN_CoExpression_Similarity.h
 * Author: husain
 *
 * Created on September 7, 2011, 2:26 PM
 */

#ifndef GCN_COEXPRESSION_SIMILARITY_H
#define	GCN_COEXPRESSION_SIMILARITY_H

#include<vector>
#include<iostream>

#include"GCN_Input_Output.h"

using namespace std;

//-----------Put OutputFile Name Here for Co Expression Similarity Measure---------------------
//#define CoExpressionSimilarity_out_FileName "/home/jisa/NetBeansProjects/GeneCoExpressionNetwork/result/Macrophages/MI_A_Module.txt"

class GCN_CoExpression_Similarity {
public:
    GCN_CoExpression_Similarity();

    double CoExpression_similarity;
    vector<double> X,Y,Z,Rx,Ry,Rz;
    bool data_ties_flag;

    vector<string> list_of_Gene;
    vector<int> degree_of_Gene;

    vector<double> GeneCoExprSimilarity[NumGene];
    vector< vector<double> > CLR_GeneCoExprSimilarity;


    //-------Filter Gene Exprssion by specific Gene List-------
    void Filter_GeneExpression_by_GeneList();
    //----parameter, gene_list_fileName,gene_expr_fileName,filtered_gene_expr_fileName--------------------
    void Filter_GeneExpression_by_GeneList(string,string,string);

    void Gene_CoExpression_Similsrity(vector<double>[], vector<string>, string);
    double Pearson_Correlation(vector<double>, vector<double>);
    double Pearson_Correlation_2(vector<double>, vector<double>);
    double Spearman_Correlation(vector<double>, vector<double>);
    //-----Filter by threshold value-------------------
    void Filter_Gene_CoExpression_Similarity(double);


    ///-----Not using now------
    void Network_info();
    void Filter_by();
    //-----------------------------
   
private:

};

#endif	/* GCN_COEXPRESSION_SIMILARITY_H */

