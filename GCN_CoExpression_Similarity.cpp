/* 
 * File:   GCN_CoExpression_Similarity.cpp
 * Author: husain
 * 
 * Created on September 7, 2011, 2:26 PM
 */

#include<iostream>
#include<vector>
#include<math.h>
#include<cstdlib>
#include<fstream>

#include "GCN_CoExpression_Similarity.h"
#include "GCN_Input_Output.h"
#include "GCN_Utility.h"
#include "GCN_Mutual_Information.h"
#include "GCN_Background_Correction.h"

using namespace std;

GCN_CoExpression_Similarity::GCN_CoExpression_Similarity() {

}

 void GCN_CoExpression_Similarity::Gene_CoExpression_Similsrity(vector<double> GeneExprValue[], vector<string> listGene, string Corre_type) {

    GCN_Utility GCN_Utility_;
    ofstream writeFile;
    //writeFile.open(CoExpressionSimilarity_out_FileName);
    writeFile.open(CoExpressionSimilarity_out_FileName);

    if(!writeFile) {
        cout << "Output File, CoExpressionSimilarity_out_FileName is Invalid !!!" << "\n";
        GCN_Utility_.Write_LogFile("Output File, CoExpressionSimilarity_out_FileName is Invalid !!!\n");
        exit(0);
    }

    if(Corre_type == "pearson") {

        cout << "Pearson Correlation have been choosen for co-expression similarity...\n";
        GCN_Utility_.Write_LogFile("Pearson Correlation have been choosen for co-expression similarity...\n");
        for(int i=0;i<NumGene;i++) {
            for(int j=i+1;j<NumGene;j++) {
                X = GeneExprValue[i];
                Y = GeneExprValue[j];
                CoExpression_similarity = Pearson_Correlation(X,Y);
                GeneCoExprSimilarity[i].push_back(CoExpression_similarity);
                writeFile << listGene.at(i) << '\t' << listGene.at(j) << '\t' << CoExpression_similarity << '\n';
                //writeFile << listGene.at(i) << '\t' << listGene.at(j) << '\t' << setprecision(3)  << CoExpression_similarity << '\n';
            }
            cout << "Gene CoExpression Similarity Done " << i+1 << " of " << NumGene << '\n';
        }
    } else if(Corre_type == "CLR_pearson") {

        cout << "CLR with Pearson Correlation have been choosen for co-expression similarity...\n";
        GCN_Utility_.Write_LogFile("CLR with Pearson Correlation have been choosen for co-expression similarity...\n");
        for(int i=0;i<NumGene;i++) {
            for(int k=0;k<i;k++) {
                GeneCoExprSimilarity[i].push_back(GeneCoExprSimilarity[k].at(i));
                writeFile << listGene.at(i) << '\t' << listGene.at(k) << '\t' << GeneCoExprSimilarity[i].at(k) << '\n';
                }
            for(int j=i;j<NumGene;j++) {
                X = GeneExprValue[i];
                Y = GeneExprValue[j];
                CoExpression_similarity = Pearson_Correlation(X,Y);
                GeneCoExprSimilarity[i].push_back(CoExpression_similarity);
                writeFile << listGene.at(i) << '\t' << listGene.at(j) << '\t' << CoExpression_similarity << '\n';
                //writeFile << listGene.at(i) << '\t' << listGene.at(j) << '\t' << setprecision(3)  << CoExpression_similarity << '\n';
            }
            cout << "Gene CoExpression Similarity Done " << i+1 << " of " << NumGene << '\n';
        }

        cout << "\nBackground Correction......\n";
        GCN_Utility_.Write_LogFile("\nBackground Correction......\n");
        GCN_Background_Correction GCN_Background_Correction_;
        CLR_GeneCoExprSimilarity = GCN_Background_Correction_.CLR(GeneCoExprSimilarity, listGene);

    } else if(Corre_type == "spearman") {

        cout << "Spearman Correlation have been choosen for co-expression similarity...\n";
        for(int i=0;i<NumGene;i++) {
                for(int j=i+1;j<NumGene;j++) {
                    X = GeneExprValue[i];
                    Y = GeneExprValue[j];
                    CoExpression_similarity=Spearman_Correlation(X,Y);
                    GeneCoExprSimilarity[i].push_back(CoExpression_similarity);
                    writeFile << listGene.at(i) << '\t' << listGene.at(j) << '\t' << CoExpression_similarity << '\n';
                    //writeFile << listGene.at(i) << '\t' << listGene.at(j) << '\t' << setprecision(3)  << CoExpression_similarity << '\n';
                }
                cout << "Gene CoExpression Similarity Done " << i+1 << " of " << NumGene << '\n';
            }
    } else if(Corre_type == "mutual_information") {
        cout << "Mutual Information have been choosen for co-expression similarity...\n";

        GCN_Mutual_Information GCN_Mutual_Information_;
        GCN_Mutual_Information_.Set_bin(NumBin);
        GCN_Mutual_Information_.Set_spline_order(NumSplineOrder);
        GCN_Mutual_Information_.Set_knot_vector();

        for(int i=0;i<NumGene;i++) {
            for(int k=0;k<i;k++) {
                GeneCoExprSimilarity[i].push_back(GeneCoExprSimilarity[k].at(i));
                writeFile << listGene.at(i) << '\t' << listGene.at(k) << '\t' << GeneCoExprSimilarity[i].at(k) << '\n';
                }
            for(int j=i;j<NumGene;j++) {
                X = GeneExprValue[i];
                Y = GeneExprValue[j];
                CoExpression_similarity = GCN_Mutual_Information_.MI(X,Y);
                GeneCoExprSimilarity[i].push_back(CoExpression_similarity);
                writeFile << listGene.at(i) << '\t' << listGene.at(j) << '\t' << CoExpression_similarity << '\n';
                //writeFile << listGene.at(i) << '\t' << listGene.at(j) << '\t' << setprecision(3)  << CoExpression_similarity << '\n';
                }
            cout << "Gene CoExpression Similarity Done " << i+1 << " of " << NumGene << '\n';
            }

        cout << "\nBackground Correction......\n";
        GCN_Background_Correction GCN_Background_Correction_;
        CLR_GeneCoExprSimilarity = GCN_Background_Correction_.CLR(GeneCoExprSimilarity, listGene);

    } else cout << "Given Correlation Type Not Found ....";

    writeFile.close();
}

double GCN_CoExpression_Similarity::Pearson_Correlation(vector<double> X, vector<double> Y) {

    //-------------For Test Purpose that data is passing in correct manner-----
    /*
    for(int i=0;i<NumCellFile;i++)
        cout << X.at(i) << '\t' << Y.at(i) << '\n';
    */

    double pearsonCorrelation;
    double standardDeviations, covariant=0, sumX = 0, sumY = 0;

    int n = X.size(); //----Total number of Elements----
    //------------Calculate Mean--------------------
    for(int i=0;i<n;i++)
    {
        sumX = sumX + X.at(i);
        sumY = sumY + Y.at(i);
    }
    //cout << "num : " << NumCellFile << '\t' << X.size() << '\n';
    double meanX = sumX/n;
    double meanY = sumY/n;

    sumX=0;sumY=0;
    for(int i=0;i<n;i++)
    {
        //-----------For Standard Deviations-----------
        sumX = sumX + pow((X.at(i)- meanX),2);
        sumY = sumY + pow((Y.at(i)- meanY),2);
        //--------------For Covariants-----------------
        covariant = covariant + ((X.at(i)- meanX)*(Y.at(i)- meanY));
    }
    standardDeviations = sqrt(sumX * sumY);

    pearsonCorrelation = covariant/standardDeviations;

    return pearsonCorrelation;
}

double GCN_CoExpression_Similarity::Pearson_Correlation_2(vector<double> X, vector<double> Y) {

    GCN_Utility GCN_Utility_;
    double pearsonCorrelation;
    double meanX,meanY,stdX,stdY;
    int N = X.size();
    meanX = GCN_Utility_.A_Mean(X);
    meanY = GCN_Utility_.A_Mean(Y);
    stdX = GCN_Utility_.Standard_dev(X,meanX);
    stdY = GCN_Utility_.Standard_dev(Y,meanY);
    
    double sum=0;
    for(int i=0;i<N;i++) {
        sum = sum + ((X.at(i)-meanX)/stdX)*((Y.at(i)-meanY)/stdY);
    }

    pearsonCorrelation = sum/N;
    return pearsonCorrelation;
}

double GCN_CoExpression_Similarity::Spearman_Correlation(vector<double> X, vector<double> Y) {

    double spearman_Correlation;
    GCN_Utility utility_;
    utility_.b_Sort_XY(X,Y);
    Z=utility_.b_Sort_X(Y);//----Z is the sorted data of Y----

    data_ties_flag=false;
    Rx=utility_.Rank_of_data(X,data_ties_flag);
    Rz=utility_.Rank_of_data(Z,data_ties_flag);//---Rz is according to Z(sorted Y)----

    //--------Ry, Rank of Y according to X------
    vector<double>::iterator it;
    Ry.clear();
    for(it=Y.begin();it!=Y.end();it++) {
        int rank_pos = utility_.b_search(*it,Z);
        if(rank_pos > -1) {
            Ry.push_back(Rz.at(rank_pos));
        } else {
            cout << "Error have been found !!!";
        }
    }

//    cout << '\n';
//    for(int i=0;i<X.size();i++) {
//        cout << X.at(i) << '\t' << Y.at(i) << '\t' << Z.at(i) << '\t' << Rx.at(i) << '\t' << Rz.at(i) << '\t' << Ry.at(i) << '\n';
//    }

    /*
     * Now two options, either the ties could be found or not.
     */
    if(data_ties_flag) {
        spearman_Correlation = Pearson_Correlation(Rx,Ry);
    } else {
        double sum=0;
        int n = Rx.size();//---Total number of Elements-----
        for(int i=0;i<n;i++) {
            sum = sum + pow((Rx.at(i)-Ry.at(i)),2);
        }
        spearman_Correlation = 1 - ((6*sum)/(n*(pow(n,2)-1)));
    }

    return spearman_Correlation;
}

void GCN_CoExpression_Similarity::Filter_Gene_CoExpression_Similarity(double filter_by_value) {
    
    ifstream readFile;
    readFile.open("result/CLR_Pearson_MV_ALL_MAC_Illumina_20h_C_0.3_unique_A_TEML.txt");
    if(!readFile) {
        cout << "Input File is not Found...." << '\n';
        exit(0);
    }

    ofstream writeFile;
    writeFile.open("result/Filtered_CLR_Pearson_MV_ALL_MAC_Illumina_20h_C_0.3_unique_A_TEML_0.40.txt");
    if(!writeFile) {
        cout << "Output File is not Found...." << '\n';
        exit(0);
    }

    string line;
    vector<string> list_of_word;
    GCN_Utility GCN_Utility_;
    while(getline(readFile,line)) {
        list_of_word.clear();
        GCN_Utility_.Split(line,'\t',list_of_word);
        char char_val[list_of_word.at(2).size()];
        for(int i=0;i<list_of_word.at(2).size();i++) {
            char_val[i] = list_of_word.at(2)[i];
        }
        double double_value = atof(char_val);
        if(fabs(double_value) > filter_by_value) {
             writeFile << line << '\n';
        }
    }

    readFile.close();
    writeFile.close();

    //GCN_Utility GCN_Utility_;
    //GCN_Utility_.Network_info();
}


void GCN_CoExpression_Similarity::Filter_GeneExpression_by_GeneList() {

    GCN_Utility GCN_Utility_;
    /*----Set file name with path in Load_Gene_List() function----
     * ---From where it will load gene list-------*/
    vector<string> gene_list;
    cout << "Loading A Module TEML List ...\n";
    GCN_Utility_.Write_LogFile("Loading A Module TEML List ...\n");
    gene_list = GCN_Utility_.Load_Gene_List();
    int counter=0;
    string line,strTemp;
    ifstream readFile;
    readFile.open("result/TF_regulated_cluster_Macro_A_Module_TEML_0.40.txt");
    if(!readFile)
    {
        cout << "Input File, CLR_Pearson_STAGE_Liver_unique.txt is not found ..." << '\n';
        GCN_Utility_.Write_LogFile("Input File, CLR_Pearson_STAGE_Liver_unique.txt is not found ...");
        exit(0);
    }

    ofstream writeFile;
    writeFile.open("result/Second_Top_Hub_TF_regulated_cluster_Macro_A_Module_TEML_0.40.txt");
    if(!writeFile)
    {
        cout << "Result File, CLR_Pearson_STAGE_Liver_A_Module_TEML.txt is not found ..." << '\n';
        exit(0);
    }

    cout << "\nFiltering List ...\n";
    GCN_Utility_.Write_LogFile("\nFiltering List ...\n");
    vector<string> list_of_word;
    string gene1,gene2;
    int gcounter=0;
    while(getline(readFile,line)) {
        list_of_word.clear();
        GCN_Utility_.Split(line,'\t',list_of_word);
        gene1 = GCN_Utility_.Trim(list_of_word.at(0));
        gene2 = GCN_Utility_.Trim(list_of_word.at(1));
        //if((GCN_Utility_.Find_Gene(gene1,gene_list)>-1) || (GCN_Utility_.Find_Gene(gene2,gene_list)>-1)) {
        if((GCN_Utility_.Find_Gene(gene1,gene_list)>-1) && (GCN_Utility_.Find_Gene(gene2,gene_list)>-1)) {
            writeFile << line << '\n';
            gcounter++;
            cout << gcounter << '\n';
            //GCN_Utility_.Write_LogFile(string(gcounter));
            //if (gcounter == 100) break;
        }

    }

    readFile.close();
    writeFile.close();

    cout << "\nFiltering Completion...\n";
    GCN_Utility_.Write_LogFile("\nFiltering Completion...\n");
}

void GCN_CoExpression_Similarity::Filter_GeneExpression_by_GeneList(string geneList_fileName,string geneExpr_fileName,string filtered_geneExpr_fileName) {

    GCN_Utility GCN_Utility_;
    /*----Set file name with path in Load_Gene_List() function----
     * ---From where it will load gene list-------*/
    vector<string> gene_list;
    cout << "Loading Specific Gene List ...\n";
    GCN_Utility_.Write_LogFile("Loading Specific Gene List ...\n");
    gene_list = GCN_Utility_.Load_Gene_List(geneList_fileName);
    
    
    int counter=0;
    string line,strTemp;
    ifstream readFile;

    string geneExpr_filePath = "result/"+geneExpr_fileName;
    readFile.open(geneExpr_filePath.c_str());
    if(!readFile)
    {
        cout << "Input File, " << geneExpr_fileName << " is not found ...\n";
        GCN_Utility_.Write_LogFile("Input File, CLR_Pearson_STAGE_Liver_unique.txt is not found ...");
        exit(0);
    }

    ofstream writeFile;
    string filtered_geneExpr_filePath = "result/"+filtered_geneExpr_fileName;
    writeFile.open(filtered_geneExpr_filePath.c_str());
    if(!writeFile)
    {
        cout << "Result File, " << filtered_geneExpr_fileName <<" is not found ...\n";
        exit(0);
    }

    cout << "\nFiltering Gene Expression List ...\n";
    GCN_Utility_.Write_LogFile("\nFiltering List ...\n");
    vector<string> list_of_word;
    string gene1,gene2;
    int gcounter=0;
    while(getline(readFile,line)) {
        list_of_word.clear();
        GCN_Utility_.Split(line,'\t',list_of_word);
        gene1 = GCN_Utility_.Trim(list_of_word.at(0));
        gene2 = GCN_Utility_.Trim(list_of_word.at(1));
        //if((GCN_Utility_.Find_Gene(gene1,gene_list)>-1) || (GCN_Utility_.Find_Gene(gene2,gene_list)>-1)) {
        if((GCN_Utility_.Find_Gene(gene1,gene_list)>-1) && (GCN_Utility_.Find_Gene(gene2,gene_list)>-1)) {
            writeFile << line << '\n';
            gcounter++;
            cout << gcounter << '\n';
            //GCN_Utility_.Write_LogFile(string(gcounter));
            //if (gcounter == 100) break;
        }

    }

    readFile.close();
    writeFile.close();

    cout << "\nFiltering Completion...\n";
    GCN_Utility_.Write_LogFile("\nFiltering Completion...\n");
}

void GCN_CoExpression_Similarity::Network_info() {
    
    ifstream readFile;
    readFile.open("/home/jisa/NetBeansProjects/GeneCoExpressionNetwork/result/Filtered75_P_CoExpressionSimilarity_gene2expression_level_77.txt");
    //readFile.open("/home/jisa/NetBeansProjects/GeneCoExpressionNetwork/result/testfile.txt");
    if(!readFile) {
        cout << "File is not Found...." << '\n';
        exit(0);
    }

//    ofstream writeFile;
//    writeFile.open("Sorted_Gene_Node_Filtered75_P_CoExpressionSimilarity_gene2expression_level_77.txt");
//    if(!writeFile) {
//        cout << "File is not Found...." << '\n';
//        exit(0);
//    }

    string line;
    //vector<string> list_of_Gene;
    vector<string> temp_list;
    //vector<int> degree_of_Gene;
    GCN_Utility GCN_Utility_;
    while(getline(readFile,line)) {
        temp_list.clear();
        GCN_Utility_.Split(line,'\t',temp_list);
        int i=0;
        while(i<2) {
            int gene_pos = GCN_Utility_.Find_Gene(temp_list.at(i),list_of_Gene);
            if(gene_pos == -1) {
                list_of_Gene.push_back(temp_list.at(i));
                degree_of_Gene.push_back(1);
            } else {
                degree_of_Gene.at(gene_pos) = degree_of_Gene.at(gene_pos) + 1;
            }
            i++;
        }
    }

//    GCN_Utility_.b_Sort_XY(degree_of_Gene,list_of_Gene);
//
//    for(int i=0;i<list_of_Gene.size();i++) {
//        writeFile << list_of_Gene.at(i) << '\t' << degree_of_Gene.at(i) << '\n';
//
//    }

    readFile.close();
    //writeFile.close();
}

void GCN_CoExpression_Similarity::Filter_by() {

    ifstream readFile;
    readFile.open("/home/jisa/NetBeansProjects/GeneCoExpressionNetwork/result/Filtered75_P_CoExpressionSimilarity_gene2expression_level_77.txt");
    if(!readFile) {
        cout << "File is not Found...." << '\n';
        exit(0);
    }

    ofstream writeFile;
    writeFile.open("Filtered75_2_P_CoExpressionSimilarity_gene2expression_level_77.txt");
    if(!writeFile) {
        cout << "File is not Found...." << '\n';
        exit(0);
    }

    string line;
    vector<string> list_of_word;
    GCN_Utility GCN_Utility_;
    while(getline(readFile,line)) {
        list_of_word.clear();
        GCN_Utility_.Split(line,'\t',list_of_word);
        int gene_pos1 = GCN_Utility_.Find_Gene(list_of_word.at(0),list_of_Gene);
        int gene_pos2 = GCN_Utility_.Find_Gene(list_of_word.at(1),list_of_Gene);
        if((degree_of_Gene.at(gene_pos1) > 25) || (degree_of_Gene.at(gene_pos2) > 25)) {
            char char_val[list_of_word.at(2).size()];
            for(int i=0;i<list_of_word.at(2).size();i++) {
                char_val[i] = list_of_word.at(2)[i];
            }
            double double_value = atof(char_val);
            if(fabs(double_value) > 0.90) {
                 writeFile << line << '\n';
            }

         
        } else {
            writeFile << line << '\n';
        }
    }

    readFile.close();
    writeFile.close();
}