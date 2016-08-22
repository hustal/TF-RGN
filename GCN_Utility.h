/* 
 * File:   GCN_Utility.h
 * Author: husain
 *
 * Created on September 7, 2011, 4:51 PM
 */

#ifndef GCN_UTILITY_H
#define	GCN_UTILITY_H

#include<vector>

using namespace std;

class GCN_Utility {
public:
    GCN_Utility();

    //------Bubble sort X and arrange Y according to sorted X------
    void b_Sort_XY(vector<double>&, vector<double>&);
    void b_Sort_XY(vector<int>&, vector<string>&);
    
    vector<double> b_Sort_X(vector<double>);//--Bubble sort---
    
    /*---Search the value into sorted vector and return the index of it, 
     * if found and return -1 if not found
     * -----------------------------------------------------------------*/
    int b_search(double, vector<double>);
    
    /*---Produce and return Rank vector for sorted vector data------
     * --It will be applicable for both ties (repeated) and non ties data--
     * --bool type return the data is repeated(ties) or not---
     *----------------------------------------------------------------*/
    vector<double> Rank_of_data(vector<double>, bool&);

    void Split(const string&, char, vector<string>&);

    int Find_Gene(string,vector<string>);
    
    void Check_Duplication();//--in a list of gene---

    vector<string> Load_Gene_List();
    vector<string> Load_Gene_List(string);//---File Name---
    vector<string> Load_Gene_List(string,string);//---File Name---
    bool CheckGene_in_List(string,vector<string>);
    void Check_Two_Gene_List();//---check gene_list_1 with gene_list_2 and find how many genes are missing in gene_list_1----
    void Common_Genes();
    void Check_Files();
    bool fexists(const char*);
    
    /*----From UG_gene_ID to UG_gene_symbol-----*/
    /*Input: Set All UG_gene_ID, all UG_gene_Symbol in two files
      and set your gene_ID into third file
      Output: you will get your gene_Symbol file*/
    void UG_gene_ID_to_Symbol(string, string, string);
    
    
    double A_Mean(vector<double>);
    /*-----paramter of standard deviation ----------------
    -------list of value and mean of values--------------*/
    double Standard_dev(vector<double>, double);

     //-----For a set of Data from Whole (like for A_Module_TEML genes)-----------------
    void Normalize_Data(vector<string>, vector<double>);
    //-----For whole Data Set-----------------
    void Normalize_Data();

    void Network_info();

    void Write_LogFile(const char*);
    
    //------------String Trim Function------------------
    string LTrim(string);
    string RTrim(string);
    string Trim(string);
    
    double string_To_double(string);
    int string_To_int(string);
    
private:

    void Load_File_Path_of_Network_info();
    const char *path_Gene_list_inputFile, *path_TEML_inputFile, *path_Network_inputFile, *path_Network_info_outputFile, *path_TF_gene_inputFile;

};

#endif	/* GCN_UTILITY_H */

