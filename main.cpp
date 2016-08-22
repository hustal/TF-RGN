/* 
 * File:   main.cpp
 * Author: husain
 *
 * Created on September 6, 2011, 3:50 PM
 */

#include<cstdlib>
#include<vector>
#include<iomanip>
#include<iostream>
#include <fstream>
#include<string.h>
#include<algorithm>

#include "GCN_Input_Output.h"
#include "GCN_CoExpression_Similarity.h"
#include "GCN_Utility.h"
#include "GCN_Mutual_Information.h"
#include "GCN_TF_Based_Cluster.h"
#include "GCN_MIGen.h"
#include "GCN_Master_Table.h"

using namespace std;

/*--------------------Before run this---------------------
 *   Set all Parameters in GCN_Input_Output header file
 *------------------------------------------------------*/

int main() {
    
    GCN_Input_Output GCN_Input_Output_;
    GCN_Input_Output_.Read_SetData();

    GCN_CoExpression_Similarity GCN_CoExpression_Similarity_;
        
    /*
     * Before using this following function you have to first set the 
     * output file name and then you have to choose correlation
     * type either "pearson" or "spearman" or "mutual_information"
     * into GCN_Input_Output header file
     *  */
    GCN_CoExpression_Similarity_.Gene_CoExpression_Similsrity(GCN_Input_Output_.GeneExprValue,GCN_Input_Output_.listGene,corre_type_name);

     /*----Filter gene co-expression similarity
      *    for a set of genes
     -----like for A_module or TEML genelist---------------------------------*/
//    GCN_CoExpression_Similarity_.Filter_GeneExpression_by_GeneList();
    
    string geneList_input_fileName = "gene_list_150.txt";
    string geneExpr_input_fileName = "CLR_Pearson_MAC_Expr_36_1000_unique.txt";
    string filtered_geneExpr_output_fileName = "MAC_150_NET.txt";
    GCN_CoExpression_Similarity_.Filter_GeneExpression_by_GeneList(geneList_input_fileName,geneExpr_input_fileName,filtered_geneExpr_output_fileName);

    //------Set the threshold value in this parameter----------------
    GCN_CoExpression_Similarity_.Filter_Gene_CoExpression_Similarity(0.40);

    GCN_Utility GCN_Utility_;
    //GCN_Utility_.Check_Duplication();
    //GCN_Utility_.Common_Genes();
    //GCN_Utility_.Check_Two_Gene_List();
    //GCN_Utility_.Check_Files();
    
    //---For Whole Data------------
    //GCN_Utility_.Normalize_Data();
    
    //---For limited data like A_Module_TEML------
    vector<string> Edge;
    vector<double> Edge_value;
    GCN_Utility_.Normalize_Data(Edge,Edge_value);
    
    //GCN_Utility_.UG_gene_ID_to_Symbol("HG-U133_Plus_2_Affy_ProbeSet.txt","HG-U133_Plus_2_Affy_Gene.txt","HGU133Plus2_Hs_UG_mapping_unique.txt");
    

    /*---*****************************************************************
       Its Completely Independent Class, to reconstruct a TF based
       Cluster based on DPI algorithm.
       Before using this class you have to first set the input output file
       path into its Load_File_Path() method.
//     *****************************************************************--*/
       GCN_TF_Based_Cluster GCN_TF_Based_Cluster_;
       GCN_TF_Based_Cluster_.Create_Cluster();
  
       GCN_Utility GCN_Utility_;
       GCN_Utility_.Network_info();

     
    GCN_MIGen GCN_MIGen_;
    //GCN_MIGen_.ReadnWriteMIGen();
    //GCN_MIGen_.writeMIGen_chormosomeBychromosome();
//    vector<string> t;
//    t=GCN_MIGen_.FindGenoFromString("09-234:40:");
    //GCN_MIGen_.AddConnection_in_new_col();
   
     
    GCN_Master_Table  GCN_Master_Table_;
    
    /*======================================================================
     * tissue no's are : AAW=6,WB=7,CL=8,IMA=9,LIVER=10,MAC=11,SM=12,SF=13,
     * VF=14 
     * 
     * Set master_table_filename, tissue_cel_list_filename, tissue_no
     */
    //GCN_Master_Table_.Retrieve_Patient_ID("master_table","cel_file_list_WB",7);
       
    
    return 0;
}

