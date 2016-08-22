/* 
 * File:   Input_Output.h
 * Author: husain
 *
 * Created on September 7, 2011, 10:35 AM
 */

#ifndef GCN_INPUT_OUTPUT_H
#define	GCN_INPUT_OUTPUT_H

#include<vector>
#include<iomanip>

using namespace std;

/********************************************************************/
//----------Parameter Settings-----------------
#define NumCellFile 36
#define NumGene 1000

/*---------Set correlation type from below------------------------
 * either "pearson" or "CLR_pearson" or "spearman" or "mutual_information"
 * -------------------------------------------------------------*/
#define corre_type_name "CLR_pearson"

//----------For only Mutual Information-----------------
#define NumBin 10
#define NumSplineOrder 3

//--------------Set File Name with Path----------------------
#define geneExpression_in_FileName "result/MAC_Expr_36_1000.txt"
#define CoExpressionSimilarity_out_FileName "result/Pearson_MAC_Expr_36_1000.txt"
//   //-----For only CLR method--------//
#define CLR_out_FileName "result/CLR_Pearson_MAC_Expr_36_1000_unique.txt"

//---Also set Norm_CLR_out_FileName in utility class in Normalize_Data()---
//#define Norm_CLR_out_FileName "/home/jisa/NetBeansProjects/GeneCoExpressionNetwork/result/Macrophages/Norm_A_Module_pearson"
/*******************************************************************/

class GCN_Input_Output {
public:
    GCN_Input_Output();//---Constructor---

    vector<double> GeneExprValue[NumGene];//---Contains all expression level----
    vector<string> listGene;//---contains list of gene----
    
    void Read_SetData();//---read and set gene expresion information----
   
private:


};

#endif	/* INPUT_OUTPUT_H */

