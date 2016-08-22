/* 
 * File:   Input_Output.cpp
 * Author: husain
 * 
 * Created on September 7, 2011, 10:35 AM
 */

#include<cstdlib>
#include<iostream>
#include<fstream>

#include"GCN_Input_Output.h"
#include"GCN_Utility.h"

using namespace std;

GCN_Input_Output::GCN_Input_Output() {
};

void GCN_Input_Output::Read_SetData() {

    GCN_Utility GCN_Utility_;
    cout << "Data Loading....\n";
    GCN_Utility_.Write_LogFile("Data Loading....\n");
    string strTemp;
    ifstream readFile;
    int flag=0,row=-1;

    readFile.open(geneExpression_in_FileName);

    if(!readFile) {
        cout << "Invalid Input File !!!";
        GCN_Utility_.Write_LogFile("Invalid Input File !!!\n");
        exit(0);
    } else {
        while(!readFile.eof()) {
            if(flag==0 && row<NumGene-1) {
                readFile >> strTemp;
                listGene.push_back(strTemp);
                row++;flag++;
//                if(row == 260) 
//                    cout << row << '\n';
            } else if(flag>0) {
                double ExpVal;
                readFile >> ExpVal;
                
                GeneExprValue[row].push_back(ExpVal);
                flag++;
                if(flag==NumCellFile+1) {
                    flag=0;
                }
            } else readFile >> strTemp;

        }//----End of While--

    }//------end of outer if-else------

    readFile.close();

    //-------------For testing purpose that expression value is stored in correct way------------------
//    for(int i=0;i<NumGene;i++)
//    {
//        cout << listGene.at(i) << '\t';
//
//        //cout <<  GeneExprValue[NumGene-1].at(89);
//
//        for(int j=0;j<NumCellFile;j++)
//        {
//            cout << GeneExprValue[i].at(j) << '\t';
//        }
//
//        cout << '\n';
//    }

  
   
};


