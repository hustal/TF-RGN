/* 
 * File:   GCN_Master_Table.cpp
 * Author: husain
 * 
 * Created on May 3, 2012, 5:02 PM
 */

#include <stdlib.h>
#include <vector>
#include<iostream>
#include<fstream>

#include "GCN_Master_Table.h"
#include "GCN_Utility.h"

using namespace std;

GCN_Master_Table::GCN_Master_Table() {
}

GCN_Master_Table::GCN_Master_Table(const GCN_Master_Table& orig) {
}

GCN_Master_Table::~GCN_Master_Table() {
}

void GCN_Master_Table::Retrieve_Patient_ID(string master_table_file_name, string tissue_cel_list_file_name, int tissue_no) {
    
    GCN_Utility GCN_Utility_;
    
    vector<string> CEL_list, Patient_ID_list, tissue_CEL_ID;
    ifstream readFile;
    ofstream writeFile;
    //---------Loading CEL list--------------------------------------------
    string file_path = "result/" + tissue_cel_list_file_name + ".txt";
    readFile.open(file_path.c_str());
    
    if(!readFile) {
        cout << "tissue CEL_list File is not found !!!";
        exit(0);
    }
    
    string line;
    vector<string> split_string;
    bool header = true;
    while(getline(readFile,line)) {
        
        line = GCN_Utility_.Trim(line);
        if(line!="") {
            if(header) header = false;
            else {
                GCN_Utility_.Split(line,'.',split_string);
                tissue_CEL_ID.push_back(split_string[0]);
            }
            
        }
    }
    
    readFile.close();
    cout << tissue_CEL_ID.size() << '\n';
    
    file_path = "result/" + master_table_file_name + ".txt";
    readFile.open(file_path.c_str());
    
    if(!readFile) {
        cout << "Master File is not found !!!";
        exit(0);
    }
    
    //-------------Loading Master file info------------------------
    header = true;
    while(getline(readFile,line)) {
        
        //line = GCN_Utility_.Trim(line);
        if(line!="") {
            if(header) header = false;
            else {
                GCN_Utility_.Split(line,'\t',split_string);
                if(split_string.size()==16) {
                    Patient_ID_list.push_back(split_string[0]);
                    CEL_list.push_back(split_string[tissue_no]);
                } else {
                    cout << "Error in file reading !!!\n";
                    exit(0);
                }
                
            }
            
        }
    }
    
    readFile.close();
    
    cout << Patient_ID_list.size() << '\t' << CEL_list.size();
    
    file_path = "result/" + tissue_cel_list_file_name + "_pat_id.txt";
    writeFile.open(file_path.c_str());
    if(!writeFile) {
        cout << "output file invalid !!!\n";
        exit(0);
    }
    
    int pos=-1;
    bool first_item=true;
    for(int i=0;i<tissue_CEL_ID.size();i++) {
        pos = GCN_Utility_.Find_Gene(tissue_CEL_ID[i],CEL_list);
        if(pos>-1) {
            if(first_item) {
                writeFile << Patient_ID_list[pos];
                first_item=false;
            } else {
                writeFile << '\n' << Patient_ID_list[pos];
            }
            
        }
    }
    
    writeFile.close();
    
}