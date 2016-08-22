/* 
 * File:   GCN_Utility.cpp
 * Author: husain
 * 
 * Created on September 7, 2011, 4:51 PM
 */

#include<vector>
#include<string>
#include<math.h>
#include<iostream>
#include<fstream>
#include <algorithm>
#include<stdio.h>

#include "GCN_Utility.h"
#include "GCN_Input_Output.h"
#include "GCN_TF_Based_Cluster.h"

using namespace std;

GCN_Utility::GCN_Utility() {
}

void GCN_Utility::b_Sort_XY(vector<double>& X, vector<double>& Y) {

    for(int i=0;i<X.size()-1;i++) {
        for(int j=i+1;j<X.size();j++) {
            if(X.at(i)>X.at(j)) {
                double tempX=X.at(i), tempY=Y.at(i);
                X.at(i)=X.at(j);Y.at(i)=Y.at(j);
                X.at(j)=tempX;Y.at(j)=tempY;
            }
        }
    }


}

void GCN_Utility::b_Sort_XY(vector<int>& X, vector<string>& Y) {

     for(int i=0;i<X.size()-1;i++) {
        for(int j=i+1;j<X.size();j++) {
            if(X.at(i)>X.at(j)) {
                int tempX=X.at(i);string tempY=Y.at(i);
                X.at(i)=X.at(j);Y.at(i)=Y.at(j);
                X.at(j)=tempX;Y.at(j)=tempY;
            }
        }
    }

}

vector<double> GCN_Utility::b_Sort_X(vector<double> X) {

    for(int i=0;i<X.size()-1;i++) {
        for(int j=i+1;j<X.size();j++) {
            if(X.at(i)>X.at(j)) {
                double tempX=X.at(i);
                X.at(i)=X.at(j);
                X.at(j)=tempX;
            }
        }
    }

    return X;
}

int GCN_Utility::b_search(double key, vector<double> X) {

    int min,max,mid;
    min=0;max=X.size()-1;
    while(min <= max) {
        mid = (min + max) / 2;
        if(key > X.at(mid)) {
            min = mid + 1;
        } else if (key < X.at(mid)) {
            max = mid - 1;
        } else {
            return mid;
        }
    }
    return -1;
}

vector<double> GCN_Utility::Rank_of_data(vector<double> X, bool& tie_flag) {
    
    vector<double> Rx;
    for(int i=0;i<X.size();i++) {
        if(i+1 < X.size()) {
            if(X.at(i)==X.at(i+1)) {
                tie_flag=true;
                double Rsum=i+1;int counter=1;
                Rsum=Rsum+(i+1+counter); counter++;
                while(i+counter<X.size()) {
                    if(X.at(i)==X.at(i+counter)) {
                        Rsum=Rsum+(i+1+counter); counter++;
                    } else {
                        break;
                    }
                }
                double new_rank = Rsum/counter;
                for(int k=0;k<counter;k++) {
                    Rx.push_back(new_rank);
                }
                i=i+(counter-1);
            } else {
                Rx.push_back(i+1);
            }
        } else if(i+1 == X.size()) {
            Rx.push_back(i+1);
        }
    }
    return Rx;
}

void GCN_Utility::Split(const string& line, char split_char, vector<string>& listString) {

    listString.clear();
    string::size_type i=0;
    string::size_type j = line.find(split_char);

    while(j!= string::npos)
    {
        listString.push_back(line.substr(i,j-i));
        i=++j;
        j=line.find(split_char,j);
        if(j==string::npos)
        {
            listString.push_back(line.substr(i,line.length()));
        }

    }
}

int GCN_Utility::Find_Gene(string key_gene, vector<string> gene_list) {
    
    for(int i=0;i<gene_list.size();i++) {
        if(key_gene==gene_list.at(i)) {
            return i;
        }
    }

    return -1;
}

double GCN_Utility::A_Mean(vector<double> list_val) {

    double sum = 0;
    int N = list_val.size();
    for(int i=0;i<N;i++) {
        sum = sum + list_val.at(i);
    }
    double a_mean = sum/N;

    return a_mean;
}

double GCN_Utility::Standard_dev(vector<double> list_val, double mean_val) {

    double sum = 0;
    int N = list_val.size();
    for(int i=0;i<N;i++) {
        sum = sum + pow((list_val.at(i)-mean_val),2);
    }
    double std_dev = sqrt(sum/N);

    return std_dev;
}

vector<string> GCN_Utility::Load_Gene_List() {

    vector<string> list_of_Gene;
    ifstream readFile;
    //----Specify the file name with path----
    readFile.open("result/Second_Highest_Top_Hub_gene_Macro.txt");
    if(!readFile)
    {
        cout << "Gene List File is not found ..." << '\n';
        return list_of_Gene;
    }

    while(!readFile.eof())
    {
        string strTemp;
        readFile >> strTemp;
        strTemp = Trim(strTemp);
        if(strTemp!="")
         list_of_Gene.push_back(strTemp);
    }

    //cout << list_of_Gene.size() << '\n';
    readFile.close();

    return list_of_Gene;
}

vector<string> GCN_Utility::Load_Gene_List(string fileName) {

    vector<string> list_of_Gene;
    ifstream readFile;
    //----Specify the file name with path----
    string filePath = "result/"+fileName;
    const char * path = filePath.c_str();
    readFile.open(path);
    if(!readFile)
    {
        cout << "Gene List File is not found ..." << '\n';
        return list_of_Gene;
    }

    while(!readFile.eof())
    {
        string strTemp;
        readFile >> strTemp;
        strTemp = Trim(strTemp);
        if(strTemp!="")
         list_of_Gene.push_back(strTemp);
    }

    //cout << list_of_Gene.size() << '\n';
    readFile.close();

    return list_of_Gene;
}

vector<string> GCN_Utility::Load_Gene_List(string fileName, string tag) {

    vector<string> list_of_Gene;
    ifstream readFile;
    //----Specify the file name with path----
    string filePath = "result/"+fileName;
    const char * path = filePath.c_str();
    readFile.open(path);
    if(!readFile)
    {
        cout << "Gene List File is not found ..." << '\n';
        return list_of_Gene;
    }

    int counter=0;
    string strTemp;
    while(getline(readFile,strTemp))
    {
        
        //readFile >> strTemp;counter++;cout << strTemp << '\n';
        //strTemp = Trim(strTemp);
        //if(strTemp!="")
         list_of_Gene.push_back(strTemp);
        //else cout << strTemp << counter << "+tag" ;
         
    }

    //cout << list_of_Gene.size() << '\n';
    readFile.close();

    return list_of_Gene;
}

bool GCN_Utility::CheckGene_in_List(string strGene, vector<string> gene_list) {

    bool flag=false;
    int N = gene_list.size();
    for(int i=0;i<N;i++)
    {
        if(gene_list.at(i)==strGene)
        {
            flag=true;
            break;
        }
    }
    return flag;
}

void GCN_Utility::Check_Two_Gene_List() {
    
    ofstream writeFile;
    
    vector<string> gene_list1,gene_list2;
       
    gene_list1.clear();
    gene_list1 = Load_Gene_List("A_Module_TEML_Gene_illumina_ver_232.txt");
    
    gene_list2.clear();
    gene_list2 = Load_Gene_List("A_Module_Gene_unique.txt");
    
    writeFile.open("result/C_missing_gene_list.txt");
    if(!writeFile) {
        cout << "Output File is Invalid !!!\n";
        exit(0);
    }
    int missing_counter=0;
    for(int i=0;i<gene_list2.size();i++) {
        if(!CheckGene_in_List(gene_list2[i],gene_list1)) {
            writeFile << gene_list2[i] << '\n';
            missing_counter++;
        }
    }
    writeFile << missing_counter;
    writeFile.close();
    
}

void GCN_Utility::Normalize_Data(vector<string> strEdge, vector<double> edge_value) {

    GCN_Utility GCN_Utility_;
    if(edge_value.size()==0) {

        cout << "\nLoading Started for Normalization...\n";
        GCN_Utility_.Write_LogFile("\nLoading Started for Normalization...\n");
        ifstream readFile;

        readFile.open("result/CLR_Pearson_STAGE_AAW_Homologs30.txt");
        if(!readFile) {
            cout << "Input File, CLR_Pearson_STAGE_AAW_Homologs30.txt is not Found...." << '\n';
            GCN_Utility_.Write_LogFile("Input File, CLR_Pearson_STAGE_AAW_Homologs30.txt is not Found....\n");
            exit(0);
        }

        ofstream writeFile2;
        writeFile2.open("result/CLR_Pearson_A_Module_TEML.txt");

        if(!writeFile2) {
            cout << "Output File, CLR_Pearson_A_Module_TEML.txt is Invalid !!!" << "\n";
            GCN_Utility_.Write_LogFile("Output File, CLR_Pearson_A_Module_TEML.txt is Invalid !!!\n");
            exit(0);
        }

        string line2;
        vector<string> A_TEML;
        while(getline(readFile,line2)) {
            A_TEML.push_back(line2);
        }

        cout << A_TEML.size() << '\n';
        readFile.close();

        readFile.open("result/CLR_Pearson_STAGE_Liver_A_Module_TEML.txt");
        if(!readFile) {
            cout << "Input File, CLR_Pearson_STAGE_Liver_A_Module_TEML.txt is not Found...." << '\n';
            GCN_Utility_.Write_LogFile("Input File, CLR_Pearson_STAGE_Liver_A_Module_TEML.txt is not Found....\n");
            exit(0);
        }

        string line, tempedge;int gcounter=0;
        vector<string> list_of_word;
        while(getline(readFile,line)) {
            list_of_word.clear();
            GCN_Utility_.Split(line,'\t',list_of_word);
            if((Find_Gene(list_of_word.at(0),A_TEML)>-1) && (Find_Gene(list_of_word.at(1),A_TEML)>-1)) {
                tempedge = list_of_word.at(0) + '\t' + list_of_word.at(1);
                strEdge.push_back(tempedge);
                char char_val[list_of_word.at(2).size()];
                for(int i=0;i<list_of_word.at(2).size();i++) {
                    char_val[i] = list_of_word.at(2)[i];
                }
                double double_value = atof(char_val);
                edge_value.push_back(double_value);
                writeFile2 << line << '\n';
                tempedge="";
                gcounter++;
                cout << gcounter << '\n';
                //if(gcounter == 100 ) break;
            }

        }

        readFile.close();
        writeFile2.close();

    }

    cout << "\nNormalization is running .....\n";
    GCN_Utility_.Write_LogFile("\nNormalization is running .....\n");
    //vector<double> norm_edge_val;
    ofstream writeFile;
    writeFile.open("result/Norm_CLR_Pearson_Liver_A_Module_TEML.txt");

    if(!writeFile) {
        cout << "Output File, Norm_CLR_Pearson_Liver_A_Module_TEML.txt is Invalid !!!" << "\n";
        GCN_Utility_.Write_LogFile("Output File, Norm_CLR_Pearson_Liver_A_Module_TEML.txt is Invalid !!!\n");
        exit(0);
    }
    //----------Find Max and Min element from array---------------
    vector<double>::const_iterator it_max,it_min;
    it_max = max_element(edge_value.begin(),edge_value.end());
    it_min = min_element(edge_value.begin(),edge_value.end());

    double max_val = *it_max, min_val = *it_min;
    
    int M = edge_value.size();
    cout << M << '\t' << strEdge.size() << '\n';
    for(int i=0;i<M;i++) {
        double dblNormValue = (edge_value.at(i)-min_val) / (max_val-min_val);
        //norm_edge_val.push_back(dblNormValue);
        writeFile << strEdge.at(i) << '\t' << dblNormValue << '\n';
    }

    writeFile.close();
    cout << "Normalization is successfully done...\n";
    GCN_Utility_.Write_LogFile("Normalization is successfully done...\n");
}

void GCN_Utility::Normalize_Data() {
    
    GCN_Utility GCN_Utility_;
    ifstream readFile;
    vector<string> strEdge;
    vector<double> edge_value;
    
    readFile.open("result/2.CLR_Pearson_Final_Monocyte_Illumina_20h_C_unique_A_TEML.txt");
    if(!readFile) {
        cout << "Input File, CLR_Pearson_STAGE_AAW_Homologs30.txt is not Found...." << '\n';
        GCN_Utility_.Write_LogFile("Input File, CLR_Pearson_STAGE_AAW_Homologs30.txt is not Found....\n");
        exit(0);
    }
    
    ofstream writeFile;
    writeFile.open("result/3.Normalized_CLR_Pearson_Final_Monocyte_Illumina_20h_C_unique_A_TEML.txt");

    if(!writeFile) {
        cout << "Output File, CLR_Pearson_STAGE_AAW_Homologs30.txt is Invalid !!!" << "\n";
        GCN_Utility_.Write_LogFile("Output File, CLR_Pearson_STAGE_AAW_Homologs30.txt is Invalid !!!\n");
        exit(0);
    }
    
    string line, tempedge;
    vector<string> list_of_word;
    while(getline(readFile,line)) {
        list_of_word.clear();
        GCN_Utility_.Split(line,'\t',list_of_word);
        tempedge = list_of_word.at(0) + '\t' + list_of_word.at(1);
        strEdge.push_back(tempedge);
        char char_val[list_of_word.at(2).size()];
        for(int i=0;i<list_of_word.at(2).size();i++) {
            char_val[i] = list_of_word.at(2)[i];
        }
        double double_value = atof(char_val);
        edge_value.push_back(double_value);
        tempedge="";
    }
    
    cout << "\nNormalization is running .....\n";
    GCN_Utility_.Write_LogFile("\nNormalization is running .....\n");
    //vector<double> norm_edge_val;
   
    //----------Find Max and Min element from array---------------
    vector<double>::const_iterator it_max,it_min;
    it_max = max_element(edge_value.begin(),edge_value.end());
    it_min = min_element(edge_value.begin(),edge_value.end());

    double max_val = *it_max, min_val = *it_min;
    
    ///---For make same scale week 30,40 and 50---------------
    //double max_val = 6.68513, min_val = 0;
    
    cout << "max " << max_val << '\n';
    cout << "min " << min_val << '\n';
    //---------------------------------------------------------
    
    int M = edge_value.size();
    //---------------
    cout << M << '\t' << strEdge.size() << '\n';
    for(int i=0;i<M;i++) {
        double dblNormValue = (edge_value.at(i)-min_val) / (max_val-min_val);
        
        //norm_edge_val.push_back(dblNormValue);
        writeFile << strEdge.at(i) << '\t' << dblNormValue << '\n';
    }

    readFile.close();
    writeFile.close();
    cout << "Normalization is successfully done...\n";
    GCN_Utility_.Write_LogFile("Normalization is successfully done...\n");
    
}

void GCN_Utility::Load_File_Path_of_Network_info() {

    path_Gene_list_inputFile = "result/human_genes_W50.txt";
    path_Network_inputFile = "result/sub_net_W50.txt";
    path_Network_info_outputFile = "result/Net_Info_sub_net_W50_2.txt";
    path_TF_gene_inputFile = "result/Homologs50_TF_262.txt";
    path_TEML_inputFile = "result/TEML_Gene_113.txt";
    
}

void GCN_Utility::Network_info() {

    //---Set File path in the following method-----
    Load_File_Path_of_Network_info();

    ifstream readFile;
    //----Specify the Gene_List file name with path----
    readFile.open(path_Gene_list_inputFile);
    if(!readFile)
    {
        cout << "path_Gene_list_inputFile, File is not found ..." << '\n';
        exit(0);
    }
    
    vector<string> Gene_List;
    vector<int> Gene_degree;
    
    string line;
    while (getline(readFile,line)) {
        Gene_List.push_back(line);
        Gene_degree.push_back(0);
    }
    readFile.close();
    
    //---------Load TEML genes-------------
    readFile.open(path_TEML_inputFile);
    vector<string> unique_TEML_gene;
    if(!readFile) {
        cout << "path_unique_TEML_inputFile, File is not found ..." << '\n';
        //exit(0);
    } else {
        while (getline(readFile,line)) {
            unique_TEML_gene.push_back(line);
        }
        readFile.close();
    }
    
    
    ifstream readFile2;
    //-----------Specify the source Filtered File----------
    readFile2.open(path_Network_inputFile);
    if(!readFile2)
    {
        cout << "path_Network_inputFile, File is not found ..." << '\n';
        exit(0);
    } else {

        vector<string> list_of_word;
        while(getline(readFile2,line)) {
            list_of_word.clear();
            Split(line,'\t',list_of_word);
            int pos = Find_Gene(list_of_word.at(0),Gene_List);
            if(pos >=0) {
                Gene_degree.at(pos) = Gene_degree.at(pos) + 1;
            } 
            pos = Find_Gene(list_of_word.at(1),Gene_List);
            if(pos >=0) {
                Gene_degree.at(pos) = Gene_degree.at(pos) + 1;
            }
        }

        readFile2.close();
    }
    
    ofstream writeFile;
    //---------------------Specify the Network Info Output File---------
    writeFile.open(path_Network_info_outputFile);
    if(!writeFile) {
        cout << "path_Network_info_outputFile, File is not Found...." << '\n';
        exit(0);
    }

    //------------------------Check TF Gene------------------
    readFile.open(path_TF_gene_inputFile);
    if(!readFile)
    {
        cout << "path_TF_gene_inputFile, File is not found ..." << '\n';
        //exit(0);
    } else {
        //vector<string> TF;
        writeFile << "TF's Connections : \n\n";
        while (getline(readFile,line)) {
            //TF.push_back(line);
            int pos = Find_Gene(line,Gene_List);
            if(pos>=0) {
                writeFile << Gene_List.at(pos) << '\t' << Gene_degree.at(pos) << '\n';
            }
        }
        readFile.close();
    }
   

    //---------------Genes are in network---------------------------
    writeFile << "\nGene List in the network: \n\n";
    for(int i=0;i<Gene_List.size();i++) {
        if(Gene_degree.at(i)>0) {
            writeFile << Gene_List.at(i) << '\t' << Gene_degree[i] << '\n';
        } 
    }
    
    
    //---------------------------Missing Gene List---------------------
    writeFile << "\nMissing Gene List: \n\n";
    vector<string> str_A_module_gene_list, str_TEML_gene_list;
    vector<int> str_TEML_gene_degree;
    str_A_module_gene_list.clear();
    str_TEML_gene_list.clear();
    str_TEML_gene_degree.clear();
    int N=Gene_List.size();
    int counter = 0,A_Module_counter = 0, TEML_counter = 0;
    for(int i=0;i<N;i++) {
        if(Gene_degree.at(i)==0) {
            writeFile << Gene_List.at(i) << '\n';
            counter++;
            
        } else if(Gene_degree.at(i)>0){
            /*if(Find_Gene(Gene_List.at(i),unique_A_Module_gene) >= 0) {
                str_A_module_gene_list.push_back(Gene_List.at(i));
                A_Module_counter++;
                
            } else*/ if(Find_Gene(Gene_List.at(i),unique_TEML_gene) >= 0) {
                str_TEML_gene_list.push_back(Gene_List.at(i));
                str_TEML_gene_degree.push_back(Gene_degree[i]);
                TEML_counter++;
               
            }
           
        } 
        
    }
    
    writeFile << "\nTEML Gene List: \n\n";
    
    for(int i=0;i<str_TEML_gene_list.size();i++) {
        writeFile << str_TEML_gene_list.at(i) << '\t' << str_TEML_gene_degree[i] << '\n';
    }
    
    writeFile << "\n\nTotal missing Gene : " << counter << '\n';
    writeFile << "Genes From TEML : " << TEML_counter << '\n';
    writeFile << "\nTotal Gene was in network : " << N << '\n';
    writeFile << "Total Genes in this Cluster : " << N-counter << '\n';


    //-----------------Highest Degree Genes upto LDB2-------------
    b_Sort_XY(Gene_degree,Gene_List);

    //-----------------Highest Hub Genes-------------
    //----------Assign hub connection in here,
    int top_hub_conn = 10, second_top_hub_conn = 5;
    int top_hub_gene_counter = 0, second_top_hub_gene_counter = 0, lowest_hub_gene_counter = 0;
    int top_hub_A_Module_gene_counter = 0, top_hub_TEML_gene_counter = 0, top_hub_Common_gene_counter = 0;
    writeFile << "\nHighest Hubs containing gene: \n\n";
    int i;
    for(i=N-1;i>-1;i--) {
        if(top_hub_conn < Gene_degree[i]) {
            writeFile << Gene_List[i] << '\t' << Gene_degree[i] << '\n';
            //-----
            /*if(Find_Gene(Gene_List[i],unique_A_Module_gene)>-1) {
                top_hub_A_Module_gene_counter++;
            } else if(Find_Gene(Gene_List[i],unique_TEML_gene)) {
                top_hub_TEML_gene_counter++;
            } else top_hub_Common_gene_counter++;*/
            //---------------
            top_hub_gene_counter++;
        } else break;
    }
    
    writeFile << "\nSecond Highest Hubs containing gene: \n\n";
    for(i=i;i>-1;i--) {
        if(second_top_hub_conn <= Gene_degree[i]) {
            writeFile << Gene_List[i] << '\t' << Gene_degree[i] << '\n';
            second_top_hub_gene_counter++;
        } else break;
    }

    writeFile << "\nLowest Hubs containing gene: \n\n";
    for(i=i;i>-1;i--) {
        if(Gene_degree[i] > 0) {
            writeFile << Gene_List[i] << '\t' << Gene_degree[i] << '\n';
            lowest_hub_gene_counter++;
        } else break;
    }
    
    writeFile << "\n\nTotal Highest Top Hub gene : " << top_hub_gene_counter << "\n";
    writeFile << "Total Second Highest Top Hub gene : " << second_top_hub_gene_counter << "\n";
    writeFile << "Total Lowest Hub gene : " << lowest_hub_gene_counter << "\n";
    
    writeFile.close();
}


void GCN_Utility::Write_LogFile(const char* message) {

    FILE* logFile = fopen("result/logfile.txt","a");
    fprintf(logFile,"%s\n", message);
    fclose(logFile);
}

string GCN_Utility::LTrim(string strLTrim) {
    
   while((*strLTrim.begin()==' ') || (*strLTrim.begin()=='\t') || (*strLTrim.begin()=='\r')) {
        strLTrim.erase(strLTrim.begin());
    }
    
    return strLTrim;
}

string GCN_Utility::RTrim(string strRTrim) {
    
    string::iterator it;
    it = strRTrim.end();
    it--;
    
    while((*it==' ') || (*it=='\t') || (*it=='\r')) {
        strRTrim.erase(it);
        it--;
    }
    
    return strRTrim;
}

string GCN_Utility::Trim(string strTrim) {
    
    if(strTrim!="") {
    strTrim = LTrim(strTrim);
    strTrim = RTrim(strTrim);
    }
    return strTrim;
}

void GCN_Utility::Check_Duplication() {
    
    ifstream readFile;
    ofstream writeFile;
    
    writeFile.open("result/duplicate.txt");
    if(!writeFile) {
        cout << "Output File is Invalid !!!\n";
        exit(0);
    }
    
    readFile.open("result/CAD_gene_symbol.txt");
    if(!readFile) {
        cout << "Input File is not Found !!!\n";
        exit(0);
    } else {
        vector<string> list_gene;
        string line;
        while(getline(readFile,line)) {
            line = Trim(line);
            list_gene.push_back(line);
        }
        readFile.close();
        
        int num_gene = list_gene.size();
        bool flag = false;
        for(int i=0;i<num_gene;i++) {
            for(int j=0;j<i;j++) {
                if(list_gene[i] == list_gene[j]) {
                    flag = true;
                    break;
                }
            }
            if(flag) {
                writeFile << "Position :" << i << "\t" << list_gene[i] << '\n';
                flag = false;
            }
        }
        writeFile.close();
    }
}

void GCN_Utility::Common_Genes() {
    
    ifstream readFile;
    ofstream writeFile;
    
    vector<string> gene_list_network1,gene_list_network2;
    readFile.open("result/unique_All_gene_final.txt");
    if(!readFile) {
        cout << "Input File is Invalid !!!\n";
        exit(0);
    } else {
        gene_list_network1.clear();
        string line;
        while(getline(readFile,line)) {
            line = Trim(line);
            gene_list_network1.push_back(line);
        }
    }
    readFile.close();
    
    readFile.open("result/unique_gene_list.txt");
    if(!readFile) {
        cout << "Input File is Invalid !!!\n";
        exit(0);
    } else {
        gene_list_network2.clear();
        string line;
        while(getline(readFile,line)) {
            line = Trim(line);
            gene_list_network2.push_back(line);
        }
    }
    readFile.close();
    
    writeFile.open("result/common_list_R_C.txt");
    if(!writeFile) {
        cout << "Output File is Invalid !!!\n";
        exit(0);
    }
    
    int len_gene_list_network1, len_gene_list_network2;
    len_gene_list_network1 = gene_list_network1.size();
    len_gene_list_network2 = gene_list_network2.size();
    
    int common_counter = 0;
    if(len_gene_list_network1 < len_gene_list_network2) {
        
        for(int i=0;i<len_gene_list_network1;i++) {
            for(int j=0;j<len_gene_list_network2;j++) {
                if(gene_list_network1[i] == gene_list_network2[j]) {
                    writeFile << gene_list_network1[i] << '\n';
                    //writeFile << j+1 << '\n';
                    common_counter++;
                    break;
                }
            }
        }
    } else {
        
        for(int i=0;i<len_gene_list_network2;i++) {
            for(int j=0;j<len_gene_list_network1;j++) {
                if(gene_list_network2[i] == gene_list_network1[j]) {
                    writeFile << gene_list_network2[i] << '\n';
                    common_counter++;
                    break;
                }
            }
        }
    }
    writeFile << common_counter;
    writeFile.close();
}

double GCN_Utility::string_To_double(string strVal) {

    char *ch[strVal.size()];
    *ch = (char*)strVal.c_str();
    double dblVal = atof(*ch);
    
    return dblVal;
}

int GCN_Utility::string_To_int(string strVal) {

    char *ch[strVal.size()];
    *ch = (char*)strVal.c_str();
    double intVal = atoi(*ch);
    
    return intVal;
}

void GCN_Utility::Check_Files() {
    
    ifstream readFile;
    readFile.open("result/cel_file_list_IMA.txt");
    ofstream writeFileFound,writeFileNotFound;
    writeFileFound.open("result/found_file.txt");
    writeFileNotFound.open("result/not_found_file.txt");
    
    string line, file_path;
    while(getline(readFile,line)) {
        line = Trim(line);
        if(line!="") {
            file_path = "/home/husain/Projects/Expression_level/STAGE_IMA_CEL_79/" + line;
            if(fexists(file_path.c_str())) {
                //cout << line << "\tfound" <<'\n';
                writeFileFound << line << '\n';
            } else {
                //cout << line << "\tNot found" <<'\n';
                writeFileNotFound << line << '\n';
            }
        }
    }
    
    readFile.close();
    writeFileFound.close();
    writeFileNotFound.close();
    
}

bool GCN_Utility::fexists(const char *filename) {
    
    ifstream ifile(filename);
    return ifile;
}

void GCN_Utility::UG_gene_ID_to_Symbol(string All_UG_ID, string All_UG_Symbol, string given_UG_ID) {
    
    vector<string> All_ID, All_symbol, given_ID;
    ofstream writeFile, writeFile2;
    
    All_ID = Load_Gene_List(All_UG_ID);
    All_symbol = Load_Gene_List(All_UG_Symbol,"tag");
    given_ID = Load_Gene_List(given_UG_ID);
    
    int gene_pos;
    writeFile.open("result/UG_given_expr_gene_Symbol.txt");
    if(!writeFile) {
        cout << "Output File is Invalid !!!\n";
        exit(0);
    }
    writeFile2.open("result/UG_missing_given_expr_gene_Symbol.txt");
    if(!writeFile2) {
        cout << "Output File is Invalid !!!\n";
        exit(0);
    }
    
       
    for(int i=0;i<given_ID.size();i++) {
        
        gene_pos = Find_Gene(given_ID[i],All_ID);
        if(gene_pos>-1) {
            writeFile << All_symbol[gene_pos] << '\n';
        } else {
            writeFile2 << given_ID[i] << '\n';
        }
        cout << given_ID.size() << "/" <<  i+1 << '\n';
    }
    
    writeFile.close();
    writeFile2.close();
    
   
}