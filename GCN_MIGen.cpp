/* 
 * File:   GCN_MIGen.cpp
 * Author: husain
 * 
 * Created on April 12, 2012, 1:43 PM
 */
#include <stdlib.h>
#include <vector>
#include<iostream>
#include<fstream>
#include <algorithm>
#include <sstream>

#include "GCN_MIGen.h"
#include "GCN_Utility.h"

using namespace std;

GCN_MIGen::GCN_MIGen() {
}

GCN_MIGen::GCN_MIGen(const GCN_MIGen& orig) {
}

GCN_MIGen::~GCN_MIGen() {
}

void GCN_MIGen::ReadnWriteMIGen() {
    
    GCN_Utility GCN_Utility_;
    ifstream readFile;
    readFile.open("result/MIGen_freeze_121007_ig_fullresults_final3_logreg.txt");
    if(!readFile)
    {
        cout << "Input File, is not found ...\n";
        exit(0);
    }

    ofstream writeFile, write_complex_file,write_complex_geno_file, write_report_file;
    writeFile.open("result/varified_MIGen.txt");
    write_complex_file.open("result/complex_MIGen.txt");
    write_complex_geno_file.open("result/complex_geno_MIGen.txt");
    write_report_file.open("result/report_MIGen.txt");
    if(!writeFile || !write_complex_file || !write_complex_geno_file || !write_report_file)
    {
        cout << "Output File, is not found ...\n";
        exit(0);
    }

    int line_counter=-1, match_line_counter=0, complex_line_counter=0, complex_geno_line_counter=0;
    
    vector<string> list_of_field;
    vector<string> GENO, GENO_A, GENO_U;
    string line, match_line;
    while(getline(readFile,line)) {
        line = GCN_Utility_.Trim(line);
        if(line!="") {
            line_counter++;//--Total line counter except first line--
            GCN_Utility_.Split(line,' ',list_of_field);
            if(list_of_field.size()==21) {//---Check total field number---
                if(line_counter==0) {
                    writeFile << "CHR SNP POS GC GENOAA GENOAB GENOBB GENO_AAA GENO_AAB GENO_ABB GENO_UAA GENO_UAB GENO_UBB IMPUTED? GENE_LIST";
                } else {
                    if(list_of_field[1].substr(0,2)=="rs") {//---Check SNP with rs or not-----
                        match_line_counter++; //---accepted line counter according to SNP name --------
                        GENO.clear();GENO_A.clear();GENO_U.clear();
                        GENO = FindGenoFromString(list_of_field[14]);
                        GENO_A = FindGenoFromString(list_of_field[15]);
                        GENO_U = FindGenoFromString(list_of_field[16]);
                        if(GENO.size()==3 && GENO_A.size()==3 && GENO_U.size()==3) {
                            match_line = list_of_field[0] + " " +list_of_field[1] + " " +list_of_field[2] + " " +list_of_field[3] + " " + GENO[0] + " " + GENO[1]+ " " + GENO[2] + " " + GENO_A[0]+ " " + GENO_A[1]+ " " + GENO_A[2]+ " " + GENO_U[0]+ " " + GENO_U[1]+ " " + GENO_U[2] + " " + list_of_field[17] + " " + list_of_field[20];
                            writeFile << '\n' << match_line;
                        } else {
                            if(complex_geno_line_counter>0) write_complex_geno_file << '\n';
                            complex_geno_line_counter++;
                            write_complex_geno_file << "Original line number : " << line_counter << " Match Line number : " << match_line_counter << "\n";
                            write_complex_geno_file << line;
                        }
                    }
                }
                
            } else {
                if(complex_line_counter>0) write_complex_file << '\n';
                complex_line_counter++;
                write_complex_file << "Original line number : " << line_counter << "\n";
                write_complex_file << line;
            }
        }

        //cout << line_counter << '\n';
    }

    write_report_file << "Total line : " << line_counter << '\n';
    write_report_file << "Total Complex line (field is != 21) : " << complex_line_counter << '\n';
    write_report_file << "Total line considered for SNP is : " << line_counter-complex_line_counter << "\n\n";
    
    write_report_file << "Total line (where rs is not present) : " << (line_counter-complex_line_counter)-match_line_counter << '\n';
    write_report_file << "Total accepted line (where rs is present) : " << match_line_counter << '\n';
    write_report_file << "Total Complex GENO line (where something is wrong with format of GENO, GENO_A or GENO_U) : " << complex_geno_line_counter << '\n';
    write_report_file << "Finally, Validated line in MIGen is : " << match_line_counter-complex_geno_line_counter;
    
    readFile.close();
    writeFile.close();
    write_complex_file.close();
    write_complex_geno_file.close();
    write_report_file.close();
}

vector<string> GCN_MIGen::FindGenoFromString(string geno_string) {
    vector<string> ind_geno_number;
    GCN_Utility GCN_Utility_;
    geno_string=GCN_Utility_.Trim(geno_string);
    char Mychar;
    int ASCII_Mychar;
    string temp_geno_string;
    for(int i=0;i<geno_string.size();i++) {
        Mychar = geno_string[i];
        ASCII_Mychar = static_cast<int>(Mychar);
        if(ASCII_Mychar>=48 && ASCII_Mychar<=57) {
            temp_geno_string = temp_geno_string + Mychar;
        } else {
            if(temp_geno_string!="") ind_geno_number.push_back(temp_geno_string);
            temp_geno_string="";
        }
    }
    if(temp_geno_string!="") ind_geno_number.push_back(temp_geno_string);        

    //for(int i=0;i<ind_geno_number.size();i++) cout << ind_geno_number[i] << "\n";
    return ind_geno_number;
}

void GCN_MIGen::AddConnection_in_new_col() {
    
    GCN_Utility GCN_Utility_;
    ifstream readFile, readConnFile;
    readFile.open("result/50 weeks.txt");
    readConnFile.open("result/connection_week_50.txt");
    if(!readFile || !readConnFile)
    {
        cout << "Input File, is not found ...\n";
        exit(0);
    }

    ofstream writeFile;
    writeFile.open("result/new_50_weeks.txt");
    if(!writeFile)
    {
        cout << "Output File, is not found ...\n";
        exit(0);
    }
    
    vector<string> gene_list,deg_conn;
    vector<string> split_line;
    string line;
    while(getline(readConnFile,line)) {
        line=GCN_Utility_.Trim(line);
        if(line!="") {
            GCN_Utility_.Split(line,'\t',split_line);
            if(split_line.size()==2) {
                gene_list.push_back(split_line[0]);
                deg_conn.push_back(split_line[1]);
            }
        }
        
    }
    
//    for(int i=0;i<gene_list.size();i++) {
//        cout << gene_list[i] << '\t' << deg_conn[i] << '\n';
//    }
    
    readConnFile.close();
    bool header=true;
    while(getline(readFile,line)) {
        line=GCN_Utility_.Trim(line);
        if(line!="") {
            GCN_Utility_.Split(line,'\t',split_line);
            if(split_line.size()==3) {
                if(header) {
                    line=line+'\t'+"Network connections";
                    header=false;
                    writeFile << line;
                } else {
                    int index = GCN_Utility_.Find_Gene(split_line[1],gene_list);
                    if(index==-1) {
                        line=line+'\t'+"";
                    } else {
                        line = line + '\t' + deg_conn[index];
                    }
                    writeFile << '\n' << line;
                }
            }
        }
        
    }
    
    readFile.close();
    writeFile.close();
    
    
}

void GCN_MIGen::writeMIGen_chormosomeBychromosome() {
    
    GCN_Utility GCN_Utility_;
    ofstream writeFile, writeInfo;
    ifstream readFile;
    
    vector<int> chr;
    vector<string> chr_info;
    
    readFile.open("result/varified_MIGen.txt");
    if(!readFile) {
        cout << "Input file not found ...\n";
    }
    
    writeInfo.open("result/chr/MIGen_chr_info.txt");
    if(!writeInfo) {
        cout << "Input file not found ...\n";
    }
    
    string line;
    int chr_num;
    bool header_tag = true;
    string header;
    vector<string> split_line;
    while(getline(readFile,line)) {
        line=GCN_Utility_.Trim(line);
        if(line!="") {
            if(header_tag) {
                header = line;
                header_tag = false;
            } else {
                GCN_Utility_.Split(line,' ',split_line);
                chr_num = GCN_Utility_.string_To_int(split_line[0]);
                chr.push_back(chr_num);
                chr_info.push_back(line);
            }
        }
    }
    readFile.close();
    
    int chr_size = chr.size();
    
    bool new_file_tag = true;
    int chr_control_num;
    vector<int> chr_list;
    vector<int> chr_counter;
    int pos, chr_counter_index;
    vector<int>::iterator it;
    string file_name;
    int total_chr = 0;
    for(int i=0;i<chr_size;i++) {
        if(new_file_tag) {
            stringstream ss(stringstream::in | stringstream::out);
            ss << chr[i];
            file_name = "result/chr/chr" + ss.str() + ".txt";
            writeFile.open(file_name.c_str());
            writeFile << header;
            new_file_tag = false;
            chr_list.push_back(chr[i]);
            chr_counter.push_back(0);
            chr_counter_index=chr_counter.size()-1;
            chr_control_num = chr[i];
        }
        
        if(!new_file_tag) {
            if(chr_control_num==chr[i]) {
                writeFile << '\n' << chr_info[i];
                chr_counter[chr_counter_index]=chr_counter[chr_counter_index]+1;
            } else {
                writeFile.close();
                it = find(chr_list.begin(),chr_list.end(),chr[i]);
                if(*it==chr[i]) {
                    pos = distance(chr_list.begin(),it);
                    stringstream ss(stringstream::in | stringstream::out);
                    ss << chr_list[pos];
                    file_name = "result/chr/chr" + ss.str() + ".txt";
                    writeFile.open(file_name.c_str(),istream::app);
                    writeFile << '\n' << chr_info[i];
                    chr_counter_index = pos;
                    chr_counter[chr_counter_index]=chr_counter[chr_counter_index]+1;
                    chr_control_num = chr[i];
                } else {
                    new_file_tag = true;
                    i--;
                }
                
                //if()
            }
        }
    }
    
    for(int i=0;i<chr_counter.size();i++) {
        writeInfo << "Total in chr" << chr_list[i] << " is :" << chr_counter[i] << '\n';
        total_chr=total_chr+chr_counter[i];
    }
    
    
    writeInfo << "\nIn total chromosome line :" << total_chr;
    writeInfo.close();
    //writeFile.open("result/")
}