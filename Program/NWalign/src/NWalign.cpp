/*
 * =============================================================
 * Needleman Wunsch alignment in C++
 * 2017 Chengxin Zhang forked from Yang Cao
 * =============================================================
 */

#include <iostream>
#include <sstream> 
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstdlib> 
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include "NWalign.hpp"

using namespace std;

int SeqIdentity(const string &seq1, const string &seq2)
{
    if(seq1.length()!=seq2.length())
    {
        cout<<"Error UnAligned\n";
        return 0.0;
    }
    int i=0, j=0, m=0, n=0, k=seq1.length();
    for(i=0; i<k; i++)
    {
        if(seq1[i] == seq2[i] && seq1[i]!='-')
        {
              j++;
        }
        if(seq1[i]!='-') m++;
        if(seq2[i]!='-') n++;
    }
    
    if(m==0) m=1;
    if(n==0) n=1;
    
    //cout<<"Identical Number = "<<j<<endl;
    
    if(m<n)
        return int(1000*float(j)/(float(m)));
    else
        return int(1000*float(j)/(float(n)));
}

int main(int argc, char **argv)
{
    if(argc<3)
    {
        cout<<"Error input!\n";
        cout<<"Usage: NWalign [Query FASTA File 1] [Query FASTA File 2]\n";
        return 0;
    }
    
    int i, j, k, score, m, n, ss, score2;
    char buf, buf1, buf2 ;
    bool Fdebug=false;
    if(argc>3 && argv[3][0]=='d') Fdebug=true;
    
    
    string cseq1, cseq2;
    vector<string> sname, qname;//Sequence Name
    vector<string> ssequence, qsequence;//Sequence data
    ReadFASTAm(argv[1], qname, qsequence);
    ReadFASTAm(argv[2], sname, ssequence);
    
    
    string stag, sref, s_tag, s_ref;

    for (int a=0;a<aa_list.length();a++) r2i[aa_list[a]]=a;
    
    for(i=0; i<qsequence.size(); i++)
    {
        for(k=0; k<ssequence.size(); k++)
        {
            stag=qsequence[i];
            sref=ssequence[k];
            s_tag.clear(); 
            s_ref.clear(); 
          
            score=0; 
            cseq1.clear();
            for(j=0; j<stag.length(); j++)
            {
                  if(stag.c_str()[j] !='-' )
                  {
                        cseq1+=stag[j];
                  }
            }
            cseq2.clear();
            for(j=0; j<sref.length(); j++)
            {
                  if(sref.c_str()[j] !='-' )
                  {
                        cseq2+=sref[j];
                  }
            }
            score=Align_NW(cseq1, cseq2, s_tag, s_ref, 11, 1, Fdebug);
            score2=SeqIdentity(s_tag, s_ref);
            //cout<<">"<<qname[i]<<"\n";
            //cout<<">"<<sname[k]<<"\n";
            //cout<<score2/10.0;
            cout<<score;
            //cout<<"Score "<<score<<" "<<"#identity(%) = "<<score2/10.0<<endl;
            //cout<<">"<<qname[i]<<"\n"<<s_tag<<endl;
            //cout<<">"<<sname[k]<<"\n"<<s_ref<<endl;
            //cout<<endl;
        }
    }  
    return 0;
}
