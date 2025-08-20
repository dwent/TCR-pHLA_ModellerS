/*
 * ========================================================================
 * Description:  header file for NWalign
 * Chengxin Zhang 2016-03-13
 * Yang Cao       2017-01-08
 *
 * Backtrace of the origin program starts from lower-right cell, which might
 * cause error. It was revised so that backtrace start from the maximum value
 * of the last row and the last column. The final concensus sequence might
 * not be the matching position of multiple-choice position. That is
 * because this value will be substituted in matrix comparison
 * 2017.1.8.
 * ========================================================================
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
#include "BLOSUM62.hpp"
#include "debug.hpp"

using namespace std;

vector<char> marker;
map<char, int> r2i;

int  ReadFASTAm(const char *file, vector<string>& sname, vector<string>& sequence);

int  max1( int, int, int, int, int, int *);
void dpm_init1( int ** F, int ** traceback, int L1, int L2, int gap , int ext);
void print_traceback2 ( int ** const, string, string );

// standard NW-align. first unaligned position is gap opening. following
// unaligned positions are gap extensions
int Align_NW( string& seq_1, string& seq_2, string& seq_1_al, string& seq_2_al, int gap, int ext, bool prm);
int StandardDynamicPrograming_NW( int **F, int **traceback, const string& seq_1, const string& seq_2, string& seq_1_al, string& seq_2_al, int gap, int ext);

//Read FASTA type file 
int ReadFASTAm(const char *file, vector<string>& sname, vector<string>& sequence)
{
    string buf, res, dat;
    ifstream fd(file, ios::in);
    if(!fd)
    {
        cout<<"Error in Opening "<<file<<endl;
        return 0;
    }
    while(fd.good())
    {
        getline(fd, dat);
        if(dat.empty()) continue;
        
        buf.clear();
        for(int i=0; i<dat.length(); i++)
        {
            if(dat[0]!='>' && dat[i] == ' ')
            {
                continue;
            }
            //else if( dat[0]!='>' && (64< int(dat.c_str()[i]) && int(dat.c_str()[i])<91 ))
            else if( dat[0]!='>' && (int('A')<= int(dat.c_str()[i]) && int(dat.c_str()[i])<=int('Z') ))
            {
                buf.append(dat.substr(i, 1));
            }
            else if( dat[0]=='>')
            {
                buf.append(dat.substr(i, 1));
            }
        }
        if(buf.empty()) continue;
        
        if( buf.find(">") == 0 )
        {
            buf.erase(buf.find(">"), 1);
            sname.push_back(buf);
            if(sname.size()>1)
            {
                sequence.push_back(res);
                res.clear();
            }
        }
        else
        {
            res.append(buf);
        }
    }
    if(sname.size()>0) sequence.push_back(res);
    fd.close();
    return 1;
}



void  print_traceback2( int ** traceback, string seq_1, string seq_2 )
{
    int  L1 = seq_1.length();
    int  L2 = seq_2.length();

    cout << "        ";
    for( int j = 0; j < L1; j++ )
    {
        cout <<setw(3)<< seq_1[ j ] << " ";
    }
    cout << "\n    ";

    for( int i = 0; i <= L2; i++ )
    {
        if( i > 0 )
        {
            cout <<setw(3)<< seq_2[ i-1 ]<< " ";
        }
        for( int j = 0; j <= L1; j++ )
        {
            cout<<setw(3)<< traceback[ i ][ j ]<< " ";
        }
        cout << endl;
    }
}

//Matrix Initialization
void  dpm_init1( int ** F, int ** traceback, int L1, int L2, int gap , int ext)
{
    F[ 0 ][ 0 ] =  0 ;
    traceback[ 0 ][ 0 ] = 0 ;

    int i=0, j=0;

    F[ 0 ][ 1 ] = gap;
    traceback[ 0 ][ 1 ] =  -1 ;
    for( j = 2; j <= L1; j++ )
    {
        F[ 0 ][ j ] = F[ 0 ][ j-1 ] + ext ;
        traceback[ 0 ][ j ] =  -1 ;
    }
    F[ 1 ][ 0 ] = gap;
    traceback[ 1 ][ 0 ] =  1 ;
    for( i = 2; i <= L2; i++ )
    {
        F[ i ][ 0 ] = F[ i-1 ][ 0 ] + ext ;
        traceback[ i ][ 0 ] =  1 ;
    }
}


//Score Comparsion
int  max1( int f1, int f2, int f3, int nU, int nL, int * ptr )
{
    int  max = 0 ;
    if( f2 > f3 && f2>f1 )              
    {
        max = f2 ;
        *ptr = 0 ;
    }
    else if( f1 >= f2 && f1 >= f3 )  
    {
        max = f1 ;
        *ptr = nU ;
    }
    else if( f3 >= f2 && f3 >= f1 )  
    {
        max = f3 ;
        *ptr = -nL ;
    }
    else  
    {
        cout<<"Error "<<f1<<" "<<f2<<" "<<f3<<endl;
    }
    return  max ;   
}


int Align_NW( string& seq_1, string& seq_2, string& seq_1_al, string& seq_2_al, int gap, int ext, bool prm )
{
    int score=0 ;

    int  L1 = seq_1.length();
    int  L2 = seq_2.length();

    int ** F = new int * [ L2+1 ];
    for( int i = 0; i <= L2; i++ )  F[ i ] = new int [ L1+1 ];

    int ** traceback = new int * [ L2+1 ];
    for( int i = 0; i <= L2; i++ )  traceback[ i ] = new int [ L1+1 ];

    dpm_init1( F, traceback, L1, L2, -gap, -ext );

    score=StandardDynamicPrograming_NW( F, traceback, seq_1, seq_2, seq_1_al, seq_2_al, gap, ext);


    if( prm )
    {
        cout << "\nDynamic programming matrix: " << "\n\n";
        print_matrix( F, seq_1, seq_2 );

        cout << "\nTraceback matrix: " << "\n\n";
        print_traceback2( traceback, seq_1, seq_2 );

        cout << endl;
    }

    for( int i = 0; i <= L2; i++ )  delete F[ i ];
    delete[] F;
    for( int i = 0; i <= L2; i++ )  delete traceback[ i ];
    delete[] traceback;

    return  score ;
}


//NW algorithm
int StandardDynamicPrograming_NW( int **F, int **traceback, const string& seq_1, const string& seq_2, string& seq_1_al, string& seq_2_al, int gap, int ext)
{
    int  k = 0, x = 0, y = 0;
    int  ptr, fU, fD, fL, UU, LL ;
    char buf1, buf2 ;
    int  i = 0, j = 0, c=0, m=0, n=0, score=0;
    
    int  L1 = seq_1.length();
    int  L2 = seq_2.length();
    int  ss, tt, gg, pl, mU, mL, nU, nL;

    for( i = 1; i <= L2; i++ )
    {
        for( j = 1; j <= L1; j++ )
        {
    	     buf1=toupper(seq_1.c_str()[j-1]);
    	     buf2=toupper(seq_2.c_str()[i-1]);
    	     if(r2i.count(buf1)>0 && r2i.count(buf2)>0)
    	     {
    	         m=r2i[buf1];
    	         n=r2i[buf2];
    	         ss=scmx[m][n];
    	     }
    	     else
    	     {
    	         ss=0;
    	     }
    	     
             fU = F[ i-1 ][ j ] - gap ;
             fD = F[ i-1 ][ j-1 ] + ss ;
             fL = F[ i ][ j-1 ] - gap ;
             
             mU = fU; nU = 1;
             for(k=0; k<i; k++)
             {
                 UU = F[ k ][ j ] - ext*(i-1-k) - gap;
                 if(mU < UU)
                 {
                     mU = UU;
                     nU = i-k;
                 }
             }
             
             mL = fL; nL = 1;
             for(k=0; k<j; k++)
             {
                 LL = F[ i ][ k ] - ext*(j-1-k) - gap;
                 if(mL < LL)
                 {
                     mL = LL;
                     nL = j-k;
                 }
             }
             
             F[ i ][ j ] = max1( mU, fD, mL, nU, nL, &ptr ) ;
             traceback[ i ][ j ]  = ptr ;
        }
    }
    i-- ; j-- ; k=0;
    
    score=F[i][j];
    ss=F[i][j]; m=i;
    tt=F[i][j]; n=j;
    
    while( i > 0 || j > 0 )
    {
        if(traceback[ i ][ j ] == 0)
        {
            seq_1_al += seq_1[ j-1 ] ; 
            seq_2_al += seq_2[ i-1 ] ; 
            i-- ;  j-- ;
        }
        else if(traceback[ i ][ j ] > 0)
        {
            fU=traceback[ i ][ j ];
            for(k=0; k<fU; k++)
            {
                seq_1_al += '-' ; 
                seq_2_al += seq_2[ i-1 ] ; 
                i-- ;
            }
        }
        else if(traceback[ i ][ j ] < 0)
        {
            fL=-traceback[ i ][ j ];
            for(k=0; k<fL; k++)
            {
                seq_1_al += seq_1[ j-1 ] ; 
                seq_2_al += '-' ; 
                j-- ;
            }
        }
        else
        {
            cout<<"E\n";
        }
    }
    
    reverse( seq_1_al.begin(), seq_1_al.end() );
    reverse( seq_2_al.begin(), seq_2_al.end() );

    return  score ;
}
