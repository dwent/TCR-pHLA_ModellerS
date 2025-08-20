#include <string>
#include <iostream>

using namespace std;

void print_matrix    ( int ** const, string, string );
void print_al        ( string&, string& );

void  print_al( string& seq_1_al, string& seq_2_al )
{
    cout << seq_1_al << endl;
    cout << seq_2_al << endl;
}

void  print_matrix( int ** F, string seq_1, string seq_2 )
{
    int  L1 = seq_1.length();
    int  L2 = seq_2.length();

    cout << "        ";
    for( int j = 0; j < L1; j++ )
    {
        cout << seq_1[ j ] << "   ";
    }
    cout << "\n  ";

    for( int i = 0; i <= L2; i++ )
    {
        if( i > 0 )
        {
                cout << seq_2[ i-1 ] << " ";
        }
        for( int j = 0; j <= L1; j++ )
        {
                cout.width( 3 );
                cout << F[ i ][ j ] << " ";
        }
        cout << endl;
    }
}
