#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
using namespace std;

typedef float coef;

class ASSIGN
{
  public:
  	vector<vector<coef> > matrix;
  	vector<vector<coef> > matrix_backup;
  	
  	struct COORD
  	{
  		int x;
  		int y;
  	} coord;
  	
  	vector<COORD> match;
  	vector<vector<COORD> > matchs;
    vector<float> sums;
  	
  	int row_n;
	  int col_n;
  	int indep_0;
  	
  
  	//function
  	void initialize()
  	{
  		match.clear();
  		matchs.clear();
  		
  		row_n=matrix.size();
		  col_n=matrix[0].size();
  		indep_0=0;
  		
  	};
    int read_matrix(char *file, ofstream &log);
  	void input_matrix(int, int, ofstream &);
  	void modify_matrix(ofstream &);
  	void label_zero(ofstream &);
  	void Match(int, ofstream &);
    void outputResults(ofstream &);

};
