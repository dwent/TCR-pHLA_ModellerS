#include "Hungarian.hpp"
/*
int main(int argc, char* argv[])
{

	if(argc<2)
    {
    	cout<<"Usage: ./Hungary [out file]\n";
        return 0;
    }
    int m, n;
    cout<<"Please enter the number of the matrix's row:"<<endl;
    cin>>m;

	cout<<"Please enter the number of the matrix's column:"<<endl;
    cin>>n;

	ofstream out_f(argv[1],ios::out);
	
	ASSIGN task;
	task.input_matrix(m, n, out_f);
	task.modify_matrix(out_f);
	task.label_zero(out_f);
	task.Match(1, out_f);
	out_f.close();
	return 1;
}
*/

int main(int argc, char* argv[])
{

	if(argc<3)
    {
    	cout<<"Usage: ./Hungary [intput file] [out file]\n";
        return 0;
    }
    

	ofstream out_f(argv[2],ios::out);
	
	ASSIGN task;
	task.read_matrix(argv[1], out_f);
	
	task.modify_matrix(out_f);
	task.label_zero(out_f);
	task.Match(1, out_f);
	task.outputResults(out_f);
	out_f.close();
	
	return 1;
}
