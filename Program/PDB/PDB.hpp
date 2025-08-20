#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <string>
#include <vector>
#include <map>

using namespace std;

class PDB
{

  public:
    struct ResUnit
    {
       char chn;
       char ins;
       int seq_ori;
       int seq;
       float Bfactor;
       vector<float> xyz;
       string atm;
       string res;
       string elem;
    } unit;
	
	string mol_name;
	vector<string> chains_name;
    
    vector<ResUnit> pep;
    vector<ResUnit> hetatms;
    vector<vector<ResUnit>> ligands;
    vector<vector<int> > res_classify;
    
    map<int, char> resiSeq;
    map<char, int> ResiNumber;
    string chains_id;
    
    
    int ReadPDB(char *);
    void OutputPDB(string);
    void OutputPDB_SortAtoms(string label);
    void OutputPDB_partChains(string path, string label, string chains);
    void OutputPDBByChain_part(char chain_id, int start, int end);
    int OutputPDBByChain();
    int SplitPep();
    
    int OutputSeqByChain();
    int OutputSeq_partChains(string label, string chains);
    int OutputSeqByChain_part(char chain_id, int start, int end);
    int OutputAASeq(string outf);
    
    int SplitLigands();
    void BindingSiteRes(char *);
    void ResiSeq();
    int CountResi();
    
    void ResClassify();//记录一个氨基酸残基是否是配体结合位点残基，是否是预测结合位点残基；
    void OutputBindingSite();//输出配体结合位点的3D结构；
    //void ResiSeq();//将蛋白pdb文件中的氨基酸编号与其单字母标示相对应,得到单字母表示的氨基酸序列；
    
    void init()
    {
    	unit.Bfactor = 1.00;
    	unit.elem = "";
    }

};
