#include "../PDB/PDB.cpp"

int main(int agrc, char* argv[])
{
	PDB pro;
	pro.ReadPDB(argv[1]);
	/*
	receptor.OutputPDB();
	receptor.OutputPDBByChain();
	receptor.OutputSeqByChain();
	receptor.OutputAASeq();
	int num = receptor.SplitLigands();
	cout<<"Ligand Number: "<<num<<endl;
	receptor.BindingSiteRes(argv[1]);
	*/
	pro.OutputSeqByChain();
	//pro.OutputPDB_partChains("part", "ABC");
}
