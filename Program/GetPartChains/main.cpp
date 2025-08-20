#include "../PDB/PDB.cpp"

int main(int agrc, char* argv[])
{
	PDB pro;
	pro.ReadPDB(argv[1]);
	
	pro.OutputPDB_partChains(argv[2], argv[3], argv[3]);
	//pro.OutputSeq_partChains(argv[2], argv[2]);
}
