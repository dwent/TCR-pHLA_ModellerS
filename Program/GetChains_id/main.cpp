#include "../PDB/PDB.cpp"

int main(int agrc, char* argv[])
{
	PDB pro;
	pro.ReadPDB(argv[1]);
	
	pro.CountResi();
	cout<<pro.chains_id;
}
