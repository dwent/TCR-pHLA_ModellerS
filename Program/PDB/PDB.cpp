#include "PDB.hpp"
#include "basicFunc.hpp"

//g++ -std=c++11 *.cpp -o main -O3
//编译时加上-std=c++11，否者C++中的一些新函数无法使用，如to_string

char AAname3to1(string name3)
{
	char name1;
	if( name3=="ALA" )	name1='A';
	else if( name3=="ARG" ) name1='R';
	else if( name3=="ASN" ) name1='N';
	else if( (name3=="ASP")||(name3=="ASX") ) name1='D';
	else if( name3=="CYS" ) name1='C';
	else if( (name3=="GLN")||(name3=="GLX") ) name1='Q';
	else if( name3=="GLU" ) name1='E';
	else if( name3=="GLY" ) name1='G';
	else if( name3=="HSD" ) name1='H';
	else if( (name3=="HIS")||(name3=="HSE")||(name3=="HSD")||(name3=="HSP") ) name1='H';
	else if( name3=="ILE" ) name1='I';
	else if( name3=="LEU" ) name1='L';
	else if( name3=="LYS" ) name1='K';
	else if( name3=="MET" ) name1='M';
	else if( name3=="PHE" ) name1='F';
	else if( name3=="PRO" ) name1='P';
	else if( name3=="SER" ) name1='S';
	else if( name3=="THR" ) name1='T';
	else if( name3=="TRP" ) name1='W';
	else if( name3=="TYR" ) name1='Y';
	else if( name3=="VAL" ) name1='V';
	else
	{
        name1='*';
		//cout << "\nError in AAname3to1(): no residue name.\t" << name3 << endl;
	}
	
	return name1;
}

/* //该函数在读取PDB文件是将atom和hetatm分开读取，无法用于判断一条链中是否存在hetatm
int PDB::ReadPDB(char *file)
{
    init();
    string buf;
    unit.xyz.assign(3, 0);
    ifstream fd(file, ios::in);
    
	if(!fd)
	{
		cout<<"ERROR in Opening"<<file<<endl;
		return -1;
	}
	
	for(int i=0; file[i] != '\0'; i++)
    {
        if(file[i]=='/')
        	mol_name.clear();
        else if(file[i]=='.' && file[i+1] != '.' && file[i+1] != '/' && file[i-1] != '.')	
        	break;
        else
        	mol_name.push_back(file[i]);
    }
 
    while(fd.good())
    {
    	getline(fd, buf);
        if(buf.length()>53 && buf.substr(0,4)=="ATOM")
        {
        	if(buf[16]!=' ' && buf[16]!='A' && buf[16]!='1') continue;
            unit.res=trim(buf.substr(17,4));//残基的三字母名称。
            //2018-8-11
            if(buf.substr(12,1) == "H"){
            	unit.atm=trim(buf.substr(12,4));
            }
            else{
            	unit.atm=trim(trim(buf.substr(13,3))+buf.substr(12,1));//如果atom名称的第一个字符是数字，则进行该项调整，将数字移至最后。
            }
            //
            if(unit.atm=="H") unit.atm="HN";
            if     (unit.res=="SER"&&unit.atm=="HG") unit.atm="HG1";
            else if(unit.res=="ILE"&&unit.atm=="CD1") unit.atm="CD";
            else if(unit.res=="ILE"&&unit.atm=="HD11") unit.atm="HD1";
            else if(unit.res=="ILE"&&unit.atm=="HD12") unit.atm="HD2";
            else if(unit.res=="ILE"&&unit.atm=="HD13") unit.atm="HD3";
            else if(unit.res=="LEU"&&unit.atm=="CD1") unit.atm="CD2";
            else if(unit.res=="LEU"&&unit.atm=="CD2") unit.atm="CD1";

            else if(unit.res=="LEU"&&unit.atm=="HD11") unit.atm="HD21";
            else if(unit.res=="LEU"&&unit.atm=="HD12") unit.atm="HD22";
            else if(unit.res=="LEU"&&unit.atm=="HD13") unit.atm="HD23";
            else if(unit.res=="LEU"&&unit.atm=="HD21") unit.atm="HD11";
            else if(unit.res=="LEU"&&unit.atm=="HD22") unit.atm="HD12";
            else if(unit.res=="LEU"&&unit.atm=="HD23") unit.atm="HD13";

            else if(unit.res=="CYS"&&unit.atm=="HG") unit.atm="HG1";
            //else if(unit.res=="HIS") unit.res="HSD";

            unit.chn=buf.substr(21,1).c_str()[0];
            unit.seq_ori=atoi(buf.substr(22,4).c_str());//残基的序列编号。
            //对于存在多条链的蛋白，为了区分不同链上相同序号的残基，将残基序号重编码
            string alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
            string s_temp;
            for(int i=0; i<alpha.size(); i++)
            {
            	if(unit.chn == alpha[i] )
            	{
            		stringstream s;
            		s << (i+1);
            		s_temp=s.str();
            		break;
            	}
            }
            unit.seq=atoi((s_temp+trim(buf.substr(22,4))).c_str());	
            
            unit.ins=buf.substr(26,1).c_str()[0];
            unit.xyz[0]=atof(buf.substr(30,8).c_str());
            unit.xyz[1]=atof(buf.substr(38,8).c_str());
            unit.xyz[2]=atof(buf.substr(46,8).c_str());
            //unit.Bfactor=atof(buf.substr(60,6).c_str());
            //unit.elem=trim(buf.substr(66,14));
            
            //该部分用于原子叠加时让甘氨酸alphaC上的一个氢原子用与同其他氨基酸的betaC原子叠加
            if(unit.res == "GLY" && unit.atm=="HA2")
            {
            	unit.atm = "FLA";
            }
            //
            
            pep.push_back(unit);
		}
		else if(buf.length()>53 && buf.substr(0,6)=="HETATM")//注意：一些不常见的氨基酸也被标注为HETATM
        {
            unit.atm=trim(buf.substr(12,4));
            unit.res=trim(buf.substr(16,4));
            unit.chn=buf.substr(21,1).c_str()[0];
            unit.seq=atoi(buf.substr(22,4).c_str());//残基的序列编号。
            unit.ins=buf.substr(26,1).c_str()[0];
            unit.xyz[0]=atof(buf.substr(30,8).c_str());
            unit.xyz[1]=atof(buf.substr(38,8).c_str());
            unit.xyz[2]=atof(buf.substr(46,8).c_str());
            //unit.Bfactor=atof(buf.substr(60,6).c_str());
            //unit.elem=trim(buf.substr(66,14));
            hetatms.push_back(unit);
        }
    }
    fd.close();
    return 1;
}
*/

int PDB::ReadPDB(char *file)
{
    init();
    string buf;
    unit.xyz.assign(3, 0);
    ifstream fd(file, ios::in);
    vector<PDB::ResUnit> pep_temp;
    
	if(!fd)
	{
		cout<<"ERROR in Opening"<<file<<endl;
		return -1;
	}
	
	for(int i=0; file[i] != '\0'; i++)
    {
        if(file[i]=='/')
        	mol_name.clear();
        else if(file[i]=='.' && file[i+1] != '.' && file[i+1] != '/' && file[i-1] != '.')	
        	break;
        else
        	mol_name.push_back(file[i]);
    }
 
    while(fd.good())
    {
    	getline(fd, buf);
        if(buf.length()>53 && (buf.substr(0,4)=="ATOM" || buf.substr(0,6)=="HETATM"))
        {
        	if(buf[16]!=' ' && buf[16]!='A' && buf[16]!='1') continue;
            unit.res=trim(buf.substr(17,4));//残基的三字母名称。
            //2018-8-11
            if(buf.substr(12,1) == "H"){
            	unit.atm=trim(buf.substr(12,4));
            }
            else{
            	unit.atm=trim(trim(buf.substr(13,3))+buf.substr(12,1));//如果atom名称的第一个字符是数字，则进行该项调整，将数字移至最后。
            }
            //
            if(unit.atm=="H") unit.atm="HN";
            if     (unit.res=="SER"&&unit.atm=="HG") unit.atm="HG1";
            else if(unit.res=="ILE"&&unit.atm=="CD1") unit.atm="CD";
            else if(unit.res=="ILE"&&unit.atm=="HD11") unit.atm="HD1";
            else if(unit.res=="ILE"&&unit.atm=="HD12") unit.atm="HD2";
            else if(unit.res=="ILE"&&unit.atm=="HD13") unit.atm="HD3";
            else if(unit.res=="LEU"&&unit.atm=="CD1") unit.atm="CD2";
            else if(unit.res=="LEU"&&unit.atm=="CD2") unit.atm="CD1";

            else if(unit.res=="LEU"&&unit.atm=="HD11") unit.atm="HD21";
            else if(unit.res=="LEU"&&unit.atm=="HD12") unit.atm="HD22";
            else if(unit.res=="LEU"&&unit.atm=="HD13") unit.atm="HD23";
            else if(unit.res=="LEU"&&unit.atm=="HD21") unit.atm="HD11";
            else if(unit.res=="LEU"&&unit.atm=="HD22") unit.atm="HD12";
            else if(unit.res=="LEU"&&unit.atm=="HD23") unit.atm="HD13";

            else if(unit.res=="CYS"&&unit.atm=="HG") unit.atm="HG1";
            //else if(unit.res=="HIS") unit.res="HSD";

            unit.chn=buf.substr(21,1).c_str()[0];
            unit.seq_ori=atoi(buf.substr(22,4).c_str());//残基的序列编号。
            //对于存在多条链的蛋白，为了区分不同链上相同序号的残基，将残基序号重编码
            string alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
            string s_temp;
            for(int i=0; i<alpha.size(); i++)
            {
            	if(unit.chn == alpha[i])
            	{
            		stringstream s;
            		s << (i+1);
            		s_temp=s.str();
            		break;
            	}
            }
            unit.seq=atoi((s_temp+trim(buf.substr(22,4))).c_str());	
            
            unit.ins=buf.substr(26,1).c_str()[0];
            unit.xyz[0]=atof(buf.substr(30,8).c_str());
            unit.xyz[1]=atof(buf.substr(38,8).c_str());
            unit.xyz[2]=atof(buf.substr(46,8).c_str());
            //unit.Bfactor=atof(buf.substr(60,6).c_str());
            //unit.elem=trim(buf.substr(66,14));
            pep_temp.push_back(unit);
		}
		if(buf.substr(0,3) == "TER")
        {
        	pep = pep_temp;//获取最后一个TER之前的所有原子信息
        }
    }
    fd.close();
    return 1;
}


//将蛋白pdb文件中的氨基酸编号与其单字母标示相对应,得到单字母表示的氨基酸序列；
void PDB::ResiSeq()
{
	for(int i=0; i<pep.size(); i++)
	{
		resiSeq[pep[i].seq] = AAname3to1(pep[i].res);
	}
}

//使用babel加Ｈ后得到的文件中Ｈ全部加在最后，且部分原子名称与PDB的一般格式不同，该输出函数旨在将每个残基与其对应的Ｈ原子进行整合，并在读取的时候将Ｈ原子的命名进行矫正。
void PDB::OutputPDB_SortAtoms(string label)
{
	//对读取的所有原子安装残基序号进行排序
	vector<PDB::ResUnit> pep_sorted;
	for(int i=1; i<=pep[pep.size()-1].seq_ori; i++)
	{
		for(int j=0; j<pep.size(); j++)
    	{
    		if(pep[j].seq_ori == i)
    		{
    			pep_sorted.push_back(pep[j]);
    		}
   		}
	}
	pep = pep_sorted;
	string outfile;
    time_t timep;
    time (&timep);  
    if(label == "")  outfile=mol_name+".pdb";
    else  outfile=mol_name+"_"+label+".pdb";   	
    ofstream pdb2out(outfile.c_str(),ios::out);
    pdb2out<<"REMARK The PDB FILE IS CREATED AT "<<ctime(&timep);  
    for(int i=0; i<pep.size(); i++)
    {
    	if(pep[i].res == "GLY" && pep[i].atm=="FLA")
        {
            pep[i].atm = "HA2";
        }
        
        pdb2out<<setiosflags(ios::left)<<setw(6)<<"ATOM"<<resetiosflags(ios::left);
        pdb2out<<setw(5)<<i+1<<" ";
        
        if(pep[i].atm.size()==4)
            pdb2out<<setiosflags(ios::left)<<setw(3)<<pep[i].atm<<resetiosflags(ios::left);
        else
            pdb2out<<" "<<setiosflags(ios::left)<<setw(3)<<pep[i].atm<<resetiosflags(ios::left);
        pdb2out<<" ";
        pdb2out<<pep[i].res;
        pdb2out<<" ";
        pdb2out<<pep[i].chn<<setw(4)<<pep[i].seq_ori;
        pdb2out<<pep[i].ins<<"   ";
        pdb2out<<setiosflags(ios::fixed);
        pdb2out<<setprecision(3);
        pdb2out<<setw(8)<<pep[i].xyz[0]<<setw(8)<<pep[i].xyz[1]<<setw(8)<<pep[i].xyz[2];
        pdb2out<<setprecision(2);
        pdb2out<<"  1.00"<<setw(6)<<pep[i].Bfactor<<setw(12)<<pep[i].elem;
        pdb2out<<endl;

	}
	pdb2out<<"TER \n";
	pdb2out.close();

}

//输出纯蛋白结构,不拆分肽链，不输出Ｈ原子
void PDB::OutputPDB(string label)
{
	string outfile;
    time_t timep;
    time (&timep);
    
    if(label == "")  outfile=mol_name+".pdb";
    else  outfile=mol_name+"_"+label+".pdb";  
    	
    ofstream pdb2out(outfile.c_str(),ios::out);
    pdb2out<<"REMARK The PDB FILE IS CREATED AT "<<ctime(&timep);
    
    for(int i=0; i<pep.size(); i++)
    {
    	if(pep[i].res == "GLY" && pep[i].atm=="FLA")
        {
            pep[i].atm = "HA2";
        }
        
        if(pep[i].atm[0] == 'H')	//不输出Ｈ原子
        {
        	continue;
        }
        
        pdb2out<<setiosflags(ios::left)<<setw(6)<<"ATOM"<<resetiosflags(ios::left);
        pdb2out<<setw(5)<<i+1<<" ";
        
        if(pep[i].atm.size()==4)
            pdb2out<<setiosflags(ios::left)<<setw(3)<<pep[i].atm<<resetiosflags(ios::left);
        else
            pdb2out<<" "<<setiosflags(ios::left)<<setw(3)<<pep[i].atm<<resetiosflags(ios::left);
        pdb2out<<" ";
        pdb2out<<pep[i].res;
        pdb2out<<" ";
        pdb2out<<pep[i].chn<<setw(4)<<pep[i].seq_ori;
        pdb2out<<pep[i].ins<<"   ";
        pdb2out<<setiosflags(ios::fixed);
        pdb2out<<setprecision(3);
        pdb2out<<setw(8)<<pep[i].xyz[0]<<setw(8)<<pep[i].xyz[1]<<setw(8)<<pep[i].xyz[2];
        pdb2out<<setprecision(2);
        pdb2out<<"  1.00"<<setw(6)<<pep[i].Bfactor<<setw(12)<<pep[i].elem;
        pdb2out<<endl;

	}
	pdb2out<<"TER \n";
	pdb2out.close();
}

//输出指定的若干条肽链,不拆分肽链
void PDB::OutputPDB_partChains(string path, string label, string chains)
{
	string outfile;
    time_t timep;
    time (&timep);
     
    if(label == "")  outfile = path + '/' + mol_name +".pdb";
    else  outfile = path + '/' + mol_name+"_"+ label +".pdb";
    	
    ofstream pdb2out(outfile.c_str(),ios::out);
    pdb2out<<"REMARK The PDB FILE IS CREATED AT "<<ctime(&timep);
    
    for(int i=0; i<pep.size(); i++)
    {
    	
    	string::size_type idex;
    	idex = chains.find(pep[i].chn);
    	if(idex == string::npos)
        	//cout << "not found\n";
        	continue;
   		else
   		{
   			if(pep[i].res == "GLY" && pep[i].atm=="FLA")
		    {
		        pep[i].atm = "HA2";
		    }
		    
		    pdb2out<<setiosflags(ios::left)<<setw(6)<<"ATOM"<<resetiosflags(ios::left);
		    pdb2out<<setw(5)<<i+1<<" ";
		    
		    if(pep[i].atm.size()==4)
		        pdb2out<<setiosflags(ios::left)<<setw(3)<<pep[i].atm<<resetiosflags(ios::left);
		    else
		        pdb2out<<" "<<setiosflags(ios::left)<<setw(3)<<pep[i].atm<<resetiosflags(ios::left);
		    pdb2out<<" ";
		    pdb2out<<pep[i].res;
		    pdb2out<<" ";
		    pdb2out<<pep[i].chn<<setw(4)<<pep[i].seq_ori;
		    pdb2out<<pep[i].ins<<"   ";
		    pdb2out<<setiosflags(ios::fixed);
		    pdb2out<<setprecision(3);
		    pdb2out<<setw(8)<<pep[i].xyz[0]<<setw(8)<<pep[i].xyz[1]<<setw(8)<<pep[i].xyz[2];
		    pdb2out<<setprecision(2);
		    pdb2out<<"  1.00"<<setw(6)<<pep[i].Bfactor<<setw(12)<<pep[i].elem;
		    pdb2out<<endl;
   		}
	}
	pdb2out<<"TER \n";
	pdb2out.close();
}


//输出纯蛋白结构,并拆分肽链
int PDB::OutputPDBByChain()
{
	string outfile;
    time_t timep;
    time (&timep);
    
    int chn=0, k=0;
    while(k<pep.size())
    {
    	chn++;
    	outfile=mol_name+"-"+pep[k].chn+".pdb";
    	chains_name.push_back(outfile);
    	
        ofstream pdb2out(outfile.c_str(),ios::out);
        pdb2out<<"REMARK The PDB FILE IS SPLITED BY CHAINS AT "<<ctime(&timep);
		for(int i=k; i<pep.size(); i++)
        {
        	pdb2out<<setiosflags(ios::left)<<setw(6)<<"ATOM"<<resetiosflags(ios::left);
            pdb2out<<setw(5)<<i+1<<" ";
            if(pep[i].atm.size()==4)
            	pdb2out<<setiosflags(ios::left)<<setw(3)<<pep[i].atm<<resetiosflags(ios::left);
            else
                pdb2out<<" "<<setiosflags(ios::left)<<setw(3)<<pep[i].atm<<resetiosflags(ios::left);
            pdb2out<<" ";
            pdb2out<<pep[i].res;
            pdb2out<<" ";
            pdb2out<<pep[i].chn<<setw(4)<<pep[i].seq_ori;
            pdb2out<<pep[i].ins<<"   ";
            pdb2out<<setiosflags(ios::fixed);
            pdb2out<<setprecision(3);
            pdb2out<<setw(8)<<pep[i].xyz[0]<<setw(8)<<pep[i].xyz[1]<<setw(8)<<pep[i].xyz[2];
            pdb2out<<setprecision(2);
            pdb2out<<"  1.00"<<setw(6)<<pep[i].Bfactor<<setw(12)<<pep[i].elem;
            pdb2out<<endl;

            if(i == pep.size()-1 || pep[i].chn != pep[i+1].chn)
            {
            	pdb2out<<"TER \n";
                pdb2out.close();

                outfile="";
                k=i+1;
                break;
            }
		}
	}
    return chn;
}


// 输出某一条肽链的一部分肽段
void PDB::OutputPDBByChain_part(char chain_id, int start, int end)
{
	string outfile;
    time_t timep;
    time (&timep);
    
    ostringstream s_start;
    s_start << start;
    ostringstream s_end;
    s_end << end;
    outfile=mol_name+"-"+chain_id+"_"+s_start.str()+"-"+s_end.str()+".pdb"; 
    	
    ofstream pdb2out(outfile.c_str(),ios::out);
    pdb2out<<"REMARK The PDB FILE IS CREATED AT "<<ctime(&timep);
    
    for(int i=0; i<pep.size(); i++)
    {
    	if(pep[i].chn == chain_id && pep[i].seq_ori >= start && pep[i].seq_ori <= end)
    	{
    		if(pep[i].res == "GLY" && pep[i].atm=="FLA")
		    {
		        pep[i].atm = "HA2";
		    }
		   
		    pdb2out<<setiosflags(ios::left)<<setw(6)<<"ATOM"<<resetiosflags(ios::left);
		    pdb2out<<setw(5)<<i+1<<" ";
		    
		    if(pep[i].atm.size()==4)
		        pdb2out<<setiosflags(ios::left)<<setw(3)<<pep[i].atm<<resetiosflags(ios::left);
		    else
		        pdb2out<<" "<<setiosflags(ios::left)<<setw(3)<<pep[i].atm<<resetiosflags(ios::left);
		    pdb2out<<" ";
		    pdb2out<<pep[i].res;
		    pdb2out<<" ";
		    pdb2out<<pep[i].chn<<setw(4)<<pep[i].seq_ori;
		    pdb2out<<pep[i].ins<<"   ";
		    pdb2out<<setiosflags(ios::fixed);
		    pdb2out<<setprecision(3);
		    pdb2out<<setw(8)<<pep[i].xyz[0]<<setw(8)<<pep[i].xyz[1]<<setw(8)<<pep[i].xyz[2];
		    pdb2out<<setprecision(2);
		    pdb2out<<"  1.00"<<setw(6)<<pep[i].Bfactor<<setw(12)<<pep[i].elem;
		    pdb2out<<endl;
		}
	}
	pdb2out<<"TER \n";
	pdb2out.close();
}

//统计每条链的氨基酸数目
int PDB::CountResi()
{
    int chn=0, k=0;
    int resi_num = 0;
    chains_id = "";
    while(k<pep.size())
    {
          chn++;
          for(int i=k; i<pep.size(); i++)
          {
              if(i==0)
                  resi_num++;
              else if(i !=0 && (pep[i].seq_ori != pep[i-1].seq_ori || pep[i].ins != pep[i-1].ins))
                  resi_num++;
                  
              if(i == pep.size()-1 || pep[i].chn != pep[i+1].chn)
              { 
                  ResiNumber[pep[i].chn] = resi_num;
                  chains_id += pep[i].chn;
                  resi_num = 0;
                  k=i+1;
                  break;
              }
         }
    }
    return chn;
}

//拆分蛋白和多肽
int PDB::SplitPep()
{
	string outfile;
    time_t timep;
    time (&timep);
    
    int chn=0, k=0;
    while(k<pep.size())
    {
    	chn++;
    	if(ResiNumber[pep[k].chn] <= 25)
    	{
    		outfile=mol_name+"_pep.pdb";
    	}
    	else
    	{
    		outfile=mol_name+"_rec.pdb";
    	}
    	
        ofstream pdb2out(outfile.c_str(),ios::app);
        pdb2out<<"REMARK The PDB FILE IS SPLITED BY CHAINS AT "<<ctime(&timep);
		for(int i=k; i<pep.size(); i++)
        {
        	pdb2out<<setiosflags(ios::left)<<setw(6)<<"ATOM"<<resetiosflags(ios::left);
            pdb2out<<setw(5)<<i+1<<" ";
            if(pep[i].atm.size()==4)
            	pdb2out<<setiosflags(ios::left)<<setw(3)<<pep[i].atm<<resetiosflags(ios::left);
            else
                pdb2out<<" "<<setiosflags(ios::left)<<setw(3)<<pep[i].atm<<resetiosflags(ios::left);
            pdb2out<<" ";
            pdb2out<<pep[i].res;
            pdb2out<<" ";
            pdb2out<<pep[i].chn<<setw(4)<<pep[i].seq_ori;
            pdb2out<<pep[i].ins<<"   ";
            pdb2out<<setiosflags(ios::fixed);
            pdb2out<<setprecision(3);
            pdb2out<<setw(8)<<pep[i].xyz[0]<<setw(8)<<pep[i].xyz[1]<<setw(8)<<pep[i].xyz[2];
            pdb2out<<setprecision(2);
            pdb2out<<"  1.00"<<setw(6)<<pep[i].Bfactor<<setw(12)<<pep[i].elem;
            pdb2out<<endl;

            if(i == pep.size()-1 || pep[i].chn != pep[i+1].chn)
            {
            	pdb2out<<"TER \n";
                pdb2out.close();

                outfile="";
                k=i+1;
                break;
            }
		}
	}
    return chn;
}

//按链输出氨基酸序列
int PDB::OutputSeqByChain()
{
    string outfile;
    time_t timep;
    time (&timep);
	
	outfile=mol_name+".fasta";  
	ofstream seq2out(outfile.c_str(),ios::out);
    int chn=0, k=0;
    while(k<pep.size())
    {
          chn++;
          //outfile=mol_name+"-"+pep[k].chn+".fasta";         
          //seq2out<<">"<<mol_name+"-"+pep[k].chn<<"  "<<ctime(&timep);
		  seq2out<<"> "<<pep[k].chn<<"  "<<ctime(&timep);
          for(int i=k; i<pep.size(); i++)
          {
              if(i==0)
                  seq2out<<AAname3to1(pep[i].res);
              else if(i !=0 && (pep[i].seq_ori != pep[i-1].seq_ori || pep[i].ins != pep[i-1].ins))           	
                  seq2out<<AAname3to1(pep[i].res);
                    
              if(i == pep.size()-1 || pep[i].chn != pep[i+1].chn)
              {
                  //seq2out.close();
                  //outfile="";
                  seq2out<<"\n";
                  k=i+1;
                  break;
              }
         }
    }
    /*   
    if(chn != 2)//按肽链数进行处理，根据需要进行修改
    {
        return 0;
    }
    */
    return chn;
}


// 输出某些肽链的氨基酸序列
int PDB::OutputSeq_partChains(string label, string chains)
{
	int res=0;
	
    for(int i=0; i<pep.size(); i++)
    {
        string::size_type idex;
    	idex = chains.find(pep[i].chn);
    	if(idex == string::npos)
        	//cout << "not found\n";
        	continue;
   		else
   		{
        	if(i==0)
		    {
		         cout<<AAname3to1(pep[i].res);
		         res++;
		    }
		    if(i!=0 && (pep[i].seq_ori != pep[i-1].seq_ori || pep[i].ins != pep[i-1].ins))
		    {
		         cout<<AAname3to1(pep[i].res);
		         res++;
		    }
        } 
    }
    return res;
}


// 输出某一条肽链的一部分肽段
int PDB::OutputSeqByChain_part(char chain_id, int start, int end)
{
	int res=0;
	
    for(int i=0; i<pep.size(); i++)
    {
        if(pep[i].chn == chain_id && pep[i].seq_ori >= start && pep[i].seq_ori <= end)
        {
        	if(i==0)
		    {
		         cout<<AAname3to1(pep[i].res);
		         res++;
		    }
		    if(i!=0 && (pep[i].seq_ori != pep[i-1].seq_ori || pep[i].ins != pep[i-1].ins))
		    {
		         cout<<AAname3to1(pep[i].res);
		         res++;
		    }
        } 
    }
    return res;
}


//输出完整的氨基酸序列
int PDB::OutputAASeq(string outf)
{
    int res=0;
    string outfile;
    time_t timep;
    time (&timep);

	if(outf == ""){
		outfile=mol_name+".fasta";
	}else{
		outfile=outf+".fasta";
	}
    
    ofstream seq2out(outfile.c_str(),ios::out);
    seq2out<<">"<<mol_name<<"  "<<ctime(&timep);
    for(int i=0; i<pep.size(); i++)
    {
        if(i==0)
        {
             seq2out<<AAname3to1(pep[i].res);
             res++;
        }
        if(i!=0 && (pep[i].seq_ori != pep[i-1].seq_ori || pep[i].ins != pep[i-1].ins))
        {
             seq2out<<AAname3to1(pep[i].res);
             res++;
             if(res%50 == 0)
             {
                 seq2out<<endl;
             }
        }
    }
    return res;
}

void PDB::ResClassify()
{
	vector<int> res;
	res.assign(3,0);//1.残基ID； 2.是否为结合位点残基（1,0）； 3.是否为预测的结合位点残基（1,0）
	for(int i=0; i<pep.size(); i++)
    {
        if(i==0)
        {
        	res[0]=pep[i].seq;
			res_classify.push_back(res);
			
		}
        if(i !=0 && (pep[i].seq_ori != pep[i-1].seq_ori || pep[i].ins != pep[i-1].ins))
        {
        	res[0]=pep[i].seq;
			res_classify.push_back(res);
		}
	}
}

void PDB::OutputBindingSite()
{

    time_t timep;
    time (&timep);

    ofstream pdb2out("bindingSite.pdb",ios::out);
    pdb2out<<"REMARK The PDB FILE IS THE LIGAND BINDING SITE CREATED AT "<<ctime(&timep);

    for(int i=0; i<res_classify.size(); i++)
    {
        if(res_classify[i][1] == 1)
        {
            for(int j=0; j<pep.size(); j++)
            {
                if(pep[j].seq == res_classify[i][0])
                {
                    pdb2out<<setiosflags(ios::left)<<setw(6)<<"ATOM"<<resetiosflags(ios::left);
                    pdb2out<<setw(5)<<j+1<<" ";
                    if(pep[j].atm.size()==4)
                        pdb2out<<setiosflags(ios::left)<<setw(3)<<pep[j].atm<<resetiosflags(ios::left);
                    else
                        pdb2out<<" "<<setiosflags(ios::left)<<setw(3)<<pep[j].atm<<resetiosflags(ios::left);
                    pdb2out<<" ";
                    pdb2out<<pep[j].res;
                    pdb2out<<" ";
                    pdb2out<<pep[j].chn<<setw(4)<<pep[j].seq_ori;
                    pdb2out<<pep[j].ins<<"   ";
                    pdb2out<<setiosflags(ios::fixed);
                    pdb2out<<setprecision(3);
                    pdb2out<<setw(8)<<pep[j].xyz[0]<<setw(8)<<pep[j].xyz[1]<<setw(8)<<pep[j].xyz[2];
                    pdb2out<<setprecision(2);
                    pdb2out<<"  1.00"<<setw(6)<<pep[j].Bfactor<<setw(12)<<pep[j].elem;
                    pdb2out<<endl;
                }
            }
        }
    }
    pdb2out<<"TER \n";
}

//读取蛋白质中所有配体信息；
int PDB::SplitLigands()
{

	time_t timep;
    time (&timep);
	int ligs=0, k=0;
	vector<ResUnit> tmp;
	
    while(k<hetatms.size())
    {
    	if(hetatms[k].res == "HOH")
		{
			k++;
			continue;
		}
		
        ligs++;
        tmp.clear();
        string lig_name = mol_name.substr(0,4) + "_" + hetatms[k].res + "-" + hetatms[k].chn + "-" + to_string(hetatms[k].seq) + ".pdb";
        ofstream pdb2out(lig_name.c_str(),ios::out);
        pdb2out<<"REMARK The PDB FILE IS SPLITED BY LIGAND AT "<<ctime(&timep);
		for(int i=k; i<hetatms.size(); i++)
		{			
			pdb2out<<setiosflags(ios::left)<<setw(6)<<"HETATM"<<resetiosflags(ios::left);
            pdb2out<<setw(5)<<i+1<<" ";
            pdb2out<<setiosflags(ios::left)<<setw(4)<<hetatms[i].atm<<resetiosflags(ios::left);
            pdb2out<<" ";
            pdb2out<<hetatms[i].res;
            pdb2out<<" ";
            pdb2out<<hetatms[i].chn<<setw(4)<<hetatms[i].seq;
            pdb2out<<hetatms[i].ins<<"   ";
            pdb2out<<setiosflags(ios::fixed);
            pdb2out<<setprecision(3);
            pdb2out<<setw(8)<<hetatms[i].xyz[0]<<setw(8)<<hetatms[i].xyz[1]<<setw(8)<<hetatms[i].xyz[2];
            pdb2out<<setprecision(2);
            pdb2out<<"  1.00"<<setw(6)<<hetatms[i].Bfactor<<setw(12)<<hetatms[i].elem;
            pdb2out<<endl;

			tmp.push_back(hetatms[i]);
			
            if(i == hetatms.size()-1 || (hetatms[i].chn+"-"+to_string(hetatms[i].seq)) != (hetatms[i+1].chn+"-"+to_string(hetatms[i+1].seq)))
            {
            	pdb2out<<"TER \n";
                pdb2out.close();

                k=i+1;
                break;
            }
		}
		ligands.push_back(tmp);
	}
	return ligs; 
}

//获取配体结合位点残基
void PDB::BindingSiteRes(char *file)
{
	map<string, float> AtomRvdw;
	
	AtomRvdw["H"]=1.09;
	AtomRvdw["C"]=1.70;
	AtomRvdw["N"]=1.55;
	AtomRvdw["O"]=1.52;
	AtomRvdw["F"]=1.47;
	AtomRvdw["P"]=1.80;
	AtomRvdw["S"]=1.80;
	AtomRvdw["Cl"]=1.75;
	AtomRvdw["Br"]=1.85;
	AtomRvdw["I"]=1.98;
	AtomRvdw["K"]=2.75;
	AtomRvdw["Na"]=2.27;
	AtomRvdw["Ca"]=2.00;
	AtomRvdw["Fe"]=2.00;
	AtomRvdw["Cu"]=1.40;
    AtomRvdw["Zn"]=1.39;
    AtomRvdw["Mg"]=1.73;
    AtomRvdw["Mn"]=2.00;
    AtomRvdw["Co"]=2.00;
    AtomRvdw["Ni"]=1.63;
    
    string outfile="bindingSiteResi.txt";
	fstream outf(outfile.c_str(), ios::out|ios::app);
    map <int, int> ID_RES;
    for(int j = 0; j<ligands.size(); j++)
    {
    	ID_RES.clear();
		for(int i=0; i<ligands[j].size(); i++)
		{
			if(ligands[j][i].res == "HOH")
			{
				continue;
			}
			
		    if(ligands[j][i].elem == "H")
				continue;
			
		    for(int k=0; k<pep.size(); k++)
		    {
		        if(pep[k].elem == "H")
		            continue;
		            
		        if( (AtomRvdw[pep[k].elem] + AtomRvdw[ligands[j][i].elem] + 0.5) >= Points2Distance(ligands[j][i].xyz, pep[k].xyz))
		        {
		        	//ID_RES[pep[k].seq] = AAname3to1(pep[k].res);
		        	ID_RES[pep[k].seq] = pep[k].seq_ori;
				}
			}	
		}
		
		map<int, int>::iterator it;
		int num=0;
		
		outf<<ligands[j][0].res+"  "+"select "<<file<<", chain "<<ligands[j][0].chn<<" & resi ";	//输出某一条链上的某个配体的结合位点
		for(it=ID_RES.begin(); it!= ID_RES.end(); it++)
		{
			num++;
			outf<<it->second;
			if(num != ID_RES.size())
				outf<<"+";
			//cout<<it->first<<"  "<<it->second<<endl;
		}
		outf<<endl;
	}
	outf.close();
}



