#include "Hungarian.hpp"

#define Label1 1e4
#define Label2 2e4
using namespace std;
//输入耗费矩阵；
void ASSIGN::input_matrix(int row_n, int col_n, ofstream &log)
{
	matrix.resize(row_n,vector<coef>(col_n));
	int i,j;
	cout<<"press ENTER to input next coef!"<<endl;
	for (i=0;i<row_n;i++)
	{
		cout<<"input the "<<i+1<<"th row:"<<endl;
		for(j=0;j<col_n;j++)
		{
			cin>>matrix[i][j];
		}
	}
	
	log<<"the matrix inputted is:"<<endl;
	for(i=0;i<row_n;i++)
	{
		for(j=0;j<col_n;j++)
			log<<setw(5)<<matrix[i][j]; //注意setw()是在输出对象的前面加空格；
		log<<endl;
	}
		
}

//读取矩阵文件
int ASSIGN::read_matrix(char *file, ofstream &log)
{
	
	string buf;
	ifstream fd(file, ios::in);
    if(!fd)
    {
		cout<<"Error in Opening "<<file<<endl;
        return -1;
    }

	while(fd.good())
	{
		buf.clear();
		getline(fd, buf);
		if(buf.length()<=1)
		{
			continue;
		}
		istringstream stream(buf);
		float temp;
		vector<float> row;

		while(stream>>temp)
	    {  
	        //cout<<temp<<endl;
	        row.push_back(temp);
	    }  
	    
	    matrix.push_back(row);
				
	}
	

}

//修改矩阵，使矩阵的每一行，每一列都出现0元素；
void ASSIGN::modify_matrix(ofstream &log)
{
	initialize();
	int i,j;
	coef row_min[row_n],column_min[col_n];
	
	matrix_backup=matrix;
	
	for (i=0;i<row_n;i++)
	{
		row_min[i]=matrix[i][0];
		for(j=0;j<col_n;j++)
			if(row_min[i]>matrix[i][j])
				row_min[i]=matrix[i][j];
		for(j=0;j<col_n;j++)
			matrix[i][j]-=row_min[i];
	}

	for (j=0;j<col_n;j++)
	{
		column_min[j]=matrix[0][j];
		for(i=0;i<row_n;i++)
			if(column_min[j]>matrix[i][j])
				column_min[j]=matrix[i][j];
		if(column_min[j]==0) continue;
		for(i=0;i<row_n;i++)
			matrix[i][j]-=column_min[j];
	}
	
	/*
	log<<"the matrix modified is:"<<endl;
	
	for(i=0;i<row_n;i++)
	{
		for(j=0;j<col_n;j++)
			log<<setw(5)<<matrix[i][j]; 
		log<<endl;
	}
	*/
	
}

void ASSIGN::label_zero(ofstream &log)
{
	
	int i, j, n, n_max, min_r_id, min_c_id, min_r, min_c;
	float min=1000.0;
	vector<int> count_r_0,count_c_0,label_r,label_c;

	n = (row_n < col_n) ?  row_n : col_n;
	n_max = (row_n > col_n) ?  row_n : col_n;

	while(indep_0<n)
	{                             
		label_r.assign(row_n,0);
		label_c.assign(col_n,0);
		do
		{
			count_r_0.assign(row_n,0);
			count_c_0.assign(col_n,0);
			for(i=0;i<row_n;i++) 
			{            
				for(j=0;j<col_n;j++)
					if(abs(matrix[i][j])<1e-6)
						count_r_0[i]++;      	 //统计每一行的0元素；	
			}
			for(j=0;j<col_n;j++)
			{
				for(i=0;i<row_n;i++)
					if(abs(matrix[i][j])<1e-6) 
						count_c_0[j]++;			//统计每一列的0元素；	
			}

			min_r=row_n+1;
			min_r_id=row_n;
			for(i=0;i<row_n;i++)
				if(min_r>count_r_0[i] && count_r_0[i]!=0)
				{
					min_r_id=i;					//找到0元素最少的一行；    
					min_r=count_r_0[i];
				} 
     		min_c=col_n+1;
			min_c_id=col_n;
			for(j=0;j<col_n;j++)
				if(min_c>count_c_0[j] && count_c_0[j]!=0)
				{
					min_c_id=j;                 //找到0元素最少的一列；       
					min_c=count_c_0[j];
				}  

			if(min_r_id==row_n && min_c_id==col_n) break;	//此时已不存在0元素，不需要再往后进行；
			
			if(min_r<=min_c) //min_r<min_c ; row_n < col_n
			{					
				for(j=0;j<col_n;j++)
				{
						if(abs(matrix[min_r_id][j]) <1e-6)
						{
							indep_0++;
							matrix[min_r_id][j]=Label1;	//独立0元素标记为Label1；
							label_c[j]=1;			    //第j列用直线标记；
							for(i=0;i<row_n;i++)
								if(abs(matrix[i][j]) <1e-6)
									matrix[i][j]=Label2;//独立0元素所在列的其余0元素标记为Label2；
							break;
						}
				}
				
				for(j=0;j<col_n;j++)
					if(abs(matrix[min_r_id][j]) <1e-6)
						matrix[min_r_id][j]=Label2;		//独立0元素所在行的其余0元素标记为Label2；
			}
			else
			{
								
				for(i=0;i<row_n;i++)
				{
						if(abs(matrix[i][min_c_id]) <1e-6)
						{
							indep_0++;
							matrix[i][min_c_id]=Label1;
							label_r[i]=1;			//第i行用直线标记；
							for(j=0;j<col_n;j++)
								if(abs(matrix[i][j]) <1e-6)
									matrix[i][j]=Label2;
							break;
						}
				}
				
				for(i=0;i<row_n;i++)
					if(abs(matrix[i][min_c_id]) <1e-6)
						matrix[i][min_c_id]=Label2;
			}		
		}while(indep_0<n);
		
		/*
		log<<"the matrix marked is:"<<endl;
		
		for(i=0;i<row_n;i++)
		{
			for(j=0;j<col_n;j++)
				log<<fixed<<setprecision(3)<<setw(8)<<matrix[i][j]; 
			log<<endl;
		}
		*/

		if(indep_0==n) 
		{
			//cout<<"Complete!"<<endl;
			for(i=0;i<row_n;i++)
				for(j=0;j<col_n;j++)
				{
					if(matrix[i][j] == Label1 || matrix[i][j] == Label2) 
						matrix[i][j]=0.0;
					else matrix[i][j]=n_max+1;   //--------------------------------------------------
				}
				
			/*
			log<<"the optimal assignment is:"<<endl;
			for(i=0;i<row_n;i++)
			{
				for(j=0;j<col_n;j++)
					log<<setw(5)<<matrix[i][j]; 
				log<<endl;
			}
			*/			
			break;
		}
		else
		{
			for(i=0;i<row_n;i++)
				for(j=0;j<col_n;j++)
					if(matrix[i][j]==Label1 || matrix[i][j]==Label2) 
						matrix[i][j]=0.0;
		
			for(i=0;i<row_n;i++)
				for(j=0;j<col_n;j++)
					if(label_r[i]==0 && label_c[j]==0)
						if(matrix[i][j]!=0.0 && matrix[i][j]<=min)  min=matrix[i][j];
	
			for(i=0;i<row_n;i++)
				for(j=0;j<col_n;j++)
				{
					if(label_r[i]==0 && label_c[j]==0)
						matrix[i][j]-=min;
					if(label_r[i]==1 && label_c[j]==1)
						matrix[i][j]+=min;
				}	
		}
		
		/*
		log<<"the matrix modified_again is:"<<endl;
		
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
				log<<setw(5)<<matrix[i][j]; 
			log<<endl;
		}
		*/
		indep_0=0;
		min=1000.0; //可根据需要对min附初始值；
	}
}
//matrix的初始元素的值只有0和n_max+1（所有非零元素全部置为n_max+1）;var初始值为1；
void ASSIGN::Match(int var, ofstream &log)
{
	int i,j,k;
	static int t=0;
	coef sum=0;
	int n = (row_n < col_n) ?  row_n : col_n;
	int n_max = (row_n > col_n) ?  row_n : col_n;
	if(var > n)
	{	
		
		//int flag;
		/*
		for(i=0;i<row_n;i++)
		{
			for(k=0;k<col_n;k++)
				if(matrix[i][k]==0)		
				{
					flag=1;
					matrix[i][k]=var;		//若最后一行存在0元素，说明找到一组最优解，下面储存该最优解即可
					break;
				}

					
			if(flag)
				break;
		}
		
		
		for(j=0;j<col_n;j++)
		{
			flag = 1;
			if(matrix[var-1][j]!=0)
			{
				flag = 0;
				continue;			
			}
			matrix[var-1][j]=var;		
			for(k=0;k<col_n;k++)			
				if(matrix[var-1][k] == 0)
					matrix[var-1][k]=-var;
			for(k=0;k<row_n;k++)
				if(matrix[k][j] == 0)
					matrix[k][j]=-var;
			break;
		}
		*/
		

		//if(flag)
		//{	
			t++;	
			match.clear();
			//log<<"##############################"<<endl
				//<<"the "<<t<<"th "<<"result"<<endl;
		
			for(i=0;i<row_n;i++)
			{
				for(j=0;j<col_n;j++)
				{
					if(matrix[i][j]>0 && matrix[i][j]<n_max+1)
					{
						sum+=matrix_backup[i][j];
						coord.x=i;
						coord.y=j;
						match.push_back(coord);
						//log<<setw(3)<<0;
					}
					else 
					{
						//log<<setw(3)<<1;
					}
				}
				//log<<endl;
			 }
			 matchs.push_back(match);
			 sums.push_back(sum);
			 //log<<"SUM = "<<sum<<endl;
			 //log<<"##############################"<<endl;
		//}
		/*
		for(i=0;i<row_n;i++)		
			for(k=0;k<col_n;k++)				
				if(matrix[i][k]==var)								
					matrix[i][k]=0;
		
		var++;
		*/
	}				
	else if(var<=n)
	{	
		for(j=0;j<(row_n<col_n?col_n:row_n);j++)
		{
			if(matrix[row_n<col_n?var-1:j][row_n<col_n?j:var-1]!=0)
			{
				continue;				//若某一行不存在0元素，说明之前的标0操作不可能得到一组最优解，所以后续迭代不再进行，返回上一级继续
			}
			
			matrix[row_n<col_n?var-1:j][row_n<col_n?j:var-1]=var;		//依次标记每一行的每一个0元素为 var，一次只标记一个
			for(k=0;k<(row_n<col_n?col_n:row_n);k++)			//该0元素所在行、列的其他0元素全部标记为 -var
				if(matrix[row_n<col_n?var-1:k][row_n<col_n?k:var-1] == 0)
					matrix[row_n<col_n?var-1:k][row_n<col_n?k:var-1] = -var;
			for(k=0;k<(row_n<col_n?row_n:col_n);k++)
				if(matrix[row_n<col_n?k:j][row_n<col_n?j:k] == 0)
					matrix[row_n<col_n?k:j][row_n<col_n?j:k] = -var;
				
			Match(var+1,log);
			

			for(i=0;i<row_n;i++)
				for(k=0;k<col_n;k++)
					if(abs(matrix[i][k])==var)
						matrix[i][k]=0;
		}		
	}
}

void ASSIGN::outputResults(ofstream &log)
{
	int i,j;
	for(i=0; i<sums.size()-1; i++)
	{	
		for(j=0; j<sums.size()-1-i; j++)
		{
			if(sums[j]> sums[j+1])
			{
				swap(sums[j], sums[j+1]);
				swap(matchs[j], matchs[j+1]);
			}
		}
	}
	for(i=0; i<sums.size(); i++)
	{	
		log<<"SUM = "<<sums[i]<<endl;
		for(j=0; j<matchs[i].size(); j++)
		{
			log<<matchs[i][j].x<<"_"<<matchs[i][j].y<<endl;
		}
	}
}

