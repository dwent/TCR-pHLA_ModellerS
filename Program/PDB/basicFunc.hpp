#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <cstdlib> 
using namespace std;

//---------------------------------------------------------------
//字符串处理函数
//用于去除字符串前后的空格
static string trim(string s)
{
    if (s.empty())
    {
        return s;
    }
    s.erase(0,s.find_first_not_of(" "));
    s.erase(s.find_last_not_of(" ") + 1);
    return s;
} 

//按字符分割字符串
static void split(const string &s, vector<string> &v, const string &c)
{
	string::size_type p1, p2;
	p2 = s.find(c);
	p1=0;
	while(s.npos != p2)
	{
		v.push_back(s.substr(p1, p2-p1));
		p1=p2+c.size();
		p2 = s.find(c, p1);
	}
	
	if(p1 != s.size())
	{
		v.push_back(s.substr(p1));
	}
}

inline string separator(char c, int n)
{
	string line="\n";
	for(int i=0; i<n; i++)
	{
		line+=c;
	}
	line+="\n";
	return line;
}

inline string removeLineBreaks(string str)
{
	string newstr="";
	for(int i=0; str[i] != '\0'; i++)
	{
		if(str[i] != '\n')
		{
			newstr+=str[i];
		}
	}
	return newstr;
}

inline float Points2Distance2(const vector<float> &c1, const vector<float> &c2)
{
   float a=(c1[0]-c2[0]);
   float b=(c1[1]-c2[1]);
   float c=(c1[2]-c2[2]);
   return a*a+b*b+c*c;
}

inline float Points2Distance(const vector<float> &c1, const vector<float> &c2)
{
   float a=(c1[0]-c2[0]);
   float b=(c1[1]-c2[1]);
   float c=(c1[2]-c2[2]);
   return sqrt(a*a+b*b+c*c);
}


