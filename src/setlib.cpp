#include "setlib.h"

bool issubset(string set1[], string set2[], short m, int n) 
{
	/*checks if a set set1 is a subset of set2*/
	/*Create a dictionary from set2. Then lookup for each element of set1.
	If all elements of set1 are found -> return true, else false.	
	*/
	
	bool flag = true;	

	std::map< string, int > set2_map;
	//std::map<char*,int>::iterator it;

	for(int i=0; i<n; i++)
		set2_map.insert(std::make_pair(set2[i],0));

	//Now look up for each element of set1
	for(int i=0; i<m; i++) 
	{
		//it = set2_map.find(set1[i]);
		if(set2_map.find(set1[i]) == set2_map.end()) {return false;}
	} 
	
	return flag;

}
