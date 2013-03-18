/*******************************R-callable function***************************************/
/*decisivatoR
Software realisation of Mareike Fischer's algorithms.
*/

#include <stdio.h>
#include <iostream>
#include <string.h>
#include "comb.h"
#include "setlib.h"
#include <vector>
#include <map>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

RcppExport SEXP IsDecisiveRooted(SEXP taxa, SEXP s, SEXP n, SEXP _k) {
    
    short N = as<int>(n); //Number of taxons
    //std::vector<string> X(taxa.size());   
    Rcpp::CharacterVector cx = Rcpp::CharacterVector(taxa);  
    string *X = new string[N];
    for (int i=0; i<cx.size(); i++) {  
      X[i] = cx[i];  
    } 
    
    //Triplet
    string triplet[3];

    //Actual data set
    Rcpp::List SS = Rcpp::List(s);   
    int k= SS.size(); //Number of taxon subsets
    string **S = NULL;
    try {
        S = new string*[k];
    } catch (bad_alloc& ba) {
        cerr << "bad_alloc caught: " << ba.what() << endl;
    }
    for (int i=0; i<k; i++) {  
        Rcpp::CharacterVector SSS = SS[i];
        S[i] = new string[SSS.size()];
        for (int ii=0; ii<SSS.size(); ii++) {  
          S[i][ii] = SSS[ii];  
        }
    }  
    
    
    //Unrooted tree case
    short m = 3; //We consider triplets.
    
    //Number of possible combinations:
    unsigned long combnum = Res(N,m);
    cout << "Number of combinations is: " << combnum << endl;
    int *C = NULL;
    C = new int[m];
    for(int i=0; i<m; i++)
    {
      C[i] = i;
      triplet[i] = X[C[i]];
      
    }
    
    //Decisiveness. Main loop.
    //First, check the very first triplet:
    bool decisive = false;
    for(int j = 0; j < k; j++)
    {
       if(issubset(triplet, S[j], 3, 4) == true) {decisive = true;}
    }
    if(decisive == false) 
    {
      cout << "The given set is not decisive." << endl;
      return(R_NilValue);
    }
    
    for(int ii=0; ii < combnum-1; ++ii)
    {
      decisive = false;
      GetNext( C, N, m );
      
      //Printing result
      for(int jj=0; jj < m; jj++)
      {
	      triplet[jj] = X[C[jj]];
      }
      
      for(int j = 0; j < k; j++)
      {
        if(issubset(triplet, S[j], 3, 4) == true) {decisive = true;}
      }
      if(decisive == false) 
      {
        cout << "The given set is not decisive." << endl;
        return(R_NilValue);
      }
    }

    free(C);
    cout << "Computation is done. The given data set is decisive.\n";
    
    return(R_NilValue);
}

RcppExport SEXP IsDecisiveUnrooted(SEXP taxa, SEXP s, SEXP n, SEXP _k) {
    
    unsigned int N = as<unsigned int>(n); //Number of taxons
    //std::vector<string> X(taxa.size());   
    Rcpp::CharacterVector cx = Rcpp::CharacterVector(taxa);  
    string *X = new string[N];
    for (int i=0; i<cx.size(); i++) {  
      X[i] = cx[i];  
    } 
    //cout << N << endl;
    
    //Actual data set
    Rcpp::List SS = Rcpp::List(s);   
    int k= SS.size(); //Number of taxon subsets
//cout << k << endl;
    string **S = NULL;
    try {
        S = new string*[k];
    } catch (bad_alloc& ba) {
        cerr << "bad_alloc caught: " << ba.what() << endl;
    }
    for (int i=0; i<k; i++) {  
        Rcpp::CharacterVector SSS = SS[i];
        S[i] = new string[SSS.size()];
        for (int ii=0; ii<SSS.size(); ii++) {  
          S[i][ii] = SSS[ii];  
        }
    }  
    
    
    //Unrooted tree case
    unsigned int m = 4; //We consider quadruples.
    //Quadroople
    string quadruple[m];

    //Number of possible combinations:
    unsigned int combnum = Res(N,m);
    cout << "Number of combinations is: " << combnum << endl;
    
    int **Q = new int*[combnum]; //array of quadrooples
    
    int *C = NULL;
    C = new int[m];
    for(int i=0; i<m; i++)
    {
      C[i] = i;
      quadruple[i] = X[C[i]];
    }
    
    
    vector<int> CQ; //Cross-quadrooples
    //Determine Cross-quadrooples
    for(unsigned int i=1; i<combnum; ++i)
    {
	Q[i-1] = new int[m];
        for(int j=0;j<m;++j) { Q[i-1][j] = C[j];}

	bool cq_flag = true;    	
	for(int j=0; j<k;++j)
    	{
		if(issubset(quadruple, S[j], 4, 4) == true) {cq_flag = false; break;}
    	}

	if(cq_flag == true) {CQ.push_back(i-1);}
    	
	GetNext( C, N, m );
	
	for(int jj=0; jj < m; jj++)
    	{
		quadruple[jj] = X[C[jj]];
    	}
	
    }
	
    //Debug
    //Pringting CQs:
    //for(int i=0; i<CQ.size(); ++i) {
    //	for(int j=0; j<4;++j) {cout << Q[CQ[i]][j] << " ";}
    //	cout << endl;
    //}

    //Main loop
    //Determine is the CQ can be resolved, i.e. find a "fixed" taxon.
    vector<int> RCQ; //Resolved cross-quadruples
    unsigned int flag_size = CQ.size();
    for(unsigned int i=0; i<CQ.size(); ++i) //For each cross-quadruple:
    {
	//Determine if the taxon is not presented in the given CQ.
        std::map< int, int > cq_map;
	cq_map.insert(std::make_pair(Q[CQ[i]][0],0));
	cq_map.insert(std::make_pair(Q[CQ[i]][1],0));
	cq_map.insert(std::make_pair(Q[CQ[i]][2],0));
	cq_map.insert(std::make_pair(Q[CQ[i]][3],0));
	
	//for(int jj=0; jj < m; jj++)
    	//{
	//	cout << X[Q[CQ[i]][jj]];
	//}	
	//cout << endl;
	for(int j=0; j<N; ++j) //For each taxon in the taxon set X
	{
		//If the taxon does not belong to the i-th cross-quadruple:
		if(cq_map.find(j) == cq_map.end()) 
		{
			//cout << j << endl;			
			int fixed_taxon_counter = 0;			
			//Consequently replace all the elements in the given cross-quadruple by that taxon and check if that taxon is a fixed taxon
			for(int kk=0; kk<4;++kk) 
			{
				for(int jj=0; jj < m; jj++)
    				{
					if(jj != kk) 
					{
						quadruple[jj] = X[Q[CQ[i]][jj]];
					} else 
					{
						quadruple[jj] = X[j];
					}
					//cout << quadroople[jj] << " ";
    				} 
				//cout << endl;
				for(int ii=0; ii<k; ++ii) //For all subsets in S
				{
					if(issubset(quadruple, S[ii], 4, 4) == true) {fixed_taxon_counter += 1;break;}
					
				}
				if(fixed_taxon_counter < 4)
				{
					//Check resolved quadruples too:
					for(unsigned int r=0; r<RCQ.size(); ++r)
					{	
						string rq[4] = { X[Q[RCQ[r]][0]], X[Q[RCQ[r]][1]], X[Q[RCQ[r]][2]], X[Q[RCQ[r]][3]]};
						
						if(issubset(quadruple, rq, 4, 4) == true) {fixed_taxon_counter += 1;break;}
					}
				} 
				if(fixed_taxon_counter == 4) {break;}
			}
				
			if(fixed_taxon_counter == 4) 
			{
				//cout << fixed_taxon_counter << endl;
				//Mark quadruple as resolved. Put its id into resolved quadruples:
				RCQ.push_back(CQ[i]);
							
				flag_size--; 
				break;
			}
			//cout << fixed_taxon_counter << endl;
			

		}
	}
	
    }

    
    if(flag_size == 0) 
    {
	cout << "Computation is done. The given data set is decisive.\n";
    }
    else
    {
	cout << "Computation is done. The given data set is not decisive.\n";
    }

    //Debug: resolved cross-quadruples:
    //for(int r=0; r<RCQ.size(); ++r)
    //{	
    //	for(int iii=0; iii<4; iii++)	
    //		cout << X[Q[RCQ[r]][iii]] << " ";
    //		cout << "\n";
    //}
    
    //int set1[] = {1,6,5,4};
    //int set2[] = {1,2,3,4};
    //vector<int> dif = diff(set1,set2,4,4);
    //for(int i=0; i<dif.size(); ++i)
    //{
    //	cout << dif[i] << " ";
    //}
    //cout << endl;

    free(C);
    
    return(R_NilValue);
}
