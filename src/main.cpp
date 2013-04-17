/*******************************R-callable function***************************************/
/*decisivatoR
C++ realization of Mareike Fischer's algorithms.
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

/*============GLOBAL VARIABLES===========*/


template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
    return std::pair<B,A>(p.second, p.first);
}

template<typename A, typename B>
std::map<B,A> flip_map(const std::map<A,B> &src)
{
    std::map<B,A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()), 
                   flip_pair<A,B>);
    return dst;
}

bool GetDecisivenessUnrooted(unsigned int N, string *X, string **S, int *sizes, int **Q, int k,vector<int> &CQ, vector<int> &RCQ)
{
    //Main loop: we have to "resolve the cross-quadruples"
    //Determine is the CQ can be resolved, i.e. find a "fixed" taxon.
    unsigned int m = 4; //We consider quadruples.
    unsigned int flag_size = CQ.size();
    string quadruple[m];
    for(unsigned int i=0; i<CQ.size(); ++i) //For each cross-quadruple:
    {
	      //Determine if the taxon is not presented in the given CQ (missing taxon).
        std::map< int, int > cq_map;
	      cq_map.insert(std::make_pair(Q[CQ[i]][0],0));
	      cq_map.insert(std::make_pair(Q[CQ[i]][1],0));
	      cq_map.insert(std::make_pair(Q[CQ[i]][2],0));
	      cq_map.insert(std::make_pair(Q[CQ[i]][3],0));
	
	    	for(int j=0; j<N; ++j) //For each taxon in the taxon set X
	      {
		      //If the taxon does not belong to the i-th cross-quadruple:
		      if(cq_map.find(j) == cq_map.end()) 
		      {
			      int fixed_taxon_counter = 0;			
			      //Consequently replace all the elements in the given cross-quadruple by that taxon and check if that taxon is a "fixed" taxon
			      for(int kk=0; kk<4;++kk) 
			      {
				      for(int jj=0; jj < m; jj++)
    				  {
					      if(jj != kk)
					      {
						      quadruple[jj] = X[Q[CQ[i]][jj]];
					      } 
                else 
					      {
						      quadruple[jj] = X[j];
					      }
					    } 
				      
				      for(int ii=0; ii<k; ++ii) //For all subsets in S
				      {
					      if(sizes[ii] >= 4) 
                {
                  if(issubset(quadruple, S[ii], 4, sizes[ii]) == true) {fixed_taxon_counter += 1;break;}
                }
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
				      //Mark quadruple as resolved. Put its id into resolved quadruples:
				      RCQ.push_back(CQ[i]);
							flag_size--; 
				      break;
			      }
			    }
	      }
        
        cq_map.clear();
	
    }

    bool res = false;
    if(flag_size == 0) 
    {
	    res = true;//"Computation is done. The given data set is decisive.";
    }
        
    return res;
}


RcppExport SEXP IsDecisiveRooted(SEXP taxa, SEXP s, SEXP n, SEXP _k, SEXP fflag) {
    vector<string> CT; //This vector structure holds Cross-Triplets
    string res = "";
    std::map< int, int > ct_map_freq; //This holds cross-triples and their "weights": <T_id , weight>
    int **T = NULL;
    short N = as<int>(n); //Number of taxons
    
    Rcpp::CharacterVector cx = Rcpp::CharacterVector(taxa);  
    string *X = new string[N];
    for (int i=0; i<cx.size(); i++) 
    {  
      X[i] = cx[i];  
    } 
    
    //Triplet
    string triplet[3];

    //Actual data set
    Rcpp::List SS = Rcpp::List(s);   
    int k= SS.size(); //Number of taxon subsets
    string **S = NULL;
    int *sizes = NULL; //a set of sizes of each of taxon subset
    
    try 
    {
        S = new string*[k];
        sizes = new int[k];
    } 
    catch (bad_alloc& ba) 
    {
        cerr << "bad_alloc caught: " << ba.what() << endl;
    }
    
    for (int i=0; i<k; i++) {  
        Rcpp::CharacterVector SSS = SS[i];
        sizes[i] = SSS.size();
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
    T = new int*[combnum]; //Array of triplets. This array keeps all quadruples created from X
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
      if(sizes[j] >= 3) 
      {
         if(issubset(triplet, S[j], 3, sizes[j]) == true) {decisive = true;}
      }
    }
    if(decisive == false) 
    {
      cout << "The given set is not decisive." << endl;
      return(R_NilValue);
    }
    
    for(int ii=0; ii < combnum; ++ii)
    {
      int freq_counter = 0;
      decisive = false;
      if (ii>0) GetNext( C, N, m );
      
      T[ii] = new int[m]; //New triple
      for(int j=0;j<m;++j) { T[ii][j] = C[j]; } //Fill this new quadruple with taxons ids
      
      
      
      for(int j = 0; j < k; j++)
      {
        if(sizes[j] >= 3) 
        {
          if(issubset(triplet, S[j], 3, sizes[j]) == true) {decisive = true;}
          
          freq_counter += diff2(triplet, S[j], 3, sizes[j]).size();
        }
      }
      if(decisive == false) 
      {
        res = "Computation is done. The given set is not decisive.";
        ct_map_freq.insert(std::make_pair(ii,freq_counter));
        //break;
      }
      
      //Create a triplet
      for(int jj=0; jj < m; jj++)
      {
        triplet[jj] = X[C[jj]];
      }
    }

    if(decisive && (ct_map_freq.size() == 0)) 
    {
      res = "Computation is done. The given data set is decisive.";
    }
    else
    {
      if(as<bool>(fflag))
      {
        cout << "The given data set is not decisive. Trying to fix it..." << endl;
    
        //Fixing the data set
        //Sorting cq_map_freq by value:
        std::map<int, int> dst = flip_map(ct_map_freq);
        std::map<int, int>::reverse_iterator iter;
        int iteration = 1;
        for (iter = dst.rbegin(); iter != dst.rend(); ++iter) 
        {
          //Get the random number from 1..k:
          int r_subset = rand() % k;
          sizes[r_subset] = sizes[r_subset] + 4; //sizes[i] = SSS.size();
          string *tmp_str = new string[sizes[r_subset]];
          for(unsigned int mm=0; mm<sizes[r_subset]-4; ++mm)
          {
            tmp_str[mm] = S[r_subset][mm];
          }
          cout << "\033[A\033[2K";
          cout << "Iteration " << iteration << "..." << endl;
          for(int iii=0; iii<4; ++iii)
          {
            tmp_str[sizes[r_subset]-4+iii] = X[T[iter->second][iii]];
          }
          S[r_subset] = new string[sizes[r_subset]];
          for(int iii=0; iii<sizes[r_subset]; ++iii)
          {
            S[r_subset][iii] = tmp_str[iii];
          }
      
          //Determine a new cross-quadruples with updated data set:
          for(unsigned int i=0; i<combnum; ++i)
          {
            decisive = false;
            //Printing result
            for(int jj=0; jj < m; jj++)
            {
              triplet[jj] = X[T[i][jj]];
            }
            
            for(int j = 0; j < k; j++)
            {
              if(sizes[j] >= 3) 
              {
                if(issubset(triplet, S[j], 3, sizes[j]) == true) {decisive = true;}
              }
            }
            if(decisive == false) 
            {
              res = "Computation is done. The given set is not decisive.";
              break;
            }
    	    }
            
          if(decisive == true) 
          {
            res = "Computation is done. The given data set is decisive.";
            break;
          }
          else
          {
            res = "Computation is done. The given data set is not decisive.";
          }
          iteration+=1;
        }
      }
      else
      {
        res = "Computation is done. The given data set is not decisive.";
      }
      
    }
    
    cout << res << endl;
    
    Rcpp::List s_out(k+1);
    for(int i=0; i<k; i++) 
    {     
       vector<string> s_i(sizes[i]);
       for(int j=0; j<sizes[i]; ++j)
       {
         s_i[j] = S[i][j];
       }
       Rcpp::CharacterVector ss_i(sizes[i]);
       std::copy( s_i.begin(), s_i.end(), ss_i.begin() ) ;
       s_out[i] = ss_i;
    }
    
    s_out[k] = res;
    
    free(C);
    free(sizes);
    //free(S);
    
    //return(Rcpp::CharacterVector(res));
    return(Rcpp::wrap(s_out));
}

RcppExport SEXP IsDecisiveUnrooted(SEXP taxa, SEXP s, SEXP n, SEXP _k, SEXP fflag) {
    unsigned int N; //Total number of taxons
    string *X = NULL;
    string **S = NULL;
    int *sizes = NULL; //a set of sizes of each of taxon subset
    int **Q = NULL;
    int k = 0; //number of subsets
    vector<int> RCQ; //Resolved cross-quadruples
    vector<int> CQ; //This holds cross-quadrooples
    std::map< int, int > cq_map_freq; //This holds cross-quadruples and their "weights": <Q_id , weight>
      
    /************************ Data receiving **************************/
    
    N = as<unsigned int>(n); //Total Number of taxons
    Rcpp::CharacterVector cx = Rcpp::CharacterVector(taxa);
    X = new string[N];
    for (int i=0; i<cx.size(); i++) {
      X[i] = cx[i];
    }
    bool fix_flag = as< bool >(fflag);
    
    //Data set
    Rcpp::List SS = Rcpp::List(s);   
    k= SS.size(); //Number of taxon subsets
    
    try
    {
        S = new string*[k];
        sizes = new int[k];
    } 
    catch (bad_alloc& ba) 
    {
        cerr << "bad_alloc caught: " << ba.what() << endl;
    }
    
    for (int i=0; i<k; i++) {  
        Rcpp::CharacterVector SSS = SS[i];
        sizes[i] = SSS.size(); 
        S[i] = new string[sizes[i]];
        for (int ii=0; ii<sizes[i]; ii++) {  
          S[i][ii] = SSS[ii];
        }
    }
    
    /*****************   End of data receiving   ***********************/
    
    //Unrooted tree case
    unsigned int m = 4; //We consider quadruples.
    
    //Number of possible combinations (i.e. how many quadruples can be created from a given taxon set X):
    unsigned int combnum = Res(N,m);
    cout << "Number of combinations is: " << combnum << endl;
    
    Q = new int*[combnum]; //Array of quadrooples. This array keeps all quadruples created from X
    
    //Quadroople: an array that holds 4 taxons. Here we create an initial quadruple: ["A", "B", "C", "D",....]. Our later task is to re-arrange the 
    //taxons in this quadruple
    string quadruple[m];
    int *C = NULL;
    C = new int[m];
    
    for(int i=0; i<m; i++)
    {
      C[i] = i;
      quadruple[i] = X[C[i]];
    }
    
    //Here we re-arrange, i.e. create new quadruples and also determine cross-quadruples
    
    //Determine Cross-quadrooples
    for(unsigned int i=0; i<combnum; ++i)
    {
        int freq_counter = 0;
        
        if (i>0) GetNext( C, N, m ); //Getting next quadruple, with different arrangement of the taxons
        
        Q[i] = new int[m]; //New quadruple
        for(int j=0;j<m;++j) { Q[i][j] = C[j]; } //Fill this new quadruple with taxons ids
        
        bool cq_flag = true;      
	      for(int j=0; j<k;++j) //For each subset S[j] of the S
    	  {
		      if(sizes[j] >= 4) 
          {
            if(issubset(quadruple, S[j], 4, sizes[j]) == true) {cq_flag = false; break;}
		      }
          
          freq_counter += diff2(quadruple, S[j], 4, sizes[j]).size();
          
    	  }
        
        //If this particular quadruple is not a subset of any of the taxon set S[j], add this quadruple to the set of CQs (Cross-Quadruples)
	      if(cq_flag == true) 
        {
          CQ.push_back(i); 
          cq_map_freq.insert(std::make_pair(i,freq_counter));
        }
    	
	      
  
        for(int jj=0; jj < m; jj++)
    	  {
		      quadruple[jj] = X[C[jj]];
    	  }
        
    }
	
    //Debug
    //Pringting Qs:
    /*for(int i=0; i<combnum; ++i) {
    	for(int j=0; j<4;++j) {cout << Q[i][j] << " ";}
    	cout << endl;
    }*/

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
    
    string res_str = "";
    
    bool res = GetDecisivenessUnrooted(N, X, S, sizes, Q, k, CQ, RCQ);
    
    if(res == true) 
    {
      res_str = "Computation is done. The given data set is decisive.";
    }
    else if((res == false) && (fix_flag == true))
    {
      
      cout << "The given data set is not decisive. Trying to fix it..." << endl;
    
      //Fixing the data set
      //Sorting cq_map_freq by value:
      std::map<int, int> dst = flip_map(cq_map_freq);
      std::map<int, int>::reverse_iterator iter;
      int iteration = 1;
      for (iter = dst.rbegin(); iter != dst.rend(); ++iter) 
      {
        //Get the random number from 1..k:
        int r_subset = rand() % k;
        sizes[r_subset] = sizes[r_subset] + 4; //sizes[i] = SSS.size();
        string *tmp_str = new string[sizes[r_subset]];
        for(unsigned int mm=0; mm<sizes[r_subset]-4; ++mm)
        {
          tmp_str[mm] = S[r_subset][mm];
        }
        cout << "\033[A\033[2K";
        cout << "Iteration " << iteration << "..." << endl;
        for(int iii=0; iii<4; ++iii)
        {
          tmp_str[sizes[r_subset]-4+iii] = X[Q[iter->second][iii]];
        }
        S[r_subset] = new string[sizes[r_subset]];
        for(int iii=0; iii<sizes[r_subset]; ++iii)
        {
          S[r_subset][iii] = tmp_str[iii];
        }
      
        //Determine a new cross-quadruples with updated data set:
        CQ.clear();  
        RCQ.clear(); 
        for(unsigned int i=0; i<combnum; ++i)
        {
          int freq_counter = 0;
          for(int jj=0; jj < m; jj++)
          {
		        quadruple[jj] = X[Q[i][jj]];
    	    }
        
          bool cq_flag = true;      
          for(int j=0; j<k;++j) //For each subset S[j] of the S
    	    {
		        if(sizes[j] >= 4) 
            {
              if(issubset(quadruple, S[j], 4, sizes[j]) == true) {cq_flag = false; break;}
		        }
            //freq_counter += diff2(quadruple, S[j], 4, sizes[j]).size();
          }
        
          //If this particular quadruple is not a subset of any of the taxon set S[j], add this quadruple to the set of CQs (Cross-Quadruples)
	        if(cq_flag == true) 
          {
            CQ.push_back(i); 
            //cout << freq_counter << endl;
            //cq_map_freq.insert(std::make_pair(i,freq_counter));
          }
    	  }
      
        res = GetDecisivenessUnrooted(N, X, S, sizes, Q, k, CQ, RCQ);
        
      
        if(res == true) 
        {
          res_str = "Computation is done. The given data set is decisive.";
          break;
        }
        else
        {
          res_str = "Computation is done. The given data set is not decisive.";
        }
        iteration+=1;
      }
    }
    else
    {
      res_str = "Computation is done. The given data set is not decisive.";
    }
    cout << res_str << endl;
    
    Rcpp::List s_out(k+1);// = Rcpp::List(S); 
    for(int i=0; i<k; i++) 
    {     
       //string *s_i = new string[sizes[i]];
       vector<string> s_i(sizes[i]);
       for(int j=0; j<sizes[i]; ++j)
       {
         s_i[j] = S[i][j];
       }
       Rcpp::CharacterVector ss_i(sizes[i]);
       std::copy( s_i.begin(), s_i.end(), ss_i.begin() ) ;
       s_out[i] = ss_i;
    }
    
    s_out[k] = res_str;
    
    free(C);
    free(Q);
    free(sizes);
    //free(S);
    //free(X);
    
    
    return(Rcpp::wrap(s_out));
}


