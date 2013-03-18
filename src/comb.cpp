#include "comb.h"

void GetNext(int *C, unsigned int N, unsigned int k )
{
  int i, j;
  i = k-1;
  
  while( (C[i] + k - i + 1) > N) 
  {
    /*Search for the next element to increase*/
    i--;
  }
  
  C[i] = C[i] + 1;
  
  for( j = i+1; j < k; j++ )
  {
    C[j] = C[j-1] + 1;
  }

  //Pringting C:
  //for(int j=0; j<k;++j) {cout << C[j] << " ";}
  //cout << endl;
  
}

unsigned int Res(unsigned int N, unsigned int k) 
{
  /*Calculates the number of possible combinations*/
  int i, j;
  int *A, *B;
  A = new int[k+1];
  B = new int[k+1];
  for(int i=0; i<k+1; i++)
  {
      A[i] = 0;
      B[i] = 0;
  }
  
  A[0] = 1; A[1] = 1;
  
  for(i=1; i<N; i++)
  {
    B[0] = 1; B[1] = 1;
    
    for(j=1; j<k+1; j++)
    {
      B[j] = A[j] + A[j-1];
    }
    
    for(int ii=0; ii<k+1; ii++)
    {
      A[ii] = B[ii];
    }
  }
  
  unsigned long cnum = A[k]; 
  
  free(A);
  free(B);
  
  return cnum;
}
