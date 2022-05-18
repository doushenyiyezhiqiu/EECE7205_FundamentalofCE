#include <iostream>
#include <vector>
#include <math.h>

using namespace std; 

void printValue(int v[], int l) {
    int i;
    std::cout<<"list: {";
    for (i=0; i<l; i++){
        std::cout<<v[i] << " ";
    }
    std::cout<<"}"<<std::endl;
}   

void radixSort(int A[], int B[], int l){
    int index[10];
    int output[l];
    for (int i=0;i<10;i++){ //set index to zero
        index[i]=0;
    }
    for (int i=0;i<l;i++){ //count how many digits
        index[B[i]]=index[B[i]]+1;
    }
    for (int i=1;i<10;i++){ //create CDF
        index[i]=index[i]+index[i-1];
    }
    for (int i=l-1;i>=0;i--){ //sort A accounting to B
        output[index[B[i]]-1] = A[i];
        index[B[i]] = index[B[i]]-1;
    }    
    for (int i=0;i<l;i++){ // copy output back to A
        A[i]=output[i];
    }
}

int main()
{
    int A[] = {329, 457, 657, 839, 436, 720, 353};
    int l = sizeof(A)/sizeof(*A);
    int index[10];
    int output[l];
    int s = 3;

    int B[l]; // single digit of the array
    for (int i=0;i<l;i++){
        B[i]=0;
    }

    for (int exp_v=0;exp_v<s;exp_v++){ // loop through # of digits
        for (int i=0;i<l;i++){ //get the single digit value
            int temp = pow(10,exp_v);
            B[i] = (A[i]/temp)%10;
        }
        radixSort(A,B,l); //sorted base on B
        printValue(A,l);

    }
}