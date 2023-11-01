/*
Encoder of Polar code. Take N = 4 for example
*/
#include <stdio.h>
#include <stdlib.h>

#define N 4
#define K 3
#define n 2

// n-kronecker product of F
const int Fn[4][4] = {
    {1, 0, 0, 0},
    {1, 1, 0, 0},
    {1, 0, 1, 0},
    {1, 1, 1, 1}
}; 

int main(void)
{
    int i, j, k;                // looping indices           
    int I[K] = {1, 2, 3};       // information set
    int B[K];                   // bit reversal
    int r[n];                   // to compute bit reversal
    int temp;                   // temporary storage
    int u[N];                   // info bits + frozen bits
    int v[N];                   // after bit reversal
    int x[N];                   // polar codeword

    // obtain B, the bit reversal for info bits only
    // B[i] = bit reversal of I[i]
    for (i = 0; i < K; i++) {
        temp = I[i];
        for (j = 0; j < n; j++) {
            // find the least significant bit
            r[j] = temp % 2;    
            // right shift by 1 bit
            temp = (temp  - r[j]) / 2;
        }
        B[i] = 0;
        temp = 1;
        for (j = 0; j < n; j++) {   // bit-reversed value
            B[i] += (r[n - j - 1] * temp);   
            temp = temp * 2;
        }
    }
    /*
    // read Fn from file
    Fn = (int **)calloc(N, sizeof(int *));
    for (i = 0; i < N; i++) {
        Fn[i] = (int *)calloc(N, sizeof(int));
        for (j = 0; j < N; j++) {
            scanf("%d", &temp);
            if (temp != 0 && temp != 1)
                printf("Illegal input!\n");
            Fn[i][j] = temp;
        }
    }
    printf("Fn init completed.\n");     // for debug
    */
    // init both u and v as all-0
    for (i = 0; i < N; i++) {
        u[i] = 0;
        v[i] = 0;
    }
    // Encoder
    /* input u by hand */
   u[1] = 1;
   u[2] = 0;
   u[3] = 1;
    // perform bit reversal
    for (j = 0; j < K; j++)
        v[B[j]] = u[I[j]];
    // reset x to all-0
    for (j = 0; j < N; j++)
        x[j] = 0;
    // encode by Fn
    for (j = 0; j < K; j++)
        if (v[B[j]] == 1)
            // add the B[j]-th row to x
            for (k = 0; k < N; k++) {
                // modulo-2 addition
                x[k] += Fn[B[j]][k];
                if (x[k] > 1) x[k] -= 2;
            }
    // Print 
    printf("u:");
    for (j = 0; j < N; j++) 
        printf(" %d", u[j]);
    printf("\nv:");
    for (j = 0; j < N; j++) 
        printf(" %d", v[j]);
    printf("\nx:");
    for (j = 0; j < N; j++) 
        printf(" %d", x[j]);

    return 0;
}
