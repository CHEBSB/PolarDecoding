/*
Simulation of successive cancellation (SC) decoder
for polar code with N = 128 and rate = 0.5
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define bSNR_dB 3.046          // Eb/N0 in dB
#define N 128
#define K 64
#define n 7

// seed for generating random number in (0, 1)
const unsigned long long SEED = 1024; 
unsigned long long RANV;
int RANI = 0;
double n1, n2;              // gaussian noise
// standard deviation of Gaussian noise
const double std = pow(10, bSNR_dB / ((double)-20)); 
// 5G polar sequence (Q[0] is most unreliable)
const int Q[N] = {
0, 1, 2, 4, 8, 16, 32, 3, 5, 64, 9, 6, 17, 10, 18, 12, 33,
65, 20, 34, 24, 36, 7, 66, 11, 40, 68, 19, 13, 48, 14, 72,
21, 35, 26, 80, 37, 25, 22, 38, 96, 67, 41, 28, 69, 42, 49,
74, 70, 44, 81, 50, 73, 15, 52, 23, 76, 82, 56, 27, 97, 39,
84, 29, 43, 98, 88, 30, 71, 45, 100, 51, 46, 75, 104, 53,
77, 54, 83, 57, 112, 78, 85, 58, 99, 86, 60, 89, 101, 31,
90, 102, 105, 92, 47, 106, 55, 113, 79, 108, 59, 114, 87,
116, 61, 91, 120, 62, 103, 93, 107, 94, 109, 115, 110, 117,
118, 121, 122, 63, 124, 95, 111, 119, 123, 125, 126, 127};
// n-kronecker product of F
int **Fn; 
// LLR of each node
double L[n + 1][N];
// whether LLR of a node has been computed
int Done[n + 1][N];

// generate a uniform random number
double Ranq1();            
// generate 2 normal random number given standard deviation
void normal();
// check node operation
double CHK(double L1, double L2);
// successive cancellation decoder
int SCdecode(double *y);

int main(void)
{
    int i, j, k;                // looping indices
    int m = 0;                  // index for PN sequence
    // stepsize of m after 1 iteration
    int step_m = K % 63;                
    // to generate PN sequence
    int U[] = {0, 0, 0, 0, 0, 0};   // previous bits
    int b;                      // current bit
    int PN[63];                 // 1 period of PN sequence
    int I[K];                   // information set
    int B[K];                   // bit reversal
    int r[n];                   // to compute bit reversal
    int errBlock = 0;           // number of block errors
    int temp;                   // temporary storage
    // int w[K];                   // information bits
    int u[N];                   // info bits + frozen bits
    int v[N];                   // after bit reversal
    int x[N];                   // polar codeword
    double y[N];                // codeword + Gaussian noise

    // generate a cycle of PN sequence
    for (i = 0; i < 63; i++) {
        if (i == 0) b = 1;
        else if (i < 6) b = 0;
        else b = (U[4])? (!U[5]): U[5]; // exculsive or
        PN[i] = b;                      // write the array
        // flip-flop shift
        U[5] = U[4];        
        U[4] = U[3];
        U[3] = U[2];
        U[2] = U[1];
        U[1] = U[0];
        U[0] = b;
    }
    // determine the information set I
    for (i = 0; i < K; i++) {
        // pick most reliable channels
        I[i] = Q[N - 1 - i];
    }
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
    // init both u and v as all-0
    for (i = 0; i < N; i++) {
        u[i] = 0;
        v[i] = 0;
    }
    // run simulation until desired error blocks
    for (i = 0; errBlock < 50; i++) {
        // Encoder
        /* use PN sequence to have K bits, and
        put them into information set in u */
        for (j = 0; j < K; j++)
            u[I[j]] = PN[(m + j) % 63];
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
        // Channel
        // y = x + noise
        for (j = 0; j < N; j += 2) {
            normal();           // generate Gaussian noise
            if (x[j] == 0) y[j] = 1 + n1;
            else y[j] = -1 + n1;
            if (j + 1 < N) {
                if (x[j + 1] == 0) y[j + 1] = 1 + n2;
                else y[j + 1] = -1 + n2;
            }   
        }
        // SC decoder
        SCdecode(y);
        m += step_m;                   // increase m
        errBlock += 1;
    }
    return 0;
}

// generate an uniform random number in (0, 1)
double Ranq1()
{
    if (RANI == 0) {
        RANV = SEED ^ 4101842887655102017LL;
        RANV ^= RANV >> 21;
        RANV ^= RANV << 35;
        RANV ^= RANV >> 4;
        RANV = RANV * 2685821657736338717LL;
        RANI++;
    }
    RANV ^= RANV >> 21;
    RANV ^= RANV << 35;
    RANV ^= RANV >> 4;
    return RANV * 2685821657736338717LL * 5.42101086242752217E-20;
}

// generate 2 normal random number given standard deviation
void normal()
{
    double x1, x2, s;       // x1 and x2 are Uniform(0, 1)

    do {
        x1 = Ranq1();
        x2 = Ranq1();
        x1 = 2 * x1 - 1;    // transform (0, 1) to (-1, 1)
        x2 = 2 * x2 - 1;
        s = x1 * x1 + x2 * x2;
    } while (s >= 1.0);     // until s is in unit circle
    n1 = std * x1 * sqrt(-2 * log(s) / s);
    n2 = std * x2 * sqrt(-2 * log(s) / s);
    return;
}

// check node operation
double CHK(double L1, double L2)
{
    int s1, s2;         // sign of L1, L2
    double A1, A2;      // A1 = |L1|, A2 = |L2|
    double delta;       // the complicated term
    double sAbs, dAbs;  // |L1 + L2|, |L1 - L2|

    sAbs = fabs(L1 + L2);
    dAbs = fabs(L1 - L2);
    if (sAbs < 0.196) delta = 0.65;
    else if (sAbs < 0.433)  delta = 0.55;
    else if (sAbs < 0.71)   delta = 0.45;
    else if (sAbs < 1.05)   delta = 0.35;
    else if (sAbs < 1.508)  delta = 0.25;
    else if (sAbs < 2.252)  delta = 0.15;
    else if (sAbs < 4.5)    delta = 0.05;
    else delta = 0;
    if (dAbs < 0.196) delta -= 0.65;
    else if (dAbs < 0.433)  delta -= 0.55;
    else if (dAbs < 0.71)   delta -= 0.45;
    else if (dAbs < 1.05)   delta -= 0.35;
    else if (dAbs < 1.508)  delta -= 0.25;
    else if (dAbs < 2.252)  delta -= 0.15;
    else if (dAbs < 4.5)    delta -= 0.05;
    A1 = fabs(L1);
    A2 = fabs(L2);
    s1 = (L1 >= 0)? 1: -1;
    s2 = (L2 >= 0)? 1: -1;
    if (A1 > A2)
        return s1 * s2 * A2 + delta;
    return s1 * s2 * A1 + delta;
}

// LLR-based successive cancellation decoder
int SCdecode(double *y)
{
    int i, j;                   // looping indices

    // init LLR; set Done to 0 (undone yet)
    for (i = 0; i < n; i++) {
        for (j = 0; j < N; j++) {
            L[i][j] = 0;
            Done[i][j] = 0;
        }
    }
    // compute channel LLR
    for (j = 0; j < N; j++) {
        L[n][j] = 2 * y[j] / std / std;;
        Done[n][j] = 1;
    }

    return 0;
}