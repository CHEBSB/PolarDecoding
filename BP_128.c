/*
Simulation of successive cancellation (SC) decoder
for polar code with N = 128 and rate = 0.5
Follow factor graph in Lee's thesis
G = F^{\otimes n} (no bit reversal)
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double bSNR_dB;             // Eb/N0 in dB
#define N 128
#define K 64
#define n 7
#define iterMax 100

typedef struct node {
    double l;       // left-propagating message (LLR)
    double r;       // right-propagating message
    /* Left & right position, i.e., when it appears at
    the left/right side of a BCB, is it the upper or
    the lower node.
    1 for upper, 0 for lower
    */
    int leftP;
    int rightP;
    struct node *pU;        // upper parent
    struct node *pL;        // lower parent
    struct node *cU;        // upper child
    struct node *cL;        // lower child
} node;

// seed for generating random number in (0, 1)
const unsigned long long SEED = 1024; 
unsigned long long RANV;
int RANI = 0;
double n1, n2;              // gaussian noise
// standard deviation of Gaussian noise
double std; 
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
int I[K];                   // information set
int inI[N];                 // whether a bit is in I
// n-kronecker product of F
int **Fn; 
// all nodes on factor graph
node ***V;
// whether a node is init
int initV[n + 1][N];
// x from decoder output
int x_hat[N];

// generate a uniform random number
double Ranq1();            
// generate 2 normal random number given standard deviation
void normal();
// 2 to the power of the given exponent
int pow2(int expo);
// check node operation
double CHK(double L1, double L2);
// connect the factor graph
void connectBCB(int i, int j);
// iterative BP decoder
void BP(double *y, int *u_hat);

int main(void)
{
    int i, j, k;                // looping indices
    int run;                    // sumulation index
    int m = 0;                  // index for PN sequence
    // stepsize of m after 1 iteration
    const int step_m = K % 63;                
    // to generate PN sequence
    int U[] = {0, 0, 0, 0, 0, 0};   // previous bits
    int b;                      // current bit
    int PN[63];                 // 1 period of PN sequence
    int errBlock = 0;           // number of block errors
    int temp;                   // temporary storage
    int u[N];                   // encoder input
    int x[N];                   // polar codeword
    double y[N];                // codeword + Gaussian noise
    int u_hat[N];               // decoder's output
    int errbit = 0;             // # of bit error

    // allocate memory for V
    V = (node ***)calloc(n + 1, sizeof(node **));
    for (i = 0; i <= n; i++) {
        V[i] = (node **)calloc(N, sizeof(node *));
        for (j = 0; j < N; j++) 
            V[i][j] = (node *)calloc(1, sizeof(node));
    }
    // mark all nodes as uninitialized
    for (i = 0; i <= n; i++)
        for (j = 0; j < N; j++)
            initV[i][j] = 0;
    // init pointers in V (= connect the factor graph)
    // the leftmost stage
    for (j = 0; j < N; j++) {
        V[0][j]->pU = NULL; // first layer has no parent
        V[0][j]->pL = NULL;
        connectBCB(0, j);
    }
    for (i = 1; i < n; i++)
        for (j = 0; j < N; j++) 
            connectBCB(i, j);
    // the rightmost stage has no child
    for (j = 0; j < N; j++) {
        V[n][j]->cU = NULL;
        V[n][j]->cL = NULL;
    }
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
    // at first, no channels are picked
    for (i = 0; i < N; i++)
        inI[i] = 0;
    // determine the information set I
    for (i = 0; i < K; i++) {
        // pick most reliable channels
        I[i] = Q[N - K + i];
        inI[I[i]] = 1;
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
    // printf("Fn init completed.\n");     // for debug
    /* set u to all-0 (frozen bits
    should keep 0 throughout simualtion) */
    for (i = 0; i < N; i++)
        u[i] = 0;
for (bSNR_dB = 1.0; bSNR_dB <= 4; bSNR_dB += 0.5) {
    errBlock = 0;
    std = pow(10, bSNR_dB / ((double)-20));
    // run simulation until desired error blocks
    for (run = 0; errBlock < 100; run++) {
        // reset vectors to all-zero
        for (i = 0; i < N; i++) {
            u_hat[i] = 0;
            x[i] = 0;
        }
        // Encoder
        // use all-zero as the first step
        for (i = 0; i < K; i++) {
            u[I[i]] = PN[(m + i) % 63];
        }
        // encode by Fn
        for (i = 0; i < N; i++)
            if (u[i] == 1) {
                // add the i-th row to x
                for (j = 0; j < N; j++) {
                    // modulo-2 addition
                    x[j] += Fn[i][j];
                    if (x[j] > 1) x[j] -= 2;
                }
            }
        // Channel
        // y = x + noise
        for (i = 0; i < N; i += 2) {
            normal();           // generate Gaussian noise
            if (x[i] == 0) y[i] = 1 + n1;
            else y[i] = -1 + n1;
            if (i + 1 < N) {
                if (x[i + 1] == 0) y[i + 1] = 1 + n2;
                else y[i + 1] = -1 + n2;
            }   
        }
        // run the BP decoder
        BP(y, u_hat);
        // check info bits for block error
        temp = 0;           // flag for loop
        for (i = 0; i < K; i++) {
            if (u[I[i]] != u_hat[I[i]]) {
                temp = 1;
                errbit += 1;
            }
        }
        errBlock += temp;
        m += step_m;                   // increase m
        if (m >= 63) m -= 63;
    }
    // final output
    printf("bSNR = %.2lf\terror block = %d\trun = %d\t",
        bSNR_dB, errBlock, run);
    printf("BLER = %lf\n", ((double)errBlock) / run);
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

// 2 to the power of the given (positive) exponent
int pow2(int expo)
{
    int result = 1;         // the result to return

    if (expo < 0)
        printf("Positive exponent required!\n");
    while (expo > 0) {
        result = result * 2;
        expo -= 1;
    }
    return result;
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
    delta = 0;      // for debug
    if (A1 > A2)
        return s1 * s2 * A2 + delta;
    return s1 * s2 * A1 + delta;
}

// assign pointers given BCB: require i < n
void connectBCB(int i, int j)
{
    if (initV[i][j] == 1) {     // already done
        V[i][j]->leftP = 0;     // it is the lower node
        return;
    }
    V[i][j]->leftP = 1;         // else, upper node
    initV[i][j] = 1;            // mark upperleft node
    initV[i][j + pow2(i)] = 1;  // mark lowerleft node
    // assign childs
    // for the upperleft bit
    V[i][j]->cU = V[i + 1][j];
    V[i][j]->cL = V[i + 1][j + pow2(i)];
    V[i + 1][j]->rightP = 1;            // set their rightP
    V[i + 1][j + pow2(i)]->rightP = 0;
    // for the lowerleft bit
    V[i][j + pow2(i)]->cU = V[i + 1][j];
    V[i][j + pow2(i)]->cL = V[i + 1][j + pow2(i)];
    // assign parents conversely
    V[i + 1][j]->pU = V[i][j];
    V[i + 1][j]->pL = V[i][j + pow2(i)];
    V[i + 1][j + pow2(i)]->pU = V[i][j];
    V[i + 1][j + pow2(i)]->pL = V[i][j + pow2(i)];
    return;
}

// conventional BP decoder for polat codes
void BP(double *y, int *u_hat)
{
    int i, j;       // looping indices
    int iter;       // iteration index
    // initialize
    // left-propagating message
    for (i = 0; i < n; i++) 
        for (j = 0; j < N; j++)
            V[i][j]->l = 0; 
    for (j = 0; j < N; j++)
        V[n][j]->l = 2 * y[j] / std / std;
    // right-propagating message
    for (i = 1; i <= n; i++) 
        for (j = 0; j < N; j++)
            V[i][j]->r = 0; 
    for (j = 0; j < N; j++) {
        if (inI[j] == 0)        // frozen bit => infinity
            V[0][j]->r = 999;
        else                    // info bit => 0
            V[0][j]->r = 0;   
    }
    for (iter = 0; iter < iterMax; iter++) {
        // R propag
        for (i = 0; i < n; i++) 
            for (j = 0; j < N; j++) {
                // only start from upperleft
                if (V[i][j]->leftP == 1) {
                    (V[i][j]->cU)->r = CHK(V[i][j]->r,
                        (V[i][j]->cL)->l + (V[i][j + pow2(i)])->r);
                    (V[i][j]->cL)->r = (V[i][j + pow2(i)])->r
                        + CHK(V[i][j]->r, (V[i][j]->cU)->l);
                }
            }
        // L propag
        for (i = n - 1; i >= 0; i--) 
            for (j = 0; j < N; j++) {
                // only start from upperleft
                if (V[i][j]->leftP == 1) {
                    V[i][j]->l = CHK((V[i][j]->cU)->l,
                        (V[i][j]->cL)->l + (V[i][j + pow2(i)])->r);
                    V[i][j + pow2(i)]->l = (V[i][j]->cL)->l
                        + CHK(V[i][j]->r, (V[i][j]->cU)->l);
                }
            }
    }
    for (j = 0; j < N; j++) {
        if (inI[j] == 0) {      // frozen bit
            u_hat[j] = 0;
        } else if (V[0][j]->l + V[0][j]->r >= 0) {
            u_hat[j] = 0;
        } else {
            u_hat[j] = 1;
        }
    }
    return;
}