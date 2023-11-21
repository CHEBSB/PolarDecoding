/*
Gaussian Approximation of density evolution
for BP decoding of polar codes where N = 128, K = 64.
Assume all zero codewords transmitted
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
// numerically compute the phi function used in DE-GA
double phi(double mean);
// the inverse function of phi
double phi_inv(double x);
// the derivative of phi, used in computing phi inverse
double derivative_phi(double x);
// Density evolution with Gaussian approximation for BP
void BPDEGA();

int main(void)
{
    int i, j, k;                // looping indices
    int temp;                   // temporary storage

    // allocate memory for V
    V = (node ***)calloc(n + 1, sizeof(node **));
    for (i = 0; i <= n; i++) {
        V[i] = (node **)calloc(N, sizeof(node *));
        for (j = 0; j < N; j++) 
            V[i][j] = (node *)calloc(1, sizeof(node));
    }
// Initialize the factor graph
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
    // at first, no channels are picked
    for (i = 0; i < N; i++)
        inI[i] = 0;
    // determine the information set I
    for (i = 0; i < K; i++) {
        // pick most reliable channels
        I[i] = Q[N - K + i];
        inI[I[i]] = 1;
    }
for (bSNR_dB = 3.0; bSNR_dB <= 3; bSNR_dB += 0.5) {
    std = pow(10, bSNR_dB / ((double)-20));
    // run DEGA
    BPDEGA();
    // final output
    printf("For debug.\n");
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
// piecewise approximation of phi function
double phi(double x)
{
    if (x < 0) {
        printf("illegal input for phi function!\n");
        return 1;
    }
	if (x <= 0.1910)
		return exp(0.1047 * x * x - 0.4992 * x);
	if (x <= 0.7420)
		return 0.9981 * exp(0.05315 * x * x - 0.4795 * x);
	if (x <= 9.2254)
		return exp(-0.4527 * pow(x, 0.86) + 0.0218);
	if (x <= 15)
		return exp(-0.2832 * x - 0.4254);
	if (x <= 25)
		return exp(-0.26725134794 * x - 0.6646297809);
	return sqrt(3.14159265 / x) * exp(-x / 4) * (1 - 10.0 / 7.0 / x);
}
// piecewise approximation of inverse of phi function
double phi_inv(double x)
{
	double x0, x1, delta, epsilon;
	
	if (x <= 1 && x >= 0.91253609394)
		return (0.4992 - sqrt(0.24920064 + 0.4188 * log(x)))
            / 0.2094;
	if (x >= 0.72005453218)
		return (0.4795 - sqrt(0.22992025 + 0.2126 * log(x / 0.9981)))
            / 0.1063;
	if (x >= 0.04792905738)
		return pow((log(x) - 0.0218) / -0.4527, 1 / 0.86);
	if (x >= 0.00934045792)
		return -(log(x) + 0.4254) / 0.2832;
	if (x >= 0.0006452237)
		return -(log(x) + 0.6646297809) / 0.26725134794;
	else {
		// x0 = 25;
		x1 = 25 - (phi(25) - x) / derivative_phi(25);
		delta = fabs(x1 - 25);
		// epsilon = 1e-3;
		while (delta >= 1e-3) {
			x0 = x1;
			x1 = x1 - (phi(x1) - x) / derivative_phi(x1);
			delta = fabs(x1 - x0);
		}
		return x1;
	}
} 
// derivative of the piecewise approximation of phi function
double derivative_phi(double x)
{
    if (x < 0) {
        printf("illegal input for phi's derivative'!\n");
        return 1;
    }
	if (x <= 0.1910)
		return (0.2094 * x - 0.4992)
            * exp(0.1047 * x * x - 0.4992 * x);
	if (x <= 0.7420)
		return 0.9981 * (0.1063 * x - 0.4795)
            * exp(0.05315 * x * x - 0.4795 * x);
	if (x <= 9.2254)
		return -0.389322 * exp(0.0218 - 0.4527 * pow(x, 0.86))
            / pow(x, 0.14);
	if (x <= 15)
		return -0.2832 * exp(-0.2832 * x - 0.4254);
	if (x <= 25)
		return -0.26725134794
            * exp(-0.26725134794 * x - 0.6646297809);	
	return exp(-x / 4) * sqrt(3.14159265 / x)
    * (-0.5 / x * (1 - 10.0/7.0 / x) - 0.25 * (1 - 10.0 / 7.0 / x)
        + 10.0 / 7.0 / x / x);
}

// DEGA for iterative BP decoder
void BPDEGA()
{
    int i, j;       // looping indices
    int iter;       // iteration index
    // initialize
    // left-propagating message
    for (i = 0; i < n; i++) 
        for (j = 0; j < N; j++)
            V[i][j]->l = 0; 
    for (j = 0; j < N; j++)
        V[n][j]->l = 2 / std / std;
    // right-propagating message
    for (i = 1; i <= n; i++) 
        for (j = 0; j < N; j++)
            V[i][j]->r = 0; 
    for (j = 0; j < N; j++) {
        if (inI[j] == 0)        // frozen bit => infinity
            V[0][j]->r = 99999;
        else                    // info bit => 0
            V[0][j]->r = 0;   
    }
    for (iter = 0; iter < iterMax; iter++) {
        // R propag
        for (i = 0; i < n; i++) 
            for (j = 0; j < N; j++) {
                // only start from upperleft
                if (V[i][j]->leftP == 1) {
                    V[i][j]->cU->r = phi_inv(
                    phi(V[i][j]->r)
                    + phi(V[i][j + pow2(i)]->r + V[i][j]->cL->l)
                    - phi(V[i][j]->r)
                        * phi(V[i][j + pow2(i)]->r + V[i][j]->cL->l)
                    );
                    V[i][j]->cL->r = phi_inv(
                    phi(V[i][j]->r) + phi(V[i][j]->cU->l)
                    - phi(V[i][j]->r) * phi(V[i][j]->cU->l)
                    ) + V[i][j + pow2(i)]->r;
                }
            }
        // L propag
        for (i = n - 1; i >= 0; i--) 
            for (j = 0; j < N; j++) {
                // only start from upperleft
                if (V[i][j]->leftP == 1) {
                    V[i][j]->l = phi_inv(
                    phi(V[i][j]->cU->l)
                    + phi(V[i][j + pow2(i)]->r + V[i][j]->cL->l)
                    - phi(V[i][j]->cU->l)
                        * phi(V[i][j + pow2(i)]->r + V[i][j]->cL->l)
                    );
                    V[i][j + pow2(i)]->l = phi_inv(
                    phi(V[i][j]->r) + phi(V[i][j]->cU->l)
                    - phi(V[i][j]->r) * phi(V[i][j]->cU->l)
                    ) + V[i][j]->cL->l;
                }
            }
    }
    return;
}