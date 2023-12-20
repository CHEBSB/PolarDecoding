/*
Gaussian Approximation of density evolution
for BPR decoding of polar codes where N = 128, K = 64.
Assume all zero codewords transmitted
Evaluate the sum of all bit error rates for each stage
(whenever it is frozen or not)
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double bSNR_dB;             // Eb/N0 in dB
#define N 128
#define K 64
#define n 7
#define iterMax 30

typedef struct node {
    double l;           // left-propagating message (LLR)
    double r;           // right-propagating message
    /* Left & right position, i.e., when it appears at
    the left/right side of a BCB, is it the upper or
    the lower node.
    1 for upper, 0 for lower
    */
    double u[n + 1];    // leftward LLR for BPR analysis
    int leftP;
    int rightP;
    struct node *pU;        // upper parent
    struct node *pL;        // lower parent
    struct node *cU;        // upper child
    struct node *cL;        // lower child
} node;

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
int bRev[N];                // bit reversal for [0, 127]
// n-kronecker product of F
int **Fn; 
// all nodes on factor graph
node ***V;
// whether a node is init
int initV[n + 1][N];
// union-bound-BLER of candidate from each stage
double E[n + 1];

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
// Density evolution with Gaussian approximation for BPR
void BPRDEGA();

int main(void)
{
    int i, j, k;                // looping indices
    int temp;                   // temporary storage
    double bler;                // block error rate

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
    // build the bit reversal
    for (i = 0; i < N; i++) {
        bRev[i] = 0;            // set to 0 at first
        temp = i;               // the integer to reverse
        for (j = n - 1; j >= 0; j--) {
            if (temp % 2 == 1)  // if this bit = 1
                bRev[i] += pow2(j);
            temp = temp / 2;    // go to the next bit
        }
    }
    // for debug
    for (i = 0; i < N; i++)
        if (bRev[bRev[i]] != i)
            printf("Bit reversal in error!\n");
    // at first, no channels are picked
    for (i = 0; i < N; i++)
        inI[i] = 0;
    // determine the information set I
    for (i = 0; i < K; i++) {
        // pick most reliable channels
        I[i] = Q[N - K + i];
        inI[I[i]] = 1;
    }
printf("iterMax = %d\n", iterMax);
for (bSNR_dB = 3.0; bSNR_dB <= 3; bSNR_dB += 0.5) {
    std = pow(10, bSNR_dB / ((double)-20));
    // run DEGA
    BPRDEGA();
    bler = 0;   // BLER from union bound
    for (i = 0; i < K; i++) 
    // Q(x) = 0.5 * erfc(x / sqrt(2))
    // BER = Q(sqrt(V[0][I[i]]->l / 2))
        bler += erfc(sqrt(V[0][bRev[I[i]]]->l) / 2.0);
    bler = bler * 0.5;
    // final output
    printf("bSNR = %.2lf\tBLER = %lf\t\tBER = %lfe-2\n",
        bSNR_dB, bler, bler * 100 / K);
    // mean of left propagating message at decision stage
    /*
    printf("E[L[0][j]] for unfrozen bits:\n");
    for (j = 0; j < K; j++) 
        printf("%lf, ", V[0][bRev[I[j]]]->l);
    printf("\n");
    */
}
    return 0;
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
    int d = n - i - 1;          // in-block seperation

    if (initV[i][j] == 1) {     // already done
        V[i][j]->leftP = 0;     // it is the lower node
        return;
    }
    V[i][j]->leftP = 1;         // else, upper node
    initV[i][j] = 1;            // mark upperleft node
    initV[i][j + pow2(d)] = 1;  // mark lowerleft node
    // assign childs
    // for the upperleft bit
    V[i][j]->cU = V[i + 1][j];
    V[i][j]->cL = V[i + 1][j + pow2(d)];
    V[i + 1][j]->rightP = 1;            // set their rightP
    V[i + 1][j + pow2(d)]->rightP = 0;
    // for the lowerleft bit
    V[i][j + pow2(d)]->cU = V[i + 1][j];
    V[i][j + pow2(d)]->cL = V[i + 1][j + pow2(d)];
    // assign parents conversely
    V[i + 1][j]->pU = V[i][j];
    V[i + 1][j]->pL = V[i][j + pow2(d)];
    V[i + 1][j + pow2(d)]->pU = V[i][j];
    V[i + 1][j + pow2(d)]->pL = V[i][j + pow2(d)];
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
	// else 
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
void BPRDEGA()
{
    int i, j, k;    // looping indices
    int iter;       // iteration index
    double tempL, tempU;  // temp storage

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
        if (inI[bRev[j]] == 0)  // frozen bit => infinity
            V[0][j]->r = 999;
        else                    // info bit => 0
            V[0][j]->r = 0;   
    }
    for (iter = 0; iter < (iterMax - 4 * bSNR_dB); iter++) {
        // R propag
        for (i = 0; i < n; i++) 
            for (j = 0; j < N; j++) {
                // only start from upperleft
                if (V[i][j]->leftP == 1) {
                    V[i][j]->cU->r = phi_inv(
                    phi(V[i][j]->r)
                    + phi(V[i][j]->cL->pL->r + V[i][j]->cL->l)
                    - phi(V[i][j]->r)
                        * phi(V[i][j]->cL->pL->r + V[i][j]->cL->l)
                    );
                    V[i][j]->cL->r = phi_inv(
                    phi(V[i][j]->r) + phi(V[i][j]->cU->l)
                    - phi(V[i][j]->r) * phi(V[i][j]->cU->l)
                    ) + V[i][j]->cL->pL->r;
                }
            }
        // L propag
        for (i = n - 1; i >= 0; i--) 
            for (j = 0; j < N; j++) {
                // only start from upperleft
                if (V[i][j]->leftP == 1) {
                    V[i][j]->l = phi_inv(
                    phi(V[i][j]->cU->l)
                    + phi(V[i][j]->cL->pL->r + V[i][j]->cL->l)
                    - phi(V[i][j]->cU->l)
                        * phi(V[i][j]->cL->pL->r + V[i][j]->cL->l)
                    );
                    V[i][j]->cL->pL->l = phi_inv(
                    phi(V[i][j]->r) + phi(V[i][j]->cU->l)
                    - phi(V[i][j]->r) * phi(V[i][j]->cU->l)
                    ) + V[i][j]->cL->l;
                }
            }
        // if (iter % 2 == 1) {   // every 2 iterations
        if (iter < 10 && iter > 1) {
            printf("%2d\t", iter + 1);
            // for each stage
            for (i = 0; i <= n; i++) {
                // compute E[i], the sum of BER
                E[i] = 0;
                for (j = 0; j < N; j++) {
                    if (erfc(sqrt(V[i][j]->l + V[i][j]->r) / 2.0) < 0)
                        printf("ERR!\n");
                    E[i] += erfc(sqrt(V[i][j]->l + V[i][j]->r) / 2.0);
                }
                E[i] = E[i] * 0.5;
                // printf("E%d = %lf\t", i, E[i]);
                printf("%lf\t", E[i]);
            }
            printf("\n");
        }
    }
    return;
}