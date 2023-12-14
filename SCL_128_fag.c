/*
Simulation of successive cancellation list (SCL) decoder
for polar code with N = 128 and rate = 0.5
Assume list size L = power of 2
Follow factor graph in Kao's thesis
G = F^{\otimes n} (no bit reversal)
Perform bit reversal at decoder side
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double bSNR_dB;         // Eb/N0 in dB
#define N 128
#define K 64
#define n 7             // n = log2(N)
#define L 8             // list size
#define BLE 200         // desired # of block error

typedef struct node {
    double l[L];    // LLR for all branches
    int b[L];       // bit value
    int lDone[L];   // whether LLR is computed
    int bDone[L];   // whether bit value decided
    /* Left & right position, i.e., when it appears at
    the left/right side of a BCB, is it the upper or
    the lower node.
    1 for upper, 0 for lower */
    int leftP;
    int rightP;
    struct node *pU;        // upper parent
    struct node *pL;        // lower parent
    struct node *cU;        // upper child
    struct node *cL;        // lower child
} node;

// seed for generating random number in (0, 1)
unsigned long long SEED; 
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
int bRev[N];                // bit reversal for [0, 127]
// n-kronecker product of F
int **Fn; 
// all nodes on factor graph
node ***V;
// path metric array
double *PM;
double *PMcand;             // the 2L's PM candidates
int surviv[L];              // L's survivor out of 2L
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
// recursively compute LLR
void getLLR(node *v, int k);
// after v is decoded, update other nodes
void updateBit(node *v, int k);
// copy a factor graph from path c to path k
void copyPath(int c, int k);
// copy factor graph from c to k w/out bDone or lDone 
void simpleCopy(int c, int k);
/* phi function for Path metric update
k: path, j: bit index, u: bit value */
double PHI(int k, int j, int u);
// Quick sort: from small to large
void QuickSort(double *list, int low, int high);
// for quick sort: put the 1st element in place
// return the position of the element
int Partition(double *list, int low, int high);
// successive cancellation list decoder
void SCLdecode(double *y, int *u_hat);

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
    int errBlock;               // number of block errors
    int temp;                   // temporary storage
    int u[N];                   // encoder input
    int x[N];                   // polar codeword
    double y[N];                // codeword + Gaussian noise
    int u_hat[N];               // decoder's output
    int errbit;                 // # of bit error

    // set random seed
    SEED = ((unsigned long long)time(NULL)) % 10000;
    printf("SEED = %ld\n", SEED);
    // allocate memory for PM
    PM = (double *)calloc(2 * L, sizeof(double));
    PMcand = (double *)calloc(2 * L, sizeof(double));
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
            printf("Bit reversal wrong!\n");
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
for (bSNR_dB = 1.0; bSNR_dB <= 3.0; bSNR_dB += 0.5) {
    errBlock = 0;
    errbit = 0;
    std = pow(10, bSNR_dB / ((double)-20));
    // run simulation until desired error blocks
    for (run = 0; errBlock < BLE; run++) {
        // reset vectors to all-zero
        for (i = 0; i < N; i++) {
            u_hat[i] = 0;
            x[i] = 0;
        }
        // Encoder
        /* use PN sequence to have K bits, and
        put them into information set in u */
        for (i = 0; i < K; i++) 
            u[I[i]] = PN[(m + i) % 63];
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
        // SCL decoder
        SCLdecode(y, u_hat);
        // check info bits for block error
        temp = 0;           // flag for loop
        for (i = 0; i < K; i++)
            if (u[I[i]] != u_hat[I[i]]) {
                temp = 1;
                errbit += 1;
            }
        errBlock += temp;
        m += step_m;                   // increase m
        if (m >= 63) m -= 63;
    }
    // final output
    printf("L = %d\tbSNR = %.2lf\terror block = %d\trun = %d\tBLER = %lf * 10^-1\n",
        L, bSNR_dB, errBlock, run, ((double)errBlock) * 10 / run);
    printf("Error bit = %d\tBER = %lf\n", errbit,
        ((double)errbit) / K / run);
}
    // debug
    /*
    printf("u\tuh\tI\n");
    for (i = 0; i < N; i++) {
        printf("%d\t%d\t", u[i], u_hat[i]);
        if (inI[i] == 1) printf("1");
        printf("\n");
    }
    */
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

// recursively compute LLR at the k-th path (decoder)
void getLLR(node *v, int k)
{
    if (v->lDone[k] == 1) return;
    getLLR(v->cU, k);
    getLLR(v->cL, k);
    // if it is upperleft bit
    if (v->leftP == 1)  
        v->l[k] = CHK(v->cU->l[k], v->cL->l[k]);
    else if (v->cU->pU->bDone[k] == 1) {
        if (v->cU->pU->b[k] == 0) 
            v->l[k] = v->cL->l[k] + v->cU->l[k];
        else
            v->l[k] = v->cL->l[k] - v->cU->l[k];
    } else
        printf("Wrong propagation order!\n");
    v->lDone[k] = 1;
    return;
}

// Set bDone to 1 and update descendants at the k-th path
void updateBit(node *v, int k)
{
    if (v->bDone[k] == 1) return;
    v->bDone[k] = 1;
    if (v->cU == NULL || v->cL == NULL)
        return;     // return when we reach the rightmost
    // if it is upperleft bit
    if (v->leftP == 1) {
        // check if the lowerleft bit is done
        if (v->cU->pL->bDone[k] == 1) {
            v->cU->b[k] = (v->b[k] + v->cU->pL->b[k]) % 2;
            updateBit(v->cU, k);
        }
    } else {        // if it is the lowerleft bit
        // check the upperleft node
        if (v->cU->pU->bDone[k] == 1) {
            v->cU->b[k] = (v->b[k] + v->cU->pU->b[k]) % 2;
            updateBit(v->cU, k);
        }
        // update the lowerright bit
        v->cL->b[k] = v->b[k];
        updateBit(v->cL, k);
    }
    return;
}

// copy a factor graph from path c to path k
void copyPath(int c, int k)
{
    int i, j;           // looping indices

    // don't need to copy the rightmost stage
    for (i = 0; i < n; i++)
        for (j = 0; j < N; j++) {
            V[i][j]->l[k] = V[i][j]->l[c];
            V[i][j]->lDone[k] = V[i][j]->lDone[c];
            V[i][j]->b[k] = V[i][j]->b[c];
            V[i][j]->bDone[k] = V[i][j]->bDone[c];
        }
    return;
}

// copy factor graph from c to k w/out bDone or lDone
void simpleCopy(int c, int k)
{
    int i, j;           // looping indices

    // don't need to copy the rightmost stage
    for (i = 0; i < n; i++)
        for (j = 0; j < N; j++) {
            V[i][j]->l[k] = V[i][j]->l[c];
            V[i][j]->b[k] = V[i][j]->b[c];
        }
    return;
}

// updated Path metric for path k with the j-th bit := u
double PHI(int k, int j, int u)
{   
    double result;      // the result to return
    double absL = fabs(V[0][j]->l[k]);

    if (u != 0 && u != 1)
        printf("Illegal input u!\n");
    // table look-up for the term ln(1 + e^-x)
    if (absL < 0.196) result = 0.65;
    else if (absL < 0.433)  result = 0.55;
    else if (absL < 0.71)   result = 0.45;
    else if (absL < 1.05)   result = 0.35;
    else if (absL < 1.508)  result = 0.25;
    else if (absL < 2.252)  result = 0.15;
    else if (absL < 4.5)    result = 0.05;
    else result = 0;
    if ((u == 0 && V[0][j]->l[k] < 0)
     || (u == 1 && V[0][j]->l[k] > 0))
        // result += -(1 - 2 * u) * V[0][j]->l[k];
        result += absL;
    return result;
}

// quick sort
void QuickSort(double *list, int low, int high)
{
	int mid;

	if (low < high) {
		// put the first element in place and let mid be that index 
		mid = Partition(list, low, high);	
		QuickSort(list, low, mid - 1);		// sort LHS of mid
		QuickSort(list, mid + 1, high);		// sort RHS of mid
	}
}
// for quick sort: put the 1st element in place
// return the position of the element
int Partition(double *list, int low, int high)
{
	int i, j;			// looping indices
	double temp, v;	    // temp storage; the first element

	v = list[low];		// the first element
	i = low + 1;		// from left hand side
	j = high;			// from right hand side
	do {
		while (list[i] < v && i < high) 
		// until we find a list[i] > v
			i += 1;
		while (list[j] > v) 
		// until we find a list[j] < v
			j -= 1;	
		// list[j] < v < list[i] 
		if (i < j) {		// swap
			temp = list[i];
			list[i] = list[j];
			list[j] = temp;
		}
	} while (i < j);
	// move v to the correct position
	list[low] = list[j];
	list[j] = v;		
	return j;
}

// successive cancellation decoder
void SCLdecode(double *y, int *u_hat)
{
    int i, j, k;                // looping indices
    int r;                      // bit-reversed index
    int actL;                   // number of active paths
    double med;                 // (L + 1)-th smallest PM
    double min;                 // minimum path metric
    int min_k;                  // the index of min

    // init path metrics
    PM[0] = 0;
    // init bDone
    for (i = 0; i <= n; i++)
        for (j = 0; j < N; j++)
            // init the first decoder only
            V[i][j]->bDone[0] = 0;  // undone yet
    // propagate frozen bits
    for (j = 0; j < N; j++)
        if (inI[bRev[j]] == 0) {    // frozen bit
            V[0][j]->b[0] = 0;
            updateBit(V[0][j], 0); 
        }
    // init lDone
    for (i = 0; i < n; i++)
        for (j = 0; j < N; j++) {
            V[i][j]->lDone[0] = 0;  // undone yet
        }
    // channel LLR for all paths
    for (j = 0; j < N; j++) {
        // only need computation in the first path
        V[n][j]->l[0] = 2 * y[bRev[j]] / std / std;
        V[n][j]->lDone[0] = 1;
        for (k = 1; k < L; k++) {
            V[n][j]->l[k] = V[n][j]->l[0];
            V[n][j]->lDone[k] = 1;
        }
    }
    // Decoding: LLR propagation: top-down computation
    // before the list is full
    actL = 1;
    for (j = 0; j < N && actL < L; j++) {
        r = bRev[j];        // bit-reversing shuffle
        // each decoder compute LLR
        for (k = 0; k < actL; k++) 
            getLLR(V[0][r], k);
        if (inI[j] == 1) {              // if is info bit
            for (k = 0; k < actL; k++)
                copyPath(k, k + actL);  // copy path
            for (k = 0; k < actL; k++) {
                // branch 0 and 1
                V[0][r]->b[k] = 0;      
                V[0][r]->b[k + actL] = 1;
                // update path metric
                PM[k + actL] = PM[k] + PHI(k, r, 1);
                PM[k] = PM[k] + PHI(k, r, 0);
                // propagate bit value
                updateBit(V[0][r], k);
                updateBit(V[0][r], k + actL);
            }
            actL = actL * 2;    // double the # of paths
        } else {                // only need to update PM
            for (k = 0; k < actL; k++)
                PM[k] += PHI(k, r, 0);
        }
    }
    for (; j < N; j++) {
        r = bRev[j];        // bit-reversing shuffle
        // each decoder compute LLR
        for (k = 0; k < L; k++) 
            getLLR(V[0][r], k);
        if (inI[j] == 1) {              // if is info bit
            // compute all possible path metrics
            for (k = 0; k < L; k++) {
                PMcand[k] = PM[k] + PHI(k, r, 0);
                PMcand[k + L] = PM[k] + PHI(k, r, 1);
                PM[k] = PMcand[k];
                PM[k + L] = PMcand[k + L];
            }
            // sort PMcand from small to large
            QuickSort(PMcand, 0, 2 * L - 1);
            med = PMcand[L];
            if (PMcand[L - 1] == med)
                printf("Oops!\n");      // for debug
            // we only want PMcand[0] to PMcand[L - 1]
            for (k = 0; k < L; k++) {
                if (PM[k] < med && PM[k + L] < med)
                    surviv[k] = 2;  // both branches survive
                else if (PM[k] >= med && PM[k + L] < med)
                    surviv[k] = 1;  // only branch1 survive
                else if (PM[k] < med && PM[k + L] >= med)
                    surviv[k] = 0;  // only branch0 survive
                else 
                    surviv[k] = -1; // neither survive
            }
            /* put both-branch survivor to
            where none-branch survivor used to be */
            i = 0;          // for unoccupied space search
            for (k = 0; k < L; k++) {
                switch (surviv[k]) {
                case 0:
                    V[0][r]->b[k] = 0;
                    updateBit(V[0][r], k);
                    break;
                case 1:
                    V[0][r]->b[k] = 1;
                    updateBit(V[0][r], k);
                    PM[k] = PM[k + L];
                    break;
                case 2:
                    // find space unoccupied
                    for (; surviv[i] != -1; i++);
                    if (i >= L) printf("Error!\n");
                    // put branch1 to i
                    simpleCopy(k, i);
                    V[0][r]->b[k] = 0;
                    updateBit(V[0][r], k);
                    V[0][r]->b[i] = 1;
                    updateBit(V[0][r], i);
                    surviv[i] = -2;     // mark as copied
                    PM[i] = PM[k + L];
                }
            }
        } else {    // frozen bit => only need to update PM
            for (k = 0; k < L; k++)
                PM[k] += PHI(k, r, 0);
        }
    }
    // Decoder terminate: find the path with smallest PM
    min = PM[0];
    min_k = 0;
    for (k = 1; k < L; k++)
        if (PM[k] < min) {
            min = PM[k];
            min_k = k;
        }
    // decoder output
    for (j = 0; j < N; j++) {
        u_hat[j] = V[0][bRev[j]]->b[min_k];
    }
    return;
}