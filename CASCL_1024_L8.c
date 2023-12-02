/*
Simulation of CRC-aided SCL decoder of N = 1024, rate 0.5
with CRC given by g(D) = D^24 + D^23 + D^21 + D^20 + D^17
+ D^15 + D^13+ D^12 + D^8 + D^4 + D^2 + D + 1
Assume list size L = power of 2
Follow factor graph in Lee's thesis
G = F^{\otimes n} (no bit reversal)

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double bSNR_dB;         // Eb/N0 in dB
#define N 1024
#define K 512
#define n 10            // n = log2(N)
#define r 24            // number of bits added by CRC
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
0,1,2,4,8,16,32,3,5,64,9,6,17,10,18,128,12,33,65,20,256,34,24,36,7,129,66,512,11,40,
68,130,19,13,48,14,72,257,21,132,35,258,26,513,80,37,25,22,136,260,264,38,514,96,67,
41,144,28,69,42,516,49,74,272,160,520,288,528,192,544,70,44,131,81,50,73,15,320,133,
52,23,134,384,76,137,82,56,27,97,39,259,84,138,145,261,29,43,98,515,88,140,30,146,71,
262,265,161,576,45,100,640,51,148,46,75,266,273,517,104,162,53,193,152,77,164,768,268,
274,518,54,83,57,521,112,135,78,289,194,85,276,522,58,168,139,99,86,60,280,89,290,529,
524,196,141,101,147,176,142,530,321,31,200,90,545,292,322,532,263,149,102,105,304,296,
163,92,47,267,385,546,324,208,386,150,153,165,106,55,328,536,577,548,113,154,79,269,
108,578,224,166,519,552,195,270,641,523,275,580,291,59,169,560,114,277,156,87,197,116,
170,61,531,525,642,281,278,526,177,293,388,91,584,769,198,172,120,201,336,62,282,143,
103,178,294,93,644,202,592,323,392,297,770,107,180,151,209,284,648,94,204,298,400,608,
352,325,533,155,210,305,547,300,109,184,534,537,115,167,225,326,306,772,157,656,329,
110,117,212,171,776,330,226,549,538,387,308,216,416,271,279,158,337,550,672,118,332,
579,540,389,173,121,553,199,784,179,228,338,312,704,390,174,554,581,393,283,122,448,
353,561,203,63,340,394,527,582,556,181,295,285,232,124,205,182,643,562,286,585,299,354,
211,401,185,396,344,586,645,593,535,240,206,95,327,564,800,402,356,307,301,417,213,568,
832,588,186,646,404,227,896,594,418,302,649,771,360,539,111,331,214,309,188,449,217,
408,609,596,551,650,229,159,420,310,541,773,610,657,333,119,600,339,218,368,652,230,
391,313,450,542,334,233,555,774,175,123,658,612,341,777,220,314,424,395,673,583,355,
287,183,234,125,557,660,616,342,316,241,778,563,345,452,397,403,207,674,558,785,432,
357,187,236,664,624,587,780,705,126,242,565,398,346,456,358,405,303,569,244,595,189,
566,676,361,706,589,215,786,647,348,419,406,464,680,801,362,590,409,570,788,597,572,
219,311,708,598,601,651,421,792,802,611,602,410,231,688,653,248,369,190,364,654,659,
335,480,315,221,370,613,422,425,451,614,543,235,412,343,372,775,317,222,426,453,237,
559,833,804,712,834,661,808,779,617,604,433,720,816,836,347,897,243,662,454,318,675,
618,898,781,376,428,665,736,567,840,625,238,359,457,399,787,591,678,434,677,349,245,
458,666,620,363,127,191,782,407,436,626,571,465,681,246,707,350,599,668,790,460,249,
682,573,411,803,789,709,365,440,628,689,374,423,466,793,250,371,481,574,413,603,366,
468,655,900,805,615,684,710,429,794,252,373,605,848,690,713,632,482,806,427,904,414,
223,663,692,835,619,472,455,796,809,714,721,837,716,864,810,606,912,722,696,377,435,
817,319,621,812,484,430,838,667,488,239,378,459,622,627,437,380,818,461,496,669,679,
724,841,629,351,467,438,737,251,462,442,441,469,247,683,842,738,899,670,783,849,820,
728,928,791,367,901,630,685,844,633,711,253,691,824,902,686,740,850,375,444,470,483,
415,485,905,795,473,634,744,852,960,865,693,797,906,715,807,474,636,694,254,717,575,
913,798,811,379,697,431,607,489,866,723,486,908,718,813,476,856,839,725,698,914,752,
868,819,814,439,929,490,623,671,739,916,463,843,381,497,930,821,726,961,872,492,631,
729,700,443,741,845,920,382,822,851,730,498,880,742,445,471,635,932,687,903,825,500,
846,745,826,732,446,962,936,475,853,867,637,907,487,695,746,828,753,854,857,504,799,
255,964,909,719,477,915,638,748,944,869,491,699,754,858,478,968,383,910,815,976,870,
917,727,493,873,701,931,756,860,499,731,823,922,874,918,502,933,743,760,881,494,702,
921,501,876,847,992,447,733,827,934,882,937,963,747,505,855,924,734,829,965,938,884,
506,749,945,966,755,859,940,830,911,871,639,888,479,946,750,969,508,861,757,970,919,
875,862,758,948,977,923,972,761,877,952,495,703,935,978,883,762,503,925,878,735,993,
885,939,994,980,926,764,941,967,886,831,947,507,889,984,751,942,996,971,890,509,949,
973,1000,892,950,863,759,1008,510,979,953,763,974,954,879,981,982,927,995,765,956,887,
985,997,986,943,891,998,766,511,988,1001,951,1002,893,975,894,1009,955,1004,1010,957,
983,958,987,1012,999,1016,767,989,1003,990,1005,959,1011,1013,895,1006,1014,1017,1018,
991,1020,1007,1015,1019,1021,1022,1023};
int I[K + r];               // information set
int inI[N];                 // whether a bit is in I
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
// check if the k-th path pass cyclic redundancy check
int CRcheck(int k);
// successive cancellation list decoder
void CASCL(double *y, int *u_hat);

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
    int w[K + r];               // CRC encoder output
    int u[N];                   // encoder input
    int x[N];                   // polar codeword
    double y[N];                // codeword + Gaussian noise
    int u_hat[N];               // decoder's output
    int errbit;                 // # of bit error

    // set random seed
    SEED = ((unsigned long long)(time(NULL))) % 10000;
    printf("SEED = %ld\n", SEED);   // for debug
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
    // at first, no channels are picked
    for (i = 0; i < N; i++)
        inI[i] = 0;
    // determine the information set I
    for (i = 0; i < K + r; i++) {
        // pick most reliable channels
        I[i] = Q[N - (K + r) + i];
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
for (bSNR_dB = 1.0; bSNR_dB <= 1.5; bSNR_dB += 0.5) {
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
        // CRC encoder with
        for (i = 0; i < K + r; i++) 
            w[i] = 0;
        /* use PN sequence to generate K bits, then encode
        by g(x) = D^24 + D^23 + D^21 + D^20 + D^17 + D^15
        + D^13+ D^12 + D^8 + D^4 + D^2 + D + 1 */
        for (i = 0; i < K; i++)
            if (PN[(m + i) % 63] == 1) {
                w[i] += 1;
                w[i + 1] += 1;
                w[i + 2] += 1;
                w[i + 4] += 1;
                w[i + 8] += 1;
                w[i + 12] += 1;
                w[i + 13] += 1;
                w[i + 15] += 1;
                w[i + 17] += 1;
                w[i + 20] += 1;
                w[i + 21] += 1;
                w[i + 23] += 1;
                w[i + 24] += 1;
            }
        // Polar Encoder
        // put them into information set in u
        for (i = 0; i < K + r; i++) {
            u[I[i]] = (w[i] % 2);
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
        // SC decoder
        CASCL(y, u_hat);
        // check info bits for block error
        temp = 0;           // flag for loop
       for (i = 0; i < K + r; i++) {
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
    printf("L = %d\tbSNR = %.2lf\terror block = %d\trun = %d\tBLER = %lfe-4\n",
        L, bSNR_dB, errBlock, run, ((double)errBlock) / (run / 10000.0));
    /* printf("Error bit = %d\tBER = %lf\n", errbit,
        ((double)errbit) / (K + r) / run); */
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

// check CRC for the k-th path; return 1 for success
int CRcheck(int k)
{
    int i;          // looping index
    int C[K + r];       // CRC codeword to check

    // copy the decoded information bits to C
    for (i = 0; i < K + r; i++) 
        C[i] = V[0][I[i]]->b[k];
    /* long division by g(D) = D^24 + D^23 + D^21 + D^20
    + D^17 + D^15 + D^13+ D^12 + D^8 + D^4 + D^2 + D + 1 */
    for (i = K + r - 1; i >= r; i--) 
        if (C[i] == 1) {
            C[i] = 0;
            C[i - 1] = (C[i - 1] + 1) % 2;
            C[i - 3] = (C[i - 3] + 1) % 2;
            C[i - 4] = (C[i - 4] + 1) % 2;
            C[i - 7] = (C[i - 7] + 1) % 2;
            C[i - 9] = (C[i - 9] + 1) % 2;
            C[i - 11] = (C[i - 11] + 1) % 2;
            C[i - 12] = (C[i - 12] + 1) % 2;
            C[i - 16] = (C[i - 16] + 1) % 2;
            C[i - 20] = (C[i - 20] + 1) % 2;
            C[i - 22] = (C[i - 22] + 1) % 2;
            C[i - 23] = (C[i - 23] + 1) % 2;
            C[i - 24] = (C[i - 24] + 1) % 2;
        }
    for (i = r - 1; i >= 0; i--) 
        if (C[i] == 1) return 0;
    return 1;
}

// CRC-aided successive cancellation decoder
void CASCL(double *y, int *u_hat)
{
    int i, j, k;                // looping indices
    int actL;                   // number of active paths
    double med;                 // (L + 1)-th smallest PM
    double min;                 // minimum path metric
    int min_k;                  // the index of min
    int pass[L];                // whether a path pass CRC
    int anyPass = 0;            /* indicate if there is any
    path that pass CRC; turn to 1 when 1 path pass CRC */
    int firstCand;              // 1st path to pass CRC

    // init path metrics
    PM[0] = 0;
    // init bDone
    for (i = 0; i <= n; i++)
        for (j = 0; j < N; j++)
            // init the first decoder only
            V[i][j]->bDone[0] = 0;  // undone yet
    // propagate frozen bits
    for (j = 0; j < N; j++)
        if (inI[j] == 0) {          // frozen bit
            V[0][j]->b[0] = 0;
            updateBit(V[0][j], 0); 
        }
    // init lDone
    for (i = 0; i < n; i++)
        for (j = 0; j < N; j++) {
            V[i][j]->lDone[0] = 0;  // undone yet
        }
    // channel LLR for all paths
    for (j = 0; j < N; j++) 
        for (k = 0; k < L; k++) {
            V[n][j]->l[k] = 2 * y[j] / std / std;
            V[n][j]->lDone[k] = 1;
        }
    // Decoding: LLR propagation: top-down computation
    // before the list is full
    actL = 1;
    for (j = 0; j < N && actL < L; j++) {
        // each decoder compute LLR
        for (k = 0; k < actL; k++) 
            getLLR(V[0][j], k);
        if (inI[j] == 1) {              // if is info bit
            for (k = 0; k < actL; k++)
                copyPath(k, k + actL);  // copy path
            for (k = 0; k < actL; k++) {
                // branch 0 and 1
                V[0][j]->b[k] = 0;      
                V[0][j]->b[k + actL] = 1;
                // update path metric
                PM[k + actL] = PM[k] + PHI(k, j, 1);
                PM[k] = PM[k] + PHI(k, j, 0);
                // propagate bit value
                updateBit(V[0][j], k);
                updateBit(V[0][j], k + actL);
            }
            actL = actL * 2;    // double the # of paths
        } else {                // only need to update PM
            for (k = 0; k < actL; k++)
                PM[k] += PHI(k, j, 0);
        }
    }
    for (; j < N; j++) {
        // each decoder compute LLR
        for (k = 0; k < L; k++) 
            getLLR(V[0][j], k);
        if (inI[j] == 1) {              // if is info bit
            // compute all possible path metrics
            for (k = 0; k < L; k++) {
                PMcand[k] = PM[k] + PHI(k, j, 0);
                PMcand[k + L] = PM[k] + PHI(k, j, 1);
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
                    V[0][j]->b[k] = 0;
                    updateBit(V[0][j], k);
                    break;
                case 1:
                    V[0][j]->b[k] = 1;
                    updateBit(V[0][j], k);
                    PM[k] = PM[k + L];
                    break;
                case 2:
                    // find space unoccupied
                    for (; surviv[i] != -1; i++);
                    if (i >= L) printf("Error!\n");
                    // put branch1 to i
                    simpleCopy(k, i);
                    V[0][j]->b[k] = 0;
                    updateBit(V[0][j], k);
                    V[0][j]->b[i] = 1;
                    updateBit(V[0][j], i);
                    surviv[i] = -2;     // mark as copied
                    PM[i] = PM[k + L];
                }
            }
        } else {    // frozen bit => only need to update PM
            for (k = 0; k < L; k++)
                PM[k] += PHI(k, j, 0);
        }
    }
    // Pick only those that pass CRC
    firstCand = -1;
    for (k = 0; k < L; k++) {
        pass[k] = CRcheck(k);
        if (anyPass == 0 && pass[k] == 1) {
            anyPass = 1;
            firstCand = k;
        }
    }
    // Find the legal path with smallest PM
    if (anyPass == 1) {
        min = PM[firstCand];
        min_k = firstCand;
        for (k = 1; k < L; k++)
            if (pass[k] == 1) {
                if (PM[k] < min) {
                    min = PM[k];
                    min_k = k;
                }
            }
        
    } else {
    // No legal path => find the path with smallest PM
        min = PM[0];
        min_k = 0;
        for (k = 1; k < L; k++)
            if (PM[k] < min) {
                min = PM[k];
                min_k = k;
            }
    }
    // decoder output
    for (j = 0; j < N; j++) {
        u_hat[j] = V[0][j]->b[min_k];
    }
    return;
}