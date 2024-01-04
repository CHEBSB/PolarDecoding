/*
Gaussian Approximation of density evolution
for BPR decoding of polar codes where N = 1024, K = 512
Assume all zero codewords transmitted
Make use of matrix to propagate LLR
and compute the unfrozen BER at stage 0
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double bSNR_dB;             // Eb/N0 in dB
#define N 1024
#define K 512
#define n 10
#define iterMax 40         // not really max iteration

typedef struct node {
    double l;           // left-propagating message (LLR)
    double r;           // right-propagating message
    /* Left & right position, i.e., when it appears at
    the left/right side of a BCB, is it the upper or
    the lower node.
    1 for upper, 0 for lower
    */
    double p;           // error probability of this bit
    double W;           // ln(Pr{v = 0} / Pr{v = 1})
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
int I[K];                   // information set
int inI[N];                 // whether a bit is in I
int bRev[N];                // bit reversal for [0, 127]
// M1[k][i]: positions of 1 in the i-th column of M_k
int ***M1; 
// weight of columns of M
int Mw[n + 1][N];
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


    // allocate memory for M1
    M1 = (int ***)calloc(n + 1, sizeof(int **));
    M1[0] = NULL;       // don't use the 0-th
    for (k = 1; k <= n; k++)        // each M has N columns
        M1[k] = (int **)calloc(N, sizeof(int *));
    // read Mw from file to decide memory use
    for (k = 1; k <= n; k++) 
        for (i = 0; i < N; i++) {
            scanf("%d", &temp);
            Mw[k][i] = temp;
            M1[k][i] = (int *)calloc(temp, sizeof(int));
        }
    for (i = 0; i < N; i++)         // set the 0-th to -1
        Mw[0][i] = -1;
    // read M from file
    for (k = 1; k <= n; k++) {
    // for stage k
        for (i = 0; i < N; i++)
            for (j = 0; j < Mw[k][i]; j++) {
                scanf("%d", &temp);
                M1[k][i][j] = temp;
            }
    }
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
for (bSNR_dB = 2.0; bSNR_dB <= 4; bSNR_dB += 0.5) {
    printf("bSNR = %.2lf\t", bSNR_dB);
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
    printf("BLER = %e\t\tBER = %e\n",
        bler, bler / K);
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
    int s;          // bit-reversed info index
    int iter;       // iteration index
    double tempL;   // temp storage

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
            V[0][j]->r = 99;
        else                    // info bit => 0
            V[0][j]->r = 0;   
    }
    printf("iterMax = %d\n", (int)(iterMax / bSNR_dB));
    for (iter = 0; iter < (int)(iterMax / bSNR_dB); iter++) {
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
        if (iter > 3) {
            printf("%2d  ", iter + 1);
            // for each stage
        // for i = 0
            // compute E[i], the approx BLER
            E[i] = 0;
            for (j = 0; j < K; j++) 
                E[i] += erfc(sqrt(V[0][bRev[I[j]]]->l) / 2.0);
            E[i] = E[i] * 0.5;
            printf("%.5e  ", E[i]);
        // for the rest: need propagation
            for (i = 1; i <= n; i++) {
                // compute the LLR for bit value
                for (j = 0; j < N; j++)
                    V[i][j]->W = V[i][j]->l + V[i][j]->r;
                // compute LLR of the j-th info bit
                for (j = 0; j < K; j++) {
                    tempL = V[i][M1[i][bRev[I[j]]][0]]->W;
                    for (k = 1; k < Mw[i][bRev[I[j]]]; k++)
                        tempL = CHK(tempL,
                        V[i][M1[i][bRev[I[j]]][k]]->W);
                    // convert LLR back to error prob
                    V[0][bRev[I[j]]]->p = erfc(sqrt(tempL) / 2.0);
                }
                // compute E[i], the approx BLER
                E[i] = 0;
                for (j = 0; j < K; j++) 
                    E[i] += V[0][bRev[I[j]]]->p;
                // printf("E%d = %lf\t", i, E[i]);
                printf("%.5e  ", E[i]);
            }
            printf("\n");
        }
    }
    return;
}