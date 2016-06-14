/*
    * Title: this is a code improved from Robert Morelos-Zaragoza ,in this code we implement two algorithm bm and ibm  to finish bch decode
    * author : Leger Liu
    * data :2016/6/14
    * address : Wuhan National Laboratory for Optoelectronics, Huazhong University of Science&Technology.
    * email :liulengend@gmail.com

     * m = order of the Galois field GF(2**m)
     * n = 2**m - 1 = size of the multiplicative group of GF(2**m)
     * length = length of the BCH code
     * t = error correcting capability (max. no. of errors the code corrects)
     * d = 2*t + 1 = designed min. distance = no. of consecutive roots of g(x) + 1
     * k = n - deg(g(x)) = dimension (no. of information bits/codeword) of the code
     * p[] = coefficients of a primitive polynomial used to generate GF(2**m)
     * g[] = coefficients of the generator polynomial, g(x)
     * alpha_to [] = log table of GF(2**m)
     * index_of[] = antilog table of GF(2**m)
     * data[] = information bits = coefficients of data polynomial, i(x)
     * bb[] = coefficients of redundancy polynomial x^(length-k) i(x) modulo g(x)
     * numerr = number of errors
     * errpos[] = error positions
     * recd[] = coefficients of the received polynomial
     * decerror = number of decoding errors (in _message_ positions)
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <windows.h>

#define TEST_COUNT 200
//#define TIME_TEST 1
#define IBM 1 // by define IBM ,you can use IBM algorithm to complete the bch decode
//#define DEBUG
#define HARDWARE_TEST

#define DIMENSION 14//14 //11// 11
#define GFSIZE  16383 ////2047 //


#define MINPOLYNOMIAL {1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1}// {1,0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1} //     ////本原多项式的系数 0->m 共15位 表2-7
#define INFOSIZE 8208//  1024 // //    //有效数据位
#define CORRECTCAPACITY 40// 4  //纠错能力
#define REDUNDANCYSIZE 560// 44//////  / //冗余位
#define CODESIZE (INFOSIZE+REDUNDANCYSIZE)  //编码后的长度

void decode_bch_ibm();




int             m = DIMENSION, n = GFSIZE, k = INFOSIZE, length = CODESIZE;
int             p[DIMENSION+1] = MINPOLYNOMIAL;
int             t = CORRECTCAPACITY;
int             d = (2*CORRECTCAPACITY+1);
int             g[REDUNDANCYSIZE+1] =
{1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1,
0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1,
0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0,
0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0,
0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1,
1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1,
1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0,
1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0,
0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1};
/*
{1,0,0,0,0,
0,1,1,1,1,
0,0,0,0,0,
1,0,1,1,1,
0,0,1,1,1,
1,1,1,0,0,
0,0,1,0,0,
1,0,1,0,0,
1,1,0,0,1};*///共45个数据，生成多项式个g（x）

// GF(2**m)
int             alpha_to[GFSIZE+1], index_of[GFSIZE+1];
int             data[INFOSIZE];
int             bb[REDUNDANCYSIZE], recd[CODESIZE];









/*
 * generate GF(2**m) from the irreducible polynomial p(X) in p[0]..p[m]
 * lookup tables:  index->polynomial form   alpha_to[] contains j=alpha**i;
 * polynomial form -> index form  index_of[j=alpha**i] = i alpha=2 is the
 * primitive element of GF(2**m)
 *计算伽罗伐域上的各个数据存在 alpha_to中
 */
void generate_gf()//alpha_to[i]用向量表示
{
	register int    i, mask;
	mask = 1;
	alpha_to[m] = 0;
	for (i = 0; i < m; i++)
		{
		alpha_to[i] = mask;
		index_of[alpha_to[i]] = i;
		if (p[i] != 0)
			alpha_to[m] ^= mask;//mask始终只会有1位为1，在p[i]=1时将该位的1加到alpha的向量中既可以计算出邋alpha_to[m]的向量值
		mask <<= 1;
		}
	index_of[alpha_to[m]] = m;
	mask >>= 1;
	for (i = m + 1; i < n; i++) {
		if (alpha_to[i - 1] >= mask)
		  alpha_to[i] = alpha_to[m] ^ ((alpha_to[i - 1] ^ mask) << 1);
		//alpha_to[i] = alpha_to[m] ^ ((alpha_to[i - 1] ^ mask) << 1);mask[10]=1,移出alpha[i-1]的第m-1位的1，并左移移位，相当于*alpha,同时将移出的最高位alpha[11]加上alpha_to[m=11]，以此完成alpha[i]的计算
		else
		  alpha_to[i] = alpha_to[i - 1] << 1;//直接左移1位，相当于*alpha。
		index_of[alpha_to[i]] = i;//所有向量的指数值。
	}
	index_of[0] = -1;//此时alpha为全0，用其指数为-1表示。
}



/*
*对data中的数据进行编码，得到冗余位，共560位冗余位
*/
void encode_bch()
{
	register int    i, j;
	register int    feedback;


	for (i = 0; i < REDUNDANCYSIZE; i++) //把冗余位bb置为0，冗余位等于余数的位数，共44位、
		bb[i] = 0;

	for (i = k - 1; i >= 0; i--) {
		feedback = data[i] ^ bb[REDUNDANCYSIZE - 1];//输入数据与余数的最高位得到输入电路的数据
		if (feedback != 0) {

			for (j = REDUNDANCYSIZE - 1; j > 0; j--)//生成多项式，1~43，在j=44时，恒连接
				if (g[j] != 0)//由生成多项式决定是否有反馈连接，feedback位输入的数据
					bb[j] = bb[j - 1] ^ feedback;
				else
					bb[j] = bb[j - 1];
			bb[0] = g[0] && feedback;//g[0]   恒等于1，即一定有反馈连接，feedback由data[i]与余数最高位异或得到，相当于输入的数据，
		} else {
			for (j = REDUNDANCYSIZE - 1; j > 0; j--)//当输入电路的数据=0时，为移位寄存器。
				bb[j] = bb[j - 1];
			bb[0] = 0;
		}
	}
}


 /*
 * Simon Rockliff's implementation of Berlekamp's algorithm.
 *
 * Assume we have received bits in recd[i], i=0..(n-1).
 *
 * Compute the 2*t syndromes by substituting alpha^i into rec(X) and
 * evaluating, storing the syndromes in s[i], i=1..2t (leave s[0] zero) .
 * Then we use the Berlekamp algorithm to find the error location polynomial
 * elp[i].
 *
 * If the degree of the elp is >t, then we cannot correct all the errors, and
 * we have detected an uncorrectable error pattern. We output the information
 * bits uncorrected.
 *
 * If the degree of elp is <=t, we substitute alpha^i , i=1..n into the elp
 * to get the roots, hence the inverse roots, the error location numbers.
 * This step is usually called "Chien's search".
 *
 * If the number of errors located is not equal the degree of the elp, then//如果错位位置数。。。。。看英文好了  不好翻译。
 * the decoder assumes that there are more than t errors and cannot correct
 * them, only detect them. We output the information bits uncorrected.
 */

void decode_bch()//未使用余数计算的方法来计算。应该手动计算其最小多项式，以便求得余数多项式。
{

	register int    i, j, u, q, t2, count = 0, syn_error = 0;

	int             elp[100][100], d[100], l[100], u_lu[100], s[100];
	int             root[200], loc[200], err[100], reg[201];


	/* first form the syndromes */

	//add by legery liu
	//register int k;
	//int temp[44]={0};
	//for(k=0;k<44;k++)
	//{
	//	temp[k]=recd[k];
	//	recd[k]=recd[44+k];
	//	recd[1067-k]=temp[k];
	//}
	//	recd[44]=1;
	//	recd[45]=1;
	//add by legery liu
		t2 = 2 * t;

//		printf("S(x) =\n ");
	for (i = 1; i <= t2; i++) {
		s[i] = 0;
		for (j = 0; j < length; j++)
			if (recd[j] != 0)//高位的数据在后。alpha 更高。
				s[i] ^= alpha_to[(i * j) % n];//(i * j) % n，GF(2048)是一个域，当i*j超过n时，alpha_to[i*j]=alpha_to[(i*j)%n].

		if (s[i] != 0)
			syn_error = 1; /* set error flag if non-zero syndrome 表示有错误*/
	/*
	* Note: If the code is used only for ERROR DETECTION, then
	*       exit program here indicating the presence of errors.
	*/
	/* convert syndrome from polynomial form to index form  */
		//printf("alpha S in 8=%3o ", s[i]);
		//printf("%3x ", s[i]);
		//printf("%3d ",s[i]);
		//printf("%3d",index_of[alpha_to[342]^alpha_to[717]^alpha_to[1573]^alpha_to[1426]]);

		s[i] = index_of[s[i]];//alpha_to是向量形式，也就是多项式形式，此前的s[i]也是向量形式，如今将S[i]转换成alpha的多少次方的指数形式。
		#ifdef DEBUG
//		printf("in 8 the s[%d] = %o\n",i,s[i]);
		printf("in 10 the s[%d] = %d\n",i,s[i]);
		#endif // DEBUG


	}
	/*printf("%o",alpha_to[704]^alpha_to[2043]);
	printf("\n");
	printf("sigma_2=%d\n",index_of[154]);
	printf("%o\n",alpha_to[980]);
	printf("%o\n",alpha_to[1960]);
	printf("%o\n",alpha_to[(980*3)%2047]);
	printf("%o\n",alpha_to[(980*4)%2047]);

	printf("sigma_test %d\n",index_of[159]);
	printf("sigma_test %d\n",index_of[393]);
	printf("sigma_test %d\n",index_of[381]);
	printf("sigma_test %d\n",index_of[357]);
	printf("sigma_test %d\n",index_of[1181]);

	printf("sigma_1*980=%o\n",alpha_to[(196+980)%2047]);
	printf("sigma_2*980=%o\n",alpha_to[(1446+980*2+0)%2047]);
	printf("sigma_3*980=%o\n",alpha_to[(1092+980*3+0)%2047]);
	printf("sigma_4*980=%o\n",alpha_to[(1690+980*4+0)%2047]);

	printf("d3=%o\n",alpha_to[(1056+936)%2047]^alpha_to[(1061+2016)%2047]^alpha_to[(1419+1188)%2047]^alpha_to[(226+817)%2047]);
	printf("d3_1=%o\n",alpha_to[(index_of[1648]+index_of[1350])]);*/



	/*S(x) = 2008 1969 1831 1891 521 1615 1240 1735
sigma(x) =   0 2008 569 1158
Roots: 1067  46  45
Succesful decoding
请按任意键继续. . .*/

	 	/*	  8  10    index
      S(x) = 751 489  2008
			2303 1219  1969
			2153 1131  1831
			1042  546  1891
			206	 134    521
			3112 1610  1615
			2646 1446  1240
			3204 1668  1735

sigma(x) =   0 2008 569 1158
Roots: 1067  46  45
Succesful decoding
请按任意键继续. . .*/
    printf("starting the BM decoding!\n");
	if (syn_error) {
		/*
		 * 如果有错误就需要求出错误的位置，Berlekamp迭代算法来做
		 * if there are errors, try to correct them
		 * Compute the error location polynomial via the Berlekamp
		 * iterative algorithm. Following the terminology of Lin and
		 * Costello's book :   d[u] is the 'mu'th discrepancy, where
		 * u='mu'+1 and 'mu' (the Greek letter!) is the step number
		 * ranging from -1 to 2*t (see L&C),  l[u] is the degree of// l[u]是错误位置多项式的幂指数。
		 * the elp at that step, and u_l[u] is the difference between
		 * the step number and the degree of the elp. // u_l[u]是计算步数和错误位置多项式最高次数幂指数的差值。
		 */
		/* initialise table entries *///初始化表格。
		d[0] = 0;			/* index form */
		d[1] = s[1];		/* index form */
		elp[0][0] = 0;		/* index form */
		elp[1][0] = 1;		/* polynomial form */
		for (i = 1; i < t2; i++) {
			elp[0][i] = -1;	/* index form */
			elp[1][i] = 0;	/* polynomial form */
		}
		l[0] = 0;
		l[1] = 0;
		u_lu[0] = -1;
		u_lu[1] = 0;
		u = 0;

		do {
			u++;//++在u后面，循环执行完一次之后再使u加1.
//			printf("u = %d",u);
			if (d[u] == -1) {//d[u]为指数形式，在d[u]=0时代表，第u步差值为0，
				l[u + 1] = l[u];//elp最高指数不变。
				for (i = 0; i <= l[u]; i++) {
					elp[u + 1][i] = elp[u][i];
					elp[u][i] = index_of[elp[u][i]];//计算完elp[u+1]之后将elp[u]转换为指数形式存储。
				}
			} else
				/*
				 * search for words with greatest u_lu[q] for
				 * which d[q]!=0
				 */
			{
				q = u - 1;
				while ((d[q] == -1) && (q > 0))
					q--;
				/* have found first non-zero d[q]  *///向前找的次数最低时！
				if (q > 0) {
				  j = q;
				  do {
				    j--;
				    if ((d[j] != -1) && (u_lu[q] < u_lu[j]))//找遍步数小于q的所有步数，求解第p步差值!=0.而且u_lu[p]具有最大值。
				      q = j;
				  } while (j > 0);
				}

				/*
				 * have now found q such that d[u]!=0 and
				 * u_lu[q] is maximum
				 */
				/* store degree of new elp polynomial */
				if (l[u] > l[q] + u - q)//求解第u步差值的最大幂指数。
					l[u + 1] = l[u];
				else
					l[u + 1] = l[q] + u - q;

				/* form new elp(x) */
				for (i = 0; i < t2; i++)
					elp[u + 1][i] = 0;

				for (i = 0; i <= l[q]; i++)
					if (elp[q][i] != -1)
						elp[u + 1][i + u - q] =
                                   alpha_to[(d[u] + n - d[q] + elp[q][i]) % n];
				for (i = 0; i <= l[u]; i++) {
					elp[u + 1][i] ^= elp[u][i];//加上第u步的错误位置多项式的系数。
					elp[u][i] = index_of[elp[u][i]];//再将第u步的错误位置多项式转换成指数形式。
				}
			}
			u_lu[u + 1] = u - l[u + 1];

			/* form (u+1)th discrepancy */
			if (u < t2) {
			/* no discrepancy computed on last iteration */
			  if (s[u + 1] != -1)
			    d[u + 1] = alpha_to[s[u + 1]];//将s[u+1]转换成多项式形式的值赋给d[u+1]
			  else
			    d[u + 1] = 0;//否则，d[u+1]的值为0，因为没有alpha_to[-1]的数组表示形式，所以需要单独列出。
			    for (i = 1; i <= l[u + 1]; i++)
			      if ((s[u + 1 - i] != -1) && (elp[u + 1][i] != 0))
			        d[u + 1] ^= alpha_to[(s[u + 1 - i] //d[u+1]=d[u+1]^alpha_to[(s[u+1-i]+indes_of[elp[u+1][i]%n)];   get   没错。两个数相乘等于其指数相加。
			                      + index_of[elp[u + 1][i]]) % n];
			  /* put d[u+1] into index form */
			  d[u + 1] = index_of[d[u + 1]];
			}
        #ifdef DEBUG
        if(u <= 22 )
        {


            printf("in 10 the index of d[%d] = %d \n",u,d[u]);
            for(i = 0;i <= l[u] ;i++)
            {
                printf("in 10 the index of elp[%d][%d] = %d \n",u,i,elp[u][i]);
            }
        }
        #endif // DEBUG
		} while ((u < t2) && (l[u + 1] <= t));

		u++;




		if (l[u] <= t) {/* Can correct errors */
			/* put elp into index form */
			for (i = 0; i <= l[u]; i++)
				elp[u][i] = index_of[elp[u][i]];//执行完的最后一个错误位置多项式，将其转换为指数形式。
            #ifdef DEBUG
                for(i = 0;i <= l[u] ;i++)
                {
                    printf("in 10 the index of elp[%d][%d] = %d \n",u,i,elp[u][i]);
                }
            #endif // DEBUG
			/*printf("du in 8 =\n");
			printf("alpha=%4o index=%d\n ",alpha_to[d[0]],d[0]);
			printf("alpha=%4o index=%d \n",alpha_to[d[1]],d[1]);
			printf("alpha=%4o index=%d \n",alpha_to[d[2]],d[2]);
			printf("alpha=%4o index=%d \n",alpha_to[d[3]],d[3]);
			printf("alpha=%4o index=%d \n",alpha_to[d[4]],d[4]);
			printf("alpha=%4o index=%d \n",alpha_to[d[5]],d[5]);
			printf("alpha=%4o index=%d \n",alpha_to[d[6]],d[6]);
			printf("alpha=%4o index=%d \n",alpha_to[d[7]],d[7]);
			printf("alpha=%4o index=%d \n",alpha_to[d[8]],d[8]);*/

//			printf("sigma(x)in 8  =\n ");
			/*for (i = 0; i <= l[u]; i++)
			{
				//printf("%3d ",elp[u][i]);
			//	printf(" u=1 alpha=%d",alpha_to[elp[1][i]]);
				//printf(" u=2 alpha=%d",alpha_to[elp[2][i]]);
				//printf(" u=3 alpha=%d",alpha_to[elp[3][i]]);
				//printf("u=1 alpha=%o index=%d ",alpha_to[elp[1][i]],elp[1][i]);
				//printf("u=2 alpha=%o index=%d ",alpha_to[elp[2][i]],elp[2][i]);
			//	printf("u=1 alpha=%o index=%d ",alpha_to[elp[1][i]],elp[1][i]);
				printf("u=4 alpha=%o index=%d ",alpha_to[elp[4][i]],elp[4][i]);
				printf("u=6 alpha=%o index=%d ",alpha_to[elp[6][i]],elp[6][i]);
				printf(" u=%d index=%3d alpha=%3o\n ",u, elp[u][i],alpha_to[elp[u][i]]);
			}
				printf("\n");*/
			printf("Roots: ");

			/* Chien search: find roots of the error location polynomial */
			for (i = 1; i <= l[u]; i++)
				reg[i] = elp[u][i];
			count = 0;
			for (i = 1; i <= n; i++) {//此处是i<=n，而不是length，需要搜索整个BCH码长度范围内的根。n时GFsize。而不是bch码的码长。
				q = 1;
				for (j = 1; j <= l[u]; j++)
					if (reg[j] != -1) {//reg【i】是指数形式，=-1时表示系数为0；
						reg[j] = (reg[j] + j) % n;
						q ^= alpha_to[reg[j]];
					}
				if (!q) {	/* store root and error
						 * location number indices */
					root[count] = i;
					loc[count] = n - i;
					count++;
					printf("%3d ", n - i);
				}
			}
			printf("\n");
			if (count == l[u])
			/* no. roots = degree of elp hence <= t errors */
				for (i = 0; i < l[u]; i++)
					recd[loc[i]] ^= 1;
			else	/* elp has degree >t hence cannot solve */
				printf("Incomplete decoding: errors detected\n");
		}
	}
}









/*using IBM algorithm to implement the decode , IBM algorithm can get free from inversion operation which is easy to implement in hardware   */
void decode_bch_ibm()
{
	register int    i, j, u, q, t2, count = 0, syn_error = 0;

	int             elp[100][100], d[100], k[100][100], delta[100], s[100],deg[100],elp_to_1[100][100];
	int             root[200], loc[200], err[100], reg[201];
	int             elp_temp1,elp_temp2;
	int h;

			t2 = 2 * t;

		//printf("S(x) =\n ");
	for (i = 1; i <= t2; i++) {//the s[i] compute from 0 to t2;
		s[i] = 0;
		for (j = 0; j < length; j++)
			if (recd[j] != 0)//高位的数据在后。alpha 更高。
				s[i] ^= alpha_to[(i * j) % n];//(i * j) % n，GF(2048)是一个域，当i*j超过n时，alpha_to[i*j]=alpha_to[(i*j)%n].

		if (s[i] != 0)
			syn_error = 1; /* set error flag if non-zero syndrome 表示有错误*/

	s[i] = index_of[s[i]];
//        printf("in 8 the s[%d] = %o\n",i,s[i]);
        #ifdef DEBUG
        printf("in 10 the s[%d] = %d\n",i,s[i]);
        #endif // DEBUG
	}

    printf("starting the ibm decoding!\n");

	if(syn_error)
	{
	    memset(elp,-1,sizeof(elp));
	    memset(elp_to_1,-1,sizeof(elp_to_1));
	    memset(k,-1,sizeof(k));
		d[0] = 0;     /*index form*/
		d[1] = s[1];  /*index form*/
		elp[0][0] = 0;/*index form*/
		k[0][0] = 0;  /*index form*/
		elp[1][0] = 0; /*index*/
		k[1][0] = 0;  /*index*/

		for(i =1 ;i < t2; i++)
		{
			elp[0][i] = -1 ;/*index*/
			elp[1][i] = -1 ;/*index*/
			k[0][i] = -1 ;/*index*/
			k[1][i] = -1 ;/*index*/
		}
		delta[0] = 0 ;/*index*/
		delta[1] = 0 ;/*index*///delta computer start from 1, during the first time delta[1] == 0 in index  ;
		deg[0] = 0 ;  /*not a ff number*/
		deg[1] = 0 ;  /*not a ff number */
		u = 0;
		/*initial finished */


		do{
			u++;
			for(i = 0;i <= deg[u]+5;i++)/*compute elp*/
			{
				if(i == 0)
					elp[u+1][i] = (delta[u]+elp[u][i])%GFSIZE;//index multiplier ,just add them together
				else
                    {
                    if((d[u] != -1)&&(k[u][i-1] != -1)) //the k[u][i-1] must do not == -1 ,if that the add op can not present the ff multiplier
                        elp_temp1 = alpha_to[(d[u]+k[u][i-1])%GFSIZE];
                    else
                        elp_temp1 = 0;//in polynomial form
                    if((delta[u] != -1)&&(elp[u][i] !=-1))
                        elp_temp2 = alpha_to[(delta[u]+elp[u][i])%GFSIZE];
                    else
                        elp_temp2 = 0;
                    elp[u+1][i] =index_of[elp_temp1^elp_temp2];
//                        elp[u+1][i] = index_of[alpha_to[(delta[u]+elp[u][i])] ^ alpha_to[(d[u]+k[u][i-1])]];
                    }
//                printf("in 8 index  the elp[%d][%d]= %o\n",u+1,i,elp[u+1][i]);
                #ifdef DEBUG
                if(elp[u+1][i] != -1)
                {


                    printf("in 10 index the elp[%d][%d]= %d",u+1,i,elp[u+1][i]);
                    if(elp[u+1][i] >= elp[u+1][0])
                        printf( "  nomal elp[%d][%d] = %d \n",u+1,i,elp[u+1][i]-elp[u+1][0]);
                    else
                        printf( "  nomal elp[%d][%d] = %d \n",u+1,i,elp[u+1][i]+GFSIZE-elp[u+1][0]);
    //                if(k[u][i] != -1)
    //                printf("in 10 index the k[%d][%d] = %d\n",u,i,k[u][i]);
                }
                #endif // DEBUG

			}



			//if(((deg[u] +1)> u)|(d[u] == -1))/*compute the delta and k */
            if(((deg[u]+1)> u)|(d[u] == -1))//the elp[u] and deg[u] start from 2 instead of 1 ,so the deg[u] +1 > u
			{
				delta[u+1] = delta[u];
				for(i = 0;i <= t2 ; i++)//change the i<=deg[u]+2 to t2 ,because the degree of k may be huge
				{
					if(i < 2)
					{
						k[u+1][i] = -1;
					}
					else
                    {
                        k[u+1][i] = k[u][i-2];

                    }

//                    printf("in 8 index k[%d][%d] = %o\n",u+1,i,k[u+1][i]);
                    #ifdef DEBUG
                    if(k[u+1][i] != -1 )
                    printf("in 10 index k[%d][%d] = %d\n",u+1,i,k[u+1][i]);
                    #endif // DEBUG
				}
			}
			else
			{
				delta[u+1] = d [u];
				for(i = 0;i <= deg[u]+2;i++)
				{
					if(i == 0)
						k[u+1][i] = -1;
					else
						k[u+1][i] = elp[u][i-1];
//                    printf("in 8 index k[%d][%d] = %o\n",u+1,i,k[u+1][i]);
                    #ifdef DEBUG
                    if(k[u+1][i] != -1)
                    printf("in 10 index k[%d][%d] = %d\n",u+1,i,k[u+1][i]);
                    #endif // DEBUG
				}

			}
//			printf("in 8 index delta[%d]= %o\n",u+1,delta[u+1]);
            #ifdef DEBUG

			printf("in 10 index delta[%d]= %d\n",u+1,delta[u+1]);
            #endif // DEBUG


			/*compute the deg of next loop*/
			deg[u+1] = 0;
			for(i =1; i<= deg[u]+3 ;i++)
			{
				if(elp[u+1][i] != -1)
					deg[u+1] = i ;
			}
//			printf("in 8 index deg[%d]= %o\n",u+1,deg[u+1]);
            #ifdef DEBUG
			printf("in 10 index deg[%d]= %d\n",u+1,deg[u+1]);
			#endif // DEBUG

			/*compute the d of next loop */
			d[u+1] = -1 ;/*index */
			for(i = 0;i <= deg[u+1];i++)
			{
			    if((elp[u+1][i] != -1)&&(s[2*u+1-i] != -1))
                {
                    h = (s[2*u+1-i]+elp[u+1][i])%GFSIZE ;
                }
                else
                    h = -1; //h is in index form ,so when s[2*u+1-i]== -1 ,the h is zero which in index form is -1
                //d[u+1] = index_of[ alpha_to[(s[2*u+1-i]+elp[u+1][i])%GFSIZE] ^ alpha_to[d[u+1]] ];
                d[u+1] = index_of[ alpha_to[h] ^ alpha_to[d[u+1]] ];


			}

//			printf("in 8 index d[%d]= %o\n",u+1,d[u+1]);
            #ifdef DEBUG
			printf("in 10 index d[%d]= %d\n",u+1,d[u+1]);
//			printf("the index of alpha_to[1102]+ alpha_to[803] = %d\n",index_of[alpha_to[1102]^alpha_to[803]]);

//            h = 785 ^ 1289 ^ 185;
//			printf("h = %d\n",h);
//			printf("index_of[h] = %d\n",index_of[h]);
            #endif // DEBUG
		}while(u < t);


		u++;

/*normalization the elp*/

        for(i = 0;i <= deg[u] ;i++)
        {
            if(elp[u][i] != -1 )
                 if(elp[u][i] < elp[u][0])
                    elp_to_1[u][i] = elp[u][i] + GFSIZE - elp[u][0];
                 else
                     elp_to_1[u][i] = elp[u][i] - elp[u][0];


          //  printf("in 10 the index of elp[%d][%d] = %d \n",u,i,elp_to_1[u][i]);
        }



printf("Roots: ");
	/* Chien search: find roots of the error location polynomial */
			for (i = 1; i <= deg[u]; i++)//start from i =1,that's is to say  ,that the elp[u][0] == 1 in polynomial
				reg[i] = elp_to_1[u][i];
			count = 0;
			for (i = 1; i <= n; i++) {//此处是i<=n，而不是length，需要搜索整个BCH码长度范围内的根。n时GFsize。而不是bch码的码长。
				q = 1;
				for (j = 1; j <= deg[u]; j++)
					if (reg[j] != -1) {//reg【i】是指数形式，=-1时表示系数为0；
						reg[j] = (reg[j] + j) % n;
						q ^= alpha_to[reg[j]];
					}
				if (!q) {	/* store root and error
						 * location number indices */
					root[count] = i;
					loc[count] = n - i;
					count++;
					printf("%3d ", n - i);
				}
			}
			printf("\n");
			if (count == deg[u])
			/* no. roots = degree of elp hence <= t errors */
				for (i = 0; i < deg[u]; i++)
					recd[loc[i]] ^= 1;
			else	/* elp has degree >t hence cannot solve */
				printf("Incomplete decoding: errors detected\n");



	}



}


int main(int argc, char* argv[])
{
	int             i;
	FILE*           datafile;
	char            ch;
	int             numerr = 40;
	int             errpos[CODESIZE] = {6604,342,5070,5540,4239,6499,243,4278,479,4846};
	// {2416,2024,7678,586,2280,5543,7239,5582,6392,5757};/*this is a error mode than the decode_bch_ibm cam not correct it */
//	{2196,4753,5770,2116,278,6095,3397,7585,1166,5745};/*this is a error mode than the decode_bch_ibm cam not correct it */

                                        /*{2,34,55,66,87,             98,99,133,234,345,        431,436,456,478,498,
                                       501,505,523,567,588,       675,777,787,875,885,       990,991,993,994,997,
                                        1011,1012,1014,1019,1020,        1029, 1030,1035,1039,1045};*/
                                      // }; 1054  };//{1,55,721,631,900   ,555,66,77,345,668    ,1088,22,379,227};//,56,876};// ={3,8,10,178,129,556.123,234,222,789,  389,489,589,782,  235,981,984,872,  345,567,333,555,888,  998,876,751,672,984, 666,568,889,1022};//,1234,1111,4456,7826,4953,5689,2058,1111,2569.3487,990,7800,6542,5560,4892};//{43,42,41,40}; //{44,1064,1066};///{45,46,1063,1066};/,43,45,47,48};
	int             decerror = 0;
    int             test_error  = TEST_COUNT;



    LARGE_INTEGER encoded_begin , encoded_end,start,end_time,nFreq;// , encoded_time;
	LARGE_INTEGER decoded_begin  , decoded_end ;// , decoded_time =0;
	double encoded_time,decoded_time,cost,encoded_time_aver = 0,decoded_time_aver = 0;

	//start = clock ();
    int exec_counter;

	/*
	*Construct the Galois Field GF(2**m)
	*/
	generate_gf();



	/*
	* 通过fgetc()从文件中获取要编码的数据，存到data[]中
	*/
	datafile = fopen("test_data_nfc.txt", "r");
	if(datafile == NULL)
	{
		printf("open file test_data_1068_bit.txt failed!\n");
		return 1;
	}
	for(i=0 ; (i<INFOSIZE)&&(feof(datafile)==0); i++)
	{
		while(feof(datafile)==0)
		{
			ch = fgetc(datafile);
			if((ch=='0')||(ch=='1'))
			{
				data[i] = ch-'0';//得到data[i]提供给錯ncode_bch使用。
				break;
			}
		}
	}
	fclose(datafile);

	for(exec_counter = 0;exec_counter < TEST_COUNT ;exec_counter ++ )
{
        srand((unsigned) (time(NULL)+exec_counter));
        for(i = 0;i < numerr ;i++ )
            errpos[i] = rand()%8028; /**/

    #ifdef HARDWARE_TEST
        data[0] = 1;
        for(i=1;i<1024;i++)
        {
            data[i] = !data[i-1];

        }
        data[1023] = 1;
    #endif // HARDWARE_TEST

	/*
	*对获取的data数据进行编码，获取编码后的冗余位
	*/


	//encoded_begin = clock();
    QueryPerformanceFrequency(&nFreq);
    QueryPerformanceCounter(&encoded_begin);

	encode_bch();

    QueryPerformanceCounter(&encoded_end);
	//encoded_end = clock();
	encoded_time = (double)(encoded_end.QuadPart-encoded_begin.QuadPart)/(double)nFreq.QuadPart;

	encoded_time_aver = encoded_time + encoded_time_aver;
	#ifdef TIME_TEST
	if(exec_counter == TEST_COUNT - 1)
    {
        printf("the average encode time is %lf\n",encoded_time_aver/TEST_COUNT);
        printf("the encode time is %lf\n",encoded_time);

    }
    #endif // TEST_COUNT


	/*
	*把冗余位和数据位都存放到recd[]中，这就是编码后的所有信息位
	*其中0-559位是冗余位，560-8768是数据位，打印输出到文件中
	*/

	datafile = fopen("test_data_recd_nfc.txt", "w");
	fprintf(datafile,"c(x) = \n");
	for (i = 0; i < REDUNDANCYSIZE; i++)
	{
		recd[i] = bb[i];//先输出的是余数多项式的低位，

		if (i && ((i % 64) == 0))
			fprintf(datafile,"\n");

		fprintf(datafile,"%1d", recd[i]);
	}
	fprintf(datafile,"\n");
	for (i = 0; i < INFOSIZE; i++)
	{
		recd[i + REDUNDANCYSIZE] = data[i];//将从test_data_1068_bit.txt中读出的数据存放到test_data_1068_recd_bit.txt中，最高位存放到最后。

		if (i && ((i % 64) == 0))
			fprintf(datafile,"\n");

		fprintf(datafile,"%1d", recd[i + REDUNDANCYSIZE]);//将数据位存放到test_data_1068_recd_bit.txt中，
	}
	fprintf(datafile,"\n");



	/*
	*根据errpos[]把编码后的recd信息中某些位置为错误
	*方法是直接取反，然后打印输出到文件中
	*/
	datafile = fopen("test_data_error_nfc.txt", "w");
	if(datafile == NULL)
        printf("open file failed!\n");
    else
        printf("the file opened succeed\n");
	if (numerr)
		for (i = 0; i < numerr; i++)
			recd[errpos[i]] ^= 1;//异或1，取反，将test_data_1068_recd_bit.txt中的某些位数取反，形成错误位数。
	fprintf(datafile,"r(x) = \n");

	for (i = 0; i < REDUNDANCYSIZE; i++)
	{
		if (i && ((i % 64) == 0))
			fprintf(datafile,"\n");
		fprintf(datafile,"%1d", recd[i]);
	}
	fprintf(datafile,"\n");
	for (i = 0; i < INFOSIZE; i++)
	{
		if (i && ((i % 64) == 0))
			fprintf(datafile,"\n");
		fprintf(datafile,"%1d", recd[i + REDUNDANCYSIZE]);//将标记为错误位的数据取反之后存入test_error_1068_recd_bit.txt中。
	}
	fprintf(datafile,"\n");




	/*
	*对recd中的数据进行解码，判断数据是否出错
	*如果数据出错，就找出出错位置，并进行纠正
	*/


	//decoded_begin = clock();
    QueryPerformanceCounter(&decoded_begin);

    #ifndef IBM
    decode_bch();
    #else
    decode_bch_ibm();
    #endif

    QueryPerformanceCounter(&decoded_end);
	//encoded_end = clock();
	decoded_time = (double)(decoded_end.QuadPart-decoded_begin.QuadPart)/(double)nFreq.QuadPart;
    decoded_time_aver = decoded_time + decoded_time_aver;

	//decoded_end = clock();
//	decoded_time =(double)( decoded_end - decoded_begin);
    #ifdef TIME_TEST
    if(exec_counter == TEST_COUNT - 1)
    {
        printf("the average decode time is %lf\n",decoded_time_aver/TEST_COUNT);
        printf("the decode time is %lf\n",decoded_time);

    }
    #endif // TEST_COUNT

	/*
	*DECODING ERRORS? we compare only the data portion
	*length是编码长度，k是信息位长度 length-k为冗余位长度
	*判断接受到的数据解码后，就是recd与data是否是一样的
	*/
	for (i = length - k; i < length; i++)
	{
		if (data[i - length + k] != recd[i])
		{
			decerror++;
			printf("Dismatch positions %2d\n", i);
		}
	}
	if (decerror)
    {
        printf("There were %d decoding errors in message positions\n", decerror);
        for(i=0;i<numerr;i++)
            printf("the error position is %4d\n",errpos[i]);
        printf("this is the %dth time test\n",exec_counter);
        system("pause");
    }

	else
    {
         test_error--;
         printf("Succesful decoding\n");
    }
}

    QueryPerformanceCounter(&start);
//    for(i=0;i<923990045;i++)
//    {
//        i= i+1;
//    }
	//end_time = clock();
	//cost = (double)(end_time - start);



    QueryPerformanceCounter(&end_time);
	//encoded_end = clock();
	cost = (double)(end_time.QuadPart-start.QuadPart)/(double)nFreq.QuadPart;
	//printf("the start time is%d\n", start);
//	printf("the end time is %d\n",end_time);
    if(exec_counter ==TEST_COUNT)
	printf("the total time is %lf \n",cost);

    printf("the unsuccessful times is %d\n",test_error);
//////    i = getchar();
////    system("pause");
//    printf("over");
	return 0;


}

 /*

 S(x) = 3208 6416 10023 12832 14470 3663 3802 9281 5233 12557 1045 7326 1516 7604
 7555 2179 750 10466 4306 8731 5442 2090 7954 14652 15141 3032 10207 15208 11330
 15110 14479 4358 12816 1500 8304 4549 9929 8612 9384 1079 14948 10884 14074 418
0 1586 15908 10978 12921 3674 13899 8231 6064 5643 4031 3829 14033 14036 6277 62
49 13837 573 12575 1356 8716 11187 9249 11021 3000 11009 225 5369 9098 7835 3475
 3171 841 7642 2385 9818 2158
sigma(x) =   0 3208 15133 3695
Roots: 3543  97  55
Succesful decoding
Press any key to continue

 */





