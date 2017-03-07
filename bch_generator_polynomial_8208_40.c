/*
 * File:    bch_4148_4096_9_v100.cpp
 * Title:   Encoder/decoder for binary BCH codes in C (Version 3.1)
 * Author:  Robert Morelos-Zaragoza
 * Date:    August 1994
 * Revised: June 13, 1997
 *
 * ===============  Encoder/Decoder for binary BCH codes in C =================
 *
 * Version 1:   Original program. The user provides the generator polynomial
 *              of the code (cumbersome!).
 * Version 2:   Computes the generator polynomial of the code.
 * Version 3:   No need to input the coefficients of a primitive polynomial of
 *              degree m, used to construct the Galois Field GF(2**m). The
 *              program now works for any binary BCH code of length such that:
 *              2**(m-1) - 1 < length <= 2**m - 1
 *
 * Note:        You may have to change the size of the arrays to make it work.
 *
 * The encoding and decoding methods used in this program are based on the
 * book "Error Control Coding: Fundamentals and Applications", by Lin and
 * Costello, Prentice Hall, 1983.
 *
 * Thanks to Patrick Boyle (pboyle@era.com) for his observation that 'bch2.c'
 * did not work for lengths other than 2**m-1 which led to this new version.
 * Portions of this program are from 'rs.c', a Reed-Solomon encoder/decoder
 * in C, written by Simon Rockliff (simon@augean.ua.oz.au) on 21/9/89. The
 * previous version of the BCH encoder/decoder in C, 'bch2.c', was written by
 * Robert Morelos-Zaragoza (robert@spectra.eng.hawaii.edu) on 5/19/92.
 *
 * NOTE:
 *          The author is not responsible for any malfunctioning of
 *          this program, nor for any damage caused by it. Please include the
 *          original program along with these comments in any redistribution.
 *
 *  For more information, suggestions, or other ideas on implementing error
 *  correcting codes, please contact me at:
 *
 *                           Robert Morelos-Zaragoza
 *                           5120 Woodway, Suite 7036
 *                           Houston, Texas 77056
 *
 *                    email: r.morelos-zaragoza@ieee.org
 *
 * COPYRIGHT NOTICE: This computer program is free for non-commercial purposes.
 * You may implement this program for any non-commercial application. You may
 * also implement this program for commercial purposes, provided that you
 * obtain my written permission. Any modification of this program is covered
 * by this copyright.
 *
 * == Copyright (c) 1994-7,  Robert Morelos-Zaragoza. All rights reserved.  ==
 *
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
 *
 */





#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define DIMENSION 14//4//3  //
#define GFSIZE 16383//15 //7 //

//本原多项式的系数 0->m 共15位 表2-7
//1+x+x^6+x^10+x^14
#define MINPOLYNOMIAL {1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1}//{ 1,1,0,0,1 }//{1,1,0,1} //

#define INFOSIZE  8192 //7//4  //8208    //有效数据位
#define CORRECTCAPACITY 32 //2// 1//  40  //纠错能力
#define REDUNDANCYSIZE 448 //8//3  //560  //冗余位
#define CODESIZE (INFOSIZE+REDUNDANCYSIZE)  //编码后的长度




//生成多项式的系数
//#define GENPOLYNOMIAL

// code parameter.
int             m = DIMENSION, n = GFSIZE, k = INFOSIZE, length = CODESIZE;
int             p[DIMENSION+1] = MINPOLYNOMIAL;
int             t = CORRECTCAPACITY;
int             d = (2*CORRECTCAPACITY+1);
//int             g[REDUNDANCYSIZE+1] = GENPOLYNOMIAL;
int             g[REDUNDANCYSIZE+1];

// GF(2**m)
int             alpha_to[GFSIZE+1], index_of[GFSIZE+1];
int             data[INFOSIZE];
int             bb[REDUNDANCYSIZE], recd[CODESIZE];








//计算伽罗伐域上的各个数据存在 alpha_to中
 void generate_gf()
/*
 * generate GF(2**m) from the irreducible polynomial p(X) in p[0]..p[m]
 * lookup tables:  index->polynomial form   alpha_to[] contains j=alpha**i;
 * polynomial form -> index form  index_of[j=alpha**i] = i alpha=2 is the
 * primitive element of GF(2**m)
 */
{
	register int    i, mask;
	mask = 1;
	alpha_to[m] = 0;
	for (i = 0; i < m; i++) {
		alpha_to[i] = mask;
		index_of[alpha_to[i]] = i;
		if (p[i] != 0)
			alpha_to[m] ^= mask;
		mask <<= 1;
	}
	index_of[alpha_to[m]] = m;
	mask >>= 1;
	for (i = m + 1; i < n; i++) {
		if (alpha_to[i - 1] >= mask)
		  alpha_to[i] = alpha_to[m] ^ ((alpha_to[i - 1] ^ mask) << 1);
		else
		  alpha_to[i] = alpha_to[i - 1] << 1;
		index_of[alpha_to[i]] = i;
	}
	index_of[0] = -1;
}


int * mult_poly(int *poly1,int *poly2,int size1,int size2)
 {
	// poly1,poly2 is in alpha form 
	 int *poly;// in alpha form 
	 int polys[10000] = {0};
	 int i,j;
	 poly = (int*)malloc(sizeof(int)*(size1+size2-1));
	 if(poly == NULL)
	 {
		 printf("malloc failded!\n");
		 return NULL;
	 }
	 else
	 {
		 memset(poly,0,sizeof(int)*(size1+size2-1));
	 }
	// int poly1[] = {0,1,0,4,2};
	//int size1 = 5;
	//int poly2[] = {0,0,1,5,0,1};
	//int size2 = 6;
	 for(i=0;i<size1;i++)
	 {
		 /*printf("\n");*/
		for(j=0;j<size2;j++)
		{
			if(poly1[i]&&poly2[j])
				poly[i+j] = poly[i+j] ^ alpha_to[(index_of[poly1[i]]+index_of[poly2[j]])%n]; 
			else
				poly[i+j] = poly[i+j];
			//printf("poly[%d] = %4d;   ",i+j,poly[i+j]);
		}
	 }
	// printf("\n the poly is \n");
	//for(i=0;i<size1+size2-1;i++)
	//{
	//	//poly[i] = index_of[poly[i]];
	//	printf("%5d",poly[i]);
	//}
		

	 for(i=0;i<size1;i++)
	 {
		for(j=0;j<size2;j++)
		{
			 polys[i+j] = polys[i+j] + poly1[i]*poly2[j]; 
		}
	 }
	 //printf("the polys is \n");
	 ////for(i=0;i<size1+size2-1;i++)
		//// printf("%5d",polys[i]);

	 return poly;
		
 }

 void gen_min_poly(int adjoint[][16],int ad_size[])
 {
	 //adjoint in index form ;
	register int alpha,i,j;
	int min_poly[200][1000]={0};// in alpha form ,

	int *poly1=NULL,*poly=NULL;
	int size1 = 5;
	int poly2[2] = {0};
	int size2 = 6;
	int size =0 ;
	int *min_polys[100];
	int min_polys_size[100];
	FILE *fp = NULL;
	for(alpha=0;alpha<(d-1)/2;alpha++)
	{
		poly1 = (int*)malloc(sizeof(int)*(2));
		 if(poly1 == NULL)
		 {
			 printf("malloc failded!\n");
		 }
		 else
		 {
			 memset(poly1,0,sizeof(int)*(2));
		 }

		 poly1[1] = alpha_to[adjoint[alpha][0]];
		 poly1[0] = 1 ;
		 size  = 2 ;
		for(i=1;i<ad_size[alpha];i++)
		{
			poly2[1] = alpha_to[adjoint[alpha][i]];
			poly2[0] = 1 ;
			poly = mult_poly(poly1,poly2,size,2);
			size = size + 2 - 1;
			free(poly1);
			poly1=NULL;
			if(i<ad_size[alpha]-1)
				poly1 = poly ;
			
			for(j=0;j<size;j++)
				printf("%5d",poly[j]);
			printf("\n");

		}
		min_polys[alpha] = poly;
		min_polys_size[alpha]=size;
		printf("\n the min poly [%d]",alpha);
		for(i=0;i<size;i++)
			printf("%5d",poly[i]);
		printf("\n");
	}

	fp = fopen("poly.txt","w");

	fprintf(fp,"the min poly is :            #the first is the highest! \n");
	for(alpha=0;alpha<(d-1)/2;alpha++)
	{
		fprintf(fp,"alpha[%d]'s min_poly[%d] = ",adjoint[alpha][0],alpha);
		printf("min_poly[%d] = ",alpha);
		for(i=0;i<min_polys_size[alpha];i++)
		{
			fprintf(fp,"%d",min_polys[alpha][i]);
			printf("%d",min_polys[alpha][i]);
		}
		fprintf(fp,"\n");
		printf("\n");
	}


	poly1 = min_polys[0];
	size = min_polys_size[0];
	for(alpha=1;alpha<(d-1)/2;alpha++)
	{
		poly = mult_poly(poly1,min_polys[alpha],size,min_polys_size[alpha]);
		poly1 = poly;
		size = size + min_polys_size[alpha]-1;
	}

	printf("THE GEN_POLY IS \n");
	fprintf(fp,"THE GEN_POLY IS :               #the first is the highest ");
	for(i=0;i<size;i++)
	{
		if(i%31 == 0)
		fprintf(fp,"\n");
		printf("%5d",poly[i]);
		fprintf(fp,"%d",poly[i]);

	}
	fclose(fp);




	for(alpha=0;alpha<d-1;alpha++)
	{
		min_poly[alpha][1] = 1;
		min_poly[alpha][0] = alpha_to[adjoint[alpha][0]];
		for(j=1;j<ad_size[alpha];j++)
		{
			printf("\n");
			for(i=DIMENSION;i>=0;i--)
			{
				if(i>0)
				{
					if(min_poly[alpha][i] == 0)
						min_poly[alpha][i] = min_poly[alpha][i-1];
					else
						min_poly[alpha][i] = min_poly[alpha][i-1] ^alpha_to[(index_of[min_poly[alpha][i]] + adjoint[alpha][j])%n] ; 
				}
					
				else
				{
					min_poly[alpha][i] = ((index_of[min_poly[alpha][i]] + adjoint[alpha][j])%n);
					min_poly[alpha][i] = alpha_to[min_poly[alpha][i]];
				}
				if((min_poly[alpha][i] !=0)&&(j==ad_size[alpha]-1))
				{
					printf("the %dth min_poly[%d][%d] in alpha form  is: %d index form :%d\n",j,alpha,i,min_poly[alpha][i],index_of[min_poly[alpha][i]]);
					
				}
			}

		}

	}
	//use gen to verify the min_poly ;

	

 
 }


 /*
 * Compute generator polynomial of BCH code of length = 4096+52 = 4148, redundancy = 52 m=13
 * (OK, this is not very efficient, but we only do it once, right? :)
 * 计算生成多项式 g[];
 */

void gen_poly()
{
	register int    ii, jj, ll, kaux,kk;
	int             test, aux, nocycles, root, noterms, rdncy;
	FILE *fp;

	int             cycle[2000][15], size[2000], min[2000], zeros[2000],adjoint[100][16]= {0},ad_size[100];
	//int cycle[13][7], size[13], min[13], zeros[13];
	/* Generate cycle sets modulo 63 */
	cycle[0][0] = 0; size[0] = 1;
	cycle[1][0] = 1; size[1] = 1;
	jj = 1;			/* cycle set index */
	do {
		/* Generate the jj-th cycle set */
		ii = 0;
		do {
			ii++;
			cycle[jj][ii] = (cycle[jj][ii - 1] * 2) % n;
			size[jj]++;
			aux = (cycle[jj][ii] * 2) % n;
		} while (aux != cycle[jj][0]);
		/* Next cycle set representative */
		ll = 0;
		do {
			ll++;
			test = 0;
			for (ii = 1; ((ii <= jj) && (!test)); ii++)
			/* Examine previous cycle sets */
			  for (kaux = 0; ((kaux < size[ii]) && (!test)); kaux++)
					if (ll == cycle[ii][kaux])
						test = 1;
		}
		while ((test) && (ll < (n - 1)));

		if (!(test)) {
			jj++;	/* next cycle set index */
			cycle[jj][0] = ll;
			size[jj] = 1;
		}
	} while (ll < (n - 1));

	//printf("the compute answer is :\n");
	//for(ii=0;ii<jj;ii++)
	//{
	//	printf("size[%d]=%d\n",ii,size[ii]);
	//	printf("cycle[%d]=",ii);
	//	for(kk=0;kk<size[ii];kk++)
	//	{
	//		printf(" %5d",cycle[ii][kk]);
	//	}
	//	printf("\n");
	//}

	nocycles = jj;		/* number of cycle sets modulo n */
	/* Search for roots 1, 2, ..., d-1 in cycle sets */
	kaux = 0;
	rdncy = 0;
	for (ii = 1; ii <= nocycles; ii++) {
		min[kaux] = 0;
		for (jj = 0; jj < size[ii]; jj++)
			for (root = 1; root < d; root++)
				if (root == cycle[ii][jj])
				{
					min[kaux] = ii;
				}
					
		if (min[kaux]) {
			rdncy += size[min[kaux]];
			kaux++;
		}
	}
	noterms = kaux;
	kaux = 1;
	printf("the min poly is :\n");
	for (ii = 0; ii < noterms; ii++)
		for (jj = 0; jj < size[min[ii]]; jj++)
		{
			zeros[kaux] = cycle[min[ii]][jj];
			adjoint[ii][jj] = cycle[min[ii]][jj];
			ad_size[ii] = size[min[ii]];
			//printf("%o\n",zeros[kaux]);
			kaux++;
		}
		printf("the min ploy num is %d\n",kaux);

		gen_min_poly(adjoint,ad_size);

	printf("This is a (%d, %d, %d) binary BCH code\n", length, k, d);
	printf("g(x) 的最高为为1，在最后面");
	/* Compute generator polynomial */
	g[0] = alpha_to[zeros[1]];
	g[1] = 1;		/* g(x) = (X + zeros[1]) initially */
	for (ii = 2; ii <= rdncy; ii++) {
	  g[ii] = 1;
	  for (jj = ii - 1; jj > 0; jj--)
	    if (g[jj] != 0)
	      g[jj] = g[jj - 1] ^ alpha_to[(index_of[g[jj]] + zeros[ii]) % n];
	    else
	      g[jj] = g[jj - 1];
	  g[0] = alpha_to[(index_of[g[0]] + zeros[ii]) % n];
	}
	printf("g(x) = ");
	fp=fopen("poly.txt","a");
	fprintf(fp,"\n");
	fprintf(fp,"This is a (%d, %d, %d) binary BCH code",length,k,d);
	fprintf(fp,"g(x) =                    #the last is the highest\n");
	for (ii = 0; ii <= rdncy; ii++)
	{
	  if ((ii % 31) == 0)
	  {
	     printf("\n");
		 fprintf(fp,"\n");
	  }
	  printf("%d", g[ii]);
	  fprintf(fp,"%d",g[ii]);

	}
	fprintf(fp,"\n");
	printf("\n");
	fclose(fp);
}





 main()
{
	int             i;

	generate_gf();			/* generate the Galois Field GF(2**m) 算伽罗伐域各个元素 */
	gen_poly();				/* Compute the generator polynomial of BCH code 计算生成多项式 */




	for (i = 0; i < REDUNDANCYSIZE+1; i++) {
		printf("%d, ", g[i]);
		if (i && ((i % 64) == 0))
			printf("\n");
	}
	printf("\n");
	getchar();
	system("pause");
}



/*   总共561位
This is a (8768, 8208, 81) binary BCH code
g(x) =
10000011101101001111011100111001010000000000010100100100111011010100011
0111000001010010110110010100011001111101100100010010010110111010111111
0111111011101000111111101010000100001000101111110100011011110011010010
1000111001100111111000111011110111110111000100000110000000001000101001
0101010010010011011001100110100000101011010100010011001110010101011111
1000011011111011100101000111101111110010010101111110010011111001010111
0010111101101101111001010101110011010010010000010110100110011111001100
0011001101001101110100011101100000101100011011001000111111110000000111

1000001110110100111101110011100101000000000001010010010011101101
0100011011100000101001011011001010001100111110110010001001001011
0111010111111011111101110100011111110101000010000100010111111010
0011011110011010010100011100110011111100011101111011111011100010
0000110000000001000101001010101001001001101100110011010000010101
1010100010011001110010101011111100001101111101110010100011110111
1110010010101111110010011111001010111001011110110110111100101010
1110011010010010000010110100110011111001100001100110100110111010
0011101100000101100011011001000111111110000000111
Press any key to continue



This is a (8640, 8192, 65) binary BCH code
g(x) = 10001001111110001111100101001000010101011101000111111000011111101011100
1000110110100100011111110101111010011100010001111011011011011001001111
1001010111000111010100100100100110101011011100011111001100111101110001
1001101100011101100000001101010011011001101100101100110000001100001010
0101011010001101111010111100010000011110000110110011001111000001000010
1101011011010111111110101010011101111000111100111001111010100011010000
0001000111011010110111101001
1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1,
0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1,
0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1,
0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0,
1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0,
0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1,
1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1,



1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1,
0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1,
0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1,
0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0,
1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0,
0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1,
1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1,

This is a (15, 7, 5) binary BCH code
g(x) 的最高为为1，在最后面g(x) = 100010111
1, 0, 0, 0, 1, 0, 1, 1, 1,


*/











