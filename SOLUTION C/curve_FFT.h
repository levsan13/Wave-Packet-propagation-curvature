#ifndef CURVE_FFT_H_INCLUDED
#define CURVE_FFT_H_INCLUDED

#include<stdio.h>
#include<stdlib.h>
#include<complex.h>
#include"matriz.h"
#include<math.h>

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

typedef matrix curv(double param1,double param2,double param3,double param4,double *param5,double *param6,int nx,int ny);

/* * * * * * * * * * * * * * *
* 0 - - GAUSSIAN
* 1 - - PARABOLOIDE
* 2 - - PONTO DE CELA
* 3 - - SENO
* 4 - - COSSENO
* * * * * * * * * * * * * * * */

curv *function[5];

matrix gaussian(double xsig,double ysig,double lenght, double x_0, double y_0, double *X,double *Y,int nx,int ny){
int i,j;
matrix function=(matrix)allocmatrix(nx,ny,sizeof(double));
for(i=0;i<nx;i++)
for(j=0;j<ny;j++)
function[i][j]=lenght*(exp(-pow(X[i]-x_0,2)/(2*pow(xsig,2)))+exp(-pow(Y[j]-y_0,2)/(2*pow(ysig,2))));
return function;
}

matrix paraboloide(double C,double A_x,double A_y, double x_0, double y_0, double *X,double *Y,int nx,int ny){

int i,j;
matrix function=(matrix)allocmatrix(nx,ny,sizeof(double));
for(i=0;i<nx;i++)
for(j=0;j<ny;j++)
function[i][j]=A_x*pow(X[i]-x_0,2)+A_y*pow(Y[j]-y_0,2)+C;
return function;

}

matrix ponto_de_cela(double C,double A_x,double A_y, double x_0, double y_0, double *X,double *Y,int nx,int ny){

int i,j;
matrix function=(matrix)allocmatrix(nx,ny,sizeof(double));
for(i=0;i<nx;i++)
for(j=0;j<ny;j++)
function[i][j]=A_x*pow(X[i]-x_0,2)-A_y*pow(Y[j]-y_0,2)+C;
return function;

}


//===============================================================================================================================
//                                                                                    !FFT
//=================================================================================================================================

/* Replacesdataby itsndim-dimensional discrete Fourier transform, ifisignis input as 1.nn[1..ndim]is an integer array containing the lengths of each dimension (number of complexvalues), which MUST all be powers of 2.datais a real array of length twice the product ofthese lengths, in which the data are stored as in a multidimensional complex array: real andimaginary parts of each element are in consecutive locations, and the rightmost index of thearray increases most rapidly as one proceeds alongdata. For a two-dimensional array, this isequivalent to storing the array by rows. Ifisignis input as−1,datais replaced by its inversetransform times the product of the lengths of all dimensions.
524Chapter 12.Fast Fourier TransformSample page from NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)Copyright (C) 1988-1992 by Cambridge University Press. Programs Copyright (C) 1988-1992 by Numerical Recipes Software. Permission is granted for internet users to make one paper copy for their own personal use. Further reproduction, or any copying of machine-readable files (including this one) to any server computer, is strictly prohibited. To order Numerical Recipes books or CDROMs, visit websitehttp://www.nr.com or call 1-800-872-7423 (North America only), or send email to directcustserv@cambridge.org (outside North America).
*/

void fourn(double *data, unsigned long nn[], int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	double tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=0;idim<ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim-1;idim>=0;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
                        //printf("uuuuuu %d\n\n",i3);
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {

                        //printf("i2uuuuuu %d\n\n",i2);
						k1=i2;
						k2=k1+ifp1;
						tempr=(double)wr*data[k2]-(double)wi*data[k2+1];
						tempi=(double)wr*data[k2+1]+(double)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}

void FFT(cmatrix PSI, unsigned long nn[],int ndim, int isign,int nx, int ny){

int i,j;
double *data=(double*)malloc(2*nx*ny*sizeof(double));
if(data==NULL) exit(1);
double ainv;
int c,J,K;

for(j=0;j<ny;j++)
for(i=0;i<nx;i++){
c=2*j;
J=c*nx+2*i-1;
K=c*nx+2*i;
//printf("i=%d j=%d J=%d K=%d\n\n",i,j,J,K);

data[J] =creal(PSI[i][j]);
data[K] =cimag(PSI[i][j]);
}
nn[0]=nn[0]>>1;
nn[1]=nn[1]>>1;
//printf("ANTES\n\n");
fourn(data,nn,ndim,isign);
//printf("DEPOIS\n\n");
ainv=1.0;
if(isign==1)
ainv=1/(nx*ny);

J=0.0;
K=0.0;
c=0;

for(j=0;j<ny;j++)
for(i=0;i<nx;i++){
c=2*j;
J=c*nx+2*i-1;
K=c*nx+2*i;
//printf("i=%d j=%d J=%d K=%d\n\n",i,j,J,K);
PSI[i][j]=ainv*(data[J]+I*data[K]);
}

return;

}



#endif // CURVE_FFT_H_INCLUDED
