#include<stdio.h>
#include"curve_FFT.h"
#include"matriz.h"
#include<stdlib.h>

int main(){

int nx=1000;
int ny=1000;
int i,j;
double WX=5000;
double WY=5000;
double X[nx],Y[ny];
double dx=WX/(nx);
double dy=WY/(ny);
for(i=0;i<nx;i++)
X[i]=(double)(-WX/2)+(double)(i)*dx;
for(i=0;i<ny;i++)
Y[i]=(double)(-WY/2)+(double)(i)*dy;

matrix function=ponto_de_cela(0,10,10,0,0,X,Y,nx,ny);

FILE *graf=fopen("pc.dat","w");
for(i=0;i<nx;i++)
for(j=0;j<ny;j++)
fprintf(graf,"%0.3f \t %0.3f \t %0.3f\n",X[i],Y[j],function[i][j]);

fclose(graf);

system("python packet.py");

}
