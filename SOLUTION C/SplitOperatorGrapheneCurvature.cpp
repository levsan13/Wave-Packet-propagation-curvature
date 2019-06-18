/*********************************************************************************************************************
* Universidade Federal Do Ceará - (UFC) - Departamento de Física                                                     *
*                                                                                                                    *
*                                  *** Projeto de Iniciação Cientifica(IC)***                                        *
*         ***Propagação de pacote de onda em monocamada de grafeno com curvatura usando a tecnica Split Operator***  *
*                       ***Simulação em C da evolução temporal do pacote de onda gaussiano***                        *
*                                                                                                                    *
* Autor: Sergio Levy Nobre dos Santos - Graduando em Fisica - Bacharelado - UFC                                      *
* Orientador: João Milton Pereira Junior - Professor do Departamento de Fisica -UFC                                  *
* Co-orientador: Diego Rabelo da Costa - Professor do Departamento de Fisica -UFC                                    *
*                                                                                                                    *
**********************************************************************************************************************/



//========================================================================================================================
//                                                                                    !INCLUDES (Bibliotecas C)
//========================================================================================================================

#include<stdio.h>   //Gerenciamento de arquivos FILE*
#include<stdlib.h>  //Alocação de memoria
#include<math.h>    //funções matematicas
#include<complex.h> //Numeros complexos
#include"matriz.h"  //(by Levy) Implementação de matrizes reais e complexas
#include"curve_FFT.h"



//========================================================================================================================
//                                                                                    !DEFINIÃÃES
//========================================================================================================================

#define pi M_PI       //Definição de pi
#define TWOPI 2.0*pi  //Definição de 2*pi
#define raiz3 1.732050808
#define INITIAL 1
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

//========================================================================================================================
//                                                                                    !DefiniÃ§Ãµes de tipos de variaveis
//========================================================================================================================
typedef cmatrix wf;              // Tipo wf(wavefunction)


//========================================================================================================================
//                                                                                    !Constantes globais
//========================================================================================================================

const double a0 = 0.5292;               // 'a0' Ã© o raio de Bohr (dado em Angstron)
const double Ry = 13.6057;              // 'Ry' Ã© a energia de Rydberg (dada em eV)
const double hc = 1.0545890/1.6021890;  // 'hc' Ã© a constante h cortado (dada em eV*fs)
const double t =  2.90;                 // 't' hopping para primeiros vizinhos (eV)
const double vf = 13.9;                 // 'vf' velocidade de fermi (Angstron/fs)


const int nt=1;                       //numero de evoluções temporais
const double cdt=0.1;                  //passo do tempo (dado em fs)

const double l=1.42;                   //distancia entre os sitios A e B (dado em Angstron)
const int nx=520;                      //Numero de pontos no eixo x
const int ny=520;                      //Numero de pontos no eixo y
const int n=520;
const double WX=1000;        //tamanho do  eixo x (dado em Angstron)
const double WY=1000;    //tamanho do  eixo y (dado em Angstron)
const int ndim=2;
double dx;                              //passo no eixo x
double dy;                              //passo no eixo y
const double width=0.0;

//caracteristicas do pacote de onda inicial
const double sig=20.0;                  //sigma (Angstron)
const double phi=90.0;                   // Angulo de k (Graus)
const double ak0=(double)0.06;          //k inicial (Angstron^-1)
const int eta=1;                        // vale K ou K'

const int ipropagation=1;




/*- main
*- centermass (feito)
*- FFT (feito)
*- curv
*- evolution time
*- gaussian packet (feito)
*- POT
*- STRAIN TENSOR
*/


//========================================================================================================================
//                                                                                    !FunÃ§Ãµes/routinas
//========================================================================================================================

static matrix curv_gaussian(double sigma,double leigth, double *X,double *Y, double x0,double y0);

static matrix x2_y2(double param1, double *X,double *Y,double param2);


static void wavei(wf PSIA,wf PSIB,double *X,double *Y,double x0, double y0, double xsig, double ysig); //Criação de pacote de onda inicial

static void Norm(cmatrix PSIA,cmatrix PSIB);

static void Pot(double V0,double* X,double* Y, matrix V);                                              //Preenchimeto da matriz de potencial


static void centermass(wf CPSIA,wf CPSIB,double *X,double *Y,double *XCM,double *XCMr,double *XCMt,double *YCM,double *YCMr,double *YCMt); // Calculo de <x> e <y>

static void prob(double*X,double*Y,wf CPSIA,wf CPSIB,double*PBEF,double*PINT,double*PAFT);                                                 //Probabilidade








//========================================================================================================================
//                                                                                    !FunÃ§Ã£o Principal
//========================================================================================================================

int main(){

//====================================================Constantes, declaração de funções de onda, vetores======================================================
    int i,j,it,kx,ky;
    int nnx=nx/2,nny=ny/2;

    matrix V=(matrix)allocmatrix(nx,ny,sizeof(double)); //Potencial

    //WAVEFUNCTION
    MCarray PsiU,PsiD;
    PsiU.get_info_matrix(nx,ny);
    PsiD.get_info_matrix(nx,ny); // Psi UP and Down
    //WAVEFUNCTION

    //========== Variaveis de Evolução Temporal =====

    double ti;
    double tau=sig/vf;

    double ASUM=0.0,AINV;

    //Variaveis de posições
    double XCM; //! Calculando a média de x (<x>) no espaço real
    double XCMr;
    double XCMt;
    double YCM; //! Calculando a média de y (<y>) no espaço real
    double YCMr;
    double YCMt;

    //Probabilidade
    double PBEF;
    double PINT;
    double PAFT;

    //Arquivos
    FILE *xyavg=fopen("xyavg.dat","w");     //Valores medio da função de onda em função do tempo
    FILE *probe=fopen("Prob_end.dat","w");  //função de probalidade


    //========== Variaveis de Evolução Temporal =====

    //Caracteristicas do pacote inicial gaussiano

    double x0 = 0.0;   //PosiÃ§Ã£o x inicial do pacote de onda
    double y0 = 0.0;     //! PosiÃ§Ã£o y inicial do pacote de onda
    double xsig = sig; //largura do pacote de onda inicial
    double ysig = sig; //largura do pacote de onda inicial
    double V0=0.0;     //Potencial
    double E0=hc*vf*ak0;

    //Caracteristicas do pacote inicial gaussiano


    //=====DISCRETIZAÇÃO DO ESPAÇO REAL=======

    double X[nx],Y[ny];
    dx=WX/(nx);
    dy=WY/(ny);
    for(i=0;i<nx;i++)
    X[i]=(double)(-WX/2)+(double)(i)*dx;
    for(i=0;i<ny;i++)
    Y[i]=(double)(-WY/2)+(double)(i)*dy;

    //=====DISCRETIZAÇÃO DO ESPAÇO DE POSIÇÕES========

     //=====DISCRETIZAÇÃO DO ESPAÇO RECIPROCO=======

    double Kx[nx],Ky[ny];
    unsigned long NN[2];
    NN[0]=nx;
    NN[1]=ny;
    double NXY = nx*ny;
    double NXY2 = 2*NXY;
    double DKX = (2*pi)/WX;
    double DKY = (2*pi)/WY;
    double DKXY = DKX*DKY;
    double Q;
    double AKBX=0.0;                         //K�s OF BRILLOUIN ZONE START AT ZERO
    double AKBY=0.0;
    double AKX0 = AKBX*DKX/2.0;
    double AKY0 = AKBY*DKY/2.0;

    dx=WX/(nx);
    dy=WY/(ny);

    for(i=0;i<nx;i++){
    Q = DKX*i;
    if(i>nx/2)
     Q=DKX*(i-nx);

    Kx[i]=Q+AKX0;
    }

    Q=0.0;

    for(i=0;i<ny;i++){
    Q = DKY*i;
    if(i>ny/2)
     Q=DKY*(i-ny);
    Ky[i]=Q+AKY0;
    }

    //=====DISCRETIZAÇÃO DO ESPAÇO RECIPROCO========

//====================================================Constantes, declaração de funções de onda, vetores======================================================



    wavei(PsiU.arr,PsiD.arr,X,Y,x0,y0,xsig,ysig); //Criando pacote de onda inicial

    printf("Energy=%0.3f eV\n",E0);

    Pot(V0,X,Y,V);                        //Preenchendo matriz de Potencial V

    //================================================Começo do calculo de evolução temporal==================================================================
    for(it=0;it<nt;it++){

   //PsiU.printvalues();
    FFT(PsiU.arr,NN,ndim,-1,nx,ny);
//PsiU.printvalues();

    }
	//================================================fim do calculo de evolução temporal==================================================================

//Escrita no arquivo Prob_end.dat a densidade de probabilidade de |psi>_{t_0 + nt\Delta t}





return 0;
}
//fim da função principal


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><ROUTINAS><><><><><><><><><><><><><><><><><><><><<><><<><><><><><><><><><><><<><><><><><><><><><><><






//========================================================================================================================
//                                                                                    !Criação do pacote de onda inicial
//========================================================================================================================

void wavei(wf PSIA,wf PSIB,double *X,double *Y,double x0, double y0, double xsig, double ysig){

int i,j;
double theta,akx0,aky0,ASUM=0.0,AINV;
theta= (double)(phi*pi)/180.0;


//SENTIDO POSITIVO DE Y
akx0 =  cos(theta)*ak0;
aky0 =  sin(theta)*ak0;

for(i=1;i<nx;i++){
for(j=1;j<ny;j++){
PSIA[i][j] = cexp(-I*eta*theta/2)*cexp(I*akx0*X[i])*cexp(I*aky0*Y[j])*(com)exp(-(pow((X[i]-x0),2)/(2.0*pow(xsig,2)))-(pow((Y[j]-y0),2)/(2.0*pow(ysig,2)))); //pacote de onda inicial é uma Gaussiana
PSIB[i][j] = eta*cexp(I*eta*theta/2)*cexp(I*akx0*X[i])*cexp(I*aky0*Y[j])*(com)exp(-(pow((X[i]-x0),2)/(2.0*pow(xsig,2)))-(pow((Y[j]-y0),2)/(2.0*pow(ysig,2)))); //pacote de onda inicial é uma Gaussiana
}
}


//Normalização
for(i=1;i<nx;i++)
for(j=1;j<ny;j++)
ASUM=ASUM+pow(cabs(PSIA[i][j]),2)+pow(cabs(PSIB[i][j]),2);
ASUM=sqrt(ASUM*dx*dy);
if(ASUM==0.0)
AINV=0.0;
else
AINV=1/ASUM;
for(i=1;i<nx;i++)
for(j=1;j<ny;j++){
PSIA[i][j]=PSIA[i][j]*AINV;
PSIB[i][j]=PSIB[i][j]*AINV;
}
//Normalização


FILE* psi0=fopen("psi0.dat","w");
//Escrita no arquivo psi0.dat a densidade de probabilidade de |psi>_{t_0}
for(i=1;i<nx;i++)
for(j=1;j<ny;j++){
fprintf(psi0,"%.4f  \t  %.4f  \t  %.4f \n",X[i], Y[j],(pow(cabs(PSIA[i][j]),2)+pow(cabs(PSIB[i][j]),2)));
}

fclose(psi0);

return;
}

void Norm(cmatrix PSIA,cmatrix PSIB){

int i,j;
double ASUM,AINV;

for(i=1;i<nx;i++)
for(j=1;j<ny;j++)
ASUM=ASUM+pow(cabs(PSIA[i][j]),2)+pow(cabs(PSIB[i][j]),2);
ASUM=sqrt(ASUM*dx*dy);
if(ASUM==0.0)
AINV=0.0;
else
AINV=1/ASUM;
for(i=1;i<nx;i++)
for(j=1;j<ny;j++){
PSIA[i][j]=PSIA[i][j]*AINV;
PSIB[i][j]=PSIB[i][j]*AINV;
}



return;
}



//=================================================================================================================================
//                                                                                    !Pseudo Campo magnetico
//=================================================================================================================================
//void Pseudo_magnetic_field(



//=================================================================================================================================
//                                                                                    !Pseudo potencial escalar
//=================================================================================================================================

void Pot(double V0,double* X,double* Y, matrix V){

	double WIDTH=width/a0;
    int i,j;
    FILE *POT=fopen("POT.dat","w");//Valores do potencial

    for(i=0;i<nx;i++)
    for(j=0;j<ny;j++){
      V[i][j]=0.0;
  
    } 

	for(i=0;i<nx;i++)
    for(j=0;j<ny;j++)
    fprintf(POT,"%0.3f  \t  %0.3f  \t  %0.3f \n",X[i]*a0/sig,Y[j]*a0/sig,V[i][j]);

    fclose(POT);
return;
}

//===============================================================================================================================
//                                                                                    !<x> e <y>
//=================================================================================================================================

void centermass(wf CPSIA,wf CPSIB,double *X,double *Y,double *XCM,double *XCMr,double *XCMt,double *YCM,double *YCMr,double *YCMt){

	  *XCM = 0.0; //! Calculando a média de x (<x>) no espaço real
      *XCMr = 0.0;
      *XCMt = 0.0;
      *YCM = 0.0; //! Calculando a média de y (<y>) no espaço real
      *YCMr = 0.0;
      *YCMt = 0.0;
      int i,j;
      double dx=WX/(nx), dy=WY/(ny);

	  for(i=1;i<nx;i++)
	  for(j=1;j<ny;j++){
        *XCM = *XCM + creal(conj(CPSIA[i][j])*X[i]*CPSIA[i][j] + conj(CPSIB[i][j])*X[i]*CPSIB[i][j]);
        *YCM = *YCM + creal(conj(CPSIA[i][j])*Y[j]*CPSIA[i][j] + conj(CPSIB[i][j])*Y[j]*CPSIB[i][j]);
        if(ipropagation==0){

          if(X[i]<0.0)
            *XCMr = *XCMr + creal(conj(CPSIA[i][j])*X[i]*CPSIA[i][j] + conj(CPSIB[i][j])*X[i]*CPSIB[i][j]);

          if(X[i]>=0.0)
            *XCMt = *XCMt + creal(conj(CPSIA[i][j])*X[i]*CPSIA[i][j] + conj(CPSIB[i][j])*X[i]*CPSIB[i][j]);

         }
        if(ipropagation==1){

          if(Y[j]<0.0)
            *YCMr = *YCMr + creal(conj(CPSIA[i][j])*Y[j]*CPSIA[i][j] + conj(CPSIB[i][j])*Y[j]*CPSIB[i][j]);

          if(Y[j]>=0.0)
            *YCMt = *YCMt + creal(conj(CPSIA[i][j])*Y[j]*CPSIA[i][j] + conj(CPSIB[i][j])*Y[j]*CPSIB[i][j]);

	     }
       }

      *XCM = (double)(*XCM)*dx*dy;
      *XCMr = (double)(*XCMr)*dx*dy;
      *XCMt = (double)(*XCMt)*dx*dy;
      *YCM = (double)(*YCM)*dx*dy;
      *YCMr = (double)(*YCMr)*dx*dy;
      *YCMt = (double)(*YCMt)*dx*dy;

}

//===============================================================================================================================
//                                                                                    !função de probabilidade
//=================================================================================================================================

void prob(double*X,double*Y,wf CPSIA,wf CPSIB,double*PBEF,double*PINT,double*PAFT){

     double WIDTH=width/a0;
     double SUMB = 0.0;
     double SUMI = 0.0;
     double SUMA = 0.0;


	 int i,j;

      for(i=1;i<nx;i++)
	  for(j=1;j<ny;j++){

         if(ipropagation==0){

          if(X[i]<0.0)
   	            SUMB = SUMB + creal(conj(CPSIA[i][j])*CPSIA[i][j] + conj(CPSIB[i][j])*CPSIB[i][j]);

          if(X[i]>=0.0 && X[i]<=WIDTH)
	            SUMI = SUMI + creal(conj(CPSIA[i][j])*CPSIA[i][j] + conj(CPSIB[i][j])*CPSIB[i][j]);

          if(X[i]>=0.0) //           if(X(i).gt.WIDTH)then
 	            SUMA = SUMA + creal(conj(CPSIA[i][j])*CPSIA[i][j] + conj(CPSIB[i][j])*CPSIB[i][j]);
	    }
         if(ipropagation==1){

           if(Y[j]<0.0)
   	            SUMB = SUMB + creal(conj(CPSIA[i][j])*CPSIA[i][j] + conj(CPSIB[i][j])*CPSIB[i][j]);

           if(Y[j]>=0.0 && Y[j]<=WIDTH)
	            SUMI = SUMI + creal(conj(CPSIA[i][j])*CPSIA[i][j] + conj(CPSIB[i][j])*CPSIB[i][j]);

           if(Y[j]>=0.0) //           if(Y(j).gt.WIDTH)then
 	            SUMA = SUMA + creal(conj(CPSIA[i][j])*CPSIA[i][j] + conj(CPSIB[i][j])*CPSIB[i][j]);
        }
     }

      *PBEF = SUMB*dx*dy;
      *PINT = SUMI*dx*dy;
      *PAFT = SUMA*dx*dy;

}



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><ROUTINAS><><><><><><><><><><><><><><><><><><><><<><><<><><><><><>
