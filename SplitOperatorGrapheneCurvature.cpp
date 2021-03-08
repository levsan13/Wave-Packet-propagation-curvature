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
#include"curve_FFT.h"//(by Levy) FFT, constantes e funções de curvatura



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


const int nt=350;                       //numero de evoluções temporais
const double cdt=1.0e-15;               //passo do tempo (dado em s)

const double l=0.142;                   //distancia entre os sitios A e B (dado em Angstron)
const int nx=1024;                      //Numero de pontos no eixo x
const int ny=1024;                      //Numero de pontos no eixo y
const int n=1024;
const double WX=400;                   //tamanho do  eixo x (dado em Angstron)
const double WY=400;                   //tamanho do  eixo y (dado em Angstron)
const int ndim=2;
double dx;                              //passo no eixo x
double dy;                              //passo no eixo y
const double width=0.0;

//caracteristicas do pacote de onda inicial
const double sig=30.0;                  //sigma (nm)
const double phi=0.0;                   // Angulo de k (Graus)
const double E0=0.156;                   //Energy (eV)
const int eta=1;                       // vale K ou K'

const int ipropagation=1;
const int imovie=1;

//plot B plotV
const int plotB=0;
const int plotV=0;





//========================================================================================================================
//                                                                                    !FunÃ§Ãµes/routinas
//========================================================================================================================




static void wavei(wf PSIA,wf PSIB,double *X,double *Y,double x0, double y0, double xsig, double ysig,double ak0); //Criação de pacote de onda inicial

static void Norm(cmatrix PSIA,cmatrix PSIB);

static void Pot(double V0,double* X,double* Y, matrix V);                                              //Preenchimeto da matriz de potencial

static void centermass(wf CPSIA,wf CPSIB,double *X,double *Y,double *XCM,double *XCMr,double *XCMt,double *YCM,double *YCMr,double *YCMt); // Calculo de <x> e <y>

static void prob(double*X,double*Y,wf CPSIA,wf CPSIB,double*PBEF,double*PINT,double*PAFT);                                                 //Probabilidade

static void pseudo_pot(double g,matrix V,matrix h,double *X,double *Y,int nx,int ny); // pseudo potencial escalar

static void pseudo_pot_vec(double beta,matrix A_x,matrix A_y,matrix A_z,matrix h,double *X,double *Y,int nx,int ny); // pseudo potencial vetor






//========================================================================================================================
//                                                                                    !FunÃ§Ã£o Principal
//========================================================================================================================

int main(){

//====================================================Constantes, declaração de funções de onda, vetores======================================================
    int i,j,it,kx,ky;
    int nnx=nx/2,nny=ny/2;


    //=====DISCRETIZAÇÃO DO ESPAÇO REAL=======

    double X[nx],Y[ny];
    dx=WX/(nx);
    dy=WY/(ny);
    for(i=0;i<nx;i++)
    X[i]=(double)(-WX/(2))+(double)(i)*dx;
    for(i=0;i<ny;i++)
    Y[i]=(double)(-WY/(2))+(double)(i)*dy;

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



    //WAVEFUNCTION
    MCarray PsiU,PsiD,PsiU_aux,PsiD_aux;
    PsiU.get_info_matrix(nx,ny);
    PsiD.get_info_matrix(nx,ny); // Psi UP and Down
    
    //WAVEFUNCTION


    /* - - - - - - - Curvature - - - - - - - -*/

     Marray h,V,AX,AY,AZ; // MATRIZ REAL
     
     //function of curvature
      
     /* - - - - - GAUSSIAN - - - - -*/
  
      double Xsig=50.0;
      double Ysig=50.0;
      double lenght=10.0;
      double x_0=0.0;
      double y_0=0.0;

     /* - - - - - PONTO DE CELA - - - - -*/
/*
      double C=0.0;
      double A_x=10.0;
      double A_y=10.0;
      double x_0=0.0;
      double y_0=0.0;  
*/

     h.arr=gaussian(Xsig,Ysig,lenght,x_0,y_0,X,Y,nx,ny);
     double beta=3.0;
     double g=3.0;
        

     //PSEUDO POTENTIAL SCALAR
     V.get_info_matrix(nx,ny);

     pseudo_pot(g,V.arr,h.arr,X,Y,nx,ny);

     //PSEUDO POTENTIAL VECTOR
     AX.get_info_matrix(nx,ny);  //eV*s/nm
     AY.get_info_matrix(nx,ny);
     AZ.get_info_matrix(nx,ny);
   
     pseudo_pot_vec(beta,AX.arr,AY.arr,AZ.arr,h.arr,X,Y,nx,ny);
        

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
    
    double A_aux,k_aux;
    Varray k_x,k_y;
        
    AX=(cdt*vf*eta/(2*hc))*AX;// (nm/eV*s)*(eV*s/nm) 
    AY=(cdt*vf*eta/(2*hc))*AY;
    AZ=(cdt*vf*eta/(2*hc))*AZ;
    
    k_x.get_info_vector(nx);
    k_y.get_info_vector(ny);
    
    for(i=0;i<nx;i++){
    k_x.v[i]=(cdt*vf)*Kx[i];
    k_y.v[i]=(cdt*vf)*Ky[i];
	}
	
    

    //Arquivos
    FILE *CM=fopen("centermass.dat","w");     //Valores medio da função de onda em função do tempo
    FILE *probe=fopen("Prob_end.dat","w");    //função de probalidade


    //========== Variaveis de Evolução Temporal =====

    //Outras Caracteristicas do pacote inicial gaussiano
   
    double x0 = -150.0;               //Posição x inicial do pacote de onda
    double y0 = 0.0;                  //Posição y inicial do pacote de onda
    double xsig = sig;                //largura do pacote de onda inicial
    double ysig = sig;                //largura do pacote de onda inicial
    double k0=E0/(hc*vf);             //k do pacote de onda em (nm^-1)

    //Outras Caracteristicas do pacote inicial gaussiano


//====================================================Constantes, declaração de funções de onda, vetores======================================================



    wavei(PsiU.arr,PsiD.arr,X,Y,x0,y0,xsig,ysig,k0); //Criando pacote de onda inicial

    printf("Energy=%0.6f eV \n k0=%0.6f nm^-1\n",E0,k0);
    
    printf("\n\nTime \t XCM \t YCM\n");
          
    PsiU_aux.get_info_matrix(nx,ny);
    PsiD_aux.get_info_matrix(nx,ny); // Psi UP and Down
    char buf[100];

//================================================Começo do calculo de evolução temporal==================================================================

    for(it=0;it<nt;it++){
    

        ti = (double)it*cdt*pow(10,15);
        


        //Multiplicando pelo pseudo potencial

        for(i=0;i<nx;i++)
        for(j=0;j<ny;j++){
        PsiU.arr[i][j]=PsiU.arr[i][j]*cexp(-I*cdt*V.arr[i][j]/(2*hc));
        PsiD.arr[i][j]=PsiD.arr[i][j]*cexp(-I*cdt*V.arr[i][j]/(2*hc));
        }


        //Fazendo o calculo com o pseudo potencial vetor
        for(i=0;i<nx;i++)
        for(j=0;j<ny;j++){
        	
        A_aux=sqrt(pow(AX.arr[i][j],2)+pow(AY.arr[i][j],2)+pow(AZ.arr[i][j],2));
        
        if(A_aux==0.0)
        A_aux=1e-10;
        

        
        PsiU_aux.arr[i][j]=(cos(A_aux)-I*(sin(A_aux)/A_aux)*AZ.arr[i][j])*PsiU.arr[i][j] + I*(sin(A_aux)/A_aux)*(eta*AX.arr[i][j] - I*AY.arr[i][j])*PsiD.arr[i][j];
        PsiD_aux.arr[i][j]=(cos(A_aux)+I*(sin(A_aux)/A_aux)*AZ.arr[i][j])*PsiD.arr[i][j] + I*(sin(A_aux)/A_aux)*(eta*AX.arr[i][j] + I*AY.arr[i][j])*PsiU.arr[i][j];
		}

        
        //Fazendo o calculo no espaço de momento
        
        //Fazendo a transformada DIRETA de fourier
        FFT(PsiU_aux.arr,NN,ndim,1,nx,ny);
        FFT(PsiD_aux.arr,NN,ndim,1,nx,ny);
        

        for(i=0;i<nx;i++)
        for(j=0;j<ny;j++){
        
        k_aux=sqrt(pow(k_x.v[i],2)+pow(k_y.v[j],2));
        
        if(k_aux==0.0)
        k_aux=1e-16;
        
        
        PsiU.arr[i][j]=cos(k_aux)*PsiU_aux.arr[i][j] + I*(sin(k_aux)/k_aux)*(eta*k_x.v[i] - I*k_y.v[j])*PsiD_aux.arr[i][j];
        PsiD.arr[i][j]=cos(k_aux)*PsiD_aux.arr[i][j] + I*(sin(k_aux)/k_aux)*(eta*k_x.v[i] + I*k_y.v[j])*PsiU_aux.arr[i][j];        
		
		}
		
		//Fazendo a transformada INVERSA de fourier
        FFT(PsiU.arr,NN,ndim,-1,nx,ny);
        FFT(PsiD.arr,NN,ndim,-1,nx,ny);

        

        //Fazendo o calculo com o pseudo potencial vetor novamente
        
        for(i=0;i<nx;i++)
        for(j=0;j<ny;j++){
        	
        A_aux=sqrt(pow(AX.arr[i][j],2)+pow(AY.arr[i][j],2)+pow(AZ.arr[i][j],2));
        
        if(A_aux==0)
        A_aux=1e-10;
        
        
        PsiU_aux.arr[i][j]=(cos(A_aux)-I*(sin(A_aux)/A_aux)*AZ.arr[i][j])*PsiU.arr[i][j] + I*(sin(A_aux)/A_aux)*(eta*AX.arr[i][j] - I*AY.arr[i][j])*PsiD.arr[i][j];
        PsiD_aux.arr[i][j]=(cos(A_aux)+I*(sin(A_aux)/A_aux)*AZ.arr[i][j])*PsiD.arr[i][j] + I*(sin(A_aux)/A_aux)*(eta*AX.arr[i][j] + I*AY.arr[i][j])*PsiU.arr[i][j];
		}

        // Multiplicando pelo potencial

        for(i=0;i<nx;i++)
        for(j=0;j<ny;j++){
        PsiU.arr[i][j]=PsiU_aux.arr[i][j]*cexp(-I*cdt*V.arr[i][j]/(2*hc));
        PsiD.arr[i][j]=PsiD_aux.arr[i][j]*cexp(-I*cdt*V.arr[i][j]/(2*hc));
        }

        // Calculando o centro de massa
        centermass(PsiU.arr,PsiD.arr,X,Y,&XCM,&XCMr,&XCMt,&YCM,&YCMr,&YCMt);
        fprintf(CM,"%0.3f \t %0.3f \t %0.3f\n",ti,XCM,YCM);
        
        
        // Plotando a densidade de probabilidade
         
        if(imovie==1){
        FILE *probe=fopen("Prob_end.dat","w");  //função de probalidade

        for(i=0;i<nx;i++)
        for(j=0;j<ny;j++)
        fprintf(probe,"%0.3f \t %0.3f \t %0.9f\n",X[i],Y[j],pow(cabs(PsiU.arr[i][j]),2)+pow(cabs(PsiD.arr[i][j]),2));

        fclose(probe);
         
        sprintf(buf,"python plot.py plotgraf Prob_end.dat  x y %d",it);

        system(buf);

        remove("Prob_end.dat");
        }




	    printf("%0.3f \t %0.3f \t %0.3f\n",ti,XCM,YCM);
      

    }

//================================================fim do calculo de evolução temporal==================================================================

//Escrita no arquivo Prob_end.dat a densidade de probabilidade de |psi>_{t_0 + nt\Delta t}

if(imovie==0){

for(i=0;i<nx;i++)
for(j=0;j<ny;j++)
fprintf(probe,"%0.3f \t %0.3f \t %0.9f\n",X[i],Y[j],pow(cabs(PsiU.arr[i][j]),2)+pow(cabs(PsiD.arr[i][j]),2));

fclose(probe);

system("python plot.py plotgraf Prob_end.dat  x y density");

remove("Prob_end.dat");
}


return 0;
}
//fim da função principal


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><ROUTINAS><><><><><><><><><><><><><><><><><><><><<><><<><><><><><><><><><><><<><><><><><><><><><><><






//========================================================================================================================
//                                                                                    !Criação do pacote de onda inicial
//========================================================================================================================

void wavei(wf PSIA,wf PSIB,double *X,double *Y,double x0, double y0, double xsig, double ysig,double ak0){

int i,j;
double theta,akx0,aky0,ASUM=0.0,AINV;
theta= (double)(phi*pi)/180.0;


//SENTIDO POSITIVO DE Y
akx0 =  cos(theta)*ak0;
aky0 =  sin(theta)*ak0;

for(i=0;i<nx;i++){
for(j=0;j<ny;j++){
PSIA[i][j] = cexp(-I*eta*theta/2)*cexp(I*akx0*X[i])*cexp(I*aky0*Y[j])*(com)exp(-(pow((X[i]-x0),2)/(2.0*pow(xsig,2)))-(pow((Y[j]-y0),2)/(2.0*pow(ysig,2)))); //pacote de onda inicial é uma Gaussiana
PSIB[i][j] = eta*cexp(I*eta*theta/2)*cexp(I*akx0*X[i])*cexp(I*aky0*Y[j])*(com)exp(-(pow((X[i]-x0),2)/(2.0*pow(xsig,2)))-(pow((Y[j]-y0),2)/(2.0*pow(ysig,2)))); //pacote de onda inicial é uma Gaussiana
}
}


//Normalização
for(i=0;i<nx;i++)
for(j=0;j<ny;j++)
ASUM=ASUM+pow(cabs(PSIA[i][j]),2)+pow(cabs(PSIB[i][j]),2);
ASUM=sqrt(ASUM*dx*dy);
if(ASUM==0.0)
AINV=0.0;
else
AINV=1/ASUM;
for(i=0;i<nx;i++)
for(j=0;j<ny;j++){
PSIA[i][j]=PSIA[i][j]*AINV;
PSIB[i][j]=PSIB[i][j]*AINV;
}
//Normalização


FILE* psi0=fopen("psi0.dat","w");
//Escrita no arquivo psi0.dat a densidade de probabilidade de |psi>_{t_0}
for(i=0;i<nx;i++)
for(j=0;j<ny;j++){
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
//                                                                                    !Pseudo potencial escalar
//=================================================================================================================================

void Pot(double V0,double* X,double* Y, matrix V){

	double WIDTH=width;
    int i,j;
    FILE *POT=fopen("POT.dat","w");//Valores do potencial

    for(i=0;i<nx;i++)
    for(j=0;j<ny;j++){
      V[i][j]=0.0;

    }

	for(i=0;i<nx;i++)
    for(j=0;j<ny;j++)
    fprintf(POT,"%0.3f  \t  %0.3f  \t  %0.3f \n",X[i]/sig,Y[j]/sig,V[i][j]);

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

     double WIDTH=width;
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

//===============================================================================================================================
//                                                                           !Strain Tensor and outher functions of curvature
//=================================================================================================================================


void pseudo_pot(double g,matrix V,matrix h,double *X,double *Y,int nx,int ny){

Marray f_x;
f_x.get_info_matrix(nx,ny);
f_x.arr=d_x(h,X,nx,ny);

Marray f_y;
f_y.get_info_matrix(nx,ny);
f_y.arr=d_y(h,Y,nx,ny); 


int i,j;
for(i=0;i<nx;i++)
for(j=0;j<ny;j++){
V[i][j]=(g/2)*(pow(f_x.arr[i][j],2)+pow(f_y.arr[i][j],2));
}

if(plotV==1){

FILE* plot=fopen("pseudo_scalar_potential.dat","w");

for(i=0;i<nx;i++)
for(j=0;j<ny;j++)
fprintf(plot,"%0.3f \t %0.3f \t %0.3f\n",X[i],Y[j],V[i][j]);

fclose(plot);

system("python plot.py plotgraf pseudo_scalar_potential.dat  x(nm) y(nm) V(eV)");

remove("pseudo_scalar_potential.dat");
}

return;
}

void pseudo_pot_vec(double beta,matrix A_x,matrix A_y,matrix A_z,matrix h,double *X,double *Y,int nx,int ny){

Marray f_x;
f_x.get_info_matrix(nx,ny);
f_x.arr=d_x(h,X,nx,ny);

Marray f_y;
f_y.get_info_matrix(nx,ny);
f_y.arr=d_y(h,Y,nx,ny); 


int i,j;
for(i=0;i<nx;i++)
for(j=0;j<ny;j++){
A_x[i][j]=(hc*beta/(4*l))*(pow(f_x.arr[i][j],2)-pow(f_y.arr[i][j],2)); // eV*s/nm
A_y[i][j]=-(hc*beta/(4*l))*2*((f_x.arr[i][j])*(f_y.arr[i][j]));
A_z[i][j]=0.0;
}

if(plotB==1){

matrix B=dm(d_x(A_y,X,nx,ny),d_y(A_x,Y,nx,ny),nx,ny);

FILE* plot=fopen("pseudo_potential_vector.dat","w");

for(i=0;i<nx;i++)
for(j=0;j<ny;j++)
fprintf(plot,"%0.3f \t %0.3f \t %0.3f\n",X[i],Y[j],(1/pow(10,-18))*(eta/e_c)*B[i][j]);

fclose(plot);

system("python plot.py plotgraf pseudo_potential_vector.dat  x y B");

remove("pseudo_potential_vector.dat");
}
return;

}



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><ROUTINAS><><><><><><><><><><><><><><><><><><><><<><><<><><><><><>
