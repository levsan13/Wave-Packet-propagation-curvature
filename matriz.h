#ifndef MATRIZ_H_INCLUDED
#define MATRIZ_H_INCLUDED


/************************************************************************
 *
 * Biblioteca para Calculos de matriz no corpo dos Reais e complexos em C.
 *
 * Autor: Sergio Levy Nobre dos Santos.
 * Instuição: Universidade Federal do Ceará - UFC
 * Email para contato: levy@fisica.ufc.br  
 *
 * U L T I M A   A T U A L I Z A Ç Ã O ! ! !
 *
 *
 *
 ******************************************************************************/

//========================================================================================================================
//                                                                                    !INCLUDES (Bibliotecas C)
//========================================================================================================================
#include<stdlib.h>
#include<stdio.h>
#include<complex.h>
#include<math.h>

//========================================================================================================================
//                                                                                    !Definições
//========================================================================================================================
#define RANDOM 1
#define CONSOLE 0

//========================================================================================================================
//                                                                                    !Definições de tipos
//========================================================================================================================
typedef double _Complex com;

typedef double** matrix; //double matrix
typedef com** cmatrix;   //complex matrix

typedef double* vec;    //double vector
typedef com* cvec;      //complex vector


//========================================================================================================================
//                                                                                    !Funções
//========================================================================================================================
double det(double **m, int n);

double **matrizMenor(double **m, double **m2, int c,int l, int n);

void **allocmatrix(int lin,int col,size_t size_of_type);//Cria uma matriz de ordem n alocando seus espaÃ§os

double **inv(double **m, int s);//by Levy

double **Slinear(double **m, vec b, int s);//by Levy

void printmatrix(double **m,int lin, int col);//BY LEVY

void entermatrix(double **m, int lin, int col,int _type);//by levy

double **prodmatrix(double **m1, int lin1,int col1, double **m2, int lin2, int col2);//by levy

double **sm(double **m1, double **m2, int m, int n);// by levy

double **dm(double **m1, double **m2, int m, int n);// by levy

void **allocvetor(int s,size_t size_of_type);

double **In(int n); //matriz indentidade de n-esima ordem By levy

double **diagmatrix(double **m, int s);

double **transposta(double **p, int lin,int col);

void mrerror(const char *error_text, int _EXIT);

cmatrix csm(cmatrix m1, cmatrix m2, int m, int n);

cmatrix cdm(cmatrix m1, cmatrix m2, int m, int n);

cmatrix cprodmatrix(cmatrix m1, int lin1,int col1, cmatrix m2, int lin2, int col2); //by levy

cmatrix ctrans(cmatrix m, int lin, int col);

cmatrix mconj(cmatrix m, int lin, int col);

cmatrix cinv(cmatrix m, int dimension);

void cprintmatrix(cmatrix m,int lin, int col);

void FREE(void **pointer,int len);

double norm (double *x, int length);

double dot_product (double * x, double * y, int length); 

com max_min_value(cmatrix m,int nx,int ny,int _type);


//========================================================================================================================
//                                                                                    !Code of functions
//========================================================================================================================


//==============Menssagem de erro =====================
void mrerror(const char *error_text,int _EXIT)
{
	fprintf(stderr,"matriz.h run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	if(_EXIT==1){
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
    }
}


//================AlocaÃ§Ã£o DinÃ¢mica==================
//AlocaÃ§Ã£o de matrizes
void **allocmatrix(int lin,int col,size_t size_of_type)//Cria uma matriz de ordem n alocando seus espaÃ§os
{
int i;
void **m;
m=(void**)malloc(sizeof(void*)*lin);
if(m == NULL) mrerror("erro allocmatrix() 1: malloc()",1);
for(i=0;i<lin;i++){
*(m+i)=malloc(size_of_type*(col));
if(*(m+i) == NULL) mrerror("erro allocmatrix() 2: malloc()",1);
}
return m;
}
//AlocaÃ§Ã£o de matrizes
//AlocaÃ§Ã£o de vetores
void **allocvetor(int s,size_t size_of_type)
{
return allocmatrix(1,s,size_of_type);
}
//AlocaÃ§Ã£o de vetores
//================AlocaÃ§Ã£o DinÃ¢mica==================



//==============TRANSPOSTA=================
double **transposta(double **p, int lin,int col)
{
 if(p==NULL)
 mrerror("error transposta(): p is NULL...",1);
 int i,j;
 double **aux=(matrix)allocmatrix(col,lin,sizeof(double));
 for(j=0;j<col;j++)
 for(i=0;i<lin;i++)
 aux[j][i]=p[i][j];
 p=(matrix)allocmatrix(col,lin,sizeof(double));
 for(j=0;j<col;j++)
 for(i=0;i<lin;i++)
 p[j][i]=aux[j][i];
 free(aux);
 return p;
}
//==============TRANSPOSTA=================


//========================PARA O CALCULO DO DETERMINANTE============================
double det(double **m, int n)
{
	if(m==NULL)
	mrerror("error det(): m is NULL...",1);
	//Caso fundamental
	if(n == 1)
		return **m;

	//Criando a matriz do cofator
	int c;
        double cof, d=0, **m2 = (matrix)allocmatrix(n-1,n-1,sizeof(double));

	//Para cada coluna
	for(c=0; c<n; c++){
		//Calculando o cofator
		cof = (c%2?-1:1) * det(matrizMenor(m, m2, c,0, n), n-1);
		//Somando a parcela ao determinante final
		d += m[0][c]*cof;
		}
return d;
}
double **matrizMenor(double **m, double **m2, int c,int l, int n)
{
	if(m==NULL||m2==NULL)
	mrerror("error matrizMenor(): m or m2 is NULL...",1);

	int i, j;
        double *p2,**p=m2;
	for(i = 0; i < n; i++){//exclui a primeira linha da matriz com i=1
                if(i==l) continue;//BY LEVY
                else
                {
		p2 = *p++;//pega o endereÃ§o da matriz m2 das colunas
		for(j = 0; j < n; j++)
		if(j==c) continue;
		else *p2++ = m[i][j];//copia todos os termos da j-esima coluna da matriz com exessÃ£o da coluna c
                }
	}
	return m2;
}
//========================PARA O CALCULO DO DETERMINANTE============================


//========================CALCULO DE INVERSA========================================
double **inv(double **m, int s)
{
if(det(m,s)==0){
mrerror("error inv(): m is singular matrix...",0);
return NULL;
}
int i=0,k,j,cont=1;
double **A=(matrix)allocmatrix(s,s,sizeof(double));
double **M=(matrix)allocmatrix(s,s,sizeof(double));
double **in=(matrix)allocmatrix(s,s,sizeof(double));
double d;
d=det(m,s);
//Calculo dos cofatores
for(k=0;k<s;k++)
for(j=0;j<s;j++)
{
if((k+j)%2==0)
A[k][j]=det(matrizMenor(m,M,j,k,s),s-1);
else
A[k][j]=(-1)*det(matrizMenor(m,M,j,k,s),s-1);
}
//transposta
A=transposta(A,s,s);
for(k=0;k<s;k++)
for(j=0;j<s;j++)
{
in[k][j]=A[k][j]*(1/d);
}
free(M);
free(A);
return in;
}
//========================CALCULO DE INVERSA========================================


//========================IMPRIMIR E ENTRAR COM MATRIZ MXN=======================================
void printmatrix(double **m,int lin, int col)
{
if(m==NULL)
mrerror("error printmatrix(): m is NULL...",1);

int i,j;
for(i=0;i<lin;i++){
printf("|");
for(j=0;j<col;j++)
printf(" %0.3lf ",m[i][j]);
printf(" |\n");
}
printf("\n\n");
return;
}
void entermatrix(double **m, int lin, int col, int _type )
{
if(m==NULL)
mrerror("error entermatrix(): m is NULL...",1);
int n,i,j;
if(_type==RANDOM)
{
printf("intervalo:");
scanf("%d",&n);
for(i=0;i<lin;i++)
for(j=0;j<col;j++){
m[i][j]=(rand()%n)+1;
}
}
else{
for(i=0;i<lin;i++){
printf("\n==========linha-%d==========\n",i+1);
for(j=0;j<col;j++)
{
printf("termo[%d][%d]:",i+1,j+1);
scanf("%lf",&m[i][j]);
}
}
}
return;
}
//========================IMPRIMIR E ENTRAR COM MATRIZ MXN=======================================


//========================OPERAÃÃES COM MATRIZES MXN=============================================
//PRODUTO MATRICIAL
double **prodmatrix(double **m1, int lin1,int col1, double **m2, int lin2, int col2)
{
if(col1!=lin2){
mrerror("error prodmatrix(): col1 != lin2...",0);
return NULL;
}
if(m1==NULL||m2==NULL){
mrerror("error prodmatrix(): m1 or m2 is NULL...",1);
}
int i=0,j=0,k,n=col1;
double **m=(matrix)allocmatrix(n,n,sizeof(double)),A=0;
//erro
for(i=0;i<lin1;i++)
{
        for(j=0;j<col2;j++)
        {
                for(k=0;k<n;k++)
                 {
                 A += m1[i][k] * m2[k][j];
                 }
                 m[i][j]=A;
          A=0;
         }
}
return m;
}
//PRODUTO MATRICIAL
//resoluÃ§Ã£o de sistemas lineares
double **Slinear(double **m, vec b, int s)//sistemas lineares nxn
{
matrix aux_v=(matrix)allocvetor(s,sizeof(double));
aux_v=transposta(aux_v,s,1);
aux_v[0]=b;
aux_v=transposta(aux_v,1,s);
return prodmatrix(inv(m,s),s,s,aux_v,s,1);
}
//resoluÃ§Ã£o de sistemas lineares
//soma e diferenÃ§a
double **sm(double **m1, double **m2, int m, int n)
{
int i,j;
double **A=(matrix)allocmatrix(m,n,sizeof(double));
for(i=0;i<m;i++)
for(j=0;j<n;j++)
A[i][j]=m1[i][j]+m2[i][j];
return A;
}
double **dm(double **m1, double **m2, int m, int n)
{
int i,j;
double **A=(matrix)allocmatrix(m,n,sizeof(double));
for(i=0;i<m;i++)
for(j=0;j<n;j++)
A[i][j]=m1[i][j]-m2[i][j];
return A;
}


//========================================================================================================================
//                                                                                    !Matrizes complexas
//========================================================================================================================

//Soma e diferença

cmatrix csm(cmatrix m1, cmatrix m2, int m, int n)
{
int i,j;
cmatrix A=(cmatrix)allocmatrix(m,n,sizeof(com));
for(i=0;i<m;i++)
for(j=0;j<n;j++)
A[i][j]=m1[i][j]+m2[i][j];
return A;
}


cmatrix cdm(cmatrix m1, cmatrix m2, int m, int n)
{
int i,j;
cmatrix A=(cmatrix)allocmatrix(m,n,sizeof(com));
for(i=0;i<m;i++)
for(j=0;j<n;j++)
A[i][j]=m1[i][j]-m2[i][j];
return A;
}

//Produto de matrizes complexas

cmatrix cprodmatrix(cmatrix m1, int lin1,int col1, cmatrix m2, int lin2, int col2)
{
if(col1!=lin2){
mrerror("error cprodmatrix(): col1 != lin2...",0);
return NULL;
}
if(m1==NULL||m2==NULL){
mrerror("error cprodmatrix(): m1 or m2 is NULL...",1);
}
int i=0,j=0,k,n=col1;
cmatrix m=(cmatrix)allocmatrix(n,n,sizeof(com));
com A=0;
//erro
for(i=0;i<lin1;i++)
{
        for(j=0;j<col2;j++)
        {
                for(k=0;k<n;k++)
                 {
                 A += m1[i][k] * m2[k][j];
                 }
                 m[i][j]=A;
          A=0;
         }
}
return m;
}

//tranposta e conjugado

cmatrix ctrans(cmatrix m, int lin, int col){

 if(m==NULL)
 mrerror("error ctrans(): m is NULL...",1);
 int i,j;
 cmatrix aux=(cmatrix)allocmatrix(col,lin,sizeof(com));
 for(j=0;j<col;j++)
 for(i=0;i<lin;i++)
 aux[j][i]=m[i][j];
 m=(cmatrix)allocmatrix(col,lin,sizeof(com));
 for(j=0;j<col;j++)
 for(i=0;i<lin;i++)
 m[j][i]=aux[j][i];
 free(aux);
 return m;


}


cmatrix mconj(cmatrix m, int lin, int col){

 if(m==NULL)
 mrerror("error mconj(): m is NULL...",1);
 cmatrix aux=(cmatrix)allocmatrix(col,lin,sizeof(com));
 int i,j;
 for(i=0;i<lin;i++)
 for(j=0;j<col;j++)
 aux[i][j]=conj(m[i][j]);
 return aux;

}



void cprintmatrix(cmatrix m,int lin, int col)
{
if(m==NULL)
mrerror("error cprintmatrix(): m is NULL...",1);

printf("\n\n\n----- Parte Real --------\n");
int i,j;
for(i=0;i<lin;i++){
printf("|");
for(j=0;j<col;j++)
printf(" %0.3lf ",creal(m[i][j]));
printf(" |\n");
}

printf("\n");

printf("\n----- Parte Imaginaria --------\n");
for(i=0;i<lin;i++){
printf("|");
for(j=0;j<col;j++)
printf(" %0.3lf ",cimag(m[i][j]));
printf(" |\n");
}

printf("\n\n");

return;
}


//matriz indentidade de n-esima ordem
double **In(int n){

	matrix in=(matrix)allocmatrix(n,n,sizeof(double));
	int i,j;
	for(i=0;i<n;i++)
    for(j=0;j<n;j++){
    	if(i==j)
    	in[i][j]=1.0;
    	else
    	in[i][j]=0.0;
	}

return in;
}

void FREE(void **pointer,int len){ // free memory for array bidemensional

int i,j;
for(i=0;i<len;i++)
free(pointer[i]);

return;

}


double norm (double *x, int length) {
    int i, length5;
    double a, sum = 0;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        sum += x[i] * x[i];
    }
    for(; i < length; i += 5) {
        sum += x[i] * x[i] + x[i + 1] * x[i + 1] + x[i + 2] * x[i + 2]
                           + x[i + 3] * x[i + 3] + x[i + 4] * x[i + 4];
    }

    return sqrt(sum);
}

double dot_product (double * x, double * y, int length) {
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        sum += x[i] * y[i];
    }
    for(; i < length; i += 5) {
        sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
                           + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
    }

    return sum;
}

com max_min_value(cmatrix m,int nx,int ny,int _type) // 1 for max and 0 for min 
{
	com max,min;
	int i,j;
	max=m[0][0];
	min=m[0][0];
	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	{
		if(cabs(m[i][j])>cabs(max))
		max=m[i][j];
		if(cabs(m[i][j])<cabs(min))
		min=m[i][j];
	}
	if(_type==1)
	return max;
	if(_type==0)
	return min;
}



#ifdef __cplusplus

extern "C++"{

//========================================================================================================================
//                                                                                    !Classe de matrizes e vetores
//========================================================================================================================

class Ve{
	public:

		int dimension;
		vec v;

		void get_info_vector(int _dimension){
			dimension=_dimension;
			v=(double*)malloc(dimension*sizeof(double));
			if(v==NULL)
			mrerror("erro Ve::get_info_vector(): malloc() return's NULL...",1);
			return;
		}
		void entervalues(int _type){
			if(_type==RANDOM){
				int i,j;
				printf("intervalo:");
				scanf("%d",&j);
				for(i=0;i<dimension;i++)
				v[i]=(rand()%j)+1;
			}
			if(_type==CONSOLE){
				int i;
				for(i=0;i<dimension;i++){
					printf("%d termo:",i+1);
					scanf("%lf",&v[i]);
				}
				return;
			}
		}
		void printvalues(){
			int i;
			for(i=0;i<dimension;i++){
				if(i==0)
				printf("| ");
				if(i+1==dimension){
				printf(" %0.3f ",v[i]);
				printf(" |\n");

			    }
				else
				printf(" %0.3f ",v[i]);
		   }
		   return;
	    }
	    double norma(void){return norm(v,dimension);}
        void free_memory(){
        free((vec)v);
        return;
        }

};

class M{
	public:
		int lin,col;  //linhas e colunas
	    double** arr;   //Matriz armazenada
	    void get_info_matrix(int _l,int _c){
		    lin=_l;
		    col=_c;
		    arr=(matrix)allocmatrix(lin,col,sizeof(double));
		    return;
	    }
	    void entervalues(int _type){
		    entermatrix(arr,lin,col,_type);
		    return;
	    }
	    void printvalues(void){
		    printmatrix(arr,lin,col);
		    return;
	    }

	    Ve line(int _line){
	    	Ve aux;
	    	aux.get_info_vector(lin);
	    	aux.v=arr[_line-1];
	    	return aux;
		}

		double Det(void){

		    if(lin!=col)
		    mrerror("Error M::Det() : lin!=col.....",1);
		    else
		    return det(arr,col);
		}
        void free_memory(){
        int i,j;

        for(i=0;i<lin;i++)
        while(arr[i]!=NULL)
        free((vec)arr[i]);

        return;
        }
};

class M_complex{
	public:
		int lin,col;  //linhas e colunas
	    cmatrix arr;   //Matriz armazenada
	    void get_info_matrix(int _l,int _c){
		    lin=_l;
		    col=_c;
		    arr=(cmatrix)allocmatrix(lin,col,sizeof(com));
		    return;
	    }
	    void printvalues(void){
	    	cprintmatrix(arr,lin,col);
	    	return;
		}
        void free_memory(){
        int i,j;
        for(i=0;i<lin;i++)
        while(arr[i]!=NULL)
        free((cvec)arr[i]);

        return;
        }

};

class Ve_complex{
	public:

		int dimension;
		com* v;

		void get_info_vector(int _dimension){
			dimension=_dimension;
			v=(com*)malloc(dimension*sizeof(com));
			if(v==NULL)
			mrerror("erro Ve_complex::get_info_vector(): malloc() return's NULL...",1);
			return;
		}

		void printvalues(){
			int i;
			printf("\n\n-------Parte Real ---------\n");
			for(i=0;i<dimension;i++){
				if(i==0)
				printf("| ");
				if(i+1==dimension){
				printf(" %0.3f ",creal(v[i]));
				printf(" |\n");

			    }
				else
				printf(" %0.3f ",creal(v[i]));
		   }
		   printf("\n");
		   printf("-------Parte Imaginaria ---------\n");
		   for(i=0;i<dimension;i++){
				if(i==0)
				printf("| ");
				if(i+1==dimension){
				printf(" %0.3f ",cimag(v[i]));
				printf(" |\n");

			    }
				else
				printf(" %0.3f ",cimag(v[i]));
		   }
		   printf("\n");
		   return;
	    }
        void free_memory(){
        free((cvec)v);
        return;
        }
};

typedef class M Marray;
typedef class Ve Varray;
typedef class M_complex MCarray;
typedef class Ve_complex VCarray;


Varray e(int _i,int _dimension){Varray aux; aux.get_info_vector( _dimension);aux.v=In( _dimension)[_i-1];return aux;} //_i-esimo versor da base canonico

//===================================================================================================================================
//                                                                                    !Definição das operações com matrizes e vetores
//===================================================================================================================================

Marray operator+(const Marray m1, const Marray m2){ //soma de matrizes
	Marray aux;
	aux.get_info_matrix(m1.lin,m1.col);
	aux.arr=sm(m1.arr,m2.arr,m1.lin,m1.col);
	return aux;
}

Marray operator-(const Marray m1, const Marray m2){ //subtração de matrizes
	Marray aux;
	aux.get_info_matrix(m1.lin,m1.col);
	aux.arr=dm(m1.arr,m2.arr,m1.lin,m1.col);
	return aux;
}

Marray operator*(const Marray m1, const Marray m2){ //produto de matrizes
	Marray aux;
	aux.get_info_matrix(m1.lin,m2.col);
	aux.arr=prodmatrix(m1.arr,m1.lin,m1.col,m2.arr,m2.lin,m2.col);
	return aux;
}

Marray operator*(const double c, const Marray m2){ //produto de matriz por constante
	Marray aux;
	aux.get_info_matrix(m2.lin,m2.col);
	int i,j;
	for(i=0;i<m2.lin;i++)
	for(j=0;j<m2.col;j++)
	aux.arr[i][j]=c*m2.arr[i][j];
	return aux;
}

Marray operator!(const Marray m){                //matriz transposta
	Marray aux;
	aux.get_info_matrix(m.col,m.lin);
	aux.arr=transposta(m.arr,m.lin,m.col);
	return aux;
}

Marray operator~(const Marray m){              //inversa de matriz
	Marray aux;
	aux.get_info_matrix(m.lin,m.col);
	aux.arr=inv(m.arr,m.lin);
	return aux;
}

//========================================================================================================================
//                                                                                    !Definição das operações com vetores
//========================================================================================================================

Varray Vnulo(int _dimension){
	int i;
	Varray a;
	a.get_info_vector(_dimension);
	for(i=0;i<_dimension;i++)
	a.v[i]=0.0;
	return a;
}

Varray operator+(const Varray v1,const Varray v2){  //soma de vetores
	Varray aux;
	aux.get_info_vector(v1.dimension);
	int i;
	for(i=0;i<v1.dimension;i++)
	aux.v[i]=v1.v[i]+v2.v[i];
	return aux;
}

Varray operator-(const Varray v1,const Varray v2){  //subtração de vetores
 	Varray aux;
	aux.get_info_vector(v1.dimension);
	int i;
	for(i=0;i<v1.dimension;i++)
	aux.v[i]=v1.v[i]-v2.v[i];
	return aux;
}

double operator*(const Varray v1,const Varray v2){  //produto escalar
	return dot_product(v1.v,v2.v,v1.dimension);
}

Varray operator*(const double c, const Varray ve){  //produto por uma constante
	Varray aux;
	aux.get_info_vector(ve.dimension);
	int i;
	for(i=0;i<ve.dimension;i++)
	aux.v[i]=c*ve.v[i];
	return aux;
}

Varray operator^(const Varray v1,const Varray v2){ //produto vetorial
	Varray aux;
	aux.get_info_vector(v1.dimension);
	if(aux.dimension!=3)
	mrerror("error vetorial product: dimension not is 3...",1);
	aux=((v1.v[1]*v2.v[2])-(v2.v[1]*v1.v[2]))*e(1,aux.dimension)-((v1.v[0]*v2.v[2])-(v2.v[0]*v1.v[2]))*e(2,aux.dimension)+((v1.v[0]*v2.v[1])-(v2.v[0]*v1.v[1]))*e(3,aux.dimension);
	return aux;
}

Varray operator*(const Marray m,const Varray v){   //produto de matriz por vetor
	Varray aux1;
	matrix aux2;
	matrix aux_v=(matrix)allocvetor(v.dimension,sizeof(double));
    aux_v[0]=v.v;
    aux_v=transposta(aux_v,1,v.dimension);
    aux2=transposta(prodmatrix(m.arr,v.dimension,v.dimension,aux_v,v.dimension,1),v.dimension,1);
    aux1.get_info_vector(v.dimension);
    aux1.v=aux2[0];
    free(aux2);
    free(aux_v);
    return aux1;
}

VCarray operator*(const MCarray m,const Varray v){   //produto de matriz complexa por vetor
	VCarray aux1;
	cmatrix aux2;
	matrix aux_v=(matrix)allocvetor(v.dimension,sizeof(double));
    aux_v[0]=v.v;
    aux_v=transposta(aux_v,1,v.dimension);
    aux2=ctrans(cprodmatrix(m.arr,v.dimension,v.dimension,(cmatrix)aux_v,v.dimension,1),v.dimension,1);
    aux1.get_info_vector(v.dimension);
    aux1.v=aux2[0];
    free(aux2);
    free(aux_v);
    return aux1;
}

VCarray operator*(const MCarray m,const VCarray v){   //produto de matriz complexa por vetor complexo
	VCarray aux1;
	cmatrix aux2;
	cmatrix aux_v=(cmatrix)allocvetor(v.dimension,sizeof(com));
    aux_v[0]=v.v;
    aux_v=ctrans(aux_v,1,v.dimension);
    aux2=ctrans(cprodmatrix(m.arr,v.dimension,v.dimension,aux_v,v.dimension,1),v.dimension,1);
    aux1.get_info_vector(v.dimension);
    aux1.v=aux2[0];
    free(aux2);
    free(aux_v);
    return aux1;
}

//===================================================================================================================================
//                                                                                    !Definição das operações com matrizes complexas
//===================================================================================================================================

MCarray operator+(const MCarray m1, const MCarray m2){ //soma de matrizes
	MCarray aux;
	aux.get_info_matrix(m1.lin,m1.col);
	aux.arr=csm(m1.arr,m2.arr,m1.lin,m1.col);
	return aux;
}

MCarray operator-(const MCarray m1, const MCarray m2){ //subtração de matrizes
	MCarray aux;
	aux.get_info_matrix(m1.lin,m1.col);
	aux.arr=cdm(m1.arr,m2.arr,m1.lin,m1.col);
	return aux;
}

MCarray operator*(const com c, const MCarray m2){ //produto de matriz por constante
	MCarray aux;
	aux.get_info_matrix(m2.lin,m2.col);
	int i,j;
	for(i=0;i<m2.lin;i++)
	for(j=0;j<m2.col;j++)
	aux.arr[i][j]=c*m2.arr[i][j];
	return aux;
}
MCarray operator+(const MCarray m1, const Marray m2){ //soma de matriz complexa com matriz real
    MCarray aux;
    aux.get_info_matrix(m2.lin,m2.col);
    int i,j;
    for(i=0;i<m2.lin;i++)
	for(j=0;j<m2.col;j++)
    aux.arr=csm(m1.arr,(cmatrix)m2.arr,m1.lin,m1.col);
    return aux;
}
MCarray operator-(const MCarray m1, const Marray m2){  //diferença de matriz complexa com matriz real
    MCarray aux;
    aux.get_info_matrix(m2.lin,m2.col);
    int i,j;
    for(i=0;i<m2.lin;i++)
	for(j=0;j<m2.col;j++)
    aux.arr=cdm(m1.arr,(cmatrix)m2.arr,m1.lin,m1.col);
    return aux;
}
MCarray operator~(const MCarray m1){                   //Matriz conjugada
   MCarray aux;
   aux.get_info_matrix(m1.lin,m1.col);
   int i,j;
   aux.arr=mconj(m1.arr,m1.lin,m1.col);
   return aux;
}
MCarray operator!(const MCarray m1){                  //Transposta da matrix complexa
   MCarray aux;
   aux.get_info_matrix(m1.lin,m1.col);
   int i,j;
   aux.arr=ctrans(m1.arr,m1.lin,m1.col);
   return aux;
}
MCarray operator*(const MCarray m1, const MCarray m2){ //Produto de matrizes complexas
   MCarray aux;
   aux.get_info_matrix(m1.lin,m1.col);
   aux.arr=cprodmatrix(m1.arr,m1.lin,m1.col,m2.arr,m2.lin,m2.col);
   return aux;
}



//Matrizes de Pauli

MCarray sigma(int _ordem){
MCarray A;
int i,j;
A.get_info_matrix(2,2);
if(_ordem>3)
mrerror("error sigma(): _ordem>3",1);
if(_ordem==1){
A.arr[0][0]=0;A.arr[0][1]=1;A.arr[1][0]=1;A.arr[1][1]=0;
return A;
}

if(_ordem==2){
A.arr[0][0]=0;A.arr[0][1]=-I;A.arr[1][0]=I;A.arr[1][1]=0;
return A;
}

if(_ordem==3){
A.arr[0][0]=1;A.arr[0][1]=0;A.arr[1][0]=0;A.arr[1][1]=-1;
return A;
}
}

//Matriz nula

MCarray Nulo(int _lin,int _col){ 

MCarray A; 
int i,j;
A.get_info_matrix(_lin,_col);
for(i=0;i<_lin;i++)
for(j=0;j<_col;j++)
A.arr[i][j]=0.0;
return A;

}

//Matriz indentidade

MCarray Ind(int _dimension){
MCarray A;
int i,j;
A.get_info_matrix(_dimension,_dimension);
for(i=0;i<_dimension;i++)
for(j=0;j<_dimension;j++){
    	if(i==j)
    	A.arr[i][j]=1.0;
    	else
    	A.arr[i][j]=0.0;
}
return A;
}


//========================================================================================================================
//                                                                                    !Classe de espaço vetorial
//========================================================================================================================
//Espaço Vetorial no corpo dos Reais
class vs{
	public:
        int dimension;
	    char corpo;
		Marray base;
		void get_info_space(int _d){
			dimension=_d;
			base.get_info_matrix(_d,_d);
			return;
		}
		void enterbase(int _type){
			base.entervalues(_type);
			return;
		}
		void printbase(void){
			base.printvalues();
			return;
		}
		vec check_proj(vec v,int _PRINT){
			matrix aux_r=(matrix)allocmatrix(dimension,1,sizeof(double));
			aux_r=Slinear(base.arr,v,dimension);
			aux_r=transposta(aux_r,dimension,1);
			if(_PRINT==1)
			printmatrix(aux_r,1,dimension);
			return aux_r[0];
		}
		Marray base_change(Marray out_base,int _type,int _PRINT){   //0 for base to out_base
            Marray aux;
            aux.get_info_matrix(dimension,dimension);
		    if(_type==0){                                           //1 for out_base to base
			aux.arr = prodmatrix(inv(out_base.arr,dimension),dimension,dimension,base.arr,dimension,dimension);
		    if(_PRINT==1)
			printmatrix(aux.arr,dimension,dimension);
		    return aux;
		    }
		    if(_type==1){
		    aux.arr = prodmatrix(inv(base.arr,dimension),dimension,dimension,out_base.arr,dimension,dimension);
		    if(_PRINT==1)
			printmatrix(aux.arr,dimension,dimension);
		    return aux;
			}
	    }
};
typedef class vs space;
}
#endif

#endif // MATRIZ_H_INCLUDED
