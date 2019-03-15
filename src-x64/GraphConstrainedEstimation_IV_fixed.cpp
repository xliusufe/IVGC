// GraphConstrainedEstimation.cpp : Defines the entry point for the console application.
//

#include <string.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>

extern "C"
{
double sign(double x)
{
	if(x>0) return 1;
	else 
	{
		if(x<0) return -1;
	    else return 0;
	}
}

double max(double x,double y)
{
	if(x>y) return(x);
	else return(y);
}

double maxOFn(int n,double* x)
{
	int i;
	double re=x[0];
	for(i=1;i<n;i++) re=max(re,x[i]);
	return(re);
}

void MatrixAddition(int m,int n,double* A,double* B,double* C)
{
	int i,j;

	for(i=0;i<m;i++) for(j=0;j<n;j++) C[j*m+i]=A[j*m+i]+B[j*m+i];
}

void MatrixSubtraction(int m,int n,double* A,double* B,double* C)
{
	int i,j;

	for(i=0;i<m;i++) for(j=0;j<n;j++) C[j*m+i]=A[j*m+i]-B[j*m+i];
}

void MatrixMultiplication(int m,int n,int l,double* A,double* B,double* C)
{
	int i,j,k;
	double sum;

	for(i=0;i<m;i++) 
	{
		for(j=0;j<l;j++) 
		{
			sum=0;
			for(k=0;k<n;k++) sum+=A[k*m+i]*B[j*n+k];
			C[j*m+i]=sum;
		}
	}
}

void ScalarMatrix(double a,int m,int n,double* A,double* B)
{
	int i,j;

	for(i=0;i<m;i++) for(j=0;j<n;j++) B[j*m+i]=a*A[j*m+i];
}

double L2norm(int n,double* x)
{
	int i;
	double sum;

	sum=0;
	for(i=0;i<n;i++) sum+=x[i]*x[i];

	return (sqrt(sum));
}

void CorrX(int u,int n,int L,int J,double** Xsub,double** XsubT,double** corr_xx,double* zeroBeta)
{
	int i,j,k;

	MatrixMultiplication(L,n,L,XsubT[u],Xsub[u],corr_xx[u*J+u]);
	for(i=0;i<J;i++) 
	{
		if(zeroBeta[i]==0)
		{
			MatrixMultiplication(L,n,L,XsubT[u],Xsub[i],corr_xx[u*J+i]);
            for(j=0;j<L;j++) for(k=0;k<L;k++) corr_xx[i*J+u][k*L+j]=corr_xx[u*J+i][j*L+k];
		}
	}
}

double PlusFt(double x)
{
	if(x>0) return(x);
	else return(0);
}

void P1(int u,int n,int L,int J,double** beta,double** corr_xy,double** corr_xx,double* v1)
{
	int i;
	double* sum=new double[L];
	double* v=new double[L];

	for(i=0;i<L;i++) sum[i]=0;
	
	for(i=0;i<J;i++) 
	{
		if(L2norm(L,beta[i])>0)
		{
			MatrixMultiplication(L,L,1,corr_xx[u*J+i],beta[i],v);
			MatrixAddition(L,1,sum,v,sum);
		}
	}

	MatrixSubtraction(L,1,corr_xy[u],sum,v);
	ScalarMatrix(1.0/n,L,1,v,v);
	MatrixAddition(L,1,v,beta[u],v1);

	free(sum);
	free(v);
}

void P2(int u,int L,int J,double** beta,double* G,double lambda,double alpha,double* v3)
{
	int i;
	double* sum=new double[L];
	double* v=new double[L];

	for(i=0;i<L;i++) sum[i]=0;       

	for(i=0;i<J;i++)
	{
		if((G[u*J+i]!=0) & (L2norm(L,beta[i])>0)) 
		{
			ScalarMatrix(G[u*J+i],L,1,beta[i],v);
	        MatrixAddition(L,1,sum,v,sum);
		}
	}

	ScalarMatrix(lambda*(1-alpha),L,1,sum,v3);

	free(sum);
	free(v);
}

void SubX(int n,int L,int J,double* X,double** Xsub,double** XsubT)
{
	int i,j,k;

	for(i=0;i<J;i++) 
	{
		for(j=0;j<L;j++) 
		{
			for(k=0;k<n;k++)
			{
				Xsub[i][j*n+k]=X[n*L*i+j*n+k];
				XsubT[i][k*L+j]=Xsub[i][j*n+k];
			}
		}
	}
}

void BetaEstimationGLasso(int n,int L,int J,double* y,int lambdaNum,double alpha,double* lambdaSolution,int limit,double* G,double* dWeight,double** Xsub,double** XsubT,double** corr_xy,double** corr_xx,double* zeroBeta,double* lambda)
{
	int i,j,l,u,count,index;
	double d;
	double** beta=new double*[J];
	double** beta1=new double*[J];
	double* diff=new double[L*J];
	double* v1=new double[L];
	double* v3=new double[L];

	for(i=0;i<J;i++) 
	{
		beta[i]=new double[L];
		beta1[i]=new double[L];
	}

	for(i=0;i<J;i++) for(j=0;j<L;j++) beta[i][j]=0;

	for(l=0;l<lambdaNum;l++)
	{
		index=lambdaNum-l-1;	

		count=0;
	    do
	    {
		    count++;
	        for(i=0;i<J;i++) for(j=0;j<L;j++) beta1[i][j]=beta[i][j];
            
	        for(u=0;u<J;u++)
	        {
				P1(u,n,L,J,beta,corr_xy,corr_xx,v1);
				P2(u,L,J,beta,G,lambda[index],alpha,v3);

				MatrixAddition(L,1,v1,v3,v1);

				d=1/(1+lambda[index]*(1-alpha)*dWeight[u]);

				ScalarMatrix(d*PlusFt(1-lambda[index]*alpha/L2norm(L,v1)),L,1,v1,beta[u]);

		        if((L2norm(L,beta[u])>0) & (zeroBeta[u]==0)) 
		        {
					zeroBeta[u]=1;
				    CorrX(u,n,L,J,Xsub,XsubT,corr_xx,zeroBeta);
		        }
	        }
	        for(i=0;i<J;i++)  for(j=0;j<L;j++) diff[i*L+j]=fabs(beta[i][j]-beta1[i][j]);
	    }while(maxOFn(J*L,diff)>0.0001 && count<limit);
	    for(i=0;i<J;i++) for(j=0;j<L;j++) lambdaSolution[lambdaNum*L*i+j*lambdaNum+l]=beta[i][j];
	}

	for(i=0;i<J;i++)
	{
		free(beta[i]);
		free(beta1[i]);
	}
	free(beta);
	free(beta1);
	free(diff);
	free(v1);
	free(v3);
}

void GraphConstrainedEstimation(double* betaSolution,int* vec,double* y,double* X,double* G,double* alphaFixed,double* lambdaFixed)
{
	int i,j,k,t,n,J,L,lambdaNum,alphaNum,limit;
	
	n=vec[0];
    J=vec[1];
    L=vec[2];
    lambdaNum=vec[3];
    alphaNum=vec[4];
    limit=vec[5];
           	
	double** Xsub=new double*[J];
	double** XsubT=new double*[J];
	double** corr_xy=new double*[J];
	double** corr_xx=new double*[J*J];
	double* zeroBeta=new double[J];
	double* lambdaSolution=new double[lambdaNum*J*L];
	double* dWeight=new double[J];
	double* lambda1=new double[lambdaNum];

	for(i=0;i<J;i++) dWeight[i]=0;
	for(i=0;i<J;i++) for(j=0;j<J;j++) dWeight[i]=dWeight[i]+fabs(G[i*J+j]);
	
	for(i=0;i<J;i++) zeroBeta[i]=0;

	for(i=0;i<J;i++) corr_xy[i]=new double[L];

	for(i=0;i<J*J;i++) corr_xx[i]=new double[L*L];

	for(i=0;i<J;i++)
	{
		Xsub[i]=new double[n*L];
	    XsubT[i]=new double[L*n];
	}

	SubX(n,L,J,X,Xsub,XsubT);

	for(i=0;i<J;i++) MatrixMultiplication(L,n,1,XsubT[i],y,corr_xy[i]);

	CorrX(0,n,L,J,Xsub,XsubT,corr_xx,zeroBeta);

	for(t=0;t<alphaNum;t++)
	{
		for(i=0;i<lambdaNum;i++) lambda1[i]=lambdaFixed[i*alphaNum+t];

		BetaEstimationGLasso(n,L,J,y,lambdaNum,alphaFixed[t],lambdaSolution,limit,G,dWeight,Xsub,XsubT,corr_xy,corr_xx,zeroBeta,lambda1);
    
	    for(i=0;i<lambdaNum;i++) for(j=0;j<J;j++) for(k=0;k<L;k++) betaSolution[j*alphaNum*lambdaNum*L+k*alphaNum*lambdaNum+t*lambdaNum+i]=lambdaSolution[j*L*lambdaNum+k*lambdaNum+i];
	}

	free(Xsub);
	free(XsubT);
	free(corr_xy);
	free(corr_xx);
    free(zeroBeta);
    free(lambdaSolution);
	free(dWeight);
	free(lambda1);
}
}
