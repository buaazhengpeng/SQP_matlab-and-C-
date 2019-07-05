#include "problem.h"

//double m=1.0;
//double Fx=1.0;
//double u0=0.1;
//double FN=10.0;

double objfun(double x1, double x2,double *para)
{
	double m=para[0];
	double Fx=para[1];
	double u0=para[2];
	double FN=para[3];
	double f = 1/m*pow(x1,2)-1/m*(Fx+u0*FN)*x1+2*u0*FN*x2;//+     pow(x1,2) +  pow(x2,2);
	return f;
}

void Gradient(double x1, double x2, Vector<double>& c,double *para)
{
	double m=para[0];
	double Fx=para[1];
	double u0=para[2];
	double FN=para[3];
	double grad[2] = {0.0};
	grad[0] = 2/m*x1-1/m*(Fx+u0*FN);
	grad[1] = 2*u0*FN;//+           2*x1+2*x2;

	for(int i = 0; i < c.size(); i++)
	{
		c[i] = grad[i];
	}
}

// The Hessian Matrix of the object function
void Hessian(double x1, double x2, Matrix<double>& Q,double *para)
{
	double m=para[0];
	double Fx=para[1];
	double u0=para[2];
	double FN=para[3];
	double Hess[2][2] = {0};

	Hess[0][0] = 2/m+1/pow(10,3);//+1; 
	Hess[0][1] = 0;
	Hess[1][0] = 0;
	Hess[1][1] = 1/pow(10,3);//2;

	for(int i = 0; i < Q.nrows(); i++)
	{
		for(int j = 0; j < Q.ncols(); j++)
		{
			Q[i][j] = Hess[i][j];
		}
	}
}

//The gradient of equality constraint function
void GradientCE(double x1, double x2,Matrix<double>& CE,double *para)
{

}

//The constant term of equality constraint function
void CE0(double x1, double x2, Vector<double>& ce0,double *para)
{
	
}


//for inequality constraint 不等式约束的系数
void GradientCI(double x1, double x2, Matrix<double>& CI,double *para)  
{
	double m=para[0];
	double Fx=para[1];
	double u0=para[2];
	double FN=para[3];
	double gradie[2][4]={0};

	gradie[0][0] = 1;
	gradie[0][1] = 0;
	gradie[0][2] = 1/m;
	gradie[0][3] = -1;
	
	gradie[1][0] = 0;
	gradie[1][1] = 1;
	gradie[1][2] = 1;
	gradie[1][3] = 0;
	

	

	for(int i = 0; i < CI.nrows(); i++)
	{
		for(int j = 0; j < CI.ncols(); j++)
		{
			CI[i][j] = gradie[i][j];
		}
	}
}

//不等式约束的常数项
void CI0(double x1, double x2, Vector<double>& ci0,double *para)
{
	double m=para[0];
	double Fx=para[1];
	double u0=para[2];
	double FN=para[3];
	
	double ci0x[4] = {0};
	
	// The inequality constraint 1 不等式约束 1
	double ine1 = x1;
	
	// The inequality constraint 2 不等式约束 2
	double ine2 = x2;
	
	// The inequality constraint 3 不等式约束 3
	double ine3 = 1/m*x1+x2-1/m*(Fx+u0*FN);
	
	// The inequality constraint 4 不等式约束 4
	double ine4 = 2*u0*FN-x1;
	
	

	ci0x[0] = ine1;
    ci0x[1] = ine2;
    ci0x[2] = ine3;
    ci0x[3] = ine4;

	for(int i = 0; i< ci0.size(); i++)
	{
		ci0[i] = ci0x[i];
	}
}

double Quadfun(Matrix<double> Q, Vector<double> c, Vector<double> d)
{
	//f(x)=1/2*(d^T)*Q*d+(c^T)*d
	double sum = 0.0;
	for (int i = 0; i < Q.nrows(); i++)
	{
		for (int j = 0; j < Q.ncols(); j++)
		{
			sum += d[i] * Q[i][j] * d[j];
		}
	}

	sum *= 0.5;	

	for (int i = 0; i < c.size(); i++)
	{
		sum += c[i] * d[i];
	}

	return sum;
}
