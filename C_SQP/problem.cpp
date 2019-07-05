#include "problem.h"

double objfun(double x1, double x2)
{
	double f = pow(x2-x1*x1,2)+pow(1-x1,2);
	return f;
}

void Gradient(double x1, double x2, Vector<double>& c)
{
	double grad[2] = {0.0};
	grad[0] = 4*pow(x1,3)-4*x1*x2+2*x1-2;
	grad[1] = -2*pow(x1,2)+2*x2;

	for(int i = 0; i < c.size(); i++)
	{
		c[i] = grad[i];
	}
}

// The Hessian Matrix of the object function
void Hessian(double x1, double x2, Matrix<double>& Q)
{
	double Hess[2][2] = {0};

	Hess[0][0] = 12*pow(x1,2)-4*x2+2;  
	Hess[0][1] = -4*x1;
	Hess[1][0] = -4*x1;
	Hess[1][1] = 2;

	for(int i = 0; i < Q.nrows(); i++)
	{
		for(int j = 0; j < Q.ncols(); j++)
		{
			Q[i][j] = Hess[i][j];
		}
	}
}

//The gradient of equality constraint function
void GradientCE(double x1, double x2,Matrix<double>& CE)
{

}

//The constant term of equality constraint function
void CE0(double x1, double x2, Vector<double>& ce0)
{
	
}


//for inequality constraint 不等式约束的系数
void GradientCI(double x1, double x2, Matrix<double>& CI)  
{
	double gradie[2][2]={0};

	gradie[0][0] = x2;
	gradie[0][1] = -0.5*x1;
	gradie[1][0] = x1;
	gradie[1][1] = 1;

	for(int i = 0; i < CI.nrows(); i++)
	{
		for(int j = 0; j < CI.ncols(); j++)
		{
			CI[i][j] = gradie[i][j];
		}
	}
}

//不等式约束的常数项
void CI0(double x1, double x2, Vector<double>& ci0)
{
	double ci0x[2] = {0};
	
	// The inequality constraint 1 不等式约束 1
	double ine1 = x1*x2-4;
	
	// The inequality constraint 2 不等式约束 2
	double ine2 = -0.25*pow(x1,2)+x2-2;

	ci0x[0] = ine1;
    ci0x[1] = ine2;

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
