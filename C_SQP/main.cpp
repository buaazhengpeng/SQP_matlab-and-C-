/*********************************************************************
File main.cpp for SQP

This file contains an example on how to construct a SQP algorithm using
the solve_quadprog() function.

The test problem is a nonlinear programming problem as follow:

Minimize   f(x)=(x2-x1^2)^2 + (1-x1)^2
Subject to
           x1*x2-4         >= 0
		   -0.25*x1^2+x2-2 >= 0
The solution is x^T = [1.54813 2.59431] and f(x) = 0.3407


Given a solution estamate xk, and a small step d.

f(xk+d)=f(xk)+[¨Œf(xk)]^T*d + 1/2*(d^T)[¨Œ^2f(xk)]*d+....
h(xk+d)=h(xk)+[¨Œh(xk)]^T*d + 1/2*(d^T)[¨Œ^2h(xk)]*d+....  = 0
g(xk+d)=g(xk)+[¨Œg(xk)]^T*d + 1/2*(d^T)[¨Œ^2g(xk)]*d+.... >= 0

Form the linearly-constrainted/quadratic minimization problem:

minimize: f(xk)+[¨Œf(xk)]^T*d + 1/2*(d^T)[¨Œ^2f(xk)]*d
subject to:
          h(xk)+[¨Œh(xk)]^T*d = 0;
		  g(xk)+[¨Œg(xk)]^T*d >=0;
   

Author: zhchshen
Advanced Dynamic Inc.
Date:   2009-11-25 20:13
*********************************************************************/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#include "QuadProg.h"
#include "problem.h"
using namespace QuadProg;

int main (int argc, char *const argv[]) 
{

	Matrix<double> Q, CE, CI;  

	Vector<double> c, ce0, ci0;

	Vector<double> d;  //The optimal solution of a Quadratic Programming(QP)

	//The initial start point:
	double x1 = 2;
	double x2 = 2;

	int n;   // The number of variables in objective function
	int m;   // The number of equality constraints
	int p;   // The number of inequality constraints

	//The value of quadratic function.
	double sum = 0.0;  

	n = 2;
	d.resize(n);

	// The Hessian Matrix of Objective function
	Q.resize(n, n);
	c.resize(n);

	//The number of equality constraints
	m = 0;
	//coefficient matrix of equality constraints 
	CE.resize(n, m);      
	ce0.resize(m);

	//The number of inequality constraint
	p = 2;
	//coefficient matrix of inequality constraints 
	CI.resize(n, p);
	ci0.resize(p);

	double epsilon = 1e-6;
	double err = 10;

	// The loop of the SQP
	while(err > epsilon)
	{
		Hessian(x1, x2, Q);
		Gradient(x1, x2, c);
		
		//The equality constraint
		GradientCE(x1, x2, CE); 
		CE0(x1, x2, ce0);

		//The inequality constraint
		GradientCI(x1, x2, CI); 
		CI0(x1, x2, ci0);
		
		// Call a QP-solver in the SQP loop.
		solve_quadprog(Q, c, CE, ce0, CI, ci0, d);

		std::cout << "x = [" << x1  <<" "<< x2 <<"]"<< std::endl;
		std::cout << "d = [" << d[0]<<" "<<d[1]<<"]"<< std::endl;
		
		/* FOR DOUBLE CHECKING COST since in the solve_quadprog routine the matrix Q is modified */
		Hessian(x1, x2, Q);

		double xa = x1 + d[0];
		double xb = x2 + d[1];

		/*Calculate the value of the Quadratic Function in the approximating QP problem.*/
		sum = Quadfun(Q, c, d);
		cout<<"Quadratic Fun("<<d[0]<<","<<d[1]<<") ="<< sum <<endl;

		double temp = 0.0;
		for(int i=0; i<d.size(); i++)
		{
			temp =+ d[i]*d[i];
		}

		err = pow(temp,0.5);

		cout << "x_new = ["<<xa<<" "<<xb<<"]"<<endl;
		
		double objf = objfun(xa, xb);
	    cout<< "f("<<xa<<","<<xb<<") = "<< objf <<endl;
		
		cout<< "err= "<< err <<endl;
		cout<< endl;

		//system("pause");

		x1 = xa;
		x2 = xb;

	}

	system("pause");  // This line is valid only in Visual Studio IDE. Delete it, if in Linux environment.

	return 0;
}
