//(弹簧连滑块)滑块在外激励和弹簧力作用下在非光滑水平面上运动的动力学仿真程序 

 #include <cstdio>
 #include <assert.h>
 #include <iostream>   
 #include <cstdlib> 
 #include <time.h>
 #include "dos.h"
 #include "windows.h"
 #include "math.h"
 #include "stddef.h"
 #include "SQP.h"

//以下变量名不要改动
 #define   MAX_OF_LANK	9   //整个过程中多个LCP求解所用到的最大的矩阵阶数                                                                                                                                                                                                                                                                                                                                                                           
 int LANK;                  //单个LCP调用过程中矩阵的阶数，LCP_solution函数中的n
 #define   COUNT	100     //LCP求解用到的循环数
 double  LCP_y[9],LCP_x[9];  //LCP_y即为未知向量y；LCP_x即为未知向量x，对应y=Ax+b
void LCP_solution(int n,double* LCP_y,double* LCP_x,double *LCP_a,double*LCP_b);   //先声明LCP函数
void gerkt( double t, double h,double* y,int n,double eps);


double
g=9.8,  //重力加速度
m=1.0,  //滑块质量
K=15.0, //弹簧刚度
ww=1.0, //外激励频率
A=10.0, //外激励振幅
u=0.3,  //动滑动摩擦系数
u0=0.4, //静滑动摩擦系数
ddx,    //滑块加速度
FN,     //支撑力
FS,     //摩擦力
ddxM,   //负向加速度
ddxP,   //正向加速度
FP,     //正向摩擦余量
FM,     //负向摩擦余量
Fx,    //主动力水平力的和
veps=0.000001;  //判断速度为零的系数


double aaa[4];//对应LCP公式y=Ax+b中的矩阵A(n,n)//行排列 n*n 
double bbb[2];//对应y=Ax+b中的向量b

double param[4];//m,Fx,u0,FN  for SQP method
double *para;


int main()      //主程序
{
 long int i,mm=0;
 double t,h,eps;
 double *y;
 FILE *pf1, *pf2, *fopen();
 pf1=std::fopen("sliderLCP1.dat", "w+"); //数值结果输出到指定的数据（sliderLCP1-3.dat)文件
 pf2=std::fopen("sliderLCP2.dat", "w+"); 

t=0.0;         //初始时间
h=0.005;        //计算步长
eps=0.000001;  //常微分方程的计算精度
double yy[2];
yy[0]=0.0;      //初始位置
yy[1]=0.0;      //初始速度
y=&yy[0];
FN=m*g;       //支撑力等于重力

 for(i=0;i<=6000;i++)   //循环语句开始
 {
  gerkt(t,h,y,2,eps);   //调用变步长RK数值计算程序，求解该系统的动力学方程
  if(t>=0.0)
      {
      if(mm==0)
      {
        fprintf(pf1,"%f %f %f\n",t,y[0],y[1]);   //数值结果输出到指定文件（画图用）
        fprintf(pf2,"%f %f %f\n",t,y[1],FS);
        
       // printf("%f %f %f %f\n",t,y[0],y[1],FS);
        printf("t= %e",t); //输出到屏幕上的信息
      }
      }
      if(mm>=0)mm=mm+1;
      if(mm==10)mm=0; //每计算10个步长输出一组计算结果
      t=t+h;
 }  // 循环语句结束
}  // 主程序结束


//以下是动力学方程数值求解子程序
void gerktf( double t,double* y,int n,double *d)    //用RK计算方法求解该系统的动力学方程子程序
     //int n;
     //double t,y[],d[];
     {
       Fx=A*sin(ww*t)-K*y[0];   //作用在滑块上的主动力在x轴投影的代数和
	  
	 if(fabs(y[1])>=veps)   //滑块速度不为零时的动力学方程
       {
	    if(y[1]>0.0)  FS=-FN*u;  //判断滑动摩擦力方向
	    if(y[1]<0.0)  FS=FN*u;  //判断滑动摩擦力方向
	    d[0]=y[1];            
	    d[1]=(Fx+FS)/m;
	    ddx=d[1];          //滑块的加速度
       }
     else       //滑块速度为零时的动力学方程（用LCP求该时刻滑块的加速度)
       {
/*
//LCP method
        //矩阵A赋值,对应y=Ax+b中的A	
        //aaa[0][0]=1.0/m;  aaa[0][1]=1.0;
        //aaa[1][0]=-1.0;   aaa[1][1]=0.0;
       aaa[0]=1.0/m; aaa[1]=1.0; aaa[2]=-1.0; aaa[3]=0.0;
       //向量b赋值，对应y=Ax+b中的b	
       bbb[0]=-1.0*(Fx+u0*FN)/m;
       bbb[1]=2.0*u0*FN;

       LCP_solution(2,LCP_y,LCP_x,aaa,bbb);//调用LCP求解线性互补方程y=Ax+b 中的x和y
	   //LCP_y和LCP_x对应于y=Ax+b中的互补向量y,x
       ddxM=LCP_y[0];   //负向加速度
       FP=LCP_y[1];     //正向摩擦余量
       FM=LCP_x[0];     //负向摩擦余量
       ddxP=LCP_x[1];   //正向加速度
       FS=u0*FN-FM;     //摩擦力
       ddx=ddxP-ddxM;   //滑块加速度
	   d[0]=0.0;
       d[1]=ddx;
*/

//SQP method

		param[0]=m;
		param[1]=Fx;
		param[2]=u0;
		param[3]=FN;
		para=&param[0];
		
	   solve_SQP(FM,ddxP,para);

	   FS=u0*FN-FM;     //摩擦力
	   
	   ddxM=1/m*(-Fx-u0*FN+FM)+ddxP;
	   ddx=ddxP-ddxM; 
       //ddx=1/m*(Fx+FS);   //滑块加速度
       
    //printf("---------------------------------\n");
      //printf(" %f, %f,  %f,   %f\n",y[1],Fx,FM,ddxP);
      //printf(" %f,   %f\n",FS,ddx);
       
	   d[0]=0.0;
       d[1]=ddx;

       }
      return;
    }	
//
//积分器	
void gerkt( double t, double h,double* y,int n,double eps)
{
    int m,i,j,k;
    double hh,p,dt,x,tt,q,a[4]; 
    double *g=(double*)std::malloc(n*sizeof(double));
    double *b=(double*)std::malloc(n*sizeof(double));
    double *c=(double*)std::malloc(n*sizeof(double));
    double *d=(double*)std::malloc(n*sizeof(double));
    hh=h; m=1; p=1.0+eps; x=t;
    for (i=0; i<=n-1; i++) c[i]=y[i];
    while (p>=eps)
    {
 		a[0]=hh/2.0; a[1]=a[0]; a[2]=hh; a[3]=hh;
        for (i=0; i<=n-1; i++)
            {
 			g[i]=y[i]; y[i]=c[i];
			}
        dt=h/m; t=x;
        for (j=0; j<=m-1; j++)
          	{
 			gerktf(t,y,n,d);
            for (i=0; i<=n-1; i++) b[i]=y[i];
            for (k=0; k<=2; k++)
                {
 				for (i=0; i<=n-1; i++)
                  {
 					y[i]=y[i]+a[k]*d[i];
                    b[i]=b[i]+a[k+1]*d[i]/3.0;
					}
                tt=t+a[k];
                gerktf(tt,y,n,d);
				}
            for (i=0; i<=n-1; i++)
              y[i]=b[i]+hh*d[i]/6.0;
            t=t+dt;
			}			
        p=0.0;
        for (i=0; i<=n-1; i++)
            {
 			q=fabs(y[i]-g[i]);
            if (q>p) p=q;
			}
        hh=hh/2.0; m=m+m;
	}
    free(g); free(b); free(c); free(d);
    return;  
}

// LCP求解子程序（请不要修改，上述程序的变量名不要与下面程序中的变量名重名！！！！）
void LCP_solution(int n,double *LCP_y,double *LCP_x,double*LCP_a,double* LCP_b)//对应y=Ax+b，LCP_y为其中的向量y；LCP_x为其中的向量x，
{
	int i,j,jj,k,kk,count,note[MAX_OF_LANK],LANK;
	double a[MAX_OF_LANK][MAX_OF_LANK],marix_w[MAX_OF_LANK][MAX_OF_LANK];
	double b[MAX_OF_LANK],b1[MAX_OF_LANK],q[MAX_OF_LANK],q1[MAX_OF_LANK],min;
	double rti,max,w[MAX_OF_LANK],z[MAX_OF_LANK],z0,marix_comon[MAX_OF_LANK][2*MAX_OF_LANK];
	double marix_comon1[MAX_OF_LANK][2*MAX_OF_LANK];
//对应y=Ax+b，a[][]即为矩阵A；q[]即为向量b
//当滑动时，a[i][j]对应-mmm10[i][j]；q[i]对应q10[i]
//当粘滞时，a[i][j]对应-mmm20[i][j]；q[i]对应q20[i]
	LANK=n;

	for(i=0;i<LANK;i++)
	{
		for(j=0;j<LANK;j++)
		{
			if(i==j) marix_w[i][j]=1.0;
			else marix_w[i][j]=0.0;
		}
	}
	for(i=0;i<LANK;i++)
	{
		for(j=0;j<LANK;j++)
		{
            
				{
					a[i][j]=-LCP_a[j+i*LANK];
					//a[i][j]=-LCP_a[i][j];
				}

		}
	}
	for(i=0;i<LANK;i++)
	{

			{
				q[i]=LCP_b[i];
			}

	}
	for(i=0;i<LANK;i++)//矩阵b对应lemke方法中的初始系数矩阵
	{
		b[i]=-1;
	}

	for(i=0;i<LANK;i++)
	{
		w[i]=q[i];
		z[i]=0.0;
		z0=0.0;
		note[i]=i;
	}
	min=0.0;
	for(i=0;i<LANK;i++)
	{
		if(q[i]<min)
		{
			k=i;
			min=q[i];
		}
	}
	if(min<0.0)
	{
		for(i=0;i<LANK;i++)
		{
			for(j=0;j<LANK;j++)
			{
				marix_comon[i][j]=marix_w[i][j];
			}
		}
		for(i=0;i<LANK;i++)
		{
			for(j=LANK;j<2*LANK;j++)
			{
				marix_comon[i][j]=a[i][j-LANK];
			}
		}
		for(i=0;i<LANK;i++)
		{
			if(i==k)
			{
				for(j=0;j<2*LANK;j++)
				{
					marix_comon1[i][j]=marix_comon[i][j]/b[i];
				}
				q1[i]=q[i]/b[i];
				b1[i]=b[i]/b[i];
				note[i]=2*LANK;
			}
			else
			{
				for(j=0;j<2*LANK;j++)
				{
					marix_comon1[i][j]=marix_comon[i][j]-marix_comon[k][j]*(b[i]/b[k]);
				}
				q1[i]=q[i]-q[k]*(b[i]/b[k]);
				b1[i]=b[i]-b[k]*(b[i]/b[k]);
			}
		}
		kk=k+LANK;
		if(LANK ==4)//***********************************************************
		{
			k=k;
		}
		for(jj=0;jj<COUNT;jj++)
		{
			for(i=0;i<LANK;i++)
			{
				b[i]=b1[i];
				q[i]=q1[i];
				for(j=0;j<2*LANK;j++)
				{
					marix_comon[i][j]=marix_comon1[i][j];
				}
			}
			count=0;
			for(i=0;i<LANK;i++)
			{
				if((b1[i]+1.0) != 1.0)count=count+1;
			}
			if(count>1) break;
			max=0.0;
			for(i=0;i<LANK;i++)
			{
				if(marix_comon[i][kk]>max) max=marix_comon[i][kk];
			}
			if(max==0.0)
			{
				printf("not found the solution");
				break;
			}
			else
			{
				min=pow(10.0,61.0);
				for(i=0;i<LANK;i++)
				{
					if(marix_comon[i][kk]>0.0)
					{
						rti=q[i]/marix_comon[i][kk];
						if((min-rti)>0.0001)
						{
							min=rti;
							k=i;
						}
					}
				}
			}

			for(i=0;i<LANK;i++)
			{
				if(i==k)
				{
					q1[i]=q[i]/marix_comon[k][kk];
					b1[i]=b[i]/marix_comon[k][kk];
				}
				else
				{
					q1[i]=q[i]-q[k]*(marix_comon[i][kk]/marix_comon[k][kk]);
					b1[i]=b[i]-b[k]*(marix_comon[i][kk]/marix_comon[k][kk]);
				}
				for(j=0;j<2*LANK;j++)
				{
					if(i==k)
					{
						marix_comon1[i][j]=marix_comon[i][j]/marix_comon[k][kk];
					}
					else
					{
						marix_comon1[i][j]=marix_comon[i][j]-marix_comon[k][j]*(marix_comon[i][kk]/marix_comon[k][kk]);
					}
				}
			}
			if(note[k]<LANK)
			{
				for(i=0;i<LANK;i++)
				{
					if(note[k]==i)
					{
						note[k]=kk;
						kk=i+LANK;
						break;
					}
				}
			}
			else
			{
				for(i=LANK;i<=2*LANK;i++)
				{
					if(note[k]==i)
					{
						note[k]=kk;
						kk=i-LANK;
						break;
					}
				}
			}
		}
		if(jj==COUNT)
		{
			printf("not found the solution,%d\t",jj);
		}
		else
		{
			for(i=0;i<LANK;i++)
			{
				z[i]=0.0;
				w[i]=0.0;
			}
			for(i=0;i<LANK;i++)
			{
				if(note[i]<LANK)
				{
					for(j=0;j<LANK;j++)
					{
						if(note[i]==j)w[j]=q[i];
					}
				}
				else
				{
					for(j=LANK;j<2*LANK;j++)
					{
						if(note[i]==j)
						{
							z[j-LANK]=q[i];
						}
						else if(note[i]==2*LANK)
						{
							printf("not found the solution");
							break;
						}
					}
				}
			}
		}
		for(i=0;i<LANK;i++)
		{
			LCP_y[i]=w[i];
			LCP_x[i]=z[i];
		}
	}
	else
	{
		for(i=0;i<LANK;i++)
		{
			LCP_y[i]=w[i];
			LCP_x[i]=z[i];
		}
	}	
}





