//(����������)�������⼤���͵������������ڷǹ⻬ˮƽ�����˶��Ķ���ѧ������� 

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

//���±�������Ҫ�Ķ�
 #define   MAX_OF_LANK	9   //���������ж��LCP������õ������ľ������                                                                                                                                                                                                                                                                                                                                                                           
 int LANK;                  //����LCP���ù����о���Ľ�����LCP_solution�����е�n
 #define   COUNT	100     //LCP����õ���ѭ����
 double  LCP_y[9],LCP_x[9];  //LCP_y��Ϊδ֪����y��LCP_x��Ϊδ֪����x����Ӧy=Ax+b
void LCP_solution(int n,double* LCP_y,double* LCP_x,double *LCP_a,double*LCP_b);   //������LCP����
void gerkt( double t, double h,double* y,int n,double eps);


double
g=9.8,  //�������ٶ�
m=1.0,  //��������
K=15.0, //���ɸն�
ww=1.0, //�⼤��Ƶ��
A=10.0, //�⼤�����
u=0.3,  //������Ħ��ϵ��
u0=0.4, //������Ħ��ϵ��
ddx,    //������ٶ�
FN,     //֧����
FS,     //Ħ����
ddxM,   //������ٶ�
ddxP,   //������ٶ�
FP,     //����Ħ������
FM,     //����Ħ������
Fx,    //������ˮƽ���ĺ�
veps=0.000001;  //�ж��ٶ�Ϊ���ϵ��


double aaa[4];//��ӦLCP��ʽy=Ax+b�еľ���A(n,n)//������ n*n 
double bbb[2];//��Ӧy=Ax+b�е�����b

double param[4];//m,Fx,u0,FN  for SQP method
double *para;


int main()      //������
{
 long int i,mm=0;
 double t,h,eps;
 double *y;
 FILE *pf1, *pf2, *fopen();
 pf1=std::fopen("sliderLCP1.dat", "w+"); //��ֵ��������ָ�������ݣ�sliderLCP1-3.dat)�ļ�
 pf2=std::fopen("sliderLCP2.dat", "w+"); 

t=0.0;         //��ʼʱ��
h=0.005;        //���㲽��
eps=0.000001;  //��΢�ַ��̵ļ��㾫��
double yy[2];
yy[0]=0.0;      //��ʼλ��
yy[1]=0.0;      //��ʼ�ٶ�
y=&yy[0];
FN=m*g;       //֧������������

 for(i=0;i<=6000;i++)   //ѭ����俪ʼ
 {
  gerkt(t,h,y,2,eps);   //���ñ䲽��RK��ֵ�����������ϵͳ�Ķ���ѧ����
  if(t>=0.0)
      {
      if(mm==0)
      {
        fprintf(pf1,"%f %f %f\n",t,y[0],y[1]);   //��ֵ��������ָ���ļ�����ͼ�ã�
        fprintf(pf2,"%f %f %f\n",t,y[1],FS);
        
       // printf("%f %f %f %f\n",t,y[0],y[1],FS);
        printf("t= %e",t); //�������Ļ�ϵ���Ϣ
      }
      }
      if(mm>=0)mm=mm+1;
      if(mm==10)mm=0; //ÿ����10���������һ�������
      t=t+h;
 }  // ѭ��������
}  // ���������


//�����Ƕ���ѧ������ֵ����ӳ���
void gerktf( double t,double* y,int n,double *d)    //��RK���㷽������ϵͳ�Ķ���ѧ�����ӳ���
     //int n;
     //double t,y[],d[];
     {
       Fx=A*sin(ww*t)-K*y[0];   //�����ڻ����ϵ���������x��ͶӰ�Ĵ�����
	  
	 if(fabs(y[1])>=veps)   //�����ٶȲ�Ϊ��ʱ�Ķ���ѧ����
       {
	    if(y[1]>0.0)  FS=-FN*u;  //�жϻ���Ħ��������
	    if(y[1]<0.0)  FS=FN*u;  //�жϻ���Ħ��������
	    d[0]=y[1];            
	    d[1]=(Fx+FS)/m;
	    ddx=d[1];          //����ļ��ٶ�
       }
     else       //�����ٶ�Ϊ��ʱ�Ķ���ѧ���̣���LCP���ʱ�̻���ļ��ٶ�)
       {
/*
//LCP method
        //����A��ֵ,��Ӧy=Ax+b�е�A	
        //aaa[0][0]=1.0/m;  aaa[0][1]=1.0;
        //aaa[1][0]=-1.0;   aaa[1][1]=0.0;
       aaa[0]=1.0/m; aaa[1]=1.0; aaa[2]=-1.0; aaa[3]=0.0;
       //����b��ֵ����Ӧy=Ax+b�е�b	
       bbb[0]=-1.0*(Fx+u0*FN)/m;
       bbb[1]=2.0*u0*FN;

       LCP_solution(2,LCP_y,LCP_x,aaa,bbb);//����LCP������Ի�������y=Ax+b �е�x��y
	   //LCP_y��LCP_x��Ӧ��y=Ax+b�еĻ�������y,x
       ddxM=LCP_y[0];   //������ٶ�
       FP=LCP_y[1];     //����Ħ������
       FM=LCP_x[0];     //����Ħ������
       ddxP=LCP_x[1];   //������ٶ�
       FS=u0*FN-FM;     //Ħ����
       ddx=ddxP-ddxM;   //������ٶ�
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

	   FS=u0*FN-FM;     //Ħ����
	   
	   ddxM=1/m*(-Fx-u0*FN+FM)+ddxP;
	   ddx=ddxP-ddxM; 
       //ddx=1/m*(Fx+FS);   //������ٶ�
       
    //printf("---------------------------------\n");
      //printf(" %f, %f,  %f,   %f\n",y[1],Fx,FM,ddxP);
      //printf(" %f,   %f\n",FS,ddx);
       
	   d[0]=0.0;
       d[1]=ddx;

       }
      return;
    }	
//
//������	
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

// LCP����ӳ����벻Ҫ�޸ģ���������ı�������Ҫ����������еı�������������������
void LCP_solution(int n,double *LCP_y,double *LCP_x,double*LCP_a,double* LCP_b)//��Ӧy=Ax+b��LCP_yΪ���е�����y��LCP_xΪ���е�����x��
{
	int i,j,jj,k,kk,count,note[MAX_OF_LANK],LANK;
	double a[MAX_OF_LANK][MAX_OF_LANK],marix_w[MAX_OF_LANK][MAX_OF_LANK];
	double b[MAX_OF_LANK],b1[MAX_OF_LANK],q[MAX_OF_LANK],q1[MAX_OF_LANK],min;
	double rti,max,w[MAX_OF_LANK],z[MAX_OF_LANK],z0,marix_comon[MAX_OF_LANK][2*MAX_OF_LANK];
	double marix_comon1[MAX_OF_LANK][2*MAX_OF_LANK];
//��Ӧy=Ax+b��a[][]��Ϊ����A��q[]��Ϊ����b
//������ʱ��a[i][j]��Ӧ-mmm10[i][j]��q[i]��Ӧq10[i]
//��ճ��ʱ��a[i][j]��Ӧ-mmm20[i][j]��q[i]��Ӧq20[i]
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
	for(i=0;i<LANK;i++)//����b��Ӧlemke�����еĳ�ʼϵ������
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





