function [x,mu,lam,val,k]=sqpm(x0,mu0,lam0)
% ����: �û����������պ���Hesse���SQP�������Լ���Ż�����:  
%           min f(x)     s.t. h_i(x)=0, i=1,..., l.
%����:  x0�ǳ�ʼ��, mu0�ǳ��������ĳ�ʼֵ
%���:  x, mu �ֱ��ǽ������ŵ㼰��Ӧ�ĳ���, 
%  val������ֵ, mh��Լ��������ģ, k�ǵ�������.
maxk=100;   %����������
n=length(x0); l=length(mu0); m=length(lam0);
rho=0.5; eta=0.1;  B0=eye(n);
x=x0; mu=mu0;  lam=lam0; 
Bk=B0; sigma=0.8; 
epsilon1=1e-6; epsilon2=1e-5;
[hk,gk]=cons(x);  dfk=df1(x);
[Ae,Ai]=dcons(x); Ak=[Ae; Ai]; 
k=0;  
while(k<maxk)
    [dk,mu,lam]=qpsubp(dfk,Bk,Ae,hk,Ai,gk);  %���������
    mp1=norm(hk,1)+norm(max(-gk,0),1);
    if(norm(dk,1)<epsilon1)&(mp1<epsilon2)
        break; 
    end  %������ֹ׼��
    deta=0.05;  %����������
    tau=max(norm(mu,inf),norm(lam,inf));
    if(sigma*(tau+deta)<1)
        sigma=sigma;
    else
        sigma=1.0/(tau+2*deta);
    end
    im=0;   
    while(im<=20)   % Armijo����
        if(phi1(x+rho^im*dk,sigma)-phi1(x,sigma)<eta*rho^im*dphi1(x,sigma,dk))
            mk=im;
            break;
        end
        im=im+1;
        if(im==20),   mk=10;   end
    end
    alpha=rho^mk; x1=x+alpha*dk; 
    [hk,gk]=cons(x1);  dfk=df1(x1);
    [Ae,Ai]=dcons(x1);  Ak=[Ae; Ai];
    lamu=pinv(Ak)'*dfk;  %������С���˳���
    if(l>0&m>0)
        mu=lamu(1:l); lam=lamu(l+1:l+m);
    end
    if(l==0), mu=[]; lam=lamu;  end
    if(m==0), mu=lamu; lam=[];  end
    sk=alpha*dk;  %���¾���Bk
    yk=dlax(x1,mu,lam)-dlax(x,mu,lam);
    if(sk'*yk>0.2*sk'*Bk*sk)
        theta=1; 
    else
        theta=0.8*sk'*Bk*sk/(sk'*Bk*sk-sk'*yk);
    end
    zk=theta*yk+(1-theta)*Bk*sk;
    Bk=Bk+zk*zk'/(sk'*zk)-(Bk*sk)*(Bk*sk)'/(sk'*Bk*sk);
    x=x1;
    k=k+1;
end
val=f1(x);
%p=phi1(x,sigma)
%dd=norm(dk)
%%%%%%%% l1��ȷ��ֵ���� %%%%%%%
function p=phi1(x,sigma)
f=f1(x); [h,g]=cons(x); gn=max(-g,0);
l0=length(h);  m0=length(g);
if(l0==0), p=f+1.0/sigma*norm(gn,1); end
if(m0==0),  p=f+1.0/sigma*norm(h,1); end
if(l0>0&m0>0)
    p=f+1.0/sigma*(norm(h,1)+norm(gn,1)); 
end

%%%%% ��ֵ�����ķ�����%%%%%
function dp=dphi1(x,sigma,d)
df=df1(x); [h,g]=cons(x);  gn=max(-g,0);
l0=length(h);  m0=length(g);
if(l0==0),  dp=df'*d-1.0/sigma*norm(gn,1); end
if(m0==0), dp=df'*d-1.0/sigma*norm(h,1); end
if(l0>0&m0>0)
        dp=df'*d-1.0/sigma*(norm(h,1)+norm(gn,1));
end
%%%%%%%%% �������պ��� L(x,mu) %%%%%%%%%%%%%
function l=la(x,mu,lam)
f=f1(x);   %����Ŀ�꺯���ļ�
[h,g]=cons(x);  %����Լ�������ļ�
l0=lemgth(h);  m0=length(g);
if(l0==0), l=f-lam*g;  end
if(m0==0),  l=f-mu'*h;  end
if(l0>0&m0>0)
        l=f-mu'*h-lam'*g;  
end
%%%%%%%%% �������պ������ݶ� %%%%%%%%%%%%%
function dl=dlax(x,mu,lam)
df=df1(x); %����Ŀ�꺯���ݶ��ļ�
[Ae,Ai]=dcons(x);  %����Լ������Jacobi�����ļ�
[m1,m2]=size(Ai); [l1,l2]=size(Ae);
if(l1==0),  dl=df-Ai'*lam;  end
if(m1==0), dl=df-Ae'*mu;  end
if(l1>0&m1>0),  dl=df-Ae'*mu-Ai'*lam;  end
%%% Ŀ�꺯�� f(x) %%%%%%%%%%%
function f=f1(x)
f=-pi*x(1)^2*x(2);
%f=x(1)^2+x(2)^2-16*x(1)-10*x(2);
%f=-5*x(1)-5*x(2)-4*x(3)-x(1)*x(3)-6*x(4)-5*x(5)/(1+x(5));
%f=f-8*x(6)/(1+x(6))-10*(1-2*exp(-x(7))+exp(-2*x(7)));
%%%% Ŀ�꺯�� f(x) ���ݶ�%%%%%
function df=df1(x)
df=[-2*pi*x(1)*x(2), -pi*x(1)^2]';
%df=[2*x(1)-16; 2*x(2)-10];
%df=[-5-x(3); -5; -4-x(1); -6; -5/(1+x(5))^2; -8/(1+x(6))^2; ...
%       -20*exp(-x(7))+20*exp(-2*x(7))];
%%%% Լ������ %%%%%%
function [h,g]=cons(x)
h=[pi*x(1)*x(2)+pi*x(1)^2-150];
g=[x(1);x(2)];
%h=[ ];
%g=[-x(1)^2+6*x(1)-4*x(2)+11; x(1)*x(2)-3*x(2)-exp(x(1)-3)+1; x(1); x(2)];
%h=[2*x(4)+x(5)+08*x(6)+x(7)-5; x(2)^2+x(3)^2+x(5)^2+x(6)^2-5];
%g=[10-x(1)-x(2)-x(3)-x(4)-x(5)-x(6)-x(7); 5-x(1)-x(2)-x(3)-x(4);...
%      5-x(1)-x(3)-x(5)-x(6)^2+x(7)^2; x(1); x(2); x(3); x(4); x(5); x(6); x(7)];
%% Լ������ Jacobi����%%%%
function [dh,dg]=dcons(x)
dh=[pi*x(2)+2*pi*x(1), pi*x(1)];
dg=[1 0; 0 1];
%dh=[ ];
%dg=[-2*x(1)+6, -4; x(2)-exp(x(1)-3), x(1)-3; 1, 0; 0, 1];
%dh=[0, 0, 0, 2, 1, 0.8, 1; 0, 2*x(2), 2*x(3), 0, 2*x(5), 2*x(6), 0];
%dg=[-1 -1 -1 -1 -1 -1 -1; -1 -1 -1 -1 0 0 0; -1 0 -1 0 -1 -2*x(6) 2*x(7); ...
 %       1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 1 0 0 0; ...
 %       0 0 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1];





