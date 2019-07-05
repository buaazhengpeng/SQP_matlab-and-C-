function [x,mu,val,mh,k]=newtlagr(x0,mu0)
% 功能: 用牛顿-拉格朗日法求解约束优化问题:  
%           min f(x)     s.t. h_i(x)=0, i=1,..., l.
%输入:  x0是初始点, mu0是乘子向量的初始值
%输出:  x, mu 分别是近似最优点及相应的乘子, 
%  val是最优值, mh是约束函数的模, k是迭代次数.
maxk=200;   %最大迭代次数
n=length(x0); l=length(mu0);
rho=0.6; gamma=0.2;
x=x0; mu=mu0;
k=0;  epsilon=1e-12;
while(k<maxk)
    dl=dla(x,mu);  %计算乘子函数的梯度
    if(norm(dl)^2<epsilon), break; end  %检验终止准则
    N=N1(x,mu);    % 计算拉格朗日矩阵
    dz=-N\dl;    %解方程组得搜索方向
    dx=dz(1:n); du=dz(n+1:n+l);
    m=0; mk=0;
    while(m<20)   % Armijo搜索
        if(norm(dla(x+rho^m*dx,mu+rho^m*du))^2<=(1-gamma*rho^m)*norm(dl)^2)
            mk=m; break;
        end
        m=m+1;
    end
    x=x+rho^mk*dx; mu=mu+rho^mk*du;
    k=k+1;
end
val=f1(x);
mh=norm(h1(x));
%%%%%%%%% 拉格朗日函数 L(x,mu) %%%%%%%%%%%%%
function l=la(x,mu)
f=f1(x);   %调用目标函数文件
h=h1(x);  %调用约束函数文件
l=f-mu'*h; % 计算乘子函数
%%%%%%%%% 拉格朗日函数的梯度 %%%%%%%%%%%%%
function dl=dla(x,mu)
df=df1(x); %调用目标函数梯度文件
h=h1(x);  %调用约束函数文件
dh=dh1(x);  %调用约束函数Jacobi矩阵文件
dl=[df-dh'*mu; -h];  %计算乘子函数梯度文件
%%%%%%%%% 拉格朗日函数的Hesse阵 %%%%%%%%%%%
function d2l=d2la(x,mu)
d2f=d2f1(x);   %调用目标函数Hesse阵文件
[d2h1,d2h2,d2h3]=d2h(x); %调用约束函数二阶导数文件
d2l=d2f-mu(1)*d2h1-mu(2)*d2h2-mu(3)*d2h3; %计算乘子函数的Hesse阵
%%%%%%%%% 系数矩阵N(x,mu) %%%%%%%%%%%
function N=N1(x,mu)  %计算拉格朗日矩阵
l=length(mu);
d2l=d2la(x,mu);   dh=dh1(x);
N=[d2l, -dh'; -dh, zeros(l,l)];
%%%%%%%%% 目标函数 f(x) %%%%%%%%%%%%%
function f=f1(x)
s=x(1)*x(2)*x(3)*x(4)*x(5);
f=exp(s)-0.5*(x(1)^3+x(2)^3+1)^2;
%%%%%%%%% 约束函数 h(x) %%%%%%%%%%%%%
function h=h1(x)
h=[x(1)^2+x(2)^2+x(3)^2+x(4)^2+x(5)^2-10; x(2)*x(3)-5*x(4)*x(5); x(1)^3+x(2)^3+1];
%%%%%%%%% 目标函数 f(x) 的梯度%%%%%%%%%%%%%
function df=df1(x)
s=x(1)*x(2)*x(3)*x(4)*x(5);
df(1)=s/(x(1))*exp(s)-3*(x(1)^3+x(2)^3+1)*x(1)^2;
df(2)=s/(x(2))*exp(s)-3*(x(1)^3+x(2)^3+1)*x(2)^2;
df(3)=s/(x(3))*exp(s); 
df(4)=s/(x(4))*exp(s);
df(5)=s/(x(5))*exp(s);
df=df(:);
%%%%%%%%% 约束函数 h(x) 的Jacobi矩阵A(x)%%%%%%%%%%%%%
function dh=dh1(x)
dh=[2*x(1),2*x(2),2*x(3),2*x(4),2*x(5); 0,x(3),x(2),-5*x(5),-5*x(4); 3*x(1)^2,3*x(2)^2,0,0,0];
%%%%%%%%% 目标函数 f(x) 的Hesse阵%%%%%%%%%%%%%
function d2f=d2f1(x)
s=x(1)*x(2)*x(3)*x(4)*x(5);
d2f=[(s/(x(1)))^2*exp(s)-6*x(1)*(x(1)^3+x(2)^3+1)-9*x(1)^4, ...
        (1+s)*x(3)*x(4)*x(5)*exp(s)-9*x(1)^2*x(2)^2,  (1+s)*x(2)*x(4)*x(5)*exp(s), ...
        (1+s)*x(2)*x(3)*x(5)*exp(s),  (1+s)*x(2)*x(3)*x(4)*exp(s); ...
        (1+s)*x(3)*x(4)*x(5)*exp(s)-9*x(1)^2*x(2)^2,  ...
        (s/(x(2)))^2*exp(s)-6*x(2)*(x(1)^3+x(2)^3+1)-9*x(2)^4, ...
        (1+s)*x(1)*x(4)*x(5)*exp(s),  (1+s)*x(1)*x(3)*x(5)*exp(s),  (1+s)*x(1)*x(3)*x(4)*exp(s);  ...
        (1+s)*x(2)*x(4)*x(5)*exp(s),  (1+s)*x(1)*x(4)*x(5)*exp(s),  s^2/(x(3))*exp(s), ...
        (1+s)*x(1)*x(2)*x(5)*exp(s),  (1+s)*x(1)*x(2)*x(4)*exp(s);  ...
        (1+s)*x(2)*x(3)*x(5)*exp(s),  (1+s)*x(1)*x(3)*x(5)*exp(s), ...
        (1+s)*x(1)*x(2)*x(5)*exp(s),  s^2/(x(4))*exp(s),  (1+s)*x(1)*x(2)*x(3)*exp(s); ...
        (1+s)*x(2)*x(3)*x(4)*exp(s),  (1+s)*x(1)*x(3)*x(4)*exp(s),  ...
        (1+s)*x(1)*x(2)*x(4)*exp(s),   (1+s)*x(1)*x(2)*x(3)*exp(s),  s^2/(x(5))*exp(s)]'; 
%%%%%%%%% 约束函数 h(x) 的Hesse阵%%%%%%%%%%%%%
function [d2h1,d2h2,d2h3]=d2h(x)
d2h1=[2 0 0 0 0; 0 2 0 0 0; 0 0 2 0 0; 0 0 0 2 0; 0 0 0 0 2]';
d2h2=[0 0 0 0 0; 0 0 1 0 0; 0 1 0 0 0; 0 0 0 0 -5; 0 0 0 -5 0]';
d2h3=[6*x(1) 0 0 0 0; 0 6*x(2) 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]';


