function [x,P1,val,k]=lagsqp(x0,mu0)
% ����: �û����������պ���Hesse���SQP�������Լ���Ż�����:  
%           min f(x)     s.t. h_i(x)=0, i=1,..., l.
%����:  x0�ǳ�ʼ��, mu0�ǳ��������ĳ�ʼֵ
%���:  x, mu �ֱ��ǽ������ŵ㼰��Ӧ�ĳ���, 
%  val������ֵ, mh��Լ��������ģ, k�ǵ�������.
maxk=500;   %����������
n=length(x0); l=length(mu0);
rho=0.6; gamma=0.2;
x=x0; mu=mu0;  tau=1.55;
k=0;  epsilon=1e-12;
while(k<maxk)
    P1=P(x,mu);  %���㷣������ֵ
    if(P1<epsilon), break; end  %������ֹ׼��
    H=B(x,mu,tau);    % ����KT����
    c=df1(x);    %����Ŀ�꺯���ݶ�
    be=-h1(x); %����Լ������
    Ae=dh1(x); %����Լ��������Jacobi����
    [dx,lam]=qsubp(H,c,Ae,be);
    du=lam-mu-1.0/(2*tau)*dh1(x)*dx;
    m=0; mk=0;
    while(m<20)   % Armijo����
        if(P(x+rho^m*dx,mu+rho^m*du)<=(1-gamma*rho^m)*P1)
            mk=m;
            break;
        end
        m=m+1;
    end
    x=x+rho^mk*dx; mu=mu+rho^mk*du;
    k=k+1;
end
val=f1(x);
P1=P(x,mu);
%%%%%%%% ��������� %%%%%%%%%%%%%%%
function [x,mu1]=qsubp(H,c,Ae,be)
ginvH=pinv(H);
[m,n]=size(Ae);
if (m>0)
    rb = Ae*ginvH*c + be;
    mu1 = pinv(Ae*ginvH*Ae')*rb;
    x = ginvH*(Ae'*mu1-c);
else
    x = -ginvH*c;
    mu1 = zeros(m,1);
end
%%%%%%%%% �������պ��� L(x,mu) %%%%%%%%%%%%%
function l=la(x,mu)
f=f1(x);   %����Ŀ�꺯���ļ�
h=h1(x);  %����Լ�������ļ�
l=f-mu'*h; % ������Ӻ���
%%%%%%%%% �������պ������ݶ� %%%%%%%%%%%%%
function dl=dla(x,mu)
df=df1(x); %����Ŀ�꺯���ݶ��ļ�
h=h1(x);  %����Լ�������ļ�
dh=dh1(x);  %����Լ������Jacobi�����ļ�
dl=[df-dh'*mu; -h];  %������Ӻ����ݶ��ļ�
%%%%%%%%% ������P(x,mu) %%%%%%%%%%%%%%%
function s=P(x,mu)
dl=dla(x,mu);
s=norm(dl)^2;
%%%%%%%%% �������պ�����Hesse�� %%%%%%%%%%%
function d2l=d2la(x,mu)
d2f=d2f1(x);   %����Ŀ�꺯��Hesse���ļ�
[d2h1,d2h2,d2h3]=d2h(x); %����Լ���������׵����ļ�
d2l=d2f-mu(1)*d2h1-mu(2)*d2h2-mu(3)*d2h3;  
%%%%%%%%% KT����B(x,mu) %%%%%%%%%%%
function H=B(x,mu,tau)  %����KT����
d2l=d2la(x,mu);   %����Hesse��
dh=dh1(x);  %Լ��������Jacobi����
H=d2l+1.0/(2*tau)*dh'*dh;
%%%%%%%%% Ŀ�꺯�� f(x) %%%%%%%%%%%%%
function f=f1(x)
s=x(1)*x(2)*x(3)*x(4)*x(5);
f=exp(s)-0.5*(x(1)^3+x(2)^3+1)^2;
%%%%%%%%% Լ������ h(x) %%%%%%%%%%%%%
function h=h1(x)
h=[x(1)^2+x(2)^2+x(3)^2+x(4)^2+x(5)^2-10; x(2)*x(3)-5*x(4)*x(5); x(1)^3+x(2)^3+1];
%%%%%%%%% Ŀ�꺯�� f(x) ���ݶ�%%%%%%%%%%%%%
function df=df1(x)
s=x(1)*x(2)*x(3)*x(4)*x(5);
df(1)=s/(x(1))*exp(s)-3*(x(1)^3+x(2)^3+1)*x(1)^2;
df(2)=s/(x(2))*exp(s)-3*(x(1)^3+x(2)^3+1)*x(2)^2;
df(3)=s/(x(3))*exp(s); 
df(4)=s/(x(4))*exp(s);
df(5)=s/(x(5))*exp(s);
df=df(:);
%%%%%%%%% Լ������ h(x) ��Jacobi����A(x)%%%%%%%%%%%%%
function dh=dh1(x)
dh=[2*x(1),2*x(2),2*x(3),2*x(4),2*x(5); 0,x(3),x(2),-5*x(5),-5*x(4); 3*x(1)^2,3*x(2)^2,0,0,0];
%%%%%%%%% Ŀ�꺯�� f(x) ��Hesse��%%%%%%%%%%%%%
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
%%%%%%%%% Լ������ h(x) ��Hesse��%%%%%%%%%%%%%
function [d2h1,d2h2,d2h3]=d2h(x)
d2h1=[2 0 0 0 0; 0 2 0 0 0; 0 0 2 0 0; 0 0 0 2 0; 0 0 0 0 2]';
d2h2=[0 0 0 0 0; 0 0 1 0 0; 0 1 0 0 0; 0 0 0 0 -5; 0 0 0 -5 0]';
d2h3=[6*x(1) 0 0 0 0; 0 6*x(2) 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




