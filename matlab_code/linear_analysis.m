Dum=0.3;
Duc=0.3;
Dvm=0.012;
sigmaratio=1;
sigma=log(2)/50*sigmaratio*1;
sigma1=log(2)/50*sigmaratio*1;
sigma2=1e-3*sigma1;
ep=3*log(2)/50*sigmaratio;
gamma=3.6;
alpha=0.5;
beta=1.5*1e-4;
beta1=0.01*beta;
L=6;
Dtrialsize=50;
eigen=zeros(1,10);
Reciprocalstatistic=zeros(1,Dtrialsize);
Diffusionratiod=zeros(1,Dtrialsize);
%for j=1:1:Dtrialsize
%   Dum=0.06+0.01*j;
for i=1:30
    beta1=(0.65-0.02*i)*beta;
z=solve_fixed_point(sigma,sigma1,ep,beta,gamma,alpha,beta1,sigma2);
Um0=z(1);
Vm0=z(2);
Uc0=z(3);
for n=1:1:20
k=(n-1)*pi/L;
A=[-k^2*Dum-alpha-beta*Vm0^2-sigma, -2*beta*Um0*Vm0+gamma+3*beta1*Vm0^2, ep;
    alpha+beta*Vm0^2, -k^2*Dvm+2*beta*Um0*Vm0-gamma-sigma1-3*beta1*Vm0^2, sigma2;
    sigma, sigma1, -ep-k^2*Duc-sigma2];
e=eig(A);
e=sort(e);
eigen(1,n)=e(end);
end
E=max(eigen);
if(E>0)
    break;
end
end
beta1/beta
%Reciprocalstatistic(1,j)=beta/beta1 * gamma/alpha;
%Diffusionratiod(1,j)=Dum/Dvm;
%end
%kc=sqrt((Dum*(2*beta*Um0*Vm0-3*beta1*Vm0^2-gamma)+(-beta*Vm0^2-alpha)*Dvm)/(2*Dum*Dvm));
%lambda=2*pi/kc;
% figure(4)
% N=0:1:19;
% plot(N,eigen,'.')
% eigen2=zeros(1,10);
% for n=1:1:20
% k=(n-1)*pi/L;
% B=[-k^2*Dum-alpha-beta*Vm0^2-sigma, -2*beta*Um0*Vm0+gamma+3*beta1*Vm0^2, ep,0,-beta*Um0+3*beta1*Vm0,0,-2*beta*Vm0,0,0;
%      alpha+beta*Vm0^2, -k^2*Dvm+2*beta*Um0*Vm0-gamma-sigma1-3*beta1*Vm0^2, sigma2,0,beta*Um0-3*beta1*Vm0,0,2*beta*Vm0,0,0;
%      sigma, sigma1, -ep-k^2*Duc-sigma2,0,0,0,0,0,0;
%      0,0,0,2*(-k^2*Dum-alpha-beta*Vm0^2-sigma),0,0,2*( -2*beta*Um0*Vm0+gamma+3*beta1*Vm0^2),0,2*ep;
%      0,0,0,0,2*(-k^2*Dvm+2*beta*Um0*Vm0-gamma-sigma1-3*beta1*Vm0^2),0,2*(alpha+beta*Vm0^2),0,2*sigma2;
%      0,0,0,0,0,2*(-ep-k^2*Duc-sigma2),0,2*sigma,2*sigma1;
%      0,0,0,alpha+beta*Vm0^2,-2*beta*Um0*Vm0+gamma+3*beta1*Vm0^2,0, -k^2*Dvm+2*beta*Um0*Vm0-gamma-sigma1-3*beta1*Vm0^2-k^2*Dum-alpha-beta*Vm0^2-sigma,sigma2,ep;
%      0,0,0,sigma,0,ep,sigma1,-ep-k^2*Duc-sigma2-k^2*Dum-alpha-beta*Vm0^2-sigma,-2*beta*Um0*Vm0+gamma+3*beta1*Vm0^2;
%      0,0,0,0,sigma1,sigma2,sigma,alpha+beta*Vm0^2,-k^2*Dvm+2*beta*Um0*Vm0-gamma-sigma1-3*beta1*Vm0^2-ep-k^2*Duc-sigma2;];
%  e=eig(B);
%  e=real(e);
% e=sort(e);
% eigen2(1,n)=e(end);
% end;
% figure(5)
% N=0:1:19;
% plot(N,eigen2)
% 
% A11=-alpha-beta*Vm0^2-sigma;
% A12=-2*beta*Um0*Vm0+gamma+3*beta1*Vm0^2;
% A21=alpha+beta*Vm0^2;
% A22=2*beta*Um0*Vm0-gamma-sigma1-3*beta1*Vm0^2;
% Delta0=A11*A22-A12*A21;
% Dratio=(sqrt(Delta0) -sqrt(Delta0-A11*A22))/(A11)