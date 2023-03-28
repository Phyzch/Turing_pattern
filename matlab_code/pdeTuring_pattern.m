% 采用pdepe编写解一维偏微分方程的程序
% [Reciprocalstatistic,Diffusionratiod,Wcost,Vmmaxstatistic]=pdeTuring_pattern1;
% function [Reciprocalstatistic,Diffusionratiod,Wcost,Vmmaxstatistic]= pdeTuring_pattern1
function pdeTuring_pattern
cd 'C:\学习\课件\本研\ouyangzu\Turing_pattern\self-positioned Turing pattern';
global Dum Duc Dvm  sigma sigma1 ep gamma alpha beta beta1 dx sigma2 L mark um vm uc   
mark =false;
if mark==true
    load('data5.mat');
end
Dum=0.3*1;
Duc=0.3*1;
Dvm=0.012*1;
sigmaratio=1;
sigma=log(2)/50*sigmaratio*1;
sigma1=log(2)/50*sigmaratio*1;
sigma2=1e-3*sigma1*1;
ep=3*log(2)/50*sigmaratio*1;
gamma=3.6;
alpha=0.5;
beta=1.5*1e-4;
beta1=0.35*beta;
dx=0.05;
L=6;
Dtrialsize=25;
Gammatrialsize=12;
Reciprocalstatistic=zeros(Gammatrialsize,1);
Diffusionratiod=zeros(Dtrialsize,1);
Wcost=zeros(Gammatrialsize,Dtrialsize);
Vmmaxstatistic=zeros(Gammatrialsize,Dtrialsize);

for j=1:1:Dtrialsize
    Dum=0.03+0.02*j;
for i=1:1:Gammatrialsize
    beta1=(0.65-0.05*i)*beta;
     Diffusionratiod(j,1)=Dum/Dvm;
     Reciprocalstatistic(i,1)=beta/beta1 * gamma/alpha;
% z=solve_fixed_point(sigma,sigma1,ep,beta,gamma,alpha,beta1,sigma2,gamma2,alpha2);
z=solve_fixed_point(sigma,sigma1,ep,beta,gamma,alpha,beta1,sigma2);
m=0;
x=0:dx:L;
dt=1;
t=0:dt:2000;
sol=pdepe(m,@pdex1pde,@(x)pde1lic(z,x),@pde1bc,x,t);
um=sol(:,:,1);
vm=sol(:,:,2);
uc=sol(:,:,3);
%Vmmax=max(vm(end,:));
meanvm=mean(vm(end,:));
Vmvariancestatistic(i,j)= sum((vm(end,:)-meanvm).^2);
% Freeenergy=sum(um(end,:).*log(um(end,:)./z(1))+vm(end,:).*log(vm(end,:)./z(2))+uc(end,:).*log(uc(end,:)./z(3)))
% dxum=Deriv1(um);
% dxxum=Deriv1(dxum);
% dxvm=Deriv1(vm);
% dxxvm=Deriv1(dxvm);
% dxuc=Deriv1(uc);
% dxxuc=Deriv1(dxuc);
Num=sum(um(:,:),2);
Nvm=sum(vm(:,:),2);
Nuc=sum(uc(:,:),2);
totalJ2=(gamma*Nvm(end,1)-alpha*Num(end,1));
totalJ4=(ep*Nuc(end,1)-sigma*Num(end,1));
W1=totalJ2*log(gamma*beta/(alpha*beta1))+totalJ4*log(beta*sigma1*ep/(beta1*sigma2*sigma));
Wcost(i,j)=W1;
end
end
% J1=vm(end,:).^2.*(beta.*um(end,:)-beta1.*vm(end,:));
% J2=gamma* vm(end,:) -alpha*um(end,:);
% J4=-sigma*um(end,:)+ep*uc(end,:);
% J3=-sigma2*uc(end,:)+sigma1*vm(end,:);
% Jum=-Dum.*dxum(end,:);
% Jvm=-Dvm.*dxvm(end,:);
% Juc=-Duc.*dxuc(end,:);
% dxJuc=Deriv1(Juc);
% dxJum=Deriv1(Jum);
% dxJvm=Deriv1(Jvm);
% 
% J1sum=sum(J1);

% figure(5)
% imagesc(t,x,vm');
% figure(5)
% hold on;
% xlabel('x (μ m)');
% ylabel('Time (s)');
% plot(x,um(end,:));
% 
%  figure(6)
%  hold on;
% plot(x,vm(end,:));
% xlabel('x (μ m)');
% ylabel('Time (s)');
%  figure(7)
%  hold on;
% plot(x,uc(end,:));
% xlabel('x (μ m)');
% ylabel('Time (s)');

% Potential1=log(beta.*um(end,:)./(beta1.*vm(end,:)));
% Potential2=log(gamma.*vm(end,:)./(alpha.*um(end,:)));
% Potential4=-log(sigma.*um(end,:)./(ep.*uc(end,:)));
% Potential3=-log(sigma2.*uc(end,:)./(sigma1.*vm(end,:)));
% Wdiffusion=Dum.*dxum(end,:).^2./(um(end,:))*dx+Dvm.*dxvm(end,:).^2./(vm(end,:))*dx+Duc.*dxuc(end,:).^2./(uc(end,:))*dx;
% Wreaction=J1.*Potential1+J2.*Potential2+J3.*Potential3+J4.*Potential4;
% Wreaction2=sum(J3.*Potential3+J4.*Potential4+Duc.*dxuc(end,:).^2./(uc(end,:))*dx) ;
% W2=sum(Wdiffusion+Wreaction);
% Wdiffusion1=sum(Wdiffusion);
% Wdiffusionum=sum(Dum.*dxum(end,:).^2./(um(end,:))*dx);
% relative_error=abs(W1-W2)/W1;

% calculate energy dissipation without diffusion
%
% um0,vm0,uc0代表最基本的浓度值
% um0=z(1);
% vm0=z(2);
% uc0=z(3);
%Wstable=L/dx*((gamma*vmc-alpha*umc)*log(gamma*beta/(alpha*beta1))+(ep*ucc-sigma*umc)*log(beta*sigma1*ep/(beta1*sigma2*sigma)));
save('orderparameter.mat','Reciprocalstatistic','Diffusionratiod','Wcost','Vmvariancestatistic');
% vm1=vm(:,1:1:end/3);
% vm2=vm(:,end/3:1:end*2/3);
% vm3=vm(:,end*2/3:1:end);
% Time=0;
% for i=1:1:size(vm,1)
%     if( (abs(find(vm1(i,:)==max(vm1(i,:)),1)-20)<=3) && (abs(find(vm2(i,:)==max(vm2(i,:)),1)-20)<=3) && (abs(find(vm3(i,:)==max(vm3(i,:)),1)-20)<=3) )
%         Time=i;
%         break;
%     end
% end
% Time=Time*dt; %定位时间
end

function [c,f,s]=pdex1pde(x,t,u,DuDx)
global Dum Dvm Duc alpha beta gamma sigma sigma1 ep beta1 sigma2 
c=[1;1;1];
f=[Dum;Dvm;Duc].*DuDx;
Um=u(1);
Vm=u(2);
Uc=u(3);

x0=-alpha.*Um-beta.*Um.*Vm.^2+gamma.*Vm+ep.*Uc-sigma.*Um+beta1*Vm^3;
x1=alpha.*Um+beta.*Um.*Vm.^2-gamma.*Vm-sigma1.*Vm-beta1*Vm^3 +sigma2*Uc;
x2=-ep.*Uc+sigma.*Um+sigma1.*Vm-sigma2*Uc;

% x0=-alpha2.*Um-beta.*Um.*Vm.^2+gamma2.*Wm+ep.*Uc-sigma.*Um+beta1*Vm^3;
% x1=alpha.*Wm+beta.*Um.*Vm.^2-gamma.*Vm-sigma1.*Vm-beta1*Vm^3 +sigma2*Uc;
% x2=-ep.*Uc+sigma.*Um+sigma1.*Vm-sigma2*Uc;
% x3= alpha2*Um+gamma*Vm-gamma2*Wm-alpha*Wm;

s=[x0;x1;x2];
% s=[x0;x1;x2;x3];
end

function u0= pde1lic(z,x)
global um vm uc  mark dx
std=0.01;
if mark==false
u1=std*rand(1,1)*z(1)+z(1);
u2=std*rand(1,1)*z(2)+z(2);
u3=std*rand(1,1)*z(3)+z(3);
% u1=std*z(1)+z(1);
% u2=std*z(2)+z(2);
% u3=std*z(3)+z(3);

u0=[u1;u2;u3];
% u0=[u1;u2;u3;u4];
else 
    u1=um(end,floor(x/dx)+1);
    u2=vm(end,floor(x/dx)+1);
    u3=uc(end,floor(x/dx)+1);
    u0=[u1;u2;u3];
    
end
end

function [pl,ql,pr,qr]=pde1bc(xl,ul,xr,ur,t)
pl=[0;0;0];
pr=[0;0;0];
ql=[1;1;1];
qr=[1;1;1];

% pl=[0;0;0;0];
% pr=[0;0;0;0];
% ql=[1;1;1;1];
% qr=[1;1;1;1];
end

function [dxF]=Deriv1(F)
Nx=size(F,2);
global dx
dxF = zeros(size(F));
dxF(:,1) = (-1.5*F(:,1) + 2.0*F(:,2) - 0.5*F(:,3)) / dx;
dxF(:,2:Nx-1) = (F(:,3:Nx) - F(:,1:Nx-2) ) / (2*dx);
dxF(:,Nx) = (1.5*F(:,Nx) - 2.0*F(:,Nx-1) + 0.5*F(:,Nx-2)) / dx;
end