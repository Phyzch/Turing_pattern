% [Reciprocal_statistic,Diffusion_ratio,Wcost,Vmmaxstatistic]=pdeTuring_pattern1;
% function [Reciprocal_statistic,Diffusion_ratio,Wcost,Vmmaxstatistic]= pdeTuring_pattern1
function pdeTuring_pattern
cd 'path';
global DX1 DX3 DX2  k13 k23 k31 k21 k12 beta_12 beta_21 dx k32 L mark X1 X2 X3   

DX1 = 0.3*1;
DX3 = 0.3*1;
DX2 = 0.012*1;
k13 = log(2)/50;
k23 = log(2)/50;
k32 = 1e-3*k23*1;
k31 = 3*log(2)/50;
k21 = 3.6;
k12 = 0.5;
beta_12 = 1.5*1e-4;
beta_21 = 0.35*beta_12;
dx = 0.05;
L = 6;
D_trial_size = 25;
k21_trial_size = 12;
Reciprocal_statistic = zeros(k21_trial_size,1);
Diffusion_ratio = zeros(D_trial_size,1);
Wcost = zeros(k21_trial_size,D_trial_size);
Vmmaxstatistic = zeros(k21_trial_size,D_trial_size);

for j=1:1:D_trial_size
    DX1 = 0.03 + 0.02 * j;
	
	for i = 1:1:k21_trial_size
	
		beta_21 = (0.65-0.05*i)*beta_12;
		 Diffusion_ratio(j,1) = DX1/DX2;
		 Reciprocal_statistic(i,1) = beta_12/beta_21 * k21/k12;

		z=solve_fixed_point(k13,k23,k31,beta_12,k21,k12,beta_21,k32);
		m=0;
		x=0:dx:L;
		dt=1;
		t=0:dt:2000;
		sol=pdepe(m,@pdex1pde,@(x)pde1lic(z,x),@pde1bc,x,t);  % https://www.mathworks.com/help/matlab/ref/pdepe.html
		% pdex1pde: the equation being solved 
		X1=sol(:,:,1);
		X2=sol(:,:,2);
		X3=sol(:,:,3);

		mean_vm=mean(X2(end,:));
		Vm_variance_statistic(i,j)= sum((X2(end,:)-mean_vm).^2);

		NX1 = sum(X1(:,:),2);
		NX2 = sum(X2(:,:),2);
		NX3 = sum(X3(:,:),2);
		totalJ2 = (k21 * NX2(end,1) - k12 * NX1(end,1));
		totalJ4 = (k31 * NX3(end,1) - k13 * NX1(end,1));
		W1 = totalJ2 * log(k21 * beta_12 /(k12 * beta_21)) + totalJ4 * log(beta_12 * k23 * k31 / (beta_21 * k32 * k13));
		Wcost(i,j) = W1;
		
	end
	
end
save('orderparameter.mat','Reciprocal_statistic','Diffusion_ratio','Wcost','Vm_variance_statistic');

end

function [c,f,s]=pdex1pde(x,t,u,DuDx)  % dX/dt 
global DX1 DX2 DX3 k12 beta_12 k21 k13 k23 k31 beta_21 k32 
c=[1;1;1];
f=[DX1;DX2;DX3].*DuDx;
X1=u(1);
X2=u(2);
X3=u(3);

x0=-k12.*X1-beta_12.*X1.*X2.^2+k21.*X2+k31.*X3-k13.*X1+beta_21*X2^3;
x1=k12.*X1+beta_12.*X1.*X2.^2-k21.*X2-k23.*X2-beta_21*X2^3 +k32*X3;
x2=-k31.*X3+k13.*X1+k23.*X2-k32*X3;

s=[x0;x1;x2];
% s=[x0;x1;x2;x3];
end

function u0= pde1lic(z,x)
	global X1 X2 X3  mark dx
	std=0.01;
	u1=std*rand(1,1)*z(1)+z(1);
	u2=std*rand(1,1)*z(2)+z(2);
	u3=std*rand(1,1)*z(3)+z(3);

	u0=[u1;u2;u3];
end

function [pl,ql,pr,qr]=pde1bc(xl,ul,xr,ur,t)  % boundary condition 
pl=[0;0;0];
pr=[0;0;0];
ql=[1;1;1];
qr=[1;1;1];

end

function [dxF]=Deriv1(F)
Nx=size(F,2);
global dx
dxF = zeros(size(F));
dxF(:,1) = (-1.5*F(:,1) + 2.0*F(:,2) - 0.5*F(:,3)) / dx;
dxF(:,2:Nx-1) = (F(:,3:Nx) - F(:,1:Nx-2) ) / (2*dx);
dxF(:,Nx) = (1.5*F(:,Nx) - 2.0*F(:,Nx-1) + 0.5*F(:,Nx-2)) / dx;
end
