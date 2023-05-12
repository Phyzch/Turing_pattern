DX1 = 0.3;
DX3 = 0.3;
DX2 = 0.012;
k13 = log(2)/50*1;
k23 = log(2)/50*1;
k32 = 1e-3 * k23;
k31 = 3*log(2)/50;
k21 = 3.6;
k12 = 0.5;
beta_12 = 1.5*1e-4;
beta_21 = 0.01*beta_12;
L = 6;

D_trialsize = 50;
eigen = zeros(1,20);


for i = 1:30
    beta_21 = (0.65 - 0.02*i) * beta_12;
	z = solve_fixed_point(k13, k23, k31, beta_12, k21, k12, beta_21, k32);
	X1_0 = z(1);
	X2_0 = z(2);
	X3_0 = z(3);
	
    for n=1:1:20
		k =(n-1)*pi/L;
		A =[-k^2 * DX1 - k12 - beta_12 * X2_0^2 - k13, -2 * beta_12 * X1_0 * X2_0 + k21 + 3 * beta_21 * X2_0^2, k31;
			k12 + beta_12*X2_0^2, -k^2 * DX2 + 2 * beta_12 * X1_0 * X2_0 - k21 - k23 - 3 * beta_21 * X2_0^2, k32;
			k13, k23, -k31 - k^2 * DX3 - k32];
		e =eig(A);
		e =sort(e);
		eigen(1,n) =e(end);
	end
	
	E =max(eigen);
	
	if(E >0)
		break;
	end
	
end
beta_21/beta_12
