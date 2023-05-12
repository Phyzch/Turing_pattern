function z=solve_fixed_point(k13,k23,k31,beta_12,k21,k12,beta_21,k32)

    function F=root3(z)
    X1=z(1);
    X2=z(2);
    X3=z(3);
    F(1) = -k12 * X1 - beta_12 * X1 * X2^2 + k21 * X2 + k31 * X3 - k13 *X1 +beta_21 * X2^3;
    F(2) = k12 * X1 + beta_12 * X1 * X2^2 - k21 * X2 - k23 * X2 - beta_21 * X2^3 + k32 * X3;
    F(3) = - k31 * X3 + k13 * X1 + k23 * X2 - k32 * X3;
    F(4) = X1+X2+X3-400;
    end
z0=[100,190,110];
ff=optimset('display','off');
options=optimoptions('fsolve','Display','none');
z=fsolve(@root3,z0,ff);
end