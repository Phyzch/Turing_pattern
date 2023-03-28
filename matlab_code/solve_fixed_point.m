function z=solve_fixed_point(sigma,sigma1,ep,beta,gamma,alpha,beta1,sigma2)

    function F=root3(z)
    um=z(1);
    vm=z(2);
    uc=z(3);
    F(1)=-alpha*um-beta*um*vm^2+gamma*vm+ep*uc- sigma *um+beta1*vm^3;
    F(2)=alpha*um+beta*um*vm^2-gamma*vm- sigma1 *vm-beta1*vm^3+sigma2*uc;
    F(3)=- ep *uc+sigma*um+sigma1*vm-sigma2*uc;
    F(4)=um+vm+uc-400;
    end
z0=[100,190,110];
ff=optimset('display','off');
options=optimoptions('fsolve','Display','none');
z=fsolve(@root3,z0,ff);
end