function u_steady = findSteadySolution( mutarget, numSolitons, omega, delta_x, xgrid)

%SHO solutions, first four quantum states;
if numSolitons==0,
    B=exp(-(omega*xgrid.^2)/2);
elseif numSolitons==1,
    B=sqrt(2)*xgrid.*exp(-(omega*xgrid.^2)/2);
elseif numSolitons==2,
    B=(2*xgrid.^2-1)/sqrt(2).*exp(-(omega*xgrid.^2)/2);
elseif numSolitons==3,
    B=(2*xgrid.^3-3*xgrid)/sqrt(3).*exp((-omega*xgrid.^2)/2);
end

%u_quantum harmonic oscillator
uqho=(omega/pi)^(1/4)*B;         %wave function initial guess for newton
u = uqho;

%perhaps adaptive method for traversing mu
setofMus = linspace((numSolitons+1/2)*omega, mutarget, 20);

s = length(xgrid);
finiteDifference=-diag(ones(1,s))*2+diag(ones(1,s-1),1)+diag(ones(1,s-1),-1);       
finiteDifference=finiteDifference/(delta_x.^2);
MPotential=((omega^2)*diag(xgrid.^2))/2.0;       %potential



for mu = setofMus
    
    muMatrix=diag(ones(1,s))*mu;
    npse_without_nonlinearity=-.5*finiteDifference-muMatrix+MPotential;


    tolerance=1;    %c represents the error.  of course.
    while tolerance>1e-10,
        f0=(npse_without_nonlinearity)*u'+diag((3*abs(u).^2+2)./(2*((1+abs(u).^2).^.5)))*(u');
        F=(npse_without_nonlinearity)+diag((6*u.^4+9*abs(u).^2+2)./(2*(1+abs(u).^2).^1.5));  %jacobian
        du=F\f0;
        u1=u-du';
        tolerance=norm(u1-u);
        u=u1;
    end
    
    if (max(abs(u))<0.07)
            u=uqho;
    end

end   

u_steady = u;
    
end

