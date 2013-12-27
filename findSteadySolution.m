function u_steady = findSteadySolution( mu, numSolitons, omega, delta_x, xgrid )
s = length(xgrid);
finiteDifference=-diag(ones(1,s))*2+diag(ones(1,s-1),1)+diag(ones(1,s-1),-1);
finiteDifference(1,end)=1;
finiteDifference(end,1)=1;
finiteDifference=finiteDifference/(delta_x.^2);
MPotential=((omega^2)*diag(xgrid.^2))/2.0;       %potential

%constant in front of each solution
A=(omega/pi)^(1/4);
 
    
    %SHO solutions, first four quantum states;
    
    if numSolitons==0,
        %B=exp(-omega*(x.^2)/2);
        B=exp(-((omega*xgrid).^2)/2);
    elseif numSolitons==1,
        %B=sqrt(2)*x.*exp(-omega*(x.^2)/2);
        B=sqrt(2)*xgrid.*exp(-((omega*xgrid).^2)/2);
    elseif numSolitons==2,
        %B=(2*x.^2-1)/sqrt(2).*exp(-omega*(x.^2)/2);
        B=(2*xgrid.^2-1)/sqrt(2).*exp(-((omega*xgrid).^2)/2);
    elseif numSolitons==3,
        B=(2*xgrid.^3-3*xgrid)/sqrt(3).*exp(-omega*(xgrid.^2)/2);
    end
    
    u=A*B;         %wave function initial condition
    

        muMatrix=diag(ones(1,s))*mu;
        Mc=-.5*finiteDifference-muMatrix+MPotential;         %finite difference method - mu + potential
        
        tolerance=1;    %c represents the error.  of course.
        
        while tolerance>1e-10,
            
            Mnonlinear=diag((3*u.^2+2)./(2*((1+u.^2).^.5)));         %matrix of nonlinearity
            dfdu = diag((6*u.^4+9*u.^2+2)./(2*(1+u.^2).^1.5));   %derivative of nonlinear term
            
            Mv=Mc+Mnonlinear;
            Mj=Mc+dfdu;        %Mj is the jacobian
            
            f0=Mv*(u');
            F=Mj;
            du=F\f0;
            u1=u-du';
            tolerance=norm(u1-u);
            u=u1;
        end
        
        
        
        
        
        
        
        u_steady = u;
end

