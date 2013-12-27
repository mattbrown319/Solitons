function [u_perturbed, w] = perturbSolution(u,xgrid, delta_x, omega, mu, epsilon)        
        s = length(xgrid);
        finiteDifference=-diag(ones(1,s))*2+diag(ones(1,s-1),1)+diag(ones(1,s-1),-1);
        finiteDifference(1,end)=1;
        finiteDifference(end,1)=1;
        finiteDifference=finiteDifference/(delta_x.^2);
        MPotential=((omega^2)*diag(xgrid.^2))/2.0;       %potential
        muMatrix=diag(ones(1,s))*mu;


        f_u = diag((9*abs(u).^4 + 14*abs(u).^2 + 4)./(4*(1+abs(u).^2).^1.5));
        g_u = diag((3*abs(u).^4+4*abs(u).^2)./(4*(1+abs(u).^2).^1.5));
        
        H(1:s,1:s) = -.5*finiteDifference + MPotential - muMatrix + f_u;
        H(1:s,(1:s)+s) = g_u;
        H((1:s)+s,(1:s)+s) = -(-.5*finiteDifference + MPotential - muMatrix + f_u);
        H((1:s)+s,1:s) = -g_u;
        
        [eigenfunctions, w]=eig(H);
        
        maxrealomega=max(real(diag(w)));
        sortedrealomega = sort(abs(real(diag(w))));
        minnonzerorealomega=sortedrealomega(3);
        
        

        w_vector = abs(diag(w));
        w_vector(w_vector <= 0.001) = nan;
        [min_w,min_w_index] = min(w_vector);
        theEigenfunction = eigenfunctions(:,min_w_index);
        p = theEigenfunction(1:length(theEigenfunction)/2);         %p is u.
        q = theEigenfunction((length(theEigenfunction)/2)+1:end);   %q is v.
        u_perturbed = (u + (epsilon*(p + conj(q)))');
end

