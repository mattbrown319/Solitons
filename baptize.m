function u_baptized = baptize( omega, x_start, xgrid, delta_x)

%omega=.1 and mu = 1 work
%omega=.25 and mu = 2 work
numSolitons = 0;
mu = 2;


    u_steady = findSteadySolution(mu,numSolitons,omega,delta_x, xgrid);       %find steady state solution for the ground state
    %plot(xgrid,u_pdf(u_steady));
    x_start_vector=repmat(x_start,1,length(xgrid));
    u_baptized = u_steady .* (tanh(xgrid - x_start_vector));

end