clear,clc

%set initial parameters
numSolitons = 3;
if numSolitons==0,
    mu = 2;
elseif numSolitons==1,
    mu = 3;
elseif numSolitons==2,
    mu = 3;
elseif numSolitons==3,
    mu = 4;
end

epsilon = 2;      %perturbation strength
delta_t = 0.01;     %change in t.
maxtime = 20;       %maxTime
delta_x=0.1;       %change in x.  'h'.
omega = 0.25;        %potential trap constant
x_length = 30;        %total length across xgrid
x_center = 0;
%derived parameters
numSteps = maxtime/delta_t;     %number of steps forward in RK4
xgrid=x_center - (x_length/2):delta_x:x_center + (x_length/2);         %xgrid

u_steady = findSteadySolution(mu,numSolitons,omega,delta_x, xgrid);       %newton to find steady state solution
%plot(xgrid, abs(u_steady));
[u_perturbed, w] = perturbSolution(u_steady, xgrid, delta_x, omega, mu, epsilon);
%plot(xgrid, abs(u_perturbed),xgrid,abs(u_steady));
u_xt = npse_rk4(u_perturbed, xgrid, delta_x, delta_t, maxtime, omega );
imagesc((abs(u_xt)').^2)






%build linear operator matrix
%calc eigenvalues (w) and eigenfunctions for linearization
%choose an eigenvalue and eigenfunction and build perturbation 
%add perturbation to steady state solution
%rk4 to step new solution forward in time
%plot this - heat map - should have frequency w.
%calc center of mass with respect ot time x(t)
%plot this
%fft on x(t) and take strongest frequency.  this should be w.

%which of this could be functions and which could be the driver code?



%newton function that takes in mu and numSolitons and potential trap constant and outputs a steady state
%perturb function that takes in steady state and outputs a perturbed state.
    %involves building H and finding smallest eigenvalue/eigenfunction pair
%rk4 function takes in an arbitrary initial solution (with perturbation or
    %without) and outputs a matrix of that solution stepped T steps through
    %time.  
%center of mass function takes in matrix of solutions and output x(t)
    %values for t up to T.
%fft function to take in an x(t) function and output its strongest
    %frequency.
    
%so final program looks like:

%parameters
%u_steady = findSteadySolution(mu,numSolitons)
%[u_perturbed, w] = perturbSolution(u)
%u_xt = npse_rk4(u_perturbed)
%x_t = center_of_mass(u_xt)
%w_experimental = npse_fft(x_t)
%
%test:
%does w =  w_experimental?
%OR IT AT LEAST ON THE SAME ORDER OF MAGNITUDE?