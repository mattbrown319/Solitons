clear,clc

% set initial parameters
%numSolitons = 1;
delta_t = 0.01;     %change in t.
maxtime = 500;       %maxTime
delta_x = 0.1;       %change in x.  'h'.
x_length = 15;        %total length across xgrid
x_center = 0;          
x_start = [0.1];
omega = .5;        %potential trap constant
mu = 2.0;

u_xt = zeros((x_length/delta_x)+1,maxtime/delta_t);
E_t = zeros(1,maxtime/delta_t);
x_t = zeros(1,maxtime/delta_t);

numSolitons = length(x_start);
plotlive = 0;
saveplot = 0;
numSteps = maxtime/delta_t;     %number of steps forward in RK4
trackProgress = 1;
showProgressPercentInterval = 1;
modForProgress = showProgressPercentInterval*numSteps*0.01;
xgrid=x_center - (x_length/2):delta_x:x_center + (x_length/2);         %xgrid



%omega_array = linspace(0.05,.20,16);
%mu_array = linspace(1.8,2.2,10);
%omega_array = [0.22632];
%mu_array = [2.0];
%for omega_index = 1:length(omega_array)
%      for mu_index = 1:length(mu_array)
 %         omega = omega_array(omega_index);
 %         mu = mu_array(mu_index);

           


        u_steady = findSteadySolution(mu,0,omega,delta_x, xgrid);       %newton to find steady state solution
        u_baptized = baptize(u_steady, omega, mu, x_start, xgrid, delta_x);
        u_xt(:,1) = u_baptized';
        E_t(1) = trapz(xgrid,u_pdf(u_xt(:,1)));
        x_t(1) = centerofmass( u_xt(:,1), xgrid, numSolitons );
        
        index = 1;
        while index < numSteps
            u_xt(:,index+1) = onestep_rk4( u_xt(:,index), xgrid', delta_x, index*delta_t, delta_t, omega);
            E_t(index+1) = trapz(xgrid,u_pdf(u_xt(:,index)));
            x_t(index+1) = centerofmass( u_xt(:,index), xgrid, numSolitons );

            if(trackProgress == 1 && (mod(index+1, modForProgress) == 0 ))
                [num2str(100*((index+1)/numSteps)) '%']
            end

            
            if(plotlive == 1)
                drawnow;
                subplot(2,2,1);
                plot(xgrid,0.5*(omega^2)*(xgrid.^2),xgrid,u_pdf(u_xt(:,index+1))),'-o';
                drawnow;
                subplot(2,2,2);
                plot(E_t(2:index+1));
                drawnow;
                subplot(2,2,3);
                plot(x_t(2:index+1));
                drawnow;
                subplot(2,2,4);
                imagesc(u_pdf(u_xt(:,2:index+1)));

            end
            index = index + 1;
        end


                if(saveplot == 1)
                       imagesc(u_pdf(u_xt));
                       print( '-djpeg', ['crazyplots/','mu=',num2str(mu),',','omega=',num2str(omega),'.jpg']);
                       %for windows the command will probably be the
                       %following instead:
                       %print( '-djpeg', ['crazyplots\','mu=',num2str(mu),',','omega=',num2str(omega),'.jpg'])
                       %note that all i did was change the direction of the
                       %slash... that's a big difference between mac and
                       %windows file systems.


                end

%     end
      
%end





%N = length(x_t);
%fspace = linspace(1/(N*delta_t),(N-1)/(N*delta_t),N);
%plot(fspace,abs(fft(x_t)));


%imagesc(u_pdf(u_xt));
%plot(x_t)
%max(abs(fft(x_t')))


guessWavelength = floor((1/((.21)*(1/sqrt(2))*(omega)))/delta_t);
[S, F, T, P] = spectrogram(x_t,guessWavelength, floor(guessWavelength/2), length(x_t), 1/delta_t);
surf(T,F,P,'edgecolor','none');
axis tight;
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');
[maxP, indexMaxP] = max(P);
plot(x_t(1:floor((1/F(indexMaxP))/delta_t)));
trapz(x_t(1:floor((1/F(indexMaxP))/delta_t)))


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