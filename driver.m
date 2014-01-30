clear,clc

% set initial parameters
%numSolitons = 1;
delta_t = 0.01;     %change in t.
maxtime = 300;       %maxTime
delta_x = 0.2;       %change in x.  'h'.
x_length = 100;        %total length across xgrid
x_center = 0;          
omega = .05;       %potential trap constant
mu = 2;
x_start = [2.0];
numSolitons = length(x_start);
if(mu < omega/2)
    print('THIS ISN''T A VALID MU VALUE');
end


plotlive = 1;
saveplot = 1;
numSteps = maxtime/delta_t;     %number of steps forward in RK4
trackProgress = 1;
showProgressPercentInterval = .5;
modForProgress = showProgressPercentInterval*numSteps*0.01;
xgrid=x_center - (x_length/2):delta_x:x_center + (x_length/2);         %xgrid



u_xt = zeros((x_length/delta_x)+1,maxtime/delta_t);
%E_t = zeros(1,maxtime/delta_t);
x_t = zeros(1,maxtime/delta_t);
    

u_steady = findSteadySolution(mu,0,omega,delta_x, xgrid);       %newton to find steady state solution
u_baptized = baptize(u_steady, omega, mu, x_start, xgrid, delta_x);
u_xt(:,1) = u_baptized';
%E_t(1) = trapz(xgrid,u_pdf(u_xt(:,1)));
x_t(1) = centerofmass( u_steady, u_xt(:,1), xgrid );

starttime = now;            %records the starting time as a float for end estimate
index = 1;
while index < numSteps
    averageTimePerIteration = (now - starttime)/index;
    
    u_xt(:,index+1) = onestep_rk4( u_xt(:,index), xgrid', delta_x, index*delta_t, delta_t, omega);
    %E_t(index+1) = trapz(xgrid,u_pdf(u_xt(:,index)));
    x_t(index+1) = centerofmass( u_steady, u_xt(:,index), xgrid );

    if(trackProgress == 1 && (mod(index+1, modForProgress) == 0 ))
        disp([num2str(100*((index+1)/numSteps)) '%']);
        disp(['Start time:  ' datestr(starttime,'mm/dd HH:MM PM')]);
        disp(['Elapsed time:  ' datestr(now - starttime,'HH:MM:SS')]);
        disp(['Estimated time left:  ' datestr((numSteps -  index)*averageTimePerIteration,'HH:MM:SS')]);
        disp(['Estimated end time: ' datestr((numSteps -  index)*averageTimePerIteration+starttime,'mm/dd HH:MM PM')]);
    end


    if(plotlive == 1)
        drawnow;
        subplot(2,2,1);
        plot(xgrid,0.5*(omega^2)*(xgrid.^2),xgrid,u_pdf(u_xt(:,index+1))),'-o';
        %drawnow;
        %subplot(2,2,2);
        %plot(E_t(2:index+1));
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
end

%guessWavelength = floor((1/((.21)*(1/sqrt(2))*(omega)))/delta_t);
%[S, F, T, P] = spectrogram(x_t(1:index),guessWavelength, floor(guessWavelength/2), length(x_t(1:index)), 1/delta_t);
%surf(T,F,P,'edgecolor','none');
%axis tight;
%view(0,90);
%xlabel('Time (Seconds)'); ylabel('Hz');
%[maxP, indexMaxP] = max(P);
%plot(x_t(1:floor((1/F(indexMaxP))/delta_t)));








%N = length(x_t(1:index));
%fspace = linspace(1/(N*delta_t),(N-1)/(N*delta_t),N);
%plot(fspace,abs(fft(x_t(1:index))));


%imagesc(u_pdf(u_xt));
%plot(x_t)
%max(abs(fft(x_t')))



%need to set an omega, mu value, and number of solitons.
%then start at mu = omega/sqrt(2) (i think?) and slowly walk forward
%through mu values using each steady state as the basis for the next newton
%calculation until you hit the mu value you actually want. THIS is the
%steady state you should use for moving forward.



