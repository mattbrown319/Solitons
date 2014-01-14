function [u_xt, E_t, x_t] = step_forward( u, xgrid, delta_x, delta_t, maxtime, omega, numSolitons, plotlive )
u_xt = u';
numSteps = maxtime/delta_t;
E_t = [];
x_t = [];
index = 1;
while index < numSteps
    u_xt = [u_xt, onestep_rk4( u_xt(:,end), xgrid', delta_x, index*delta_t, delta_t, omega)];
    E_t = [E_t trapz(xgrid,u_pdf(u_xt(:,index)))];
    x_t = [x_t centerofmass( u_xt(:,index), xgrid, numSolitons )];

    if(plotlive == 1)
        drawnow;
        subplot(2,2,1);
        plot(xgrid,0.5*(omega^2)*(xgrid.^2),xgrid,u_pdf(u_xt(:,index))),'-o';
        %plot(xgrid,u_pdf(u_xt(:,index))),'-o';
        drawnow;
        subplot(2,2,2);
        plot(E_t);
        drawnow;
        subplot(2,2,3);
        plot(x_t);
        drawnow;
        subplot(2,2,4);
        imagesc(u_pdf(u_xt));
        %axis([0 numSteps 7 8]); %this needs to be dynamically set
        
    end
    index = index + 1;
end

end