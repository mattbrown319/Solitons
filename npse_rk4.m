function u_xt = npse_rk4( u, xgrid, delta_x, delta_t, maxtime, omega, plotlive )
u_xt = u';
numSteps = maxtime/delta_t;

index = 1;
while index < numSteps
    u_xt = [u_xt, onestep_rk4( u_xt(:,end), xgrid', delta_x, index*delta_t, delta_t, omega)];
    if(plotlive == 1)
        plot(xgrid,u_pdf(u_xt(:,index))),'-o';
        drawnow;
    end
    index = index + 1;
end

end