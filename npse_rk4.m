function u_xt = npse_rk4( u, xgrid, delta_x, delta_t, maxtime, omega )
u_xt = u';
numSteps = maxtime/delta_t;

index = 1;
while index < numSteps
    u_xt = [u_xt, onestep_rk4( u_xt(:,end), xgrid', delta_x, index*delta_t, delta_t, omega)];
    index = index + 1;
end

end

