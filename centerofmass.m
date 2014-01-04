function x_t = centerofmass( u_xt, xgrid )

index = 1;
x_t = [];
numSteps = length(u_xt);

    while index < numSteps
      x_t_onestep = trapz(xgrid,xgrid'.*(abs(u_xt(:,index)).^2))/trapz(xgrid,(abs(u_xt(:,index))).^2); 
      x_t = [x_t x_t_onestep];
      index = index + 1;
    end

end