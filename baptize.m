function u_baptized = baptize( u, omega, mu, x_start, xgrid, delta_x, plotbaptized)

%omega=.1 and mu = 1 work
%omega=.25 and mu = 2 work



    %plot(xgrid,u_pdf(u_steady));
    x_start_vector=repmat(x_start,1,length(xgrid));
    
    %u_baptized = (u' * diag(tanh(xgrid - x_start_vector)))';
    u_baptized = u .* (tanh(xgrid - x_start_vector));
    %u_baptized = u .* 1;
    if(plotbaptized == 1)
        plot(xgrid,u_pdf(u_baptized));
    end
end