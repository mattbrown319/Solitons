function u_baptized = baptize( u, omega, mu, x_start, xgrid, delta_x, plotbaptized)

%omega=.1 and mu = 1 work
%omega=.25 and mu = 2 work
numSolitons = length(x_start);
index = 1;
while(index <= numSolitons)
   x_start_vector=repmat(x_start(index),1,length(xgrid));

   u = u .* (tanh(xgrid - x_start_vector));

        
   index = index + 1;
end
    u_baptized = u;

    %plot(xgrid,u_pdf(u_steady));
    
    
    %u_baptized = (u' * diag(tanh(xgrid - x_start_vector)))';
    %u_baptized = u .* 1;
    if(plotbaptized == 1)
        plot(xgrid,u_pdf(u_baptized));
    end
end