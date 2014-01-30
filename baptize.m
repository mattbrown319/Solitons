function u_baptized = baptize( u, omega, mu, x_start, xgrid, delta_x)


%numSolitons = length(x_start);
index = 1;
%while(index <= numSolitons)
   x_start_vector=repmat(x_start(index),1,length(xgrid));
   u = u .* (tanh(xgrid - x_start_vector));
   index = index + 1;
%end
    u_baptized = u;

end