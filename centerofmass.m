function x_t = centerofmass( u_x, xgrid, numSolitons )
%need to incorporate the number of solitons somehow in this.  maybe using
%higher moments, or the variance of u_xt?  or just cut it off at x=0, find
%the center of mass for all positive x, and mirror that across y axis.




      x_t = -trapz(xgrid,xgrid'.*(abs(u_x).^2))/trapz(xgrid,(abs(u_x)).^2); 

end