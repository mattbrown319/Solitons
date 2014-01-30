function x_t = centerofmass( u_steady, u_x, xgrid )
%need to incorporate the number of solitons somehow in this.  maybe using
%higher moments, or the variance of u_xt?  or just cut it off at x=0, find
%the center of mass for all positive x, and mirror that across y axis.

    mass = u_pdf(u_steady) - u_pdf(u_x)';

 x_t = trapz(xgrid,xgrid.*mass)/trapz(xgrid,mass); 
     % x_t = -trapz(xgrid,xgrid'.*(abs(u_x).^2))/trapz(xgrid,(abs(u_x)).^2); 

end