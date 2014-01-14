%input:
%t:  the current timestep
%u:  a 1xp vector
%h: change in x
%del_t: change in t
%x:  vector of all x values being used
function u_next = onestep_rk4( u_curr, x, del_x, t, del_t, omega  )
    imaginary = complex(0,1);
    
    ik1 = del_t*f_for_rk4(u_curr, x, del_x, t, del_t, omega);
    ik2 = del_t*f_for_rk4(u_curr+(1/2)*(-imaginary)*ik1,x,del_x,t+(1/2)*del_t,del_t, omega);
    ik3 = del_t*f_for_rk4(u_curr+(1/2)*(-imaginary)*ik2,x,del_x,t+(1/2)*del_t,del_t, omega);
    ik4 = del_t*f_for_rk4(u_curr+(-imaginary)*ik3,x,del_x,t+del_t,del_t, omega);
    
    
    
    u_next = u_curr + -imaginary*(1/6)*(ik1 + 2*ik2 + 2*ik3 + ik4);
end