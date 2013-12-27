function ik = f_for_rk4( u, x, del_x, t, del_t, gamma )
    central_diff=-diag(ones(1,length(x)))*2+diag(ones(1,length(x)-1),1)+diag(ones(1,length(x)-1),-1);
    central_diff=central_diff/(del_x^2);       %was /2*h.^2, but we're doing the division by 2 later
     u_xx=(-1/2)*central_diff*u;
    potential=(gamma^2/2)*((x).^2).*u;
    nonlinear=((3*abs(u).^2+2).*u)./(2*(1+abs(u).^2).^(0.5));
    
    ik = u_xx+potential+nonlinear;       %might be x*x' and not x'*x


end