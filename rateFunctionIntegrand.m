function out = rateFunctionIntegrand(k, t, J, h_init, h)


    theta_k_h_init = (1/2)*atan( sin(k).*(1./(h_init - cos(k))) );
    theta_k_h = (1/2)*atan( sin(k).*(1./(h - cos(k))) );
    
    phi_k = theta_k_h_init - theta_k_h;

    e_k_h = 2*J*sqrt( (h-cos(k)).^2 + sin(k).^2);
    
    out = (-1/(2*pi))*log( cos(phi_k).^2 + exp(-2*1i*t*e_k_h).*(sin(phi_k).^2));


end

