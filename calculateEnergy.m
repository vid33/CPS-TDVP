function [ energy ] = calculateEnergy(C, vl, vr, J, h)

    [~, D] = size(C);
    siteOverlap = log2(D);

    vl_C_Cbar = diag(vl)*C;
    vl_C_Cbar = transpose(vl_C_Cbar)*conj(C);
    vl_C_Cbar_vec = reshape(transpose(vl_C_Cbar), 1, D^2); 
  
    C_Cbar_vr = C*diag(vr);
    C_Cbar_vr = C_Cbar_vr*C';
    C_Cbar_vr_vec = reshape(transpose(C_Cbar_vr), D^2, 1);  

    [ C_H_Cbar  ] = calculate_H_effective( C, C, siteOverlap, D, J, h );
      
    energy= vl_C_Cbar_vec*C_H_Cbar*C_Cbar_vr_vec;
    energy = real(energy);

end

