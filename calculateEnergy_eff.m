function [ energy ] = calculateEnergy_eff(C, L, R, Hlocal)

    [~, D] = size(C);

    vl_C_Cbar = bsxfun(@times, transpose(L(1,:)), C);
    vl_C_Cbar = transpose(vl_C_Cbar)*conj(C);
    
    O1 = vl_C_Cbar.*Hlocal{1};
    O2 = transpose(C)*O1*conj(C);
    OL = O2.*eye(D);
    
    O1 = vl_C_Cbar.*eye(D);
    O2 = transpose(C)*O1*conj(C);
    OL = OL + O2.*Hlocal{4};
    
    O1 = vl_C_Cbar.*Hlocal{2};
    O2 = transpose(C)*O1*conj(C);
    OL = OL + O2.*Hlocal{3};
    
    O3 = transpose(OL)*bsxfun(@times, C, R(1,:));
    
    energy = trace(C'*O3); 
    
    energy = real(energy);

end

