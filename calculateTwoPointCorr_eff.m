function [ corr ] = calculateTwoPointCorr_eff(C, L, R, twoPointLocal, PP)

    [~, D] = size(C);

    vl_C_Cbar = bsxfun(@times, transpose(L(1,:)), C);
    vl_C_Cbar = transpose(vl_C_Cbar)*conj(C);
    
    O1 = vl_C_Cbar.*twoPointLocal{end};
    O2 = transpose(C)*O1*conj(C);
    OL = O2.*eye(D);
    
    O1 = vl_C_Cbar.*eye(D);
    O2 = transpose(C)*O1*conj(C);
    OL = OL + O2.*twoPointLocal{end};
    
    
    for kk=1:(PP-1)
        O1 = vl_C_Cbar.*twoPointLocal{kk};
        O2 = transpose(C)*O1*conj(C);
        OL = OL + O2.*twoPointLocal{PP-1-(kk-1)};
    end
    
    
    O3 = transpose(OL)*bsxfun(@times, C, R(1,:));
    
    corr = trace(C'*O3); 
    
    %corr = real(corr);

end

