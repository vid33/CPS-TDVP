function [ C_H_Cbar  ] = calculate_H_effective( C, M, siteOverlap, D, J, h )
    sx = [0,1;1,0];  sz = [1,0;0,-1]; 

    sz_vec_l = cell(siteOverlap, 1);
    sz_vec_r = cell(siteOverlap, 1);
    
    for kk=1:siteOverlap
        sz_vec_l{kk} = reshape(transpose( kron(kron(eye(2^(kk-1)), sz), eye(2^(siteOverlap-kk))) ), 1, D^2);
        sz_vec_r{kk} = reshape(transpose( kron(kron(eye(2^(kk-1)), sz), eye(2^(siteOverlap-kk))) ), D^2, 1);
    end
    
    
    if siteOverlap > 1
    sx_sx_vec_l = cell(siteOverlap-1, 1);
    sx_sx_vec_r = cell(siteOverlap-1, 1);
        for kk=1:(siteOverlap-1)
            sx_sx_vec_l{kk} = reshape(transpose( kron(kron(eye(2^(kk-1)), kron(sx,sx)), eye(2^(siteOverlap-1-kk))) ), 1, D^2);
            sx_sx_vec_r{kk} = reshape(transpose( kron(kron(eye(2^(kk-1)), kron(sx,sx)), eye(2^(siteOverlap-1-kk))) ), D^2, 1);
        end
    end
    
    sx_split_vec_l =  reshape(transpose(kron(eye(2^(siteOverlap-1)), sx)), 1, D^2);
    sx_split_vec_r =  reshape(transpose(kron(sx, eye(2^(siteOverlap-1)))), D^2, 1);

    eye_vec_l = reshape(eye(D), 1, D^2); eye_vec_r = reshape(eye(D), D^2, 1);
    
    C_kron_Cbar = kron(C, conj(M));

    C_H_Cbar_sz_l = cell(siteOverlap, 1);
    C_H_Cbar_sz_r = cell(siteOverlap, 1);
    
    if siteOverlap > 1
        C_H_Cbar_sx_sx_l = cell(siteOverlap-1, 1);
        C_H_Cbar_sx_sx_r = cell(siteOverlap-1, 1);
    end
    
    C_H_Cbar_sx_split = diag(sx_split_vec_l)*C_kron_Cbar;
    C_H_Cbar_sx_split = C_H_Cbar_sx_split*diag(sx_split_vec_r);
  
    for mm=1:siteOverlap
        C_H_Cbar_sz_l{mm} = diag(sz_vec_l{mm})*C_kron_Cbar;
        C_H_Cbar_sz_l{mm} = C_H_Cbar_sz_l{mm}*diag(eye_vec_r);
        
        C_H_Cbar_sz_r{mm} = C_kron_Cbar*diag(sz_vec_r{mm});
        C_H_Cbar_sz_r{mm} = diag(eye_vec_l)*C_H_Cbar_sz_r{mm};   
    end
    
    for mm=1:(siteOverlap-1)
        
        C_H_Cbar_sx_sx_l{mm} = diag(sx_sx_vec_l{mm})*C_kron_Cbar;
        C_H_Cbar_sx_sx_l{mm} = C_H_Cbar_sx_sx_l{mm}*diag(eye_vec_r);
        
        C_H_Cbar_sx_sx_r{mm} = C_kron_Cbar*diag(sx_sx_vec_r{mm});
        C_H_Cbar_sx_sx_r{mm} = diag(eye_vec_l)*C_H_Cbar_sx_sx_r{mm};
    
    end
    
    C_H_Cbar = -(J/siteOverlap)*C_H_Cbar_sx_split;
    
    for kk=1:(siteOverlap-1)
        C_H_Cbar = C_H_Cbar - (J/(2*siteOverlap))*(C_H_Cbar_sx_sx_l{kk} + C_H_Cbar_sx_sx_r{kk});
    end
    
    for kk=1:siteOverlap
        C_H_Cbar = C_H_Cbar + (h/(2*siteOverlap))*(C_H_Cbar_sz_l{kk} + C_H_Cbar_sz_r{kk});
    end
end

