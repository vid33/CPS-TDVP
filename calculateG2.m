function [ G ] = calculateG2( C, R, L, Rtrunc, Ltrunc )

    [~, D] = size(C);

    E = C.*conj(C);
    
    lambda = diag(L(1,:))*(1./E)*diag(R(1,:));
    
    factor_l = kron(Rtrunc', L');
    
    lambda_vec = reshape(transpose(lambda), D^2, 1);
    
    factor_l = diag(lambda_vec)*factor_l;
    
    factor_r = kron(transpose(Rtrunc), transpose(L));
    
    G = transpose(factor_l)*factor_r;
   


end

