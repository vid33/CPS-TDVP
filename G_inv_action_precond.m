function [ Ginv_vec ] = G_inv_action_precond( vec, C, L, R, Ltrunc, Rtrunc, maxiter, tolerance, vec_zero )

 %   if nargin < 8
 %       tolerance=1e-12;
 %   end
 %   if nargin < 7
 %       maxiter=10000;
 %   end

    % max num. of trials before bicgstab decides to quit trying to
    % calculate F to prescribed tolerance.
    errorCount_max = 2;
    
    EE = C.*conj(C);
    
    function out =M_action(in)
        [~, D] = size(C);

        in = transpose(reshape(in, D, (D-1)));

        Lambda = bsxfun(@times, (1./transpose(L(1,:))), EE);
        Lambda = bsxfun(@times, Lambda, (1./R(1,:)));
        
        
        T1 = transpose(conj(Ltrunc))*in*conj(R);
        T2 = T1.*Lambda;

        out = (Ltrunc)*T2*transpose(R);
        out = reshape(transpose(out), D*(D-1), 1);
    end

    function out =G_action(in)
        [~, D] = size(C);

        in = transpose(reshape(in, D, (D-1)));

        lambda = bsxfun(@times, transpose(L(1,:)), (1./EE));
        lambda = bsxfun(@times, lambda, R(1,:));
        
        
        T1 = transpose(Rtrunc)*in*L;
        T2 = T1.*lambda;

        out = conj(Rtrunc)*T2*(L');
        out = reshape(transpose(out), D*(D-1), 1);
    end
    
    [Ginv_vec ,flag, relres, iter] =bicgstab(@G_action, vec, tolerance, maxiter, @M_action, [], vec_zero);
   % relres
    errorCount = 0; relres_best = relres; Ginv_vec_best = Ginv_vec;

    while flag ~= 0
        errorCount = errorCount + 1;
        
        if errorCount >= errorCount_max
            fprintf('Giving up trying to G_inv_action_precond to prescribed tolerance\n');
            Ginv_vec = Ginv_vec_best; relres = relres_best;
            fprintf('Using best relres: %d\n', relres);
            break;
        end
        
        fprintf('in G_inv_action precond\n');
        fprintf('Flag: %d, Relres: %d, iter: %d \n', flag, relres, iter);
        fprintf('TRYING bicstab with random vec_zero with 4*maxiter!\n');
        
%        vec_zero_tmp = cpxrand(D*(D-1), 1);   
        vec_zero_tmp = zeros(D*(D-1), 1);
        
        [Ginv_vec, flag, relres, iter] =bicgstab(@G_action, vec, tolerance, 4*maxiter, @M_action, [], vec_zero_tmp);
        
        if relres < relres_best
           relres_best = relres; Ginv_vec_best = Ginv_vec; 
        end
        
        fprintf('Is this any better?\n');        
        if flag ~= 0
            fprintf('NO!\n');
            fprintf('Uncorrected error in G_inv_action precond\n');
            fprintf('Flag: %d, Relres: %d, iter: %d \n', flag, relres, iter);
        else
            fprintf('YES!\n');
            fprintf('Flag: %d, Relres: %d, iter: %d \n', flag, relres, iter);
        end 
        
        
    end
    
end

