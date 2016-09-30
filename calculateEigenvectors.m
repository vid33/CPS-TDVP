%Calculate lr eigenvectors, and normalise stuff

function [ C, Eeigvals, L, R ] = calculateEigenvectors(C)

    E = C.*conj(C); 
    [Rtmp,DE,Ltmp] = eig(E);
    Eeigvals = diag(DE);
    Ltmp = Ltmp';
    Rtmp = transpose(Rtmp);
    %normalise eigvectors
    Diag_tmp = diag(Ltmp*transpose(Rtmp));
    Rtmp = transpose(Rtmp)*diag(1./Diag_tmp);
    L = Ltmp;
    R = transpose(Rtmp);
    fprintf('1- Max eigenvalue is %d\n', 1-max(Eeigvals));
    C = C/sqrt(max(Eeigvals));
    
        E = C.*conj(C); 
    [Rtmp,DE,Ltmp] = eig(E);
        Eeigvals = diag(DE);
    Ltmp = Ltmp';
    Rtmp = transpose(Rtmp);
    %normalise eigvectors
    Diag_tmp = diag(Ltmp*transpose(Rtmp));
    Rtmp = transpose(Rtmp)*diag(1./Diag_tmp);
    L = Ltmp;
    R = transpose(Rtmp);
    
  fprintf('max(Eeigvals) is: %d\n', max(Eeigvals));
    
    if max(Eeigvals) == Eeigvals(1)
        sortEigvals = 0;
    else
        sortEigvals = 1;
    end
   
    
    Eeigvals = Eeigvals/max(Eeigvals);
 %   lambda = lambda/max(lambda);
   
    %fprintf('Max lambda is %d, min lambda is %d\n', max(lambda), lm);
    
    if sortEigvals == 0
     %   L = Ltmp';
     %   R = Rtmp';
        
     %   diag(L*R')
        
        fprintf('Not sorting the eigenvalues\n');
        
    elseif sortEigvals == 1;
    
        fprintf('Were have sorted your eigenvalues...\n');
        
        L = zeros(numel(Eeigvals), numel(Eeigvals));
        R = zeros(numel(Eeigvals), numel(Eeigvals));
    
        [Eeigvals, eigvals_index] = sort(Eeigvals,'descend'); 
        for kk=1:numel(Eeigvals)
            R(kk, :) = Rtmp(:, eigvals_index(kk)); 
            L(kk, :) = Ltmp(:, eigvals_index(kk))';
        end

    end
    
end

