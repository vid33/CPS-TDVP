 E = C.*conj(C); 
    [Rtmp, DE] = eig(E);
    Ltmp = inv(Rtmp');
    Eeigvals = diag(DE);
    
     L = Ltmp';
        R = Rtmp';
        
    for kk=1:D
        fprintf('Normtest L kk=%d: %d\n', kk, norm(L(kk,:)*E - Eeigvals(kk)*L(kk,:)));
        
        fprintf('Normtest R kk=%d: %d\n', kk, norm(E*transpose(R(kk, :)) - Eeigvals(kk)*transpose(R(kk,:))));
    end
    
    
Etest  = zeros(D);
for kk =1:D
    Etest = Etest  +Eeigvals(kk)*kron(transpose(R(kk,:)),L(kk,:));
end

fprintf('Testing eig expansion of E, norm(Etest-E)=%d\n', norm(Etest- E));