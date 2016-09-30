 E = C.*conj(C); 
 
 [Rtmp,DE,Ltmp] = eig(E);
 
 Ltmp = Ltmp';
 Rtmp = transpose(Rtmp);
 
 Diag_tmp = diag(Ltmp*transpose(Rtmp));
 
 Rtmp = transpose(Rtmp)*diag(1./Diag_tmp);
 L = Ltmp;
 R = transpose(Rtmp);
 
 Eeigvals = diag(DE);
 
        
    for kk=1:D
        fprintf('Normtest L kk=%d: %d\n', kk, norm(L(kk,:)*E - Eeigvals(kk)*L(kk,:)));
        
        fprintf('Normtest R kk=%d: %d\n', kk, norm(E*transpose(R(kk, :)) - Eeigvals(kk)*transpose(R(kk,:))));
    end
    
    
Etest  = zeros(D);
for kk =1:D
    Etest = Etest  +Eeigvals(kk)*kron(transpose(R(kk,:)),L(kk,:));
end

fprintf('Testing eig expansion of E, norm(Etest-E)=%d\n', norm(Etest- E));