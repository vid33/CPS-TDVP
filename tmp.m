[ Ctmp, Dtmp, Ltmp, Rtmp ] = calculateEigenvectors_inv(C);

EE = C.*conj(C);

norm(EE*transpose(Rtmp) - transpose(Rtmp)*diag(Dtmp))

norm(Ltmp*EE - diag(Dtmp)*Ltmp)