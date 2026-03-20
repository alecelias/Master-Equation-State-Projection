function dptilde01dt = ReduOe1SystemRHS(t,ptilde01)

global S Q QT FNorm invFNorm indxpiv

if size(ptilde01,1) == 1
    ptilde01 = ptilde01.';
end

n = length(ptilde01)/2;

ptilde0 = ptilde01(1:n);
ptilde1 = ptilde01(n+1:2*n);

rhstmp = QT*ptilde0;
rhstmp(~indxpiv) = 0;
p0 = invFNorm*rhstmp;

dptilde0dt = Q*S*p0;

rhstmp = QT*dptilde0dt;
rhstmp(~indxpiv) = 0;
dp0dt = invFNorm*rhstmp;

rhstmp = QT*ptilde1;
rhstmp2 = dp0dt - S*p0;
rhstmp(~indxpiv) = rhstmp2(~indxpiv);
p1 = invFNorm*rhstmp;

dptilde1dt = Q*S*p1;

dptilde01dt = [dptilde0dt; dptilde1dt];

return

end