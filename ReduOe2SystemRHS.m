function dptilde012dt = ReduOe2SystemRHS(t,ptilde012)

global S Q QT FNorm invFNorm indxpiv

if size(ptilde012,1) == 1
    ptilde012 = ptilde012.';
end

n = length(ptilde012)/3;

ptilde0 = ptilde012(1:n);
ptilde1 = ptilde012(n+1:2*n);
ptilde2 = ptilde012(2*n+1:3*n);

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

ddptilde0dt2 = Q*S*dp0dt;

rhstmp = QT*ddptilde0dt2;
rhstmp(~indxpiv) = 0;
ddp0dt2 = invFNorm*rhstmp;

rhstmp = QT*dptilde1dt;
rhstmp2 = ddp0dt2 - S*dp0dt;
rhstmp(~indxpiv) = rhstmp2(~indxpiv);
dp1dt = invFNorm*rhstmp;

rhstmp = QT*ptilde2;
rhstmp2 = dp1dt - S*p1;
rhstmp(~indxpiv) = rhstmp2(~indxpiv);
p2 = invFNorm*rhstmp;

dptilde2dt = Q*S*p2;

dptilde012dt = [dptilde0dt; dptilde1dt; dptilde2dt];

return

end