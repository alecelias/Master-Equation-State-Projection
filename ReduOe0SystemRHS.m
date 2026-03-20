function dptilde0dt = ReduOe0SystemRHS(t,ptilde0)

global S Q QT FNorm invFNorm indxpiv

if size(ptilde0,1) == 1
    ptilde0 = ptilde0.';
end

rhstmp = QT*ptilde0;
rhstmp(~indxpiv) = 0;
p0 = invFNorm*rhstmp;

dptilde0dt = Q*S*p0;

return

end