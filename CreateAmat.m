function outM1 = CreateAmat()

global kads kdes krxn2 Ncoord NL

kads = 8.5e+5;
kdes = 4.0e+5;
krxn2 = 1.E+06;
Ncoord = 4;
NL = 25;

for i = 0:NL
    indx = i+1;
    if i >= 1;
        outM1(indx,indx-1) = Aads(i-1);
    end
    if i <= NL-1;
        outM1(indx,indx+1) = Ades(i+1);
    end
    if i <= NL-2;
        outM1(indx,indx+2) = Arxn2(i+2);
    end
    outM1(indx,indx) = -Aads(i)-Ades(i)-Arxn2(i);    
end

outM1(end,:) = 1;

end

function out1 = Aads(n)

global kads NL

out1 = kads*(NL-n);

if out1 < 0;
    out1 = 0;
end

end

function out1 = Ades(n)

global kdes NL

out1 = kdes*n;

if out1 < 0;
    out1 = 0;
end

end

function out1 = Arxn2(n)

global krxn2 Ncoord NL

out1 = krxn2/2*Ncoord/(NL-1)*n*(n-1);
% out1 = krxn2*Ncoord/(NL-1)*n*(n-1);

if out1 < 0;
    out1 = 0;
end

end

