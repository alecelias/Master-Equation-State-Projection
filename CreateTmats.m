function [F, S] = CreateTmats(NL, kdiff, kads, kdis, krx2)
    %CREATEFMAT Creates fast reaction matrix (diffusion reactions)
    %   Detailed explanation goes here
    
    N = 2^NL;
    F = zeros(N, N);
    S = zeros(N, N);
    for i = 0:N-1
        indx = i + 1;
        state = zeros(NL, 1);
        binary = mynum2bin(i);
        l = length(binary);
        state(NL-l+1:NL) = binary;

        indx2 = checkdiff(state);
        for j = 1:length(indx2)
            F(indx2(j) + 1, indx) = F(indx2(j) + 1, indx) + kdiff;
            F(indx, indx) = F(indx, indx) - kdiff;
        end

        indx2 = checkads(state);
        for j = 1:length(indx2)
            S(indx2(j) + 1, indx) = S(indx2(j) + 1, indx) + kads;
            S(indx, indx) = S(indx, indx) - kads;
        end

        indx2 = checkdis(state);
        for j = 1:length(indx2)
            S(indx2(j) + 1, indx) = S(indx2(j) + 1, indx) + kdis;
            S(indx, indx) = S(indx, indx) - kdis;
        end

        indx2 = checkrx2(state);
        for j = 1:length(indx2)
            S(indx2(j) + 1, indx) = S(indx2(j) + 1, indx) + krx2;
            S(indx, indx) = S(indx, indx) - krx2;
        end

    end
end

function outnum = mybin2num(arr)
    N = length(arr);
    outnum = 0;
    for i = 1:N
        outnum = 2 * outnum;
        if arr(i)
            outnum = outnum + 1;
        end
    end
end

function outarr = mynum2bin(num)
    outarr = [];
    while num >= 1
        if mod(num, 2)
            outarr = [1, outarr];
            num = num - 1;
        else
            outarr = [0, outarr];
        end
        num = num / 2;
    end
end

function outnums = checkdiff(state)
    NL = length(state);
    outnums = [];
    for i = 1:NL
        if state(i) == 0
            continue
        end

        neigh = i - 1;
        if neigh == 0
            neigh = NL;
        end

        if state(neigh) == 0
            temp = state;
            temp(neigh) = 1;
            temp(i) = 0;

            indx = mybin2num(temp);
            outnums = [outnums, indx];
        end
        
        neigh = i + 1;
        if neigh > NL
            neigh = 1;
        end

        if state(neigh) == 0
            temp = state;
            temp(neigh) = 1;
            temp(i) = 0;

            indx = mybin2num(temp);
            outnums = [outnums, indx];
        end
    end
end

function outnums = checkads(state)
    NL = length(state);
    outnums = [];
    for i = 1:NL
        if state(i) == 1
            continue
        end
        temp = state;
        temp(i) = 1;

        indx = mybin2num(temp);
        outnums = [outnums, indx];
    end
end

function outnums = checkdis(state)
    NL = length(state);
    outnums = [];
    for i = 1:NL
        if state(i) == 0
            continue
        end
        temp = state;
        temp(i) = 0;

        indx = mybin2num(temp);
        outnums = [outnums, indx];
    end
end

function outnums = checkrx2(state)
    NL = length(state);
    outnums = [];
    for i = 1:NL
        if state(i) == 0
            continue
        end

        neigh = i + 1;
        if neigh > NL
            neigh = 1;
        end

        if state(neigh) == 0
            continue
        end

        temp = state;
        temp(i) = 0;
        temp(neigh) = 0;
        
        indx = mybin2num(temp);
        outnums = [outnums, indx];
    end
end


% [F, S] = CreateTmats(4, 1E+9, 8.5E+5, 4.0E+5, 1E+6);