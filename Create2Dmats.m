function [F, S] = Create2Dmats(L, kdiff, kads, kdis, krx2)
    %CREATEFMAT Creates fast reaction matrix (diffusion reactions)
    %   Detailed explanation goes here
    
    NL = L^2;
    N = 2^NL;
    F = zeros(N, N);
    S = zeros(N, N);
    for i = 0:N-1
        indx = i + 1;
        state = mynum2binmat(i, L);

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

function outnum = mybin2num(mat)
    N = numel(mat);
    arr = reshape(mat, 1, []);
    outnum = 0;
    for i = 1:N
        outnum = 2 * outnum;
        if arr(i)
            outnum = outnum + 1;
        end
    end
end

function outarr = mynum2binarr(num)
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

function outmat = mynum2binmat(num, L)
    NL = L^2;
    binary = mynum2binarr(num);
    l = length(binary);
    arr = zeros(NL, 1);
    arr(NL-l+1:NL) = binary;
    outmat = reshape(arr, [L, L]);
end

function outnums = checkdiff(state)
    [a, b] = size(state);
    outnums = [];
    for i = 1:a
        for j = 1:b
            if state(i, j) == 0
                continue
            end

            neigh = i - 1;
            if neigh == 0
                neigh = a;
            end
    
            if state(neigh, j) == 0
                temp = state;
                temp(neigh, j) = 1;
                temp(i, j) = 0;
    
                indx = mybin2num(temp);
                outnums = [outnums, indx];
            end
            
            neigh = i + 1;
            if neigh > a
                neigh = 1;
            end
    
            if state(neigh, j) == 0
                temp = state;
                temp(neigh, j) = 1;
                temp(i, j) = 0;
    
                indx = mybin2num(temp);
                outnums = [outnums, indx];
            end

            neigh = j - 1;
            if neigh == 0
                neigh = b;
            end
    
            if state(i, neigh) == 0
                temp = state;
                temp(i, neigh) = 1;
                temp(i, j) = 0;
    
                indx = mybin2num(temp);
                outnums = [outnums, indx];
            end
            
            neigh = j + 1;
            if neigh > b
                neigh = 1;
            end
    
            if state(i, neigh) == 0
                temp = state;
                temp(i, neigh) = 1;
                temp(i, j) = 0;
    
                indx = mybin2num(temp);
                outnums = [outnums, indx];
            end
        end
    end
end

function outnums = checkads(state)
    flat = reshape(state, 1, []);
    NL = length(flat);
    outnums = [];
    for i = 1:NL
        if flat(i) == 1
            continue
        end
        temp = flat;
        temp(i) = 1;

        indx = mybin2num(temp);
        outnums = [outnums, indx];
    end
end

function outnums = checkdis(state)
    flat = reshape(state, 1, []);
    NL = length(flat);
    outnums = [];
    for i = 1:NL
        if flat(i) == 0
            continue
        end
        temp = flat;
        temp(i) = 0;

        indx = mybin2num(temp);
        outnums = [outnums, indx];
    end
end

function outnums = checkrx2(state)
    [a, b] = size(state);
    outnums = [];
    for i = 1:a
        for j = 1:b
            if state(i, j) == 0
                continue
            end
    
            neigh = i + 1;
            if neigh > a
                neigh = 1;
            end
    
            if state(neigh, j) == 1
                temp = state;
                temp(i, j) = 0;
                temp(neigh, j) = 0;
                
                indx = mybin2num(temp);
                outnums = [outnums, indx];
            end
            
            neigh = j + 1;
            if neigh > b
                neigh = 1;
            end
    
            if state(i, neigh) == 1
                temp = state;
                temp(i, j) = 0;
                temp(i, neigh) = 0;
                
                indx = mybin2num(temp);
                outnums = [outnums, indx];
            end
        end
    end
end


% [F, S] = Anything(3, 1, 8.5, 4, 10);