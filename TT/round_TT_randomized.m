function [A_star, r] = round_TT_randomized(A, epsilon, step, dir, tolR, rbase)
    if ~exist('step', 'var')
        step = 1;
    end
    if ~exist('dir', 'var')
        dir = 'left';
    end
    if ~exist('tolR', 'var')
        tolR = 2;
    end
    if ~strcmpi(dir, 'right') && ~strcmpi(dir, 'left') 
        error('Unknown direction specified. Choose either LEFT or RIGHT') 
    end

    d = A.order;
    if ~exist('rbase', 'var')
        rbase = [1, repmat(step, 1, d-1), 1];
    elseif numel(rbase) == 1
        rbase = [1, repmat(rbase, 1, d-1), 1];
    elseif numel(rbase) == d + 1
        rbase([1,end]) = 1;
    else
        error('Wrong format for the guessed TT-ranks. It must either be a scalar or a vector of length (A.order + 1)') 
    end
    r0 = A.rank;
    A_star = A;
    Ys = cell(1, d);

    rold = 0;
    r = rbase;
    
    if strcmpi(dir, 'right')
        diropp = 'left';
    else %if strcmpi(dir, 'left') 
        diropp = 'right';
    end
    while ~isequal(rold, r)
        [A_prime, Ys] = truncate_TT_randomized(A, r, 0, dir, Ys);
        A_star = round_nonortho(A_prime, epsilon * (d - 1)^-.5, diropp);
        posupdate = find(A_star.rank >= r - tolR);
        posupdate = posupdate(2:end-1);
        if isempty(posupdate)
            break;
        end
        rold = r;
        r(posupdate) = r(posupdate) + step;
        posexceed = find(r > r0); 
        if ~isempty(posexceed)
            r(posexceed) = r0(posexceed);
        end
    end
end
