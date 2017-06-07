function [C_star, r] = round_hada_TT_randomized(A, B, epsilon, step, dir, tolR, rbase)
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
    if ~isequal(d, B.order)
        error('The order of the two matrices must coincide')
    end
    if ~isequal(A.size, B.size)
        error('The size of the two matrices must coincide')
    end
        
    if ~exist('rbase', 'var')
        rbase = [1, repmat(step, 1, d-1), 1];
    elseif numel(rbase) == 1
        rbase = [1, repmat(rbase, 1, d-1), 1];
    elseif numel(rbase) == d + 1
        rbase([1,end]) = 1;
    else
        error('Wrong format for rbase. It must either be a scalar or a vector of length (X.order + 1)') 
    end
    
    r0 = A.rank;
    
    C_star = A;
    Rs = cell(1, d);

    rold = 0;
    r = rbase;
    
    if strcmpi(dir, 'right')
        diropp = 'left';
    else %if strcmpi(dir, 'left') 
        diropp = 'right';
    end
    
    while ~isequal(rold, r)
        [C_prime, Rs] = truncate_hada_TT_randomized(A, B, r, 0, dir, Rs);
        C_star = round_nonortho(C_prime, epsilon * (d - 1)^-.5, diropp);
        posupdate = find(C_star.rank >= r - tolR);
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
