function [k,d] = qaryreedmullerparams(m,q)
%QARYREEDMULLERPARAMS Calculates the code-dimension and minimum distance
%for given parameters m and q (alphabet size).

% calculate n
n = q^m;

vVec = 0 : (m*(q-1));
k = zeros(length(vVec),1);
d = zeros(length(vVec),1);
for iv = 1 : length(vVec)
    v = vVec(iv);
    % calculate k
    k(iv) = 0;
    for jj = 0 : 1 : m
        if v-jj*q<0, continue; end
        k(iv) = k(iv) + (-1)^jj .* nchoosek(m,jj) .* nchoosek(m+v-jj*q, v-jj*q);
    end
    % find distance
    R=0; Q=0;
    cond = false;
    while 1
        for R = (q-1) : -1 : 0 
            if (m*(q-1)-v) == (Q*(q-1)+R)
                cond=true; break;
            end
        end
        if cond, break;
        else Q = Q+1; end
    end
    d(iv) = (R+1)*q^Q;
end
