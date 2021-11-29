function [A,w] = reedmullerweights(r,m)
% function [A,w] = reedmullerweights(r,m)
%REEDMULLERWEIGHTS Outputs the weight distribution of a binary reed-muller 
% code with parameters r=1 (first order) or r=2 (second order), for any m.
% 
% The calculation is based on:
% for r=1, the results is well-known. See any ECC course for proof.
% for r=2, see:
%   NEIL J. A. SLOANE, AND ELWYN R. BERLEKAMP, Weight Enumerator for Second-Order
%   Reed-Muller Codes, 1970.

switch r
    case 1
        w = [0,2.^(m-1),2^m];
        A = [1,2.^(m+1)-2,1];
        
    case 2
        jj = 1 : 1 : floor(m/2);
        
        w_minus = zeros(size(jj));
        A_minus = zeros(size(jj));
        for idx = 1 : 1 : length(jj)
            w_minus(idx) = 2.^(m-1) -2.^(m-1-jj(idx));
            ii = 1 : jj(idx);
            A_minus(idx) = 2.^(jj(idx).*(jj(idx)+1)) .* prod( (2.^(m-2*ii+2)-1) .* (2.^(m-2*ii+1)-1) ./ (4.^ii-1) );
        end
        
        w_plus = zeros(size(jj));
        A_plus = zeros(size(jj));
        for idx = 1 : 1 : length(jj)
            w_plus(idx) = 2.^(m-1) + 2.^(m-1-jj(idx));
            ii = 1 : jj(idx);
            A_plus(idx) = 2.^(jj(idx).*(jj(idx)+1)) .* prod( (2.^(m-2*ii+2)-1) .* (2.^(m-2*ii+1)-1) ./ (4.^ii-1) );
        end
        
        % add zeros, half-size, ones words (in a sorted manner)
        w = [0, w_minus,                2^(m-1),                 flip(w_plus), 2^m];
        A = [1, A_minus, 2*(2.^((m*(m+1))/2) - sum([1,A_minus])), flip(A_plus), 1];

    otherwise
        k = 0; for rr=0:1:r, k = k+sum(nchoosek(m,rr)); end
        filename = sprintf('%d_%d_%d.txt',2^m,k,2^(m-r));
        if exist(filename,'file')~=2, error('Unknown weight enumerator for given r.'); end
        
        WeightDistTable = readtable(filename);
        
        w = WeightDistTable{:,1}.';
        A = double(WeightDistTable{:,2}.');
end
