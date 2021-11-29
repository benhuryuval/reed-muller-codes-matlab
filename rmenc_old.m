function [codeword,G,monomG] = rmenc_old(message,r,m,G,monomG)
% function [codeword,G,monomG] = rmenc_old(message,r,m,G,monomG)
% 
% RMENC encodes a binary message to an RM(r,m) codeword using a generator matrix. The generator matrix can
% be supplied as input, or calculated inside the function.
% 
% Source:   [1] Reed-Muller Codes, Sebastian Raaphorst, Carleton University, 2003.
% 
% Written by Yuval Ben-Hur, 26/02/2020
% 


message = message(:).';

% Code length and dimension
n = 2^m;
k = 0; for rr = 0:1:r, k = k+sum(nchoosek(m,rr)); end

if nargin<5
    % Create generator matrix
    G = ones(1,2^m);
    monomG = zeros(1,m+1);
    for nvars = 1 : 1 : r
        monompows = nchoosek(1:m,nvars);
        for powvecidx = 1 : 1 : nchoosek(m,nvars)
            monom = zeros(m+1,1);
            monom(monompows(powvecidx,:)) = 1;
            G = [G; phitransform(monom)];
            monomG = [monomG; monom.']; % for use in decoder
        end
    end
end

codeword = mod(message*G,2);


function [ phi_poly ] = phitransform(poly)
% function [ phi_poly ] = phitransform(poly)
%PHITRANSFORM implements the \phi transform for converting a polynom to a 
%   binary vector.
%   size of poly is [m+1,n] where each column represents a monom and n is 
%   the number of monoms in poly.
% 
%   poly = m_1 + m_2 + ... + m_n 
%       where each m_i is a monom in R_m, that
%       equals (x_0)^p_{1,i} + ... * (x_{m-1})^p_{m,i} * 0^p_{0,i}
% 
% for example, 
%   let p = x_0 + x_1*x_2 in R_4. the corresponding
%   matrix is [ [0;0;0;0;1], [0;1;1;0;0] ].
% 
%   let p = 1 + x_1 + x_0*x_2 + x_0*x_1*x_2 in R_3. the corresponding
%   matrix is [ [0;0;0;0], [0;1;0;0], [1;0;1;0], [1;1;1;0] ].
% 

m = size(poly,1)-1;
n = size(poly,2);
phi_poly = zeros(1,2^m);
for ii = 1 : 1 : n
    phi_poly = mod(phi_poly + phitransform_monom(poly(:,ii).'),2);
end


function [ phi_monom ] = phitransform_monom(monom)
% function [ phi_monom ] = phitransform(monom)
%PHITRANSFORM_MONOM implements the \phi transform for converting a monom to a binary vector.
% 
%   Input: binary vector of size [1,m+1] containing reduced monom powers.
%       monom = [r_0,r_1,...,r_{m-1},r_m] correspond to x_0^r_0 * x_1^r_1 * ...
%       * x_{m-1}^r_{m-1} * 0^{r_m}.
% 
%   Output: binary vector of size [1,m] containing elements of phi(monom).
% 
% for example, 
%   let m = 1 in R_3. the corresponsing vector is [0;0;0;0]
%   let m = 0 in R_4. the corresponsing vector is [0;0;0;0;1]
%   let m = x_0*x_2 in R_3. the corresponsing vector is [1;0;1;0]
% 
% The notations and calculations are based on:
%   Reed-Muller Codes, Sebastian Raaphorst, Carleton University, 2003.
% 

m = size(monom,2)-1;

if monom(m+1)==1 % the power of 0 is 1, 0^1=0
    phi_monom = zeros(1,2^m); % phi(0) = 0...0
else % the power of 0 is 0, 0^0=1
    phi_monom = ones(1,2^m); % phi(1) = 1...1
end

% Multiply phi transoforms of single variable monoms
for idx = 0 : 1 : m-1
    mat = repmat(repmat([1,0],[1,2^idx]),[2^(m-idx-1),1]); % phi(x_idx)
    if monom(idx+1)==1
        phi_monom = phi_monom.*mat(:).';
    end
end
