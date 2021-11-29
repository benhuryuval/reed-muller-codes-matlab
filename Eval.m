function [ c ] = Eval(f,z)
% function [ c ] = Eval(f,z)
%EVAL outputs the evaluation of a polynomial f in vector z. In case the 
% function is called only with the polynomial, it outputs all evaluations
% of the polynomial.
% 
% Example: In case f = m_1 + m_2 + ... + m_n 
%       where each m_i = (x_0)^p_{1,i} + ... * (x_{m-1})^p_{m,i} is a monom
% 
% for example, 
%   let f = x_3 + x_1*x_2 and m=4. the corresponding
%   matrix is [ [0;0;0;1], [0;1;1;0] ].
% 
%   let f = 1 + x_0*x_1*x_2 and m=3. the corresponding
%   matrix is [ [0;0;0], [1;1;1] ].
% 

m = size(f,1); % # of variables
n = size(f,2); % # of monoms

if nargin<2 % {f(z) : z \in F_2 ^m}
    Z = de2bi(0 : 1 : 2^m-1).';
else % f(z)
    Z = z;
end
c = zeros(1,size(Z,2));
for zz = 1 : 1 : size(Z,2)
    z = Z(:,zz);
    c(zz) = sum(prod(repmat(z,[1,n]) .^ f,1),2);
end
c = mod(c,2);
