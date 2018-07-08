function [B,Q,R,T] = LLL_LR(A,delta)
% ==============================================================================
% Function:  LLL reduction of complex valued input matrix A with LLL 
%            parameter delta
%
% Interface: [B,Q,R,T] = LLL(A,delta)
%
% Input:     A     : input basis
%            delta : LLL reduction parameter
%
% Output:    B : LLL-reduced basis, B = A*T with QR decomposition B=Q*R
%            Q : matrix with orthogonal columns of length one
%            R : upper triangular matrix
%            T : unimodular transformation matrix
%
% Sources:   D. Wuebben, R. Boehnke, V. Kuehn, and K.-D. Kammeyer:
%            MMSE-Based Lattice-Reduction for Near-ML Detection of MIMO Systems
%            International ITG/IEEE Workshop on Smart Antennas (WSA 2004) 
%            Munich, Germany, March 2004
%
%            D Wuebben, D. Seethaler, J. Jaldéen, and G. Matz:
%            Lattice Reduction - A Survey with Applications in
%            Wireless Communications
%            IEEE Signal Processing Magazine, March 2011
%
% ==============================================================================
% Author:    Dirk Wuebben (University of Bremen, Bremen, Germany)
% EMail:     wuebben@ant.uni-bremen.de
% www:       www.ant.uni-bremen.de/staff/wuebben
% ==============================================================================

% Input data -------------------------------------------------------------------
if nargin <2                           % delta not specified
   fprintf('\nLLL parameter delta not specified!')
   fprintf('\nDefault value delta = 3/4 is used.\n');
   delta = 3/4;
elseif delta > 1 | delta < 0.25        % invalid value for delta
   fprintf('\nInvalid input value for LLL prameter delta!')
   fprintf('\nDefault value delta = 3/4 is used.\n');
   delta = 3/4;
end

% Initialization ---------------------------------------------------------------
B     = A;                             % reduced matrix B
[Q,R] = qr(B,0);                       % QR decomposition of matrix B
[n,m] = size(B);                       % matrix dimensions
T     = eye(m);                        % unimodular m x m matrix

% ==============================================================================
% LLL - reduction
% ==============================================================================
l  = 2;                                % 

while l <= m
  
   % Size reduction of column vector B(:,l) ------------------------------------
   for k=l-1:-1:1
      mu = round(R(k,l)/R(k,k));       % abs(R(k,l))>0.5*abs(R(k,k))
      if abs(mu)
         B(:,l)   = B(:,l)   - mu * B(:,k);
         R(1:k,l) = R(1:k,l) - mu * R(1:k,k);
         T(:,l)   = T(:,l)   - mu * T(:,k);
      end
   end

   % Lovasz condition ----------------------------------------------------------
   len = norm(R(l-1:l,l));
   if delta*abs(R(l-1,l-1))^2 > len^2
      
      % swapping of columns l-1 and l in B, T and R
      B(:,[l-1 l])   = B(:,[l l-1]);
      T(:,[l-1 l])   = T(:,[l l-1]);
      R(1:l,[l-1 l]) = R(1:l,[l l-1]);

      % reconstruction of upper triangular structure by Givens rotation 
      % mutliplication with matrix Theta achieves R(l,l-1) = 0
      c     = R(l-1,l-1) / len;        % len = ||R(l-1:l,l-1)|| after swapping
      s     = R(l,l-1)   / len;
      Theta = [c' s; -s c];
      
      R(l-1:l,l-1:end) = Theta * R(l-1:l,l-1:end);
      Q(:,l-1:l)       = Q(:,l-1:l) * Theta' ;
      l                = max(l-1,2);
   else
      l = l+1;
   end
end