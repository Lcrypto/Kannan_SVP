function [B, C, Precoding_mat, Precoding_mat_red, num_it] = Brun(H)
% ==============================================================================
% This function implements Brun's algorithm for lattice-reduction 
% assisted precoding (see [1,2]) 
% ------------------------------------------------------------------------------
% 
% INPUT: 
% 
% H                 ... complex-valued (channel) matrix, K x M with K <= M  
% ------------------------------------------------------------------------------
%
% OUTPUT:  
%
% B                 ... unimodular transformation matrix reducing pinv(H) 
% C                 ... inverse of B
% Precoding_mat     ... precoding matrix H'*inv(H*H') 
% Precoding_mat_red ... corresponding reduced precoding matrix Precoding_mat*B  
% num_it            ... number of iterations (number of basis updates)  
% ------------------------------------------------------------------------------
% 
% References: 
%  
% [1] D. Seethaler and G. Matz, "Efficient Vector Perturbation in Multi-Antenna 
% Multi-User Systems Based on Approximate Integer Relations, " EUSIPCO 2006
%
% [2] A. Burg, D. Seethaler, and G. Matz, "VLSI Implementation of a Lattice-
% Reduction Algorithm for Multi-Antenna Broadcast Precoding," IEEE ISCAS 2007  
% ------------------------------------------------------------------------------
% 
% Written by Dominik Seethaler, June 2009, as supplementary Matlab code 
% for paper: 
% 
% D. Wuebben, D. Seethaler, J. Jalden, and G. Matz, 
% "Lattice Reduction - A Survey with Applications in Wireless Communications" 
% IEEE Signal Processing Magazine, March 2011
% ------------------------------------------------------------------------------

% get size ---------------------------------------------------------------------
[K, M] = size(H); 

% check for non-empty input and full row rank condition ------------------------
if isempty(H) 
   error(['Input basis is empty']); 
elseif rank(H) < K
   error(['Input basis has not full row rank']); 
end 

% pre-computations -------------------------------------------------------------
Gram_mat        = H*H';             % Gram matrix 
Gram_mat_inv    = inv(Gram_mat);    % inverse Gram matrix 
Precoding_mat   = H'*Gram_mat_inv;  % Precoding matrix 
metric_vec      = sqrt(sum(abs(Precoding_mat).^2)); % column norms of precoding matrix 

% take an arbitrary column of the inverse Gram matrix as an approximation 
% to the left singular vector of H associated to the smallest singular 
% value (Brun's algorithm tries to find a good approximate integer relation 
% to this vector) 
% ------------------------------------------------------------------------------
mu  = Gram_mat_inv(:,1); 

% initialize -------------------------------------------------------------------
do_Brun             = 1;        % Brun termination flag 
num_it              = 0;        % number of iterations 
B                   = eye(K);   % unimodular transformation matrix 
C                   = eye(K);   % inverse transformation matrix 
Precoding_mat_red   = Precoding_mat; 

% do Brun's algorithm ----------------------------------------------------------
while do_Brun

   % increase iterations counter by one ---------------------------------------
   num_it = num_it + 1;
	    
   % Calculation of Brun's update value ---------------------------------------
   [zw s] = max(mu);		
   zw = mu; zw(s) = 0; [zw t] = max(zw); 
   r = round(mu(s)/mu(t));

   % basis update -------------------------------------------------------------
   updated_precoding_vec = Precoding_mat_red(:,s)-r'*Precoding_mat_red(:,t); 
    
   % if metric can be improved ------------------------------------------------
   updated_metric = norm(updated_precoding_vec);
   if updated_metric < metric_vec(s)          
    
      mu(s) = mu(s)-r*mu(t);
    	
      % update unimodular transformation matrix and its inverse --------------
      B(:,s) = B(:,s) - r'*B(:,t);    
      C(t,:) = C(t,:) + r'*C(s,:);   
 
      % copy updated precoding vector and corresponding metric ---------------
      Precoding_mat_red(:,s)  = updated_precoding_vec;
      metric_vec(s)           = updated_metric;
	
   % else terminate -----------------------------------------------------------
   else
      do_Brun = 0;    
   end
    	
end % - end do Brun's algorithm ------------------------------------------------