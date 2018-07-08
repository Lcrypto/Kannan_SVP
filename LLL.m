function [B,T] = LLL(B)

% function: LLL algorithm; 
% input: lattice basis B
% output: reduced basis B, and the unimodular transform T

% author: Shanxiang Lyu, Imperial College London
% IEEE Tran. Signal Process., in press, https://arxiv.org/abs/1703.03303

[Q,R]=qr(B);
[m,n]=size(B);
T=eye(n);
delta=0.9728;%Lovasz constant

i=2;
V=[0 1;1 0];   

while i<=n

     for k=i-1:-1:1 
         q=round(R(k,i)/R(k,k));
              if q~=0  %do a unimodular for B, BS is not changed
                  R(:,i)=R(:,i)-q*R(:,k);
                  T(:,i)=T(:,i)-q*T(:,k);
              end    
     end
      
    if delta*norm(R(i-1,i-1))^2>abs(R(i,i))^2+abs(R(i-1,i))^2  %Lovasz fails

           %Givens matrix, don't put these before the def of R
           tempsum=sqrt(R(i-1,i)^2+R(i,i)^2);
           alpha=R(i-1,i)/tempsum;
           beta=R(i,i)/tempsum;
           
          R(:,i-1:i)=R(:,i-1:i)*V;%SWAP
          T(:,i-1:i)=T(:,i-1:i)*V;       
 
           G0=[alpha,beta;-beta,alpha];
           %Restore
           R(i-1:i,1:n)=G0*R(i-1:i,1:n);
           Q(1:m,i-1:i)=Q(1:m,i-1:i)*G0';

           i=max(i-1,2);

    else
        i=i+1;
    end   
end
B=Q*R;