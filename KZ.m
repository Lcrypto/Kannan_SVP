function B= KZ(B)

% function: KZ algorithm; 
% input: lattice basis B 
% output: KZ reduced basis B 

% author: Shanxiang Lyu, Imperial College London
% IEEE Tran. Signal Process., in press, https://arxiv.org/abs/1703.03303

%B=LLL(B);%LLL preprocessing, can be deleted
n=size(B,2);
for i=1:n

    block=i:n;
    [~,R]=qr(B);
    Ublock=firstU(R(block,block));%find the shortest vector and expand it to the basis
    
    U=eye(n);
    U(block,block)=Ublock;
    B=B*U;

   if i>1 %size reduction
       [~,R]=qr(B);
      for k=i-1:-1:1 
         q=round(R(k,i)/R(k,k));
              if q~=0   
                  B(:,i)=B(:,i)-q*B(:,k);
                  R(:,i)=R(:,i)-q*R(:,k);
              end    
      end

   end

end



