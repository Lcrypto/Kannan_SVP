function [ U_total ] = firstU( B )

% function: a subroutine in KZ and boosted KZ algorithms
% input: lattice basis B 
% output: unimodular matrix U that makes the first vector of B*U shortest

% author: Shanxiang Lyu, Imperial College London
% IEEE Tran. Signal Process., in press, https://arxiv.org/abs/1703.03303

[m,n]=size(B);
[B,T]=LLL(B); 

radius=norm(B(:,1));
z=SDSVP(B,radius);


if any(z)<1 %avoid problems
    z(1)=1;
end
m=min([m,n]);

U_total=eye(m);%declaration

flag1=0;%the pos of one
flag2=0;%the pos of the first non-zero value
for i=1:m
    if abs(z(i))==1 %check one
        flag1=i;
        break; %then no need to check flag2
    end
    
    if z(1)==0&&z(i)~=0
        flag2=i;
        break;
    end
end


if flag1>0 %put inside new pos and swap
        U_total(:,flag1)=z;
        Uswap=eye(m);
        Uswap(1,1)=0;Uswap(flag1,flag1)=0;
        Uswap(flag1,1)=1;Uswap(1,flag1)=1;
        U_total=U_total*Uswap;
else
        if z(1)==0 %swap flag2 and 1
                Uswap=eye(m);
                Uswap(1,1)=0;Uswap(flag2,flag2)=0;
                Uswap(flag2,1)=1;Uswap(1,flag2)=1;
                U_total=U_total*Uswap;%swap matrix
                ztemp=z(1);z(1)=z(flag2);z(flag2)=ztemp;%swap z
        end
        for ind=2:m %the normal case
        U0=eye(m);
        l2=z([1,ind]);
           if   l2(2)~=0 %only cancel non-zeros
                [r2,r1,d]=extended_gcd(l2(1),l2(2));%avoid 0,0 input
                Utemp=[l2(1)/d,-r1;l2(2)/d,r2];

                z([1,ind])=Utemp\z([1,ind]);
                U0([1,ind],[1,ind])=Utemp;%put U into larger U0
           end
        U_total=U_total*U0;
        end

end
U_total=T*U_total;
end

function [old_s, old_t,old_r]=extended_gcd(a, b)
%OUTPUT: coe1 coe2 and GCD
    s = 0;    old_s =1;
    t =1;    old_t = 0;
    r=b;    old_r = a;
    while r ~=0
        quotient =floor(old_r/r);
         
        tempr=r;
        r = old_r - quotient * r;
        old_r =tempr;
        
        temps=s;
        s=old_s - quotient * s;
        old_s=temps;
        
        tempt=t;
        t=old_t - quotient * t;
        old_t=tempt;
    end      
end