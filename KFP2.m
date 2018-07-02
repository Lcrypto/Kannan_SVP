%Copyright(c) 2013, USAtyuk Vasiliy 
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met :
%*Redistributions of source code must retain the above copyright
%notice, this list of conditions and the following disclaimer.
%* Redistributions in binary form must reproduce the above copyright
%notice, this list of conditions and the following disclaimer in the
%documentation and / or other materials provided with the distribution.
%* Neither the name of the <organization> nor the
%names of its contributors may be used to endorse or promote products
%derived from this software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED.IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
%DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function  [solution, comb, normsolution]  = KFP2( A,norma2 )    

    %[B,Q,R,T] = LLL(A,0.999999);
    %[A,Q,R,T] = LLL(B,0.999999);
    
    
    %B=A; %basis after reduction
   
    [Q, R, D] = mgso(A);
    %R*Q-B must be zero
    digits(340);
    m=size(A,1);
    c=zeros(1,m);
    l=zeros(1,m+1);
    y=zeros(1,m);
    x=zeros(1,m);
    b=zeros(1,m);
    for i=1:m
    b(i)=vpa(norm(Q(i,:),2)^2);
    end
    x(1)=1;
    dx=x;
    d2x(1:m)=-1;
    d2x(1)=1;
    solution=-1;
    normsolution=-1;
    AAA=true;
    i=1;
    combination=0;
    while(AAA)
       y(i)=vpa(abs(x(i)-c(i)));
   
       l(i)=vpa(l(i+1)+vpa(vpa(y(i))^2*b(i)));
       
       if (l(i)<=norma2^2)&& (i==1)
       
             solution=x*A'
             comb=x
             normsolution=norm(solution,2)
            norm2=normsolution;
            
           break;
       end
       if (l(i)<=norma2^2)&& (i>1)
       
              i=i-1;
              for j=i+1:m
              c(i)=vpa(c(i)+vpa(x(j)*R(j,i)));
              end
           
              c(i)=-c(i);
              x(i)=round(c(i));
              dx(i)=0;
              if c(i)<x(i)
                  d2x(i)=1;
              else
                  d2x(i)=-1;
              end
       elseif (l(i)>norma2^2) &&(i==m)
       
                AAA=false;
				display('Doesn found shortest vector');
                abort;
       else
             
             i=i+1;
             d2x(i)=-d2x(i);
             dx(i)=-dx(i)+d2x(i);
             x(i)=x(i)+dx(i);
             end
          end
       end
    



