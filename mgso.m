function [Q, R,D] = mgso(A)
[d,n] = size(A);
m = min(d,n);
R = eye(m,n);
Q = zeros(d,m);
D = zeros(1,m);
for i = 1:m
    v = A(:,i);
    for j = 1:i-1
        R(j,i) = Q(:,j)'*v/D(j); %mu
        v = v-R(j,i)*Q(:,j);
    end
    Q(:,i) = v; %b ort
    D(i) = dot(Q(:,i),Q(:,i));
end
R(:,m+1:n) = bsxfun(@times,Q,1./D)'*A(:,m+1:n);
Q=Q';
R=R';
