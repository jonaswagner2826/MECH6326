C = 6;
D = 2;

X = [0:C]';
U = [0:C]';
W = [0:D]';

n = length(X);
m = length(U);
p = length(W);

p_W = [0.7;0.2;0.1]; 

P = zeros(n,n,m);
for i = 1:n
    for j = 1:m
        for k = 1:p
            x_next = min(max(0,X(i)+U(j)-W(k)),X(end));
            i_next = find(X==x_next,1);
            if U(j) < 2 - X(i)
                P(i,i_next,j) = 0;
            elseif U(j) > C - X(i)
                P(i,i_next,j) = 0;
            else
                P(i,i_next,j) = p_W(k);
            end
        end
    end
end


T = 50;

pi = zeros(n,T);
J = zeros(n,m,T+1);
Jstar = zeros(n,T+1);
figure; hold on
for t = [T-1:-1:0]
    for i = 1:n
        for j = 1:m
            J(i,j,t+1) = g(X(i),U(j),C)+P(i,:,j)*Jstar(:,t+2);
        end
        [~,indx] = min(J(i,:,t+1));
        pi(i,t+1) = U(indx);
        Jstar(i,t+1) = J(i,indx,t+1);
    end
    plot(X,Jstar(:,t+1))
end
pi

%%
function cost = g(x,u,C)
if u < 2 - x
    cost = inf;
    return
elseif u > C - x
    cost = inf;
    return
end
if u >= 1
    cost = 0.1*x+1;
else
    cost = 0.1*x;
end
end
