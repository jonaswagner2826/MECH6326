C = 6;
D = 2;

X_set = [0:C]';
U_set = [0:C]';
D_set = [0:D]';

n = length(X_set);
m = length(U_set);
p = length(D_set);

p_d = [0.7;0.2;0.1]; 

T = 100;


n_runs = 10;
figure;hold on
for j = 1:n_runs
rng(j)
x = zeros(1,T+1);
x(1) = C;
for i = 1:T
    if x(i) <= 1
        u = C-x(i);
    else
        u = 0;
    end
    r = 1-rand;
    d = find(r < cumsum(p_d),1) - 1;
    x(i+1) = min(max(0,x(i)-d+u),C);
end

stairs(0:T,x)
end


%%
P = zeros(C+1,C+1);
P(1,end:-1:end-D) = p_d;
P(2,end:-1:end-D) = p_d;
for i = 3:C+1
    P(i,i-2+[D:-1:0]) = p_d;
end

figure;hold on
for j = 1:n_runs
rng(j)
x = zeros(1,T+1);
x(1) = C;
for i = 1:T
    r = rand;
    x(i+1) = find(r < cumsum(P(x(i)+1,:)),1) - 1;
end

stairs(0:T,x)
end


%% What is the probability that we reorder at time 25 - P(x_25 \in {0,1})
d0 = [zeros(1,C),1];
d_25 = d0*P^25;

P_reorder = d_25(1)+d_25(2)

%% Confirm this with Monte Carlo

T = 25;
n_runs = 1e4;

% figure;hold on
prob = zeros(1,n_runs);
for j = 1:n_runs
rng(j)
x = zeros(1,T+1);
x(1) = C;
for i = 1:T
    r = rand;
    x(i+1) = find(r < cumsum(P(x(i)+1,:)),1) - 1;
end
if j >= 2
    prob(j) = (prob(j-1)*(j-1) + (x(end) <= 1))/j;
end
% plot(0:T,x)
end
figure;hold on
plot(prob)
plot([1 n_runs],[P_reorder P_reorder])

%% Ergodic Chain
P = [0.5 0.5 0.0 0.0;...
     0.3 0.4 0.3 0.0;...
     0.0 0.3 0.4 0.3;...
     0.0 0.0 0.5 0.5];
figure; 
for i = 1:30
    imagesc(P^i)
    pause(0.2)
end
%%
d0 = [1 0 0 0];
d20 = d0*P^25
d0 = [0 0 0 1];
d20 = d0*P^25

%% Absorbing Chain
% P = [1.0 0.0 0.0 0.0 0.0 0.0 0.0;...
%      0.4 0.3 0.3 0.0 0.0 0.0 0.0;...
%      0.0 0.4 0.3 0.3 0.0 0.0 0.0;...
%      0.0 0.0 0.4 0.3 0.3 0.0 0.0;...
%      0.0 0.0 0.0 0.4 0.3 0.3 0.0;...
%      0.0 0.0 0.0 0.0 0.4 0.3 0.3;...
%      0.0 0.0 0.0 0.0 0.0 0.0 1.0];
P = [0.3 0.3 0.0 0.0 0.0 0.0 0.4;...
     0.4 0.3 0.3 0.0 0.0 0.0 0.0;...
     0.0 0.4 0.3 0.3 0.0 0.0 0.0;...
     0.0 0.0 0.4 0.3 0.3 0.0 0.0;...
     0.0 0.0 0.0 0.4 0.3 0.3 0.0;...
     0.0 0.0 0.0 0.0 0.0 1.0 0.0;...
     0.0 0.0 0.0 0.0 0.0 0.0 1.0];
figure; 
for i = 1:50
    imagesc(P^i)
    pause(0.01)
end 
 