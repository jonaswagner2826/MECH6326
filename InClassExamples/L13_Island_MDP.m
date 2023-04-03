%% problem data
L = 30; W = 20; % grid size
I = 2; % island size

% plot state space
X = zeros(W, L);
X(W/2-I:W/2+I, L/2-I:L/2+I) = 1; % island in the middle
X(:,1) = -1; % waterfall on left boundary
figure;
imagesc(X);

n = W*L; % number of states
m = 4; % number of inputs
T = 30; % time horizon

% stage costs
G = zeros(W, L);
G(W/2-I:W/2+I, L/2-I:L/2+I) = 1;

%% Probability transition matrix
% m heading actions (right = 1, left = 2, up = 3, down = 4)
% m+1 transitions (right = 1, left = 2, up = 3, down = 4, stay = 5)
% 6 scenarios (1: island or waterfall, 2: interior states, 3: top, 4: bottom,
% 5: right, 6: top right, and 7: top left)
P = zeros(m,m+1,7); 
P(:,end,1) = 1; % if at the island or waterfall, for all heading actions, boat stays
% Interior states
P(:,:,2) = [0.6 0.2 0.1 0.1 0;
            0.0 0.8 0.1 0.1 0;
            0.0 0.3 0.6 0.1 0;
            0.0 0.3 0.1 0.6 0];
% Top states
P(:,:,3) = [0.7 0.2 0.0 0.1 0;
            0.0 0.9 0.0 0.1 0;
            0.0 0.0 0.0 0.0 0;
            0.0 0.3 0.0 0.7 0];
% Bottom states
P(:,:,4) = [0.7 0.2 0.1 0.0 0;
            0.0 0.9 0.1 0.0 0;
            0.0 0.3 0.7 0.0 0;
            0.0 0.0 0.0 0.0 0];
% Right states
P(:,:,5) = [0.0 0.0 0.0 0.0 0;
            0.0 0.8 0.1 0.1 0;
            0.0 0.3 0.6 0.1 0;
            0.0 0.3 0.1 0.6 0];
% Top Right states
P(:,:,6) = [0.0 0.0 0.0 0.0 0;
            0.0 0.9 0.0 0.1 0;
            0.0 0.0 0.0 0.0 0;
            0.0 0.3 0.0 0.7 0];
% Bottom Right states
P(:,:,7) = [0.0 0.0 0.0 0.0 0;
            0.0 0.9 0.1 0.0 0;
            0.0 0.3 0.7 0.0 0;
            0.0 0.0 0.0 0.0 0];

%% DP
J = zeros(W,L,T+1);
pistar = zeros(W,L,T);

J(:,:,T+1) = G; 

for t=T:-1:1 % backwards time recursion
    for i=1:W % river width
        for j=1:L % river length
            % Determine scenario
            if (G(i,j) == 1) || (j == 1) % Island or waterfall
                scenario = 1;
            elseif (1 < i) && (i < W) && (1 < j) && (j < L) % Interior
                scenario = 2;
            elseif (i == 1) && (1 < j) && (j < L) % Top
                scenario = 3;
            elseif (i == W) && (1 < j) && (j < L) % Bottom
                scenario = 4;
            elseif (j == L) && (1 < i) && (i < W) % Right
                scenario = 5;
            elseif (i == 1) && (j == L) % Top right
                scenario = 6;
            elseif (i == W) && (j == L) % Bottom right
                scenario = 7;
            end
            Ju = zeros(m,1);
            for k = 1:m % heading actions (right = 1, left = 2, up = 3, down = 4)
                Ju(k) = G(i,j) + P(k,:,scenario)*[J(i,min(L,j+1),t+1) J(i,max(1,j-1),t+1) J(max(1,i-1),j,t+1) J(min(W,i+1),j,t+1) J(i,j,t+1)]';
            end
            [J(i,j,t), pistar(i,j,t)] = max(Ju);
        end
    end
end
% plot value function and policy at initial time
figure;
surf(J(:,:,1));
shading interp
xlabel('length');
ylabel('width');

figure;
imagesc(pistar(:,:,1));
xlabel('length');
ylabel('width');

figure; hold on
for i = 1:T
    imagesc(pistar(:,:,i));
    xlabel('length');
    ylabel('width');
    drawnow
    pause(0.2)
end
%% Simulation
figure;
imagesc(pistar(:,:,1));
hold on 
xlabel('length');
ylabel('width');
colors = {'r','b','k','g','m'};
simT = 30;
nSims = 10;
for n = 1:nSims
rng(n);
x = zeros(2,simT+1);
% x(:,1) = [1;L-1];
x(:,1) = [randi([1 W]);randi([2 L])];
plot(x(2,1),x(1,1),'ro','MarkerSize',15);
for t = 1:simT
    if (G(x(1,t),x(2,t)) == 1) || (x(2,t) == 1) % Island or waterfall
        scenario = 1;
    elseif (1 < x(1,t)) && (x(1,t) < W) && (1 < x(2,t)) && (x(2,t) < L) % Interior
        scenario = 2;
    elseif (x(1,t) == 1) && (1 < x(2,t)) && (x(2,t) < L) % Top
        scenario = 3;
    elseif (x(1,t) == W) && (1 < x(2,t)) && (x(2,t) < L) % Bottom
        scenario = 4;
    elseif (x(2,t) == L) && (1 < x(1,t)) && (x(1,t) < W) % Right
        scenario = 5;
    elseif (x(1,t) == 1) && (x(2,t) == L) % Top right
        scenario = 6;
    elseif (x(1,t) == W) && (x(2,t) == L) % Bottom right
        scenario = 7;
    end
    
    u = pistar(x(1,t),x(2,t),1); % Optimal input
    p_u = P(u,:,scenario);
    r = rand;
    dir = find(r < cumsum(p_u),1);
    
    if dir == 1
        x(:,t+1) = x(:,t) + [0;1];
    elseif dir == 2
        x(:,t+1) = x(:,t) + [0;-1];
    elseif dir == 3
        x(:,t+1) = x(:,t) + [-1;0];
    elseif dir == 4
        x(:,t+1) = x(:,t) + [1;0];
    elseif dir == 5
        x(:,t+1) = x(:,t);
    end
    
    plot(x(2,t+1),x(1,t+1),'o','MarkerSize',15,'Color',colors{mod(t,length(colors))+1});
    drawnow
    pause(0.1)
end
delete(findobj(gca,'type','line'));
end