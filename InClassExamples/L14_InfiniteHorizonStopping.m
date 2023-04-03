%% Infinite Horizon Stopping Problem
k = 20; % Grid size
m = 2;  % Cardinality of input set

% Stage cost
G = zeros(k,k,m);
G(:,:,1) = 1; % Cost for waiting
G(5,5,2) = -120; % Costs for stopping
G(17,10,2) = -70;
G(10,15,2) = -150; 

%%
J = zeros(k,k); % Optimal costs
pi_Star = zeros(k,k); % Optimal policy
J_Plus = ones(k,k); % Optimal costs at next iteration

t = 0; % value iteration counter

%[down, up, right, left]

while norm(J - J_Plus) >= 1e-8 % Value iteration
    t = t + 1 % Count iterations
    J = J_Plus; % Update costs
    for i = 1:k
        for j = 1:k
            Ju = zeros(m,1);
            if i == 1
                if j == 1 % Top left corner
                    Ju(1) = G(i,j,1) + [1/2 1/2]*[J(i+1,j) J(i,j+1)]';
                elseif j == k % Top right corner
                    Ju(1) = G(i,j,1) + [1/2 1/2]*[J(i+1,j) J(i,j-1)]';
                else % Top boundary
                    Ju(1) = G(i,j,1) + [1/3 1/3 1/3]*[J(i+1,j) J(i,j+1) J(i,j-1)]';
                end
            elseif i == k 
                if j == 1 % Bottom left corner
                    Ju(1) = G(i,j,1) + [1/2 1/2]*[J(i-1,j) J(i,j+1)]';
                elseif j == k % Bottom right corner
                    Ju(1) = G(i,j,1) + [1/2 1/2]*[J(i-1,j) J(i,j-1)]';
                else % Bottom boundary
                    Ju(1) = G(i,j,1) + [1/3 1/3 1/3]*[J(i-1,j) J(i,j+1) J(i,j-1)]';
                end
            elseif j == 1 % Left boundary
                Ju(1) = G(i,j,1) + [1/3 1/3 1/3]*[J(i+1,j) J(i-1,j) J(i,j+1)]';
            elseif j == k % Right boundary
                Ju(1) = G(i,j,1) + [1/3 1/3 1/3]*[J(i+1,j) J(i-1,j) J(i,j-1)]';
            else
                Ju(1) = G(i,j,1) + [1/4 1/4 1/4 1/4]*[J(i+1,j) J(i-1,j) J(i,j+1) J(i,j-1)]';
            end
            Ju(2) = G(i,j,2);
            [J_Plus(i,j), pi_Star(i,j)] = min(Ju);
        end
    end
end

%% Plotting
figure;
imagesc(pi_Star)
figure;
imagesc(J_Plus)

%% Simulation
figure;
imagesc(pi_Star);
hold on 
colors = {'r','b','k','g','m'};
simT = 1000;
nSims = 10000;
Costs = zeros(nSims,1);
for n = 1:nSims
    Done = 0;
    rng(n);
    x = zeros(2,simT+1);
    x(:,1) = [10;12];
%     plot(x(2,1),x(1,1),'ro','MarkerSize',15);
    for t = 1:simT
        
        if Done == 1
            x(:,t+1) = x(:,t);
        else
            u = pi_Star(x(1,t),x(2,t)); % Optimal input
            if u == 1 % Wait
                r = rand;
                if x(1,t) == 1
                    if x(2,t) == 1 % Top left corner
                        dir = find(r < [1/2 1/2 1 1],1);
                    elseif x(2,t) == k % Top right corner
                        dir = find(r < [1/2 1/2 1/2 1],1);
                    else % Top boundary
                        dir = find(r < [1/3 1/3 2/3 1],1);
                    end
                elseif x(1,t) == k
                    if x(2,t) == 1 % Bottom left corner
                        dir = find(r < [0 1/2 1 1],1);
                    elseif x(2,t) == k % Bottom right corner
                        dir = find(r < [0 1/2 1/2 1],1);
                    else % Bottom boundary
                        dir = find(r < [0 1/3 2/3 1],1);
                    end
                elseif x(2,t) == 1 % Left boundary
                    dir = find(r < [1/3 2/3 1 1],1);
                elseif x(2,t) == k % Right boundary
                    dir = find(r < [1/3 2/3 2/3 1],1);
                else
                    dir = find(r < [0.25 0.5 0.75 1],1);
                end
%                 dir
                if dir == 1
                    x(:,t+1) = x(:,t) + [1;0];
                elseif dir == 2
                    x(:,t+1) = x(:,t) + [-1;0];
                elseif dir == 3
                    x(:,t+1) = x(:,t) + [0;1];
                elseif dir == 4
                    x(:,t+1) = x(:,t) + [0;-1];
                end
                Costs(n) = Costs(n) + 1;
            else % Stop
                x(:,t+1) = x(:,t);
                Costs(n) = Costs(n) + G(x(1,t),x(2,t),2);
                Done = 1;
                break
            end
        end
%         plot(x(2,t+1),x(1,t+1),'o','MarkerSize',15,'Color',colors{mod(t,length(colors))+1});
%         drawnow
%         pause(0.1)
    end
    delete(findobj(gca,'type','line'));
end

figure;hold on
plot(cumsum(Costs)./[1:nSims]')
plot([1 nSims],J_Plus(x(1,1),x(2,1))*[1 1])

