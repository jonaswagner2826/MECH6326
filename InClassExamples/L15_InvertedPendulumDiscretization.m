%% System Parameters
g = 9.81;
L = 1;
dt = 0.05;

%% Discretize State and Input Spaces
xMin = [-pi; -3*pi];
xMax = [ pi;  3*pi];

n = 51; % Number of state grid points in each dimension (must be odd)

xPoints = [linspace(xMin(1),xMax(1),n);linspace(xMin(2),xMax(2),n)];

uMax = 5;

m = 5; % Number of input grid points (must be odd)

uPoints = linspace(-uMax,uMax,m);

% Stage cost
G = ones(n,n);
G(ceil(n/2),ceil(n/2)) = 0; % Cost-free terminal state.

%% Inifinite Horizon Dynamic Programming
rng(1)
J = zeros(n,n); % Optimal costs
pi_Star = zeros(n,n); % Optimal policy
J_Plus = rand(n,n); % Optimal costs at next iteration
figure;
t = 0; % value iteration counter
diff = 1;
while diff >= 1e-8 % Value iteration
    t = t + 1 % Count iterations
    for i = 1:ceil(n/2)
        for j = 1:n
            Ju = zeros(m,1);
            for k = 1:m
                xNext = pendulumUpdate([xPoints(1,i);xPoints(2,j)],uPoints(k),dt,g,L);
                [~,iNext] = min(abs(xNext(1)-xPoints(1,:)));
                [~,jNext] = min(abs(xNext(2)-xPoints(2,:)));
                Ju(k) = G(i,j) + 0*abs(uPoints(k)) + J(iNext,jNext);
            end
            [J_Plus(i,j),uIndx] = min(Ju);
            pi_Star(i,j) = uPoints(uIndx);
            % Exploit symmetry
            if i ~= ceil(n/2)
                J_Plus(n-i+1,n-j+1) = J_Plus(i,j);
                pi_Star(n-i+1,n-j+1) = -pi_Star(i,j);
            end
        end
    end
%     J
%     pi_Star
    imagesc(J)
    drawnow
    diff = norm(J - J_Plus)
    J = J_Plus; % Update costs
end

%% Plotting
figure;
imagesc(pi_Star)
figure;
imagesc(J_Plus)


%% Open-loop Simulation
Tsim = 200;
x = zeros(2,Tsim+1);
x(:,1) = [0.1;0];
% x(:,1) = [pi;0];
u = 0*ones(Tsim,1);
figure; hold on
h1 = plot(L*sin(x(1)),L*cos(x(1)),'ok');
h2 = plot([0;L*sin(x(1))],[0;L*cos(x(1))],'k');
xlim([-L L])
ylim([-L L])
axis square
for t = 1:Tsim
    delete(h1); delete(h2);
    [x(:,t+1)] = pendulumUpdate(x(:,t),u(t),dt,g,L);
    h1 = plot(L*sin(x(1,t+1)),L*cos(x(1,t+1)),'ok');
    h2 = plot([0;L*sin(x(1,t+1))],[0;L*cos(x(1,t+1))],'k');
    xlim([-L L])
    ylim([-L L])
    drawnow
%     pause(dt)
end
 figure; hold on
 plot(x(1,:),x(2,:))
 
 %% Closed-loop Simulation
Tsim = 200;
x = zeros(2,Tsim+1);
% x(:,1) = [0.1;0];
x(:,1) = [pi;0];
figure; hold on
h1 = plot(L*sin(x(1)),L*cos(x(1)),'ok');
h2 = plot([0;L*sin(x(1))],[0;L*cos(x(1))],'k');
xlim([-L L])
ylim([-L L])
axis square
for t = 1:Tsim
    delete(h1); delete(h2);
    [~,i] = min(abs(x(1,t)-xPoints(1,:)));
    [~,j] = min(abs(x(2,t)-xPoints(2,:)));
    pi_Star(i,j)
    [x(:,t+1)] = pendulumUpdate(x(:,t),pi_Star(i,j),dt,g,L);
    %%%%% Include to simulate on discretized system
    [~,iNext] = min(abs(x(1,t+1)-xPoints(1,:)));
    [~,jNext] = min(abs(x(2,t+1)-xPoints(2,:)));
    x(:,t+1) = [xPoints(1,iNext);xPoints(2,jNext)];
    %%%%%
    h1 = plot(L*sin(x(1,t+1)),L*cos(x(1,t+1)),'ok');
    h2 = plot([0;L*sin(x(1,t+1))],[0;L*cos(x(1,t+1))],'k');
    xlim([-L L])
    ylim([-L L])
    drawnow
%     pause(10*dt)
end
 figure; hold on
 plot(x(1,:),x(2,:))


%% System Dynamics
function  [xPlus] = pendulumUpdate(x,u,dt,g,L)
%     xPlus(1) = x(1) + dt*x(2);
%     xPlus(2) = x(2) + dt*(g/L*sin(x(1)) + u);
    if [x;u] == [0;0;0]
        xPlus = x;
    else
        [t,xAll] = ode23(@pendulum,[0 dt],x);
        xPlus = xAll(end,:)';
        xPlus(1) = pi*sawtooth(xPlus(1)-pi);
    end
    % Continuous-time pendulum model
    function dxdt = pendulum(t,x)
        dxdt(1,1) = x(2);
        dxdt(2,1) = (g/L*sin(x(1)) + u);
    end
end