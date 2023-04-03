%% Problem Definition

% Probabilities of a win, lose, draw for Aggressive play
A.w = 0.5;
A.l = 0.4;
A.d = 0.1;
% Probabilities of a win, lose, draw for Defensive play
D.w = 0.1;
D.l = 0.2;
D.d = 0.7;

%% Open-loop play

% All options
OL = ['A' 'A' 'D' 'D';...
      'A' 'D' 'A' 'D'];
% Theoretical probabilities
p_OL = [0.555, 0.515, 0.515, 0.415];

nOptions = size(OL,2);
nSamples = 1e5;
OLRecord = zeros(nSamples,nOptions);
runningAve = zeros(nSamples,1);
rng(1)
tStart = tic;
for i = 1:nOptions
    figure; hold on
    for j = 1:nSamples
        % Game 1
        [points1] = winLoseDraw(OL(1,i),A,D);
        % Game 2
        [points2] = winLoseDraw(OL(2,i),A,D);
        % Total points
        pointsTot = points1 + points2;
        if pointsTot >= 3
            OLRecord(j,i) = 1; % Win
        elseif pointsTot == 2
            OLRecord(j,i) = randi([0 1],1,1); % Penalty kicks
        else
            OLRecord(j,i) = 0; % Lose
        end
        % Record running averages
        runningAve(j) = sum(OLRecord(1:j,i))/j;
    end
    % Plot running averages
    plot(1:nSamples,runningAve,'o')
    plot([1 nSamples],[p_OL(i) p_OL(i)],'--k')
end
toc(tStart)

% Compute final sample averages
sample_ave = sum(OLRecord)/nSamples

%% Closed-loop play

% Theoretical probability
p_CL = 0.605;

CLRecord = zeros(nSamples,1);
runningAve = zeros(nSamples,1);
rng(1)

figure; hold on
tStart = tic;
for j = 1:nSamples
    % Game 1
    [points1] = winLoseDraw('A',A,D);
    % Game 2
    if points1 == 2
        [points2] = winLoseDraw('D',A,D);
    else
        [points2] = winLoseDraw('A',A,D);
    end
    
    % Total points
    pointsTot = points1 + points2;
    if pointsTot >= 3
        CLRecord(j,1) = 1; % Win
    elseif pointsTot == 2
        CLRecord(j,1) = randi([0 1],1,1); % Penalty kicks
    else
        CLRecord(j,1) = 0; % Lose
    end
    % Record running averages
    runningAve(j) = sum(CLRecord(1:j,1))/j;
end
toc(tStart)
% Plot running averages
plot(1:nSamples,runningAve,'o')
plot([1 nSamples],[p_CL p_CL],'--k')

% Compute final sample averages
sample_ave = sum(CLRecord)/nSamples



%% Oracle function for winning or losing

function [points] = winLoseDraw(playStyle,A,D)
    randVal = rand;
    if playStyle == 'A'
        if randVal <= A.w
            points = 2;
        elseif randVal <= A.w + A.l
            points = 0;
        else
            points = 1;
        end
    else
        if randVal <= D.w
            points = 2;
        elseif randVal <= D.w + D.l
            points = 0;
        else
            points = 1;
        end
    end
end