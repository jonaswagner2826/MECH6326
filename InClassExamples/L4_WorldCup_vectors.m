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
nSamples = 1e6;
OLRecord = zeros(nSamples,nOptions);
runningAve = zeros(nSamples,1);
rng(1)
tStart = tic;
for i = 1:nOptions
    figure; hold on
    % Game 1
    [points1] = winLoseDraw(OL(1,i),A,D,nSamples);
    % Game 2
    [points2] = winLoseDraw(OL(2,i),A,D,nSamples);
    % Total points
    pointsTot = points1 + points2;
    OLRecord(:,i) = 1*(pointsTot >= 3) + randi([0 1],nSamples,1).*(pointsTot == 2);
    % Record running averages
    runningAve = cumsum(OLRecord(:,i))./[1:nSamples]';
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

tStart = tic;
% Game 1
[points1] = winLoseDraw('A',A,D,nSamples);
% Game 2
[points2A] = winLoseDraw('A',A,D,nSamples);
[points2D] = winLoseDraw('D',A,D,nSamples);
% Total points
pointsTot = points1 + points2A.*(points1 ~= 2) + points2D.*(points1 == 2);
CLRecord(:,1) = 1*(pointsTot >= 3) + randi([0 1],nSamples,1).*(pointsTot == 2);
% Record running averages
runningAve = cumsum(CLRecord(:,1))./[1:nSamples]';
toc(tStart)

% Plot running averages
figure; hold on
plot(1:nSamples,runningAve,'o')
plot([1 nSamples],[p_CL p_CL],'--k')

% Compute final sample averages
sample_ave = sum(CLRecord)/nSamples



%% Oracle function for winning or losing

function [points] = winLoseDraw(playStyle,A,D,nSamples)
    randVals = rand(nSamples,1);
    if playStyle == 'A'
        points = 2*(randVals <= A.w) + 1*(randVals>=A.w + A.l);
    else
        points = 2*(randVals <= D.w) + 1*(randVals>=D.w + D.l);
    end
end