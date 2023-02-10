%% World Cup Sim
clear
clc


%% Open Loop








%% Closed Loop








%% functions
function [seriesPoints] = PlayGame(strategy)
    randVal = rand;
    if randVal <= strategy.w
            seriesPoints = 2;
    elseif randVal <= strategy.w + strategy.l
        seriesPoints = 0;
    else
        seriesPoints = 1;
    end
end

function [winloss] = PlaySeries(strategies)
    % Game 1
    points = PlayGame(strategies[1]);
    % Game 2
    points = points + PlayGame(strategies[2]);
    % Win-Loss
    if points > 2
        winloss = 1;
    elseif points < 2
        winloss = 0;
    else
        winloss = randi([0,1]);
    end
end

    