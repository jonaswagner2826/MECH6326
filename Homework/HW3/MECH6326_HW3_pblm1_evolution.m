function X = MECH6326_HW3_pblm1_evolution(X,P_1,P_2)
    %MECH6326_HW3_pblm1_evolution 
    n = size(X,1);

    % Evolution of each state (dependent on lots of things)
    for i = 1:n
        for j = 1:n
            testvals = combvec(i-1:i+1,j-1:j+1);
            testvals = testvals(:,all([ ...
                all(testvals<=n); ...
                all(testvals>=1); ...
                any(testvals~=[i;j])...
                ]));
            P = P_2;
            for k = 1:size(testvals,2)
                if X(testvals(1,k),testvals(2,k),2) == 1; P = P_1;end
            end
            newState = randsample(4,1,true,reshape(X(i,j,:),[1,4])*P);
            X(i,j,:) = zeros(4,1);
            X(i,j,newState) = 1;
        end
    end
end