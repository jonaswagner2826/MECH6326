function X = MECH6326_HW3_pblm1_evolution(X,P_1,P_2)
    %MECH6326_HW3_pblm1_evolution 
    n = size(X,1);

%     X_update = zeros(size(X));
    % Evolution of each state (dependent on lots of things)
    for i = 1:n
        for j = 1:n
            testvals = combvec(i-1:i+1,j-1:j+1);
            testvals = testvals(:,all([ ...
                all(testvals<=n); ...
                all(testvals>=1); ...
                any(testvals~=[i;j])...
                ]));
            if any(X(testvals(1,:),testvals(2,:),2)==1)
                P=P_2;
                p = 'p_2';
            else
                P=P_1;
                p = 'p_1';
            end
            currentX = reshape(X(i,j,:),[1,4]);
            X(i,j,:) = zeros(4,1);
            X(i,j,randsample(4,1,true,currentX*P)) = 1;
%             if temp < 
% %             if temp < P(X)
%             X(i,j,:) = reshape(X(i,j,:),[1,4])*P;
        end
    end
end