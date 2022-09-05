%%% Lasso Logistic Regression Method
clear
warning off
clc

%% get data in X and Y matrices
[X, Xt, Yt, SID, SID1, variableNames] = fcnGetData;
u = unique(SID1);
%% Fit model
% Leave one out:
for i=1:length(u)
    ind = find(SID1~=u(i));
    [B,FitInfo] = lassoglm(Xt(ind,:),Yt(ind,:),'binomial','Link','logit','CV',10);
    idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
    B0 = FitInfo.Intercept(idxLambdaMinDeviance);
    b(:,i) = [B0; B(:,idxLambdaMinDeviance)];
    
end
[B{i},~] = glmfit(Xt(ind,:),Yt(ind,:),'binomial','Link','logit');

%% Evaluate model for t = 1:100, for each patient
figure(1); clf;
u = unique(SID);
for i = 1:length(u)
    ind = find(SID==u(i));
    x = X(ind,:);
    p = glmval(b(:,i),x,'logit');
    plot(p,'k'); hold on
    yfit{i} = p;
end
save('LR_prob_allp_Lasso_LOO','yfit')