function [X, Xt, Yt, SID, SID_1, variableNames] = fcnGetData

%% OUTPUTS: 
% X, Y
% X : 

A = readtable('AllBWH_5days_updated01172022.xlsx');% read data from .xlsx file
B = [A.MRN A.meanscore A.Age A.CRP A.Temp A.Ferritin A.Procalcitonin A.IL_6 A.WBC];

nan_idx = find(~isnan(B(:,2)));% remove all missing data
data_final = B(nan_idx,:);
data_final(1546:end,4) = data_final(1546:end,4)/10;
data_final(1546:end,5) = (data_final(1546:end,5)-32)/1.8;


%% try very simple forecasting method - forecast 2 days before

% fill up data structure
sid = unique(data_final(:,1)); % number of patients
onset_time = NaN(1,length(sid));
CRP  = NaN(length(sid),100);
Temp = NaN(length(sid),100);
IL_6 = NaN(length(sid),100);
WBC  = NaN(length(sid),100);
Ferr = NaN(length(sid),100);
Proc = NaN(length(sid),100);

for i = 1:length(sid)
    ind = find(data_final(:,1)==sid(i));
    L(i,1) = length(ind);
    y = zeros(1,100);
    y(1:length(ind)) = data_final(ind,2);
    CRP(i,1:L(i)) = data_final(ind,4);
    Temp(i,1:L(i))= data_final(ind,5);
    Ferr(i,1:L(i))= data_final(ind,6);
    Proc(i,1:L(i))= data_final(ind,7);
    IL_6(i,1:L(i))= data_final(ind,8);
    WBC(i,1:L(i))= data_final(ind,9);
    Age(i,1) = data_final(ind(1),3);
    
    % Calculate the onset NT:
    s = ones(1,100);
    yc = cumsum(ceil(y));
    i1 = find(yc==0);
    s(i1) = 0;
    if (~isempty(find(s,1)))
        onset_time(i) = find(s,1);
    end
end
outlier = find(onset_time>20);
load data_new

output = zeros(size(Y));
Ya = Ya(1:size(Y,1),:);
Z  = Z(1:size(Y,1),:);
Age(outlier)   = [];
L(outlier)     = [];
CRP(outlier,:) = [];
Temp(outlier,:)= [];
IL_6(outlier,:)= [];
WBC(outlier,:) = [];
Ferr(outlier,:)= [];
Proc(outlier,:)= [];
onset_time(outlier) = [];
ind = find(isnan(onset_time));

% t=1:100;output(:,t>20) = 0;
load time_LR  % s(t) = a function of time, peaks at day 9

%% Build feature matrix
% [Age, MaxT, CRP(:,k), Ferr(:,k), Proc(:,k) IL_6(:,k), WBC(:,k),t_LR]
X = []; Xt = []; Yt = []; SID = [];SID_1 = [];
variableNames = {'maxT' 'Age' 'CRP' 'Ferr' 'Proc' 'IL6' 'WBC' 'Curve'}; 
for j = 1:size(Y,1); % loop over patients
    
    % get data for training model
    Nt = min(onset_time(j), 100);
    x = [];
    yt = [];
    sid = [];
    for t=1:Nt % loop over time
        % get max of each feature for patient j up through time t
        x(t,1) = max(Temp(j,1:t));
        x(t,2) = Age(j); 
        x(t,3) = max(CRP(j,1:t)); 
        x(t,4) = max(Ferr(j,1:t)); 
        x(t,5) = max(Proc(j,1:t)); 
        x(t,6) = max(IL_6(j,1:t)); 
        x(t,7) = max(WBC(j,1:t)); 
        x(t,8) = s(t); 
        x(t,9) = t; 
        x(t,10) = t.^2;
        sid(t,1) = j;

        yt(t,1) = 0; 
    end
    if ~isnan(onset_time(j)); yt(end) = 1; end
    Xt = [Xt; x]; 
    Yt = [Yt; yt];
    SID_1 = [SID_1; sid];
    % get data for evaluating model
    x = [];
    sid = [];
    for t=1:L(j) % loop over time
        % get max of each feature for patient j up through time t
        x(t,1) = max(Temp(j,1:t));
        x(t,2) = Age(j); 
        x(t,3) = max(CRP(j,1:t)); 
        x(t,4) = max(Ferr(j,1:t)); 
        x(t,5) = max(Proc(j,1:t)); 
        x(t,6) = max(IL_6(j,1:t)); 
        x(t,7) = max(WBC(j,1:t)); 
        x(t,8) = s(t); 
        x(t,9) = t;
        x(t,10) = t.^2; 
        sid(t,1) = j; 
    end
    X = [X; x]; 
    SID = [SID; sid]; 
end

%% normalize columns
for i = 1:size(Xt,2); 
    x = Xt(:,i); 
    ind = find(isnan(x)); 
    x(ind) = nanmedian(x); 
    med(i) = nanmedian(x); 
    m(i) = mean(x); 
    sig(i) = std(x); 
    x = (x-m(i))/sig(i);
    Xt(:,i) = x; 
end

for i = 1:size(X,2); 
    x = X(:,i); 
    ind = find(isnan(x)); 
    x(ind) = med(i); 
    x = (x-m(i))/sig(i);
    X(:,i) = x; 
end