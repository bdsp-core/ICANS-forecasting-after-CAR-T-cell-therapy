close all;clear all; clc; format compact;
rng(10)

load LR_prob_allp_Lasso_LOO

A = readtable('AllBWH_5days_updated01172022.xlsx');% read data from .xlsx file
B = [A.MRN A.meanscore];
nan_idx = find(~isnan(B(:,2)));% remove all missing data
data_final = B(nan_idx,:);
sid = unique(data_final(:,1)); %
oneset_time = NaN(1,length(sid));
for i = 1:length(sid)
    ind = find(data_final(:,1)==sid(i));
    L(i,1) = length(ind);
    y = zeros(1,100);
    y(1:length(ind)) = data_final(ind,2);
    % Calculate the onset NT:
    s = ones(1,100);
    yc = cumsum(ceil(y));
    i1 = find(yc==0);
    s(i1) = 0;
    if (~isempty(find(s,1)))
        oneset_time(i) = find(s,1);
    end
    Y(i,:) = ceil(y);
end
outlier = find(oneset_time>20);
% Remove the outlier data
Y(outlier,:) = [];
L(outlier)   = [];

%% augment Y -- to get smoother transition matrices
Nt = 5000;
Ya = Y;
for i = 1:Nt
    ind = randsample(size(Y,1),1); % get random sample from Y
    y = Y(ind,:);
    ind = find(y>0);
    s = randsample([-1 1],length(ind),1);clc
    n = (rand(size(ind))<0.2).*s;
    y(ind) = y(ind)+n;
    y = max(y,0);
    y = min(y,4);
    Ya = [Ya; y];
end
Y = Y+1;Ya = Ya+1;
%% get time-dependent transmission and emission matrices
[row,Z,q,pzz] = fcnEstimateTransitionAndEmissionMx((Ya-1));
[pzzt,pxzt] = fcnGetMatrices_BaseCase(pzz,q);
%% Forecasting
t=1:2:7;
for i=1:length(t)
    for h = 1:28-t(i)
        for idx =1:size(Y,1)
            y = Y(idx,:);
            pz0x0 = zeros(20,1); pz0x0(1) = 1;
            pxxht = fcnForwardPredict(t(i),h,y,pzzt,pxzt,pz0x0,yfit{idx});
            pxxht = normc((1:5)'.*pxxht);
            F{idx}(:,t(i)+h) = pxxht;
            P_Severe{i}(idx,h) = sum(pxxht(3:5));
            P_NT{i} (idx,h) = 1-pxxht(1);
        end
    end
    F_t{i}=F;
    F={};
end
%%% P_NT: NT probability from 1:2:7 days
%%% P_Severe: severe ICANS probability from 1:2:7 days
  


