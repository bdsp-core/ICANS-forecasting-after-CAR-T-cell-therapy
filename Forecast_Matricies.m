close all;clear; clc; format compact;
rng(10)

%% Load data:
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

% pzzt: [20,20,100]
% pxzt: [5,20,100]

%% initialize
fig=figure;
for idx =1:size(Y,1)

    for h=1:7
        
        y = Y(idx,:);
        F = []; F(:,1:h) = repmat([1 0 0 0 0]',1,h);
        ye(1:h) = y(1:h);
        ym(1:h) = y(1:h);
        for t = 1:28
            pz0x0 = zeros(20,1); pz0x0(1) = 1;
            pxxht = fcnForwardPredict(t,h,y,pzzt,pxzt,pz0x0,yfit{idx});
            F(:,t+h) = pxxht;
            ye(t+h) = (1:5)*pxxht;
            [~,ym(t+h)] =max(pxxht);
        end
        F = F(:,1:28);
        T = 28;%L(idx);
        ye = (ye(1:T));
        ym = ym(1:T);
        Baias(h,idx) = sum(y(h+1:T)-ye(h+1:end))/(L(idx)-h);
        MAD(h,idx) = sum(abs(y(h+1:T)-ye(h+1:end)))/(L(idx)-h);
        wMAPE(h,idx)= sum(abs(y(h+1:T)-ye(h+1:end)))/sum(y(h+1:T));

  
    
    end

end
%% Plot:
% % create table for boxplots
% T1 = fcnPutMSE_Into_Table(Baias,'BIAS');
% T2 = fcnPutMSE_Into_Table(MAD,'MAD');
% T3 = fcnPutMSE_Into_Table(wMAPE,'WAPE');
% 
% T = [T1; T2; T3 ];
% methodOrder = { 'BIAS', 'MAD', 'WAPE'};
% T.method = categorical(T.method,methodOrder);
% 
% figure(); clf;
% b =boxchart(T.days,T.mse,'GroupByColor',T.method);
% xlabel('Forecast days ahead');
% ylabel('Metrics');
% legend
% box off
% set(gcf,'color','w');
% grid on
% set(gca,'xtick',[1,2,3,4,5,6,7]);
% 
% hold on; SP_ = 1.5;
% line([SP_ SP_],[-1.5 1.5],'LineStyle','--','Color','k')
% hold on; SP_ = 3.5;
% line([SP_ SP_],[-1.5 1.5],'LineStyle','--','Color','k')
% hold on; SP_ = 2.5;
% line([SP_ SP_],[-1.5 1.5],'LineStyle','--','Color','k')
% hold on; SP_ = 4.5;
% line([SP_ SP_],[-1.5 1.5],'LineStyle','--','Color','k')
% hold on; SP_ = 5.5;
% line([SP_ SP_],[-1.5 1.5],'LineStyle','--','Color','k')
% hold on; SP_ = 6.5;
% line([SP_ SP_],[-1.5 1.5],'LineStyle','--','Color','k')
% % ylim([0,3.5])
