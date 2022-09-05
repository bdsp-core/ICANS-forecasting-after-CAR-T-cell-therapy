close all;clear all; clc; format compact;
rng(10)

% load LR probability
load LR_prob_allp_Lasso_LOO

A = readtable('AllBWH_5days_updated01172022.xlsx');% read data from .xlsx file
B = [A.MRN A.meanscore];
nan_idx = find(~isnan(B(:,2)));% remove all missing data
data_final = B(nan_idx,:);
sid = unique(data_final(:,1)); % number of patient
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
outlier = find(oneset_time>20);% find outlier data
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

%% HMM+PLR
fig=figure;
T = 28; %time interval: [1,28]
for idx = 185
    for h=1:7
        y = Y(idx,:);
        F = []; F(:,1:h) = repmat([1 0 0 0 0]',1,h);
        
        for t = 1:T
            pz0x0 = zeros(20,1); pz0x0(1) = 1;
            pxxht = fcnForwardPredict(t,h,y,pzzt,pxzt,pz0x0,yfit{idx});
            F(:,t+h) = pxxht;
        end
        F = F(:,1:T);F = normc((1:5)'.*F);
        F(:,1:h) = NaN;
        % plot:
        subplot(8,1,h+1);
        ff=heatmap(F,'MissingDataColor',[.5, .5, .5],'CellLabelColor','none');
        colormap(flipud(hot));colorbar off
        ff.YDisplayData = flipud(ff.YDisplayData);
        ff.YDisplayLabels = {'4','3','2','1','0'};
        grid on
        title(['ICANS probability for ',num2str(h),'day(s) ahead forecasting']);
        ax = gca;
        ax.FontSize = 10;
        xlim([1,T])
        subplot(8,1,1);
        plot(1:T,y(1:T)-1,'linewidth',2);
        xticks(1:T)
        xticklabels(1:T);xtickangle(90)
        grid on
        xlim([1,T])
        ylim([0,4])
        box off
        title('Observed ICANS');
        ax = gca;
        ax.FontSize = 10;
        
    end
    han=axes(fig,'visible','off');
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'ICANS');
    xlabel(han,'Time[days]');
end
