close all; clear; clc; format compact; rng(10)
A = readtable('AllBWH_5days_deID.xlsx');
B = [A.SID A.meanscore];
nan_idx = find(~isnan(B(:,2))); data_final = B(nan_idx,:);
sid = unique(data_final(:,1)); oneset_time = NaN(1,length(sid)); Y=[];
for i = 1:length(sid)
    ind = find(data_final(:,1)==sid(i)); L(i,1)=length(ind);
    y = zeros(1,100); y(1:length(ind)) = data_final(ind,2);
    s = ones(1,100); yc = cumsum(ceil(y)); s(find(yc==0))=0;
    if ~isempty(find(s,1)); oneset_time(i)=find(s,1); end
    Y(i,:) = ceil(y);
end
outlier = find(oneset_time>20); Y(outlier,:)=[]; L(outlier)=[];
Nt = 5000; Ya = Y;
for i = 1:Nt
    ind = randsample(size(Y,1),1); y = Y(ind,:); ind = find(y>0);
    s = randsample([-1 1],length(ind),1); n = (rand(size(ind))<0.2).*s;
    y(ind) = y(ind)+n; y = max(y,0); y = min(y,4); Ya = [Ya; y];
end
[row,Z,q,pzz] = fcnEstimateTransitionAndEmissionMx(Ya);
Y=uint8(Y); Ya=uint8(Ya); Z=uint8(Z);
save('data_new.mat','Y','Ya','Z');
fprintf('data_new rebuilt: Y %dx%d, Ya %dx%d, Z %dx%d\n', size(Y,1),size(Y,2),size(Ya,1),size(Ya,2),size(Z,1),size(Z,2));
