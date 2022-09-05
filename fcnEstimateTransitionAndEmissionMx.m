function [row,Z,q,pzz] = fcnEstimateTransitionAndEmissionMx(Y)

%% get "hidden" states
S = ones(size(Y));
for i = 1:size(Y,1)
    y = ceil(Y(i,:));
    s = ones(size(y));
      
    % state before ICANS
    yc = cumsum(y); 
    i1 = find(yc==0);
    s(i1) = 1; % state 1
        
    % state after ICANS
    % find first time from end above 0 (working back)
 
    if(max(yc)~=0);yc=yc/max(yc);i4 = find(yc==1);s(i4)=4;else;i4=[];end
       
    % rising and falling ICANS states
    a = max(i1);
    b = min(i4);
    c = max(find(y==max(y)));
    if a~=length(y)
        i2 = a:c;
        i3 = c:b;
        s(i2) = 2;
        s(i3) = 3;
    end
    S(i,:) = s;
end

%% map (s,y) -> row numbers (helpful for simulation)
row = (reshape(1:20,5,4))';
%% get states z[t] = (s[t],y[t-1])
Z = ones(size(Y,1),size(Y,2));
for i = 1:size(Y,1)
    for j = 2:size(Y,2)
        y0 = Y(i,j-1)+1;
        s1 = S(i,j);
        z = row(s1,y0);
        Z(i,j) = z;
    end
end


%% estimate emission and state transition probabilities
Nss = zeros(4,4);
Nysy = zeros(20,5);
for i = 1:size(Y,1)
    for j = 2:size(Y,2)
        s0 = S(i,j-1);   s1 = S(i,j);
        y0 = Y(i,j-1)+1; y1 = Y(i,j)+1;
        r = row(s1,y0);
        Nss(s0,s1) = Nss(s0,s1)+1;
        Nysy(r,y1) = Nysy(r,y1)+1;
    end
end

for i = 1:4; p(i,:) = Nss(i,:)/sum(Nss(i,:)); end

q = [];
for i = 1:20
    q(i,:) = (Nysy(i,:))/sum(Nysy(i,:));
end


%% state space is actually 20 dimensions; consists of z = (s[t],y[t-1])

for i = 1:4 % s[t]
    for j = 1:5 % y[t-1]
        for k = 1:4 % s[t-1]
            for l = 1:5 % y[t-2]
                z1 = row(i,j);
                z0 = row(k,l);
                pzz(z0,z1) = p(k,i)*q(z0,j);
            end
        end
    end
end

% being paranoid here...
for i = 1:20
    pzz(i,:) = pzz(i,:)/sum(pzz(i,:));
end
