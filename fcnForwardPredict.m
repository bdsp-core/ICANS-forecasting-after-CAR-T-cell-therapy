function pxxht = fcnForwardPredict(t,h,y,pzzt,pxzt,pz0x0,LR_p)%yfit_old,yfit_new)

% predict-update cycle
for u = 1:t
    % predict
    if(u<length(LR_p)+1)
        pzzt(1,6,u) = LR_p(u);
        for i = 1:20
            pzzt(i,:,u) = pzzt(i,:,u)/sum(pzzt(i,:,u));
        end
    end
    pz1x0 = pzzt(:,:,u)'*pz0x0; % [20,20]
    % update
    pz1x1 = pxzt(y(u),:,u)'.*pz1x0/sum(pxzt(y(u),:,u)'.*pz1x0); % [20,1] = [20,1].*[20,1]
    pz0x0 = pz1x1;
end

% forecasting - state
pzhxt = pz1x1;
for hh = 0:(h-1)
    pzhxt = pzzt(:,:,min(100,t+hh))'*pzhxt; % [20,1] = [20,20] x [20,1]
end

% forecasting - observation / ICANS - new
pxxht = pxzt(:,:,min(100,t+h))*pzhxt; % [5,1] = [5,20]x[20,1]
