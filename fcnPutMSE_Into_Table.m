function T = fcnPutMSE_Into_Table(MSE,methodName)

ct=0; 
for d = 1:size(MSE,1)
    for j = 1:size(MSE,2)
        ct=ct+1; 
        mse(ct,1) = MSE(d,j); 
        days(ct,1) = d; 
        method{ct,1} = methodName;         
    end
end

T = table(method,days,mse); 
