function [pzzt,pxzt] = fcnGetMatrices_BaseCase(pzz,q)

% create: p(z[t]|z[t-1]) = f(z1,z0,t), p(x[t]|z[t]) = f(x1,z1,t)
for t = 1:100
  for z0 = 1:20
      for z1 = 1:20
          for x1 = 1:5
              pzzt(z1,z0,t) = pzz(z1,z0); 
              pxzt(x1,z1,t) = q(z1,x1); 
          end
      end
  end    
end