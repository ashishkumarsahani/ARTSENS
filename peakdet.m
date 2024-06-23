function [maxtab, mintab]=peakdet(v, delta)
%PEAKDET Detect peaks in a vector
%        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
%        maxima and minima ("peaks") in the vector V.
%        MAXTAB and MINTAB consists of two columns. Column 1
%        contains indices in V, and column 2 the found values.
%      
%        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
%        in MAXTAB and MINTAB are replaced with the corresponding
%        X-values.
%
%        A point is considered a maximum peak if it has the maximal
%        value, and was preceded (to the left) by a value lower by
%        DELTA.

% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
% This function is released to the public domain; Any use is allowed.
global maxtab_Mem;
global mintab_Mem;

S = size(v);
if(S(1)*S(2)==0)
    maxtab=[0,0];
    mintab =[0,0];
    return
end
max_Count =0;
min_Count =0;

mn = v(1,1); mx= v(1,1); mxpos=1;mnpos=1;
lookformax = 1;

for i=1:length(v)
  if lookformax
    if v(i) > mx, mx = v(i); mxpos =i;
    elseif v(i) < mx-delta
        if(mxpos ~=1)
          max_Count = max_Count+1;
          maxtab_Mem(max_Count,:) = [mxpos mx];
        end
      mn = v(i); mnpos = i;
      lookformax = 0;
    end  
  else
    if v(i) < mn, mn = v(i); mnpos =i;
    elseif v(i) > mn+delta
        if(mnpos~=S(1))
          min_Count = min_Count+1;
          mintab_Mem(min_Count,:) = [mnpos mn];
        end
      mx = v(i); mxpos = i;
      lookformax = 1;
    end
  end
end
maxtab = maxtab_Mem([1:max_Count],:);
mintab = mintab_Mem([1:min_Count],:);
end