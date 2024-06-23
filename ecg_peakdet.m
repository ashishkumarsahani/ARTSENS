function [Rpeak]=ecg_peakdet(v, delta)
% Detect R peak (QRS) in ECG waveform
%
%
global Rpeak_mem;

d = diff(v);
abs_d= abs(d)- mean(d);
norm_d = d /(max (abs_d));
subplot (2,1,2);
plot (norm_d);
%drawnow;
[p, n]=peakdet(norm_d, delta);

    peak_count =0; %p = floor p[]
    
for (i = 1:(length(p)))
    for (j = 0:1)
        if (length (n) < (i+j))
            break
        end
        if (((n(i+j) - p(i))<=20) && ((n(i+j) - p(i))>=0) )
            Rpeakind = floor((n(i+j) + p(i))/2);
            peak_count = peak_count + 1;
            Rpeak_mem(peak_count, :) = [Rpeakind v(Rpeakind)];
            break
        end
        
    end
end
if (peak_count>7 || peak_count ==0 )
    Rpeak_mem = [0, 0];
    Rpeak = Rpeak_mem;
    return
else
    Rpeak = Rpeak_mem([1:peak_count],:);
end

end


    

    
    %         j=1;
%         d_amp = norm_d(floor(p(i))+j);
%         j= j+1;
%         while(d_amp > 0 && (floor(p(i))+j)<length(v))
%             
%             d_amp = norm_d(floor(p(i))+j);
%             j= j+1;
%         end
%     peak_count = peak_count + 1;
        

