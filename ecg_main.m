%close all;
%clear all;


Data_Ch = zeros( 100, 1);
tstamp_Ch = Data_Ch;

tstamp_Rpeak = zeros( 2, 1);
data_size = 14*100;
ECG_Link(1);


while (1)   %for (i = 1 : 1000)
    Rpeak_detected = 0;
    [a,b, ch2, tstamp] = ECG_Link(2);
    if (a ~=0)
%        tstamp(1)- tstamp_Ch(end);
        if ((length(Data_Ch)+a )>=1400)
            n =  (length(Data_Ch)+a)-1400;
            Data_Ch   = [Data_Ch(n+1:end); ch2()];
            tstamp_Ch = [tstamp_Ch(n+1:end); tstamp()];
        else
            Data_Ch  (1:end+a) = [Data_Ch(1:end); ch2()];
            tstamp_Ch(1:end+a) = [tstamp_Ch(1:end); tstamp()];
        end
        if (length(Data_Ch) >=200)       
            ecg_data = Data_Ch /(max (Data_Ch));
            [Rpeak]=ecg_peakdet(ecg_data, 0.6);
            pospeak = Rpeak(:, 1);
            subplot (2,1,1);
            plot (ecg_data); 
            for j = 1:(length(pospeak))
                line([pospeak(j) pospeak(j)],[1.2*min(ecg_data) 1.2*max(ecg_data)],'Color',[1 0 0])
                
            end
            if (length(pospeak)>=2)
                
                
                n = floor(pospeak(2))-42;
                if (n>50)
                Rpeak_detected = 1;
                tstamp_Rpeak(1:end+1) = [tstamp_Rpeak(1:end); tstamp_Ch(floor(pospeak(2)))];
                Data_Ch   = Data_Ch(n+1:end);
                tstamp_Ch = tstamp_Ch(n+1:end);
                length(tstamp_Rpeak)
                end
            end

        end
    %plot(0:1400,Data_Ch,pospeak,signal(pospeak),'b^',negpeak,signal(negpeak),'bv')
    drawnow ;

    end 
end


% if (pospeak(j)>= (data_size -a))
%                     line([pospeak(j) pospeak(j)],[1.2*min(ecg_data) 1.2*max(ecg_data)],'Color',[0 1 0])
%                     tstamp_Rpeak = [tstamp_Rpeak(2:end); tstamp_Ch(floor(pospeak(j)))];
%                     Rpeak_detected = 1;
%                 else
