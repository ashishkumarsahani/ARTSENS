function [xi,yi,N] = fftinterpol(x,y,interval_req)
%
%  interval_req is interval required for x 
%  so xi is not fix here
%  xi(0) = x(0), xi(end)=x(end) all points in xi spaced about interval_req
%  linearly interpolated such that all 
%  (x(k),y(k)) point corresponds to same (xi(k*N),yi(k*N) point in
%  interpolated  signal
%  here N = no. of interpolated points or N= length(xi)/length(x)
%

len_x= length(x);
len_y= length(y);
if (len_x == len_y)
    diff_x=diff(x);
    N= round(median(diff_x));
    N=N/interval_req;
    
    
    temp1= (0:N:((len_x -1)*N))';
    temp2= (0:((len_x -1) *N))';
    inter_x = interp1(temp1,x,temp2,'linear');
    xi = inter_x;
    
%     [yi2,b]=resample(y,length(xi),len_y,5);
    yi3 = spline(x,y,xi);
%     fft_y=fft(y);
%     if (mod(len_y,2)==0) %even
%         center_y=len_y/2;
%         fft_y_inter=[fft_y(1:( center_y+1)) ; zeros((N -1)*(len_y -1)-1,1);fft_y((center_y +1):end)];
%     else                 %odd
%         center_y =(len_y+1)/2;
%         fft_y_inter=[fft_y(1:(center_y)) ; zeros((N -1)*(len_y-1),1);fft_y((center_y +1):end)];
%     end
%     inter_y=ifft(fft_y_inter);
%     yi = inter_y * 12;
    yi=yi3;
%     figure; plot (x,y,'b-o'); hold on;  plot (xi,yi,'r-+');
%     diff_dist = diff(yi,2);
%     diff_dist = [0;0;(diff_dist / max(diff_dist))];
%     plot (xi,diff_dist,'g');
    
    
    
if length(xi)~= length(yi)
    disp( 'output lengths are not equal ');
    xi = x;
    yi = y;
    N= 1;
end
if isreal(yi)==0
    disp( 'output yi is not real ');
    xi = x;
    yi = y;
    N= 1;
end

else
    disp( 'length x should be equal to length y');
    xi = [0];
    yi = [0];
    N= 0;
end

