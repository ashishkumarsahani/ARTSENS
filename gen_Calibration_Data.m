function [ calibration_Done, cut_Point, calibration_Frame, Array_for_File] = gen_Calibration_Data(current_Frame, sampling_Rate, samples_Per_Frame, total_No_Of_Frames, reset)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global average_Frame;
global total_Frames;
global current_Frame_Num;
global cut_Location;
global ZCR_Frame;
cut_Point =0;
calibration_Frame = zeros(samples_Per_Frame,1);
Array_for_File = zeros(samples_Per_Frame+2,1);
calibration_Done = 0;
ZCR_Win_Size = 50;

if (reset ==1)
   average_Frame = zeros(samples_Per_Frame,1);
   total_Frames = total_No_Of_Frames;
   cut_Location = 0;
   current_Frame_Num = 0;
   ZCR_Frame = zeros(floor(samples_Per_Frame/ZCR_Win_Size),total_Frames);
end

current_Frame_Num = current_Frame_Num+1;

if(current_Frame_Num < total_Frames)
    average_Frame = average_Frame + (current_Frame-mean(current_Frame))/total_Frames;
else
    cut_Point = find((current_Frame>=(max(average_Frame)/2)),1,'last');
    calibration_Frame = average_Frame;
    calibration_Done = 1;
    Array_for_File = [sampling_Rate;samples_Per_Frame;cut_Point;calibration_Frame];
end
end

