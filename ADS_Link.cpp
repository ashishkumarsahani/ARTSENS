#include "mex.h"
#include <Windows.h>
#include <stdio.h>

HINSTANCE USB_Dll;

__int32 (__cdecl* MSP430_USB_Connection_Status_Ptr)(void);
void(__cdecl* ADS1292R_Start_Live_Data_Packet_from_ADC_Ptr)(void);
void(__cdecl* ADS1292R_Stop_Live_Data_Packet_from_ADC_Ptr)(void);
__int16 (__cdecl* ADS1292R_Get_Live_Data_Packet_from_ADC_Ptr)(__int16 *);

__int16 ecg_Live_Data_Packet[31];

//######################################################
//##############--main mex function--########################
//######################################################
void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, mxArray *prhs[])
{
    int i=0, j=0;
    double *Ch1_Data_Ptr;
    double *Ch2_Data_Ptr;
    double command = mxGetScalar(prhs[0]);
    
    if (command ==1)
    {
        USB_Dll = LoadLibrary("ADS1292R_USB_lib.dll");
        if(USB_Dll != NULL)
        {
            mexPrintf("ADS DLL is loaded.");
            MSP430_USB_Connection_Status_Ptr = (__int32 (__cdecl*)(void))GetProcAddress(USB_Dll, "MSP430_USB_Connection_Status");
            ADS1292R_Start_Live_Data_Packet_from_ADC_Ptr = (void(__cdecl*)(void))GetProcAddress(USB_Dll, "ADS1292R_Start_Live_Data_Packet_from_ADC");
            ADS1292R_Stop_Live_Data_Packet_from_ADC_Ptr = (void(__cdecl*)(void))GetProcAddress(USB_Dll, "ADS1292R_Stop_Live_Data_Packet_from_ADC");
            ADS1292R_Get_Live_Data_Packet_from_ADC_Ptr = (__int16(__cdecl *)(__int16 *))GetProcAddress(USB_Dll, "ADS1292R_Get_Live_Data_Packet_from_ADC");
            if(MSP430_USB_Connection_Status_Ptr!=NULL && ADS1292R_Start_Live_Data_Packet_from_ADC_Ptr!=NULL && ADS1292R_Stop_Live_Data_Packet_from_ADC_Ptr!=NULL && ADS1292R_Get_Live_Data_Packet_from_ADC_Ptr !=NULL)
               mexPrintf("ADS DLL all functions are found."); 
            
            mexPrintf("\nWaiting for USB to connect ...");
            while(0 != MSP430_USB_Connection_Status_Ptr()){}
            mexPrintf("\nUSB device is now connected.");
            ADS1292R_Start_Live_Data_Packet_from_ADC_Ptr();
            Sleep(1000);
            mexPrintf("\nLive streamimg started.");
        }
        else
        {
            mexPrintf("DLL couldn't be loaded");
        }
    }
    else if (command ==2)
    {
        plhs[0] = mxCreateDoubleMatrix(14, 1, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(14, 1, mxREAL);
        plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
        plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
        plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
        Ch1_Data_Ptr = mxGetPr(plhs[0]);
        Ch2_Data_Ptr = mxGetPr(plhs[1]);
        while(0!=ADS1292R_Get_Live_Data_Packet_from_ADC_Ptr(ecg_Live_Data_Packet)){}
        plhs[2] = mxCreateDoubleScalar(double(ecg_Live_Data_Packet[0]));
        plhs[3] = mxCreateDoubleScalar(double(ecg_Live_Data_Packet[1]));
        plhs[4] = mxCreateDoubleScalar(double(ecg_Live_Data_Packet[2]));
        for(j = 0;j<14;j++)
        {
            Ch1_Data_Ptr[j] = ecg_Live_Data_Packet[j*2+3];
            Ch2_Data_Ptr[j] = ecg_Live_Data_Packet[j*2+4];
        }
    }
    else if (command ==3)
    {
       ADS1292R_Stop_Live_Data_Packet_from_ADC_Ptr();
       mexPrintf("ECG device is closed ...");
    }
    else
    {
        mexPrintf("Invalid Command ...");
    }
}