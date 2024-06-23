#include <windows.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <conio.h>
#include <process.h> //used for thread
#include "mex.h"
#include <time.h>
#include <Winbase.h>

HINSTANCE USB_Dll;

__int32 (__cdecl* MSP430_USB_Connection_Status_Ptr)(void);
void(__cdecl* MSP430_USB_Connection_Close_Ptr )(void); //void  MSP430_USB_Connection_Close(void );
__int8 (__cdecl* ADS1292R_Select_Filter_Ptr)(__int8);
void(__cdecl* ADS1292R_Start_Live_Data_Packet_from_ADC_Ptr)(void);
void(__cdecl* ADS1292R_Stop_Live_Data_Packet_from_ADC_Ptr)(void);
__int16 (__cdecl* ADS1292R_Get_Live_Data_Packet_from_ADC_Ptr)(__int16 *);
void cdecl ECGAcqThread(LPVOID pVoid);
void cdecl saveThread(LPVOID pVoid);


__int16 ecg_Live_Data_Packet[31];

const int samples_per_frame = 14 ; //
const int max_Store_Size = 100*samples_per_frame; // the storage size is a multiple of 14 for ease in design (each frame as 14 bytes)
float ecg_Live_Data_Ch1_Store[max_Store_Size];
float ecg_Live_Data_Ch2_Store[max_Store_Size];
float time_Stamp_Store[max_Store_Size];
long double start_time ;
__int16 ecg_Lead_Status;
__int16 ecg_HR;
__int16 ecg_Resp_Rate;

int sampling_Rate = 500;
int write_Index = 0;
int stop_Acq=0;
int thread_Started=0;
HANDLE ECGMutex; 
struct timeval tv;

char *file_Name;
HANDLE saveMutex; 
float save_frame_Store [2 * max_Store_Size] ; // saving for 28 sec data [1000 * 2 * samples_per_frame]
FILE *fs ;
int save_Index=0;
int start_save;
int total_frames;
int save_thread_Started =0;
int file_open = 0;
int saveECG = 0;


#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif
 
struct timezone 
{
  int  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};
 
int gettimeofday(struct timeval *tv)
{
  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag;
 
  if (NULL != tv)
  {
    GetSystemTimeAsFileTime(&ft);
 
    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;
 
    /*converting file time to unix epoch*/
    tmpres -= DELTA_EPOCH_IN_MICROSECS; 
    tmpres /= 10;  /*convert into microseconds*/
    tv->tv_sec = (long)(tmpres / 1000000UL);
    tv->tv_usec = (long)(tmpres % 1000000UL);
  }
 
  return 0;
}

//
// conversion big to little endian
// 
float ReverseFloat( const float inFloat )
{
   float retVal;
   char *floatToConvert = ( char* ) & inFloat;
   char *returnFloat = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];

   return retVal;
}

//######################################################
//##############--main mex function--########################
//######################################################
void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, mxArray *prhs[])
{
    
    int ThreadNr; //needed for the thread
    int i=0, j=0;
    int n=0;
    double *Ch1_Data_Ptr;
    double *Ch2_Data_Ptr;
    double command = mxGetScalar(prhs[0]);
    double *time_Stamp_Ptr;
    
    float init_data[2* 14] = {0} ;//[2* 14];
    
    if (command ==1 && thread_Started ==0)
    {
        USB_Dll = LoadLibrary("ADS1292R_USB_lib.dll");
        if(USB_Dll != NULL)
        {
            mexPrintf("ADS DLL is loaded.");
            MSP430_USB_Connection_Status_Ptr = (__int32 (__cdecl*)(void))GetProcAddress(USB_Dll, "MSP430_USB_Connection_Status");
            MSP430_USB_Connection_Close_Ptr = (void(__cdecl*)(void))GetProcAddress(USB_Dll, "MSP430_USB_Connection_Close");
            ADS1292R_Select_Filter_Ptr = (__int8 (__cdecl*)(__int8))GetProcAddress(USB_Dll, "ADS1292R_Select_Filter");
            ADS1292R_Start_Live_Data_Packet_from_ADC_Ptr = (void(__cdecl*)(void))GetProcAddress(USB_Dll, "ADS1292R_Start_Live_Data_Packet_from_ADC");
            ADS1292R_Stop_Live_Data_Packet_from_ADC_Ptr = (void(__cdecl*)(void))GetProcAddress(USB_Dll, "ADS1292R_Stop_Live_Data_Packet_from_ADC");
            ADS1292R_Get_Live_Data_Packet_from_ADC_Ptr = (__int16(__cdecl *)(__int16 *))GetProcAddress(USB_Dll, "ADS1292R_Get_Live_Data_Packet_from_ADC");
            if(MSP430_USB_Connection_Status_Ptr!=NULL && MSP430_USB_Connection_Close_Ptr!=NULL && ADS1292R_Select_Filter_Ptr !=NULL && ADS1292R_Start_Live_Data_Packet_from_ADC_Ptr!=NULL && ADS1292R_Stop_Live_Data_Packet_from_ADC_Ptr!=NULL && ADS1292R_Get_Live_Data_Packet_from_ADC_Ptr !=NULL)
               mexPrintf("ADS DLL all functions are found."); 
            
            mexPrintf("\nWaiting for USB to connect ...");
            while(0 != MSP430_USB_Connection_Status_Ptr()){}
            mexPrintf("\nUSB device is now connected.");
            ADS1292R_Select_Filter_Ptr(0);
            ADS1292R_Start_Live_Data_Packet_from_ADC_Ptr();
            Sleep(1000);
            mexPrintf("\nLive streamimg started.");
            write_Index=0;
            
            
//             save_frame_Store = (float*) malloc (sizeof(float) * 1000 * 2 * samples_per_frame);
            save_Index=0;
            saveECG = 0;
            
            ECGMutex = CreateMutex( 
            NULL,              // default security attributes
            FALSE,             // initially not owned
            NULL);             // unnamed mutex
            write_Index= 0;
            stop_Acq=0;
            gettimeofday(&tv);
            start_time = 1000 * tv.tv_sec + (tv.tv_usec/1000);
            _beginthread( ECGAcqThread, 0, &ThreadNr );
        }
        else
        {
            FreeLibrary(USB_Dll); 
            mexPrintf("DLL couldn't be loaded");
        }
    }
    else if (command ==2)
    {
        if(write_Index ==0)
        {
            plhs[0] = mxCreateDoubleScalar(0);
            plhs[1] = mxCreateDoubleScalar(0);
            plhs[2] = mxCreateDoubleScalar(0);
            plhs[3] = mxCreateDoubleScalar(0);
            plhs[4] = mxCreateDoubleScalar(0);
            plhs[5] = mxCreateDoubleScalar(0);
            plhs[6] = mxCreateDoubleScalar(0);
            return;
        }
        else
        {
            WaitForSingleObject( 
                ECGMutex,    // handle to mutex
                INFINITE);  // no time-out interval
            plhs[0] = mxCreateDoubleScalar(write_Index);
            plhs[1] = mxCreateNumericMatrix(write_Index, 1, mxSINGLE_CLASS, mxREAL);
            plhs[2] = mxCreateNumericMatrix(write_Index, 1, mxSINGLE_CLASS, mxREAL);
            Ch1_Data_Ptr =  mxGetPr(plhs[1]);
            Ch2_Data_Ptr =  mxGetPr(plhs[2]);
            memcpy(Ch1_Data_Ptr,((float*)ecg_Live_Data_Ch1_Store),write_Index*sizeof(float));
            memcpy(Ch2_Data_Ptr,((float*)ecg_Live_Data_Ch2_Store),write_Index*sizeof(float));
            plhs[3] = mxCreateNumericMatrix(write_Index, 1, mxSINGLE_CLASS, mxREAL);
            time_Stamp_Ptr =  mxGetPr(plhs[3]);
            memcpy(time_Stamp_Ptr,((float*)time_Stamp_Store),write_Index*sizeof(float));
            plhs[4] = mxCreateDoubleScalar((double)ecg_Lead_Status);
            plhs[5] = mxCreateDoubleScalar((double)ecg_HR);
            plhs[6] = mxCreateDoubleScalar((double)ecg_Resp_Rate);
//             n = write_Index;
             write_Index = 0;
             saveECG = mxGetScalar(prhs[1]);
            ReleaseMutex(ECGMutex);
            
//             if(saveECG ==1)
//             {
//                 for (i =0 ; i < (n/samples_per_frame) ; i++)
//                 {
// //                     memcpy(&(save_frame_Store[save_Index * 2 * samples_per_frame ]),&(Ch2_Data_Ptr[i*samples_per_frame]),samples_per_frame*sizeof(float));
// //                     memcpy(&(save_frame_Store[(save_Index * 2 * samples_per_frame) +samples_per_frame]),&(time_Stamp_Ptr[i*samples_per_frame]),samples_per_frame*sizeof(float));
//                     for (j=0; j <samples_per_frame ; j++)
//                     {
//                         save_frame_Store[(save_Index * 2 * samples_per_frame) + j] = Ch2_Data_Ptr[j +(i*samples_per_frame)]; 
//                         save_frame_Store[(save_Index * 2 * samples_per_frame) +samples_per_frame +j]=  time_Stamp_Ptr[j + (i*samples_per_frame) ];
//                     }
//                     save_Index = save_Index +1;
//                     if (save_Index == 500 && save_thread_Started==0)
//                     {
//                         start_save = 0;
//                         total_frames = 500;
//                         _beginthread( saveThread, 0, &ThreadNr );
//                     }
//                     else if (save_Index == 1000 && save_thread_Started==0)
//                     {
//                         start_save = 500;
//                         total_frames = 500;
//                         save_Index = 0;
//                         _beginthread( saveThread, 0, &ThreadNr );              
//                     }
//                     else if (save_thread_Started==1 && (save_Index == 1000 || save_Index == 500 ))
//                     {
//                         mexPrintf("file writing is not over yet...!\n");
//                         if (save_Index == 1000)  save_Index = 0;
//                     }
// 
//                 }
//             }  // saving statement
           
           
            
        }
    }
    else if (command ==3)
    {
        stop_Acq=1;
        mexPrintf("\nECG closing ...");
        saveECG = mxGetScalar(prhs[1]);
         if(saveECG ==1 || file_open == 1)
         {
             saveECG = 0;
             Sleep(10);
             mexPrintf("saved ECG file closing ...");
             while (save_thread_Started !=0) 
             {
                 Sleep(10);
             }
             fclose(fs);
             file_open = 0;
             mexPrintf("done\n");
         }
    }
    
    
    else if (command == 4)
    {   
        if (thread_Started ==1)
        {
            WaitForSingleObject( 
                ECGMutex,    // handle to mutex
                INFINITE);  // no time-out interval
            start_time =mxGetScalar(prhs[1]);
            write_Index=0;
            ReleaseMutex(ECGMutex);
        }
        else
        {
            write_Index=0;
            ECGMutex = CreateMutex( 
            NULL,              // default security attributes
            FALSE,             // initially not owned
            NULL);             // unnamed mutex
            write_Index= 0;
            stop_Acq=0;
            
            start_time =mxGetScalar(prhs[1]);
            
            _beginthread( ECGAcqThread, 0, &ThreadNr );
            
        }
        
            
    }   
    else if (command ==5  && thread_Started ==1)
    {
        
        if (file_open == 1) fclose(fs);
        
        WaitForSingleObject( 
                    ECGMutex,    // handle to mutex
                    INFINITE);  // no time-out interval
        start_time =mxGetScalar(prhs[1]);
        write_Index=0;
        
        
        
        file_Name = (char *) mxCalloc(mxGetN(prhs[2])+1, sizeof(char));
        mxGetString(prhs[2], file_Name, mxGetN(prhs[2])+1);
        fs = fopen(file_Name, "wb");
        file_open = 1;
        save_Index=0;
        ReleaseMutex(ECGMutex);
        
//         init_data = 0;//(float*) mxCreateNumericMatrix((2 * 14), 1, mxSINGLE_CLASS, mxREAL);
        init_data[0] = 28;
        init_data[1] = 1;
        init_data[2] = sampling_Rate;
        init_data[3] = samples_per_frame;
        init_data[4] = samples_per_frame;
        init_data[5] = start_time;
//         memcpy(&(save_frame_Store[save_Index * 2 * samples_per_frame ]),init_data,(2*14*sizeof(float)));
//         save_Index=1;
        for (i=0; i <28 ; i++)
        {
           init_data[i] = ReverseFloat(init_data[i]);
        }
        if(fwrite(((float*)(&init_data)), sizeof(float),(28), fs) != ( 28))
        { 
            printf("File write error.");
        }
        saveECG = 1;
        mexPrintf("ECG is running & file created ...\n");
        
    }
    
    else if (command ==6 && thread_Started ==1)
    {
        
       
        
//          if(saveECG ==1)
//          { 
//              while(save_thread_Started !=0)
//              {
//                  Sleep(1);
//              }
// //              while (save_thread_Started ==1) {Sleep(10);}
//          if (save_Index >= 500)      start_save = 500;
//          else if (save_Index < 500)  start_save = 0;
// 
//          total_frames = (save_Index) - start_save;
//          _beginthread( saveThread, 0, &ThreadNr );  
//          save_Index = 0;
//              
//          }
         saveECG = 0;
         Sleep(10);
         mexPrintf("saved ECG file closing ...");
         while (save_thread_Started !=0) 
         {
             Sleep(10);
         }
         fclose(fs);
         file_open = 0;
         mexPrintf("done");
         
    }
    else if (command ==1  && thread_Started ==1) mexPrintf("ECG is already running ...\n");
    else
    {
        mexPrintf("\nInvalid Command / ECG off ...");
    }
}


void cdecl ECGAcqThread(LPVOID pVoid)
{
    int j =0;
    int TV_sec_prev =0;
    thread_Started=1;
    mexPrintf("\nECG is now running ...");
    
    while(1)
    {
        if(stop_Acq==1) 
        {
            stop_Acq=0;
            break;
        }
        gettimeofday(&tv);
        TV_sec_prev = tv.tv_sec;
        while(0!=ADS1292R_Get_Live_Data_Packet_from_ADC_Ptr(ecg_Live_Data_Packet))
        {
            gettimeofday(&tv);
            if ((tv.tv_sec - TV_sec_prev) >= 2) 
            {
             mexPrintf("\n error ECG Acq");
             break;
            }
        }
        
        WaitForSingleObject( 
                ECGMutex,    // handle to mutex
                INFINITE);  // no time-out interval
        ecg_Resp_Rate = ecg_Live_Data_Packet[0];
        ecg_HR = ecg_Live_Data_Packet[1];
        ecg_Lead_Status = ecg_Live_Data_Packet[2];
        
        for(j = 0;j<samples_per_frame;j++)
        {
            ecg_Live_Data_Ch1_Store[write_Index+j] = (float)ecg_Live_Data_Packet[j*2+3];
            ecg_Live_Data_Ch2_Store[write_Index+j] = (float)ecg_Live_Data_Packet[j*2+4];
            time_Stamp_Store[write_Index+j] = (((1000 * tv.tv_sec) + (tv.tv_usec/1000) - (1000*(13-j)/sampling_Rate)) -start_time);
            save_frame_Store[j*2] = ecg_Live_Data_Ch2_Store[write_Index+j];
            save_frame_Store[(j*2)+1] = (float)time_Stamp_Store[write_Index+j];
        }
        
        write_Index=(write_Index + samples_per_frame)%max_Store_Size; 
        if (write_Index == 0) saveECG =0;
        ReleaseMutex(ECGMutex);
        if(saveECG ==1 && file_open == 1)
        {
            save_thread_Started=1;
            for (j=0; j <(2 * samples_per_frame) ; j++)
            {
               save_frame_Store[j] = ReverseFloat(save_frame_Store[j]);
            }
            if(fwrite(((float*)(&save_frame_Store)), sizeof(float),(2 * samples_per_frame), fs) != ( 2 * samples_per_frame))
            { 
                printf("File write error.");
            }
            save_thread_Started=0;
        }
    }
    
    ADS1292R_Stop_Live_Data_Packet_from_ADC_Ptr();
    thread_Started=0;
    Sleep(100);
    MSP430_USB_Connection_Close_Ptr();
//     FreeLibrary(USB_Dll); 
    mexPrintf("\nECG device is closed ...");
    
}


//################################################
//##########--save Thread Code--##################
//################################################
void cdecl saveThread(LPVOID pVoid)
{
    long k;
    
save_thread_Started=1;

mexPrintf("\nsaving data..");
for ( k = 0 ; k <(total_frames* 2 *(int)samples_per_frame) ; k++)
{
    save_frame_Store[(start_save * 2*(int)samples_per_frame) + k] = ReverseFloat(  save_frame_Store[(start_save *2 *(int)samples_per_frame) + k]);
}

if(fwrite(((float*)(&(save_frame_Store[start_save* 2 *(int)samples_per_frame]))), sizeof(float),( total_frames* 2* (int)samples_per_frame), fs) != ( total_frames* 2* (int)samples_per_frame))
{ 
    printf("File write error.");
}
save_thread_Started=0;
}

