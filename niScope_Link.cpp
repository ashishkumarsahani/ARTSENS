///Compilation command:   mex -IC:\Progra~2\IVIFou~1\VISA\WinNT\include -IC:\Progra~2\IVIFou~1\IVI\Bin -LC:\Progra~2\IVIFou~1\IVI\Lib\msc -LC:\Progra~2\IVIFou~1\VISA\WinNT\lib\msc -LC:\Progra~2\IVIFou~1\IVI\Lib\msc -lniScope niScope_Link.cpp
///run command: frame = niScope_Link(2,5195);tic;for i = 1:1000
///frame = niScope_Link(2,5195);plot(frame); drawnow;end;toc

#include <windows.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <conio.h>
#include <process.h> //used for thread
#include "mex.h"
#include <niScope.h>
#include <time.h>
#include <Winbase.h>


void cdecl AcqThread(LPVOID pVoid);
void cdecl saveThread(LPVOID pVoid);


ViSession vi;
ViInt32 numWaveform;
ViInt32 actualRecordLength;
struct niScope_wfmInfo *wfmInfoPtr = NULL;
ViReal64 *scaledWfmPtr = NULL;
void* binaryWfmPtr = NULL;
ViInt32 stop;
int frame_Count;
struct timeval tv;
float *local_frame_Store;
double *frame_new;
long *time_Stamp_Store;
int write_Frame_Index=1;
int read_Frame_Index=1;

int frame_Num_Lim = 1000;
int thread_Started =0;
double points_Per_Frame;
int stop_Acq =0;
long double start_time =0 ;
HANDLE USMutex; 

char *file_Name;
HANDLE saveMutex; 
float *save_frame_Store;
FILE *fs ;
int save_Frame_Index=1;
int start_save;
int total_frames;
int save_thread_Started =0;

#define SWAP_UINT32(x) (((x) >> 24) | (((x) & 0x00FF0000) >> 8) | (((x) & 0x0000FF00) << 8) | ((x) << 24))

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
    double command;
    float *out_Frame;
    int err_Code;
    char *dev_Name;
    double saveUS;
            
    command = mxGetScalar(prhs[0]);
    
    if(command ==1.0 && thread_Started==0)
    {
        points_Per_Frame = mxGetScalar(prhs[1]);
        dev_Name = (char *) mxCalloc(mxGetN(prhs[2])+1, sizeof(char));
        mxGetString(prhs[2], dev_Name, mxGetN(prhs[2])+1);
        
        file_Name = (char *) mxCalloc(mxGetN(prhs[3])+1, sizeof(char));
        mxGetString(prhs[3], file_Name, mxGetN(prhs[3])+1);
        fs = fopen(file_Name, "wb");
        
        err_Code = niScope_init (dev_Name, NISCOPE_VAL_FALSE, NISCOPE_VAL_FALSE, &vi);
        mexPrintf("%s",dev_Name);
        
        err_Code =  err_Code+ niScope_ConfigureVertical (vi, "0", 20, 0.0,NISCOPE_VAL_AC, 1.0, NISCOPE_VAL_TRUE);
        err_Code =  err_Code+ niScope_ConfigureHorizontalTiming (vi, 100000000.0, 5195, 0.0, 1, VI_TRUE);
        err_Code =  err_Code+ niScope_ConfigureTriggerDigital (vi, NISCOPE_VAL_PFI_1, NISCOPE_VAL_POSITIVE, 0.0, 0.0);
        err_Code =  err_Code+ niScope_ActualNumWfms (vi, "0", &numWaveform);
        err_Code =  err_Code+ niScope_ActualRecordLength (vi, &actualRecordLength);

        wfmInfoPtr = (niScope_wfmInfo*)malloc (sizeof (struct niScope_wfmInfo) * numWaveform);
        local_frame_Store = (float*) malloc (sizeof(float) * actualRecordLength * numWaveform*frame_Num_Lim);
        save_frame_Store = (float*) malloc (sizeof(float) * actualRecordLength * numWaveform*frame_Num_Lim);
        frame_new = (double*) malloc (sizeof(double) * actualRecordLength * numWaveform);
        if((local_frame_Store==NULL)&&(save_frame_Store==NULL))
        {
            mexPrintf("Unable to alloct memory.\n");
            return;
        }
        
        
        time_Stamp_Store = (long*) malloc (sizeof(long)*frame_Num_Lim);
        scaledWfmPtr = (ViReal64*) malloc (sizeof (ViReal64) * actualRecordLength * numWaveform);
        err_Code =  err_Code+ niScope_Commit (vi); 
        
        USMutex = CreateMutex( 
        NULL,              // default security attributes
        FALSE,             // initially not owned
        NULL);             // unnamed mutex
        write_Frame_Index=1;
        read_Frame_Index=1;
        save_Frame_Index=1;
        gettimeofday(&tv);
        start_time = 1000 * tv.tv_sec + (tv.tv_usec/1000);
        
        saveMutex = CreateMutex( 
        NULL,              // default security attributes
        FALSE,             // initially not owned
        NULL);             // unnamed mutex
        
        
        _beginthread( AcqThread, 0, &ThreadNr );
        
        mexPrintf("waiting.");
        while(thread_Started !=1)
        {
            Sleep(1);
            mexPrintf(".");
        }
        plhs[0] = mxCreateDoubleScalar(start_time);    
        mexPrintf(" \nScope was initialized and acquisition thread was started. Err_Code = %ld\n ", err_Code);
        mexPrintf ("start_time = %Le\n",start_time);
    }
    else if(command ==2.0)
    {
        if(thread_Started ==1)
        {
            WaitForSingleObject( 
                USMutex,    // handle to mutex
                INFINITE);  // no time-out interval
            while(read_Frame_Index == write_Frame_Index)
            {
                ReleaseMutex(USMutex);
                Sleep(1);
                WaitForSingleObject( 
                    USMutex,    // handle to mutex
                    INFINITE);  // no time-out interval
            }
            ReleaseMutex(USMutex);
             plhs[0] = mxCreateNumericMatrix(points_Per_Frame, 1, mxSINGLE_CLASS, mxREAL);
//             plhs[0] = mxCreateDoubleMatrix(points_Per_Frame, 1, mxREAL);
            out_Frame = (float *)mxGetPr(plhs[0]);
            plhs[1] = mxCreateDoubleScalar(time_Stamp_Store[read_Frame_Index-1]);
            memcpy(out_Frame,((float*)(&(local_frame_Store[(read_Frame_Index-1)*(int)points_Per_Frame]))),points_Per_Frame*sizeof(float));
            
            
            saveUS = mxGetScalar(prhs[2]);
            if(saveUS ==1)
            {
//                 WaitForSingleObject( 
//                      saveMutex,    // handle to mutex
//                      INFINITE);  // no time-out interval
                memcpy(((float*)(&(save_frame_Store[(save_Frame_Index-1)*(int)points_Per_Frame]))),((float*)(&(local_frame_Store[(read_Frame_Index-1)*(int)points_Per_Frame]))),points_Per_Frame*sizeof(float));
                save_frame_Store[(save_Frame_Index-1)*(int)points_Per_Frame] = (float) time_Stamp_Store[read_Frame_Index-1];
                
                
                if (save_Frame_Index == 500 && save_thread_Started==0)
                {
                    start_save = 0;
                    total_frames = 500;
                    _beginthread( saveThread, 0, &ThreadNr );
                }
                else if (save_Frame_Index == 1000 && save_thread_Started==0)
                {
                    start_save = 500;
                    total_frames = 500;
                    save_Frame_Index = 0;
                    _beginthread( saveThread, 0, &ThreadNr );              
                }
                else if (save_thread_Started==1 && (save_Frame_Index == 1000 || save_Frame_Index == 500 ))
                {
                    mexPrintf("file writing is not over yet...!\n");
                    if (save_Frame_Index == 1000)  save_Frame_Index = 1;
                }
                save_Frame_Index = (save_Frame_Index)+1;
//                 ReleaseMutex(saveMutex);
            }
            WaitForSingleObject( 
                    USMutex,    // handle to mutex
                    INFINITE);  // no time-out interval
            read_Frame_Index = (read_Frame_Index)%frame_Num_Lim+1;
            ReleaseMutex(USMutex);
        }
        
        else
        {
            mexPrintf("Acquisition thread is not running!");
        }
    }
    else if(command ==3.0)
    {
         WaitForSingleObject( 
            USMutex,    // handle to mutex
            INFINITE);  // no time-out interval
         if(stop_Acq!=1) stop_Acq =1;
         ReleaseMutex(USMutex);  
         
         saveUS = mxGetScalar(prhs[2]);
         if(saveUS ==1)
         { 
             while(save_thread_Started !=0)
             {
                 Sleep(1);
             }
//              while (save_thread_Started ==1) {Sleep(10);}
             if (save_Frame_Index >= 500)        start_save = 500;
             else if (save_Frame_Index < 500)  start_save = 0;
             
             total_frames = (save_Frame_Index-1) - start_save;
             _beginthread( saveThread, 0, &ThreadNr );  
             save_Frame_Index = 1;
             
         }
         
         Sleep(100);
         while (save_thread_Started !=0) 
         {
             Sleep(10);
         }
         fclose(fs);
         mexPrintf("done\n");
         
    }
    else if(command ==4.0)
    {
        if (thread_Started==1)
        {
            WaitForSingleObject( 
                    USMutex,    // handle to mutex
                    INFINITE);  // no time-out interval
            write_Frame_Index=1;
            read_Frame_Index=1;
            start_time =mxGetScalar(prhs[1]);
            ReleaseMutex(USMutex);
        }
        else
        {
            mexPrintf("Acquisition thread is not running!2");
        }
        
    }
    else
    {
        mexPrintf("Invalid command!");
    }
    //start your thread
    

//----for testing----
//Sleep(20000); //will make the mex function stay around for 20 seconds so you can see the mexPrintfs of the thread if you want
}


//################################################
//##########--Thread Code--#########################
//################################################
void cdecl AcqThread(LPVOID pVoid)
{
thread_Started=1;
mexPrintf("Scope is now running");
while(1)
{
    WaitForSingleObject( 
                USMutex,    // handle to mutex
                INFINITE);  // no time-out interval
    if(stop_Acq==1) 
    {
        stop_Acq=0;
        
        ReleaseMutex(USMutex);
        niScope_close (vi);
        break;
    }
    ReleaseMutex(USMutex);
    
    if(write_Frame_Index== read_Frame_Index)
    {
        WaitForSingleObject( 
                    USMutex,    // handle to mutex
                    INFINITE);  // no time-out interval
        gettimeofday(&tv);
        niScope_Read (vi, "0", 5, actualRecordLength, frame_new , wfmInfoPtr);
        time_Stamp_Store[write_Frame_Index-1] = (1000 * tv.tv_sec+(tv.tv_usec/1000) )- start_time ;
        
        for (int i =0 ; i <= (int)points_Per_Frame; i++)
        {
            local_frame_Store[(write_Frame_Index-1)*(int)points_Per_Frame + i] = (float)frame_new[i ];
        }
        
        
        write_Frame_Index=(write_Frame_Index)%frame_Num_Lim+1; 
        ReleaseMutex(USMutex);
    }
    else
    {
        gettimeofday(&tv);
        niScope_Read (vi, "0", 5, actualRecordLength, (frame_new), wfmInfoPtr);
        time_Stamp_Store[write_Frame_Index-1] = (1000 * tv.tv_sec+(tv.tv_usec/1000) )- start_time ;
        
        for (int i =0 ; i <= (int)points_Per_Frame; i++)
        {
            local_frame_Store[(write_Frame_Index-1)*(int)points_Per_Frame + i] = (float)frame_new[i ];
        }
        
        
        WaitForSingleObject( 
                    USMutex,    // handle to mutex
                    INFINITE);  // no time-out interval
        write_Frame_Index=(write_Frame_Index)%frame_Num_Lim+1; 
        ReleaseMutex(USMutex);
    } 
}
mexPrintf("\nScope is now stopped ...");
thread_Started=0;
}

//################################################
//##########--save Thread Code--##################
//################################################
void cdecl saveThread(LPVOID pVoid)
{
    long i;
    
save_thread_Started=1;

mexPrintf("\nUS saving data  ");
for ( i = 0 ; i <(total_frames*(int)points_Per_Frame) ; i++)
{
    save_frame_Store[(start_save *(int)points_Per_Frame) + i] = ReverseFloat(  save_frame_Store[(start_save *(int)points_Per_Frame) + i]);
}

if(fwrite(((float*)(&(save_frame_Store[start_save*(int)points_Per_Frame]))), sizeof(float),( total_frames*(int)points_Per_Frame), fs) != ( total_frames*(int)points_Per_Frame))
{ 
    printf("File write error.");
}
save_thread_Started=0;
}

