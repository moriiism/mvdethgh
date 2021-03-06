#include "mi_time.h"
#include "mvdethghlib.h"
#include "arg_mkpsf.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValMkpsf* argval = new ArgValMkpsf;
    argval->Init(argc, argv);
    argval->Print(stdout);

    FILE* fp_log = NULL;
    OpenLogfile(argval->GetOutdir(),
                argval->GetOutfileHead(),
                argval->GetProgname(),
                &fp_log);
    argval->Print(fp_log);
    
    long ntime = 0;
    double* time_arr = NULL;
    double** data_img_arr = NULL;
    MifImgInfo* img_info = NULL;
    int bitpix = 0;
    LoadData(argval->GetDataList(),
             argval->GetSubimgDat(),
             &ntime,
             &time_arr,
             &data_img_arr,
             &img_info,
             &bitpix);

    int nclip = 30;
    double significance_clip = 5.0;
    double** std_img_arr    = NULL;
    double** debias_img_arr = NULL;
    double* mean_time_arr   = NULL;
    double* stddev_time_arr = NULL;
    GenStdImgArr(ntime,
                 data_img_arr,
                 img_info,
                 nclip, significance_clip,
                 &std_img_arr,
                 &debias_img_arr,
                 &mean_time_arr,
                 &stddev_time_arr);

    double* median_img_arr = NULL;
    GenMedianImg(ntime,
                 debias_img_arr,
                 img_info,
                 &median_img_arr);
    
    double xval_peak = 0.0;
    double yval_peak = 0.0;
    GetPeakXY(img_info,
              median_img_arr,
              argval->GetNbinKernelHalf(),
              argval->GetValSmooth(),
              &xval_peak, &yval_peak);
    
    HistDataNerr2d* hd2d_psf = NULL;
    GenPsf(img_info,
           median_img_arr,
           argval->GetNbinPsfHalf(),
           xval_peak, yval_peak,
           &hd2d_psf);

    SavePsf(hd2d_psf,
            bitpix,
            argval->GetOutdir(),
            argval->GetOutfileHead());
    
    double time_ed = MiTime::GetTimeSec();
    printf("time_ed - time_st = %e\n", time_ed - time_st);

    delete argval;
    delete img_info;
    delete [] time_arr;
    for(long itime = 0; itime < ntime; itime++){
        delete [] data_img_arr[itime];
        delete [] std_img_arr[itime];
        delete [] debias_img_arr[itime];
    }
    delete [] data_img_arr;
    delete [] std_img_arr;
    delete [] debias_img_arr;
    delete [] mean_time_arr;
    delete [] stddev_time_arr;
    delete [] median_img_arr;
    delete hd2d_psf;
    fclose(fp_log);
    
    return status_prog;
}
