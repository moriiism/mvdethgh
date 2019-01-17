#include "mi_time.h"
#include "mvdethghlib.h"
#include "arg_preproc.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValPreproc* argval = new ArgValPreproc;
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

    SaveCube(ntime, img_info, data_img_arr,
             bitpix,
             argval->GetOutdir(),
             argval->GetOutfileHead(),
             "cube_org");

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

    SaveCube(ntime, img_info, std_img_arr,
             bitpix,
             argval->GetOutdir(),
             argval->GetOutfileHead(),
             "cube_std");

    double* median_img_arr = NULL;
    GenMedianImg(ntime,
                 debias_img_arr,
                 img_info,
                 &median_img_arr);
    
    printf("--- output median image ---\n");
    MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(), "median",
                           2, bitpix,
                           img_info->GetNaxesArr(),
                           median_img_arr);
    printf("=== output median image ===\n");

    double** mvobj_img_arr = NULL;
    GenMvobjImgArr(ntime,
                   img_info,
                   debias_img_arr,
                   median_img_arr,
                   &mvobj_img_arr);

    SaveCube(ntime, img_info, mvobj_img_arr,
             bitpix,
             argval->GetOutdir(),
             argval->GetOutfileHead(),
             "cube_mvobj");
    SaveImgArr(ntime, img_info, mvobj_img_arr,
               bitpix,
               argval->GetOutdir(),
               argval->GetOutfileHead(),
               "mvobj");

    double** conv_img_arr = NULL;
    GenConvPsf(ntime, img_info, mvobj_img_arr,
               argval->GetPsfDat(),
               &conv_img_arr);
    SaveImgArr(ntime, img_info, conv_img_arr,
               bitpix,
               argval->GetOutdir(),
               argval->GetOutfileHead(),
               "conv");
    
    double time_ed = MiTime::GetTimeSec();
    printf("time_ed - time_st = %e\n", time_ed - time_st);
    
    return status_prog;
}
