#include "mir_math_util.h"
#include "mir_data1d_ope.h"
#include "mir_hist2d_ope.h"
#include "mir_hist_info.h"
#include "mi_str.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "mifc_gen.h"
#include "arg_mkpsf.h"
#include "../lib/sub.h"

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
    MifImgInfo* img_info_subimg = NULL;
    int bitpix = 0;
    ReadDataList(argval->GetDataList(),
                 argval->GetSubimgDat(),
                 &ntime,
                 &time_arr,
                 &data_img_arr,
                 &img_info_subimg,
                 &bitpix);

    int nclip = 30;
    double significance_clip = 5.0;
    double** std_img_arr    = NULL;
    double** debias_img_arr = NULL;
    double* mean_time_arr   = NULL;
    double* stddev_time_arr = NULL;
    MkStdImgArr(ntime,
                data_img_arr,
                img_info_subimg,
                nclip, significance_clip,
                &std_img_arr,
                &debias_img_arr,
                &mean_time_arr,
                &stddev_time_arr);

    double* median_img_arr = NULL;
    MkMedianImg(ntime,
                debias_img_arr,
                img_info_subimg,
                &median_img_arr);
    
    double xval_peak = 0.0;
    double yval_peak = 0.0;
    GetPeakXY(img_info_subimg,
              median_img_arr,
              argval->GetFunc(),
              argval->GetParFile(),
              argval->GetNbinKernelHalf(),
              &xval_peak, &yval_peak);
    
    HistDataNerr2d* hd2d_psf = NULL;
    MkPsf(img_info_subimg,
          median_img_arr,
          argval->GetNbinPsfHalf(),
          xval_peak, yval_peak,
          &hd2d_psf);

    char file_psf[kLineSize];
    sprintf(file_psf, "%s/%s_psf.dat",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str()); 
    hd2d_psf->Save(file_psf, "x,y,z");

    MifImgInfo* img_info_psf = new MifImgInfo;
    img_info_psf->InitSetImg(1, 1,
                             hd2d_psf->GetHi2d()->GetNbinX(),
                             hd2d_psf->GetHi2d()->GetNbinY());
    MifFits::OutFitsImageD(argval->GetOutdir(),
                           argval->GetOutfileHead(),
                           "psf",
                           2, 
                           bitpix,
                           img_info_psf->GetNaxesArr(),
                           hd2d_psf->GetOvalArr()->GetVal());
    
    double time_ed = MiTime::GetTimeSec();
    printf("time_ed - time_st = %e\n", time_ed - time_st);
    
    return status_prog;
}
