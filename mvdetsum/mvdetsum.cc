#include "mi_time.h"
#include "mvdethghlib.h"
#include "arg_mvdetsum.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValMvdetsum* argval = new ArgValMvdetsum;
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


    // push image to 2d hist
    HistDataNerr2d** hd2d_arr = new HistDataNerr2d* [ntime];
    for(long itime = 0; itime < ntime; itime ++){
        hd2d_arr[itime] = new HistDataNerr2d;
        hd2d_arr[itime]->Init(img_info->GetNaxesArrElm(0), 0, img_info->GetNaxesArrElm(0),
                              img_info->GetNaxesArrElm(1), 0, img_info->GetNaxesArrElm(1));
        hd2d_arr[itime]->SetOvalArr(img_info->GetNpixelImg(), debias_img_arr[itime]);
    }

    printf("--- load parameters ---\n");
    HistInfo1d* hi1d_vel   = new HistInfo1d;
    HistInfo1d* hi1d_psi   = new HistInfo1d;
    LoadHi1dVel(argval->GetVelDat(), hi1d_vel);
    hi1d_psi->InitSetByNbin(0.0, 2.0* M_PI, argval->GetNpsi());
    hi1d_vel->Print(stdout);
    hi1d_psi->Print(stdout);
    printf("=== load parameters ===\n");

    for(long ivel = 0; ivel < hi1d_vel->GetNbin(); ivel ++){
        printf("ivel = %ld\n", ivel);
        double vel = hi1d_vel->GetBinCenter(ivel);
        double* cube_arr = new double[hi1d_psi->GetNbin() * img_info->GetNpixelTotal()];
        for(long ipsi = 0; ipsi < hi1d_psi->GetNbin(); ipsi ++){
            double psi = hi1d_psi->GetBinCenter(ipsi);
            HistDataNerr2d* hd2d_img = NULL;
            GenDetImg(hd2d_arr, time_arr, ntime,
                      vel, psi, &hd2d_img);
            for(long iarr = 0; iarr < hd2d_img->GetNbin(); iarr ++){
                long index = iarr + ipsi * img_info->GetNpixelTotal();
                cube_arr[index] = hd2d_img->GetOvalArr()->GetValElm(iarr);
            }
            delete hd2d_img;
        }

        char tag[kLineSize];
        sprintf(tag, "cube_%2.2ld", ivel);
        MifImgInfo* img_info_cube = new MifImgInfo;
        img_info_cube->InitSetCube(1, 1, 1,
                                   img_info->GetNaxesArrElm(0),
                                   img_info->GetNaxesArrElm(1),
                                   hi1d_psi->GetNbin());
        MifFits::OutFitsCubeD(argval->GetOutdir(),
                              argval->GetOutfileHead(),
                              tag,
                              3, 
                              bitpix,
                              img_info_cube->GetNaxesArr(),
                              cube_arr);
        delete img_info_cube;
        delete [] cube_arr;
    }

    double time_ed = MiTime::GetTimeSec();
    printf("time_ed - time_st = %e\n", time_ed - time_st);
    
    return status_prog;
}
