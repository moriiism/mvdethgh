#include "mi_time.h"
#include "mvdethghlib.h"
#include "arg_mkdetimg.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValMkdetimg* argval = new ArgValMkdetimg;
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
             "none",
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

    // push image to 2d hist
    HistDataNerr2d** hd2d_arr = new HistDataNerr2d* [ntime];
    for(long itime = 0; itime < ntime; itime ++){
        hd2d_arr[itime] = new HistDataNerr2d;
        hd2d_arr[itime]->Init(img_info->GetNaxesArrElm(0), 0, img_info->GetNaxesArrElm(0),
                              img_info->GetNaxesArrElm(1), 0, img_info->GetNaxesArrElm(1));
        hd2d_arr[itime]->SetOvalArr(img_info->GetNpixelImg(), debias_img_arr[itime]);
    }
    string* line_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(argval->GetParDat(), &line_arr, &nline);
    double* vel_arr = new double [nline];
    double* rho_arr = new double [nline];
    double* phi_arr = new double [nline];
    double* psi_arr = new double [nline];
    for(long iline = 0; iline < nline; iline ++){
        sscanf(line_arr[iline].c_str(), "%lf  %lf  %lf  %lf",
               &vel_arr[iline], &rho_arr[iline], &phi_arr[iline], &psi_arr[iline]);
    }
    delete [] line_arr;
    double* stddev_arr    = new double[nline];
    double* xval_mean_arr =  new double[nline];
    double* yval_mean_arr =  new double[nline];
    for(long iline = 0; iline < nline; iline ++){
        double vel = vel_arr[iline];
        double rho = rho_arr[iline];
        double phi = phi_arr[iline];
        double psi = psi_arr[iline];
        double theta = atan(vel);

        int nbin_detimg_half = 10;
        HistDataNerr2d* hd2d_img = NULL;
        GenDetImg(hd2d_arr, time_arr, ntime,
                  nbin_detimg_half,
                  theta, rho, phi, psi, &hd2d_img);
        double xval_mean = 0.0;
        double yval_mean = 0.0;
        double xval2_mean = 0.0;
        double yval2_mean = 0.0;
        double sum_weight = 0.0;
        for(long ibin = 0; ibin < hd2d_img->GetNbin(); ibin ++){
            double xval = hd2d_img->GetHi2d()->GetBinCenterXFromIbin(ibin);
            double yval = hd2d_img->GetHi2d()->GetBinCenterYFromIbin(ibin);
            double weight = hd2d_img->GetOvalElmAtXY(xval, yval);
            xval_mean += xval * weight;
            yval_mean += yval * weight;
            xval2_mean += xval * xval * weight;
            yval2_mean += yval * yval * weight;
            sum_weight += weight;
        }
        xval_mean /= sum_weight;
        yval_mean /= sum_weight;
        xval2_mean /= sum_weight;
        yval2_mean /= sum_weight;
        double xval_stddev = sqrt(xval2_mean - xval_mean * xval_mean);
        double yval_stddev = sqrt(yval2_mean - yval_mean * yval_mean);
        double stddev = sqrt( pow(xval_stddev, 2) + pow(yval_stddev, 2) ) / sqrt(2);
        stddev_arr[iline] = stddev;
        xval_mean_arr[iline] = xval_mean;
        yval_mean_arr[iline] = yval_mean;
        printf("iline = %ld, stddev = %e, (xval_mean, yval_mean) = (%e, %e)\n",
               iline, stddev, xval_mean, yval_mean);
        delete hd2d_img;
    }

    long iline_min = MirMath::GetLocMin(nline, stddev_arr);

    double vel = vel_arr[iline_min];
    double rho = rho_arr[iline_min];
    double phi = phi_arr[iline_min];
    double psi = psi_arr[iline_min];
    double theta = atan(vel);
    double xval_mean_min = xval_mean_arr[iline_min];
    double yval_mean_min = yval_mean_arr[iline_min];
    
    char outpar[kLineSize];
    sprintf(outpar, "%s/%s_par.dat",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str());
    FILE* fp_par = fopen(outpar, "w");
    fprintf(fp_par, "# xval yval   vel rho phi psi\n");
    fprintf(fp_par, "%e %e    %e %e %e %e\n",
            xval_mean_min, yval_mean_min,
            vel, rho, phi, psi);
    fclose(fp_par);
    
    HistDataNerr2d* hd2d_img = NULL;
    GenDetImg(hd2d_arr, time_arr, ntime,
              argval->GetNbinDetimgHalf(),
              theta, rho, phi, psi, &hd2d_img);
        
    char tag[kLineSize];
    sprintf(tag, "img_%3.3ld", iline_min);
    MifImgInfo* img_info_img = new MifImgInfo;
    img_info_img->InitSetImg(1, 1,
                             hd2d_img->GetHi2d()->GetNbinX(),
                             hd2d_img->GetHi2d()->GetNbinY());
    MifFits::OutFitsImageD(argval->GetOutdir(),
                           argval->GetOutfileHead(),
                           tag,
                           2, 
                           bitpix,
                           img_info_img->GetNaxesArr(),
                           hd2d_img->GetOvalArr()->GetVal());
    delete hd2d_img;
    
    
//    for(long iline = 0; iline < nline; iline ++){
//        double vel = vel_arr[iline];
//        double rho = rho_arr[iline];
//        double phi = phi_arr[iline];
//        double psi = psi_arr[iline];
//        double theta = atan(vel);
//
//        HistDataNerr2d* hd2d_img = NULL;
//        GenDetImg(hd2d_arr, time_arr, ntime,
//                  argval->GetNbinDetimgHalf(),
//                  theta, rho, phi, psi, &hd2d_img);
//        
//        char tag[kLineSize];
//        sprintf(tag, "img_%3.3ld", iline);
//        MifImgInfo* img_info_img = new MifImgInfo;
//        img_info_img->InitSetImg(1, 1,
//                                 hd2d_img->GetHi2d()->GetNbinX(),
//                                 hd2d_img->GetHi2d()->GetNbinY());
//        MifFits::OutFitsImageD(argval->GetOutdir(),
//                               argval->GetOutfileHead(),
//                               tag,
//                               2, 
//                               bitpix,
//                               img_info_img->GetNaxesArr(),
//                               hd2d_img->GetOvalArr()->GetVal());
//        delete hd2d_img;
//        
//    }
//    
    
    
    double time_ed = MiTime::GetTimeSec();
    printf("time_ed - time_st = %e\n", time_ed - time_st);

    fclose(fp_log);
    
    return status_prog;
}
