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

    char logfile[kLineSize];
    if( MiIolib::TestFileExist(argval->GetOutdir()) ){
        char cmd[kLineSize];
        sprintf(cmd, "mkdir -p %s", argval->GetOutdir().c_str());
        system(cmd);
    }
    sprintf(logfile, "%s/%s_%s.log",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str(),
            argval->GetProgname().c_str());
    FILE* fp_log = fopen(logfile, "w");
    MiIolib::Printf2(fp_log, "-----------------------------\n");
    argval->Print(fp_log);

    printf("--- read data list ---\n");
    string* line_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(argval->GetDataList(), &line_arr, &nline);
    printf("nline = %ld\n", nline);
    string* fitsfile_arr = new string[nline];
    double* time_arr     = new double[nline];
    for(long iline = 0; iline < nline; iline ++){
        int nsplit = 0;
        string* split_arr = NULL;
        MiStr::GenSplit(line_arr[iline], &nsplit, &split_arr);
        if(2 != nsplit){
            printf("Bad data_list(=%s): nsplit(%d) != 2 @ iline(%ld).\n",
                   argval->GetDataList().c_str(), nsplit, iline);
            abort();
        }
        fitsfile_arr[iline] = split_arr[0];
        time_arr[iline] = atof(split_arr[1].c_str());
        MiStr::DelSplit(split_arr);
    }
    delete [] line_arr;
    for(long iline = 0; iline < nline; iline ++){
        printf("%s  %e\n", fitsfile_arr[iline].c_str(), time_arr[iline]);
    }
    // check fits images
    double* npx_arr = new double[nline];
    double* npy_arr = new double[nline];
    for(long iline = 0; iline < nline; iline ++){
        int naxis = MifFits::GetNaxis(fitsfile_arr[iline]);
        if(2 != naxis){
            printf("fits file is not image, then abort.\n");
            abort();
        }
        npx_arr[iline] = MifFits::GetAxisSize(fitsfile_arr[iline], 0);
        npy_arr[iline] = MifFits::GetAxisSize(fitsfile_arr[iline], 1);
        printf("%ld: (npx, npy) = (%f, %f)\n", iline, npx_arr[iline], npy_arr[iline]);
    }
    double npx_amean = MirMath::GetAMean(nline, npx_arr);
    double npy_amean = MirMath::GetAMean(nline, npy_arr);
    double npx_stddev = MirMath::GetStddev(nline, npx_arr);
    double npy_stddev = MirMath::GetStddev(nline, npy_arr);
    if(npx_stddev > 1.0e-10 || npy_stddev > 1.0e-10){
        printf("bad image size, then abort\n");
        abort();
    }
    delete [] npx_arr;
    delete [] npy_arr;
    printf("=== read data list ===\n");

    printf("--- read 2d images ---\n");
    MifImgInfo* img_info_subimg = new MifImgInfo;
    if("none" == argval->GetSubimgDat()){
        img_info_subimg->InitSetImg(1, 1, npx_amean, npy_amean);
    } else {
        img_info_subimg->Load(argval->GetSubimgDat());
        img_info_subimg->PrintInfo();
    }
    long ntime = nline;
    double** data_arr = new double* [ntime];
    int bitpix = 0;
    for(int itime = 0; itime < ntime; itime ++){
        MifFits::InFitsImageD(fitsfile_arr[itime], img_info_subimg,
                              &bitpix, &data_arr[itime]);
    }
    printf("bitpix = %d\n", bitpix);
    printf("=== read 2d images ===\n");
   
    printf("--- standardize images ---\n");
    double** std_arr = new double* [ntime];
    double** debias_arr = new double* [ntime];
    for(long itime = 0; itime < ntime; itime ++){
        std_arr[itime] = new double[img_info_subimg->GetNpixelTotal()];
        debias_arr[itime] = new double[img_info_subimg->GetNpixelTotal()];
    }
    int nclip = 30;
    double significance_clip = 5.0;
    double* mean_arr = new double [ntime];
    double* stddev_arr = new double [ntime];
    for(long itime = 0; itime < ntime; itime ++){
        double mean = 0.0;
        double stddev = 0.0;
        MirMathUtil::GetMeanStddevClip(img_info_subimg->GetNpixelTotal(), data_arr[itime],
                                       nclip, significance_clip,
                                       &mean, &stddev);
        for(long iarr = 0; iarr < img_info_subimg->GetNpixelTotal(); iarr ++){
            std_arr[itime][iarr] = (data_arr[itime][iarr] - mean) / stddev;
            debias_arr[itime][iarr] = data_arr[itime][iarr] - mean;
        }
        mean_arr[itime] = mean;
        stddev_arr[itime] = stddev;
    }
    printf("=== standardize images ===\n");

    printf("--- make median image ---\n");
    double* median_arr = new double [img_info_subimg->GetNpixelImg()];
    for(long iarr = 0; iarr < img_info_subimg->GetNpixelImg(); iarr ++){
        median_arr[iarr] = 0.0;
    }
    double* time_tmp_arr = new double [ntime];
    for(long iarr = 0; iarr < img_info_subimg->GetNpixelImg(); iarr ++){
        for(long itime = 0; itime < ntime; itime ++){
            time_tmp_arr[itime] = debias_arr[itime][iarr];
        }
        double median = MirMath::GetMedian(ntime, time_tmp_arr);
        median_arr[iarr] = median;
    }
    delete [] time_tmp_arr;
    printf("=== make median image ===\n");

    printf("--- smooth image ---\n");
    // push image to 2d hist
    HistDataNerr2d* hd2d = new HistDataNerr2d;
    hd2d->Init(img_info_subimg->GetNaxesArrElm(0), 0, img_info_subimg->GetNaxesArrElm(0),
               img_info_subimg->GetNaxesArrElm(1), 0, img_info_subimg->GetNaxesArrElm(1));
    hd2d->SetOvalArr(img_info_subimg->GetNpixelImg(), median_arr);

    MirFunc* func = MifcGen::GenFunc(argval->GetFunc());
    MirFuncPar* func_par = new MirFuncPar;
    func_par->Load(argval->GetParFile());

    // kernel
    HistInfo2d* hi2d_kernel = new HistInfo2d;
    hi2d_kernel->InitSetByMidPoint(0.0, 1.0, argval->GetNbinKernelHalf(), "ceil",
                                   0.0, 1.0, argval->GetNbinKernelHalf(), "ceil");
    HistDataNerr2d* hd2d_kernel_tmp = new HistDataNerr2d;
    hd2d_kernel_tmp->Init(hi2d_kernel);
    hd2d_kernel_tmp->SetByFunc(func, func_par->GetPar());
    double sum = MirMath::GetSum(hd2d_kernel_tmp->GetNbin(), hd2d_kernel_tmp->GetOvalArr()->GetVal());
    HistDataNerr2d* hd2d_kernel = new HistDataNerr2d;
    HistData2dOpe::GetScale(hd2d_kernel_tmp, 1./sum, 0.0, hd2d_kernel);

    // smoothing
    HistDataNerr2d* hd2d_smooth = new HistDataNerr2d;
    hd2d_smooth->Init(img_info_subimg->GetNaxesArrElm(0), 0, img_info_subimg->GetNaxesArrElm(0),
                      img_info_subimg->GetNaxesArrElm(1), 0, img_info_subimg->GetNaxesArrElm(1));
    for(long ibin = 0; ibin < hd2d->GetNbin(); ibin ++){
        double oval = hd2d->GetOvalArr()->GetValElm(ibin);
        double xval = hd2d->GetHi2d()->GetBinCenterXFromIbin(ibin);
        double yval = hd2d->GetHi2d()->GetBinCenterYFromIbin(ibin);
        for(long ibin_kernel = 0; ibin_kernel < hd2d_kernel->GetNbin(); ibin_kernel++){
            double oval_kernel = hd2d_kernel->GetOvalArr()->GetValElm(ibin_kernel);
            double xval_kernel = hd2d_kernel->GetHi2d()->GetBinCenterXFromIbin(ibin_kernel);
            double yval_kernel = hd2d_kernel->GetHi2d()->GetBinCenterYFromIbin(ibin_kernel);
            double xval_add = xval + xval_kernel;
            double yval_add = yval + yval_kernel;
            if(hd2d_smooth->GetHi2d()->GetLoX() <= xval_add &&
               hd2d_smooth->GetHi2d()->GetUpX() >= xval_add &&
               hd2d_smooth->GetHi2d()->GetLoY() <= yval_add &&
               hd2d_smooth->GetHi2d()->GetUpY() >= yval_add){
                hd2d_smooth->Fill(xval_add, yval_add, oval * oval_kernel);
            }
        }
    }
    MifFits::OutFitsImageD(argval->GetOutdir(),
                           argval->GetOutfileHead(),
                           "smooth",
                           2, 
                           bitpix,
                           img_info_subimg->GetNaxesArr(),
                           hd2d_smooth->GetOvalArr()->GetVal());
    printf("--- smooth image ---\n");

    long index = MirMath::GetLocMax(hd2d_smooth->GetNbin(), hd2d_smooth->GetOvalArr()->GetVal());
    printf("%ld, %ld\n", hd2d_smooth->GetHi2d()->GetIbinX(index), hd2d_smooth->GetHi2d()->GetIbinY(index));
    
    
    HistInfo2d* hi2d_psf = new HistInfo2d;
    hi2d_psf->InitSetByMidPoint(0.0, 1.0, argval->GetNbinPsfHalf(), "ceil",
                                0.0, 1.0, argval->GetNbinPsfHalf(), "ceil");
    HistDataNerr2d* hd2d_psf = new HistDataNerr2d;
    hd2d_psf->Init(hi2d_psf);

    hd2d_psf->PrintInfo(stdout);

    double xval_max = hd2d->GetHi2d()->GetBinCenterXFromIbin(index);
    double yval_max = hd2d->GetHi2d()->GetBinCenterYFromIbin(index);
    for(long ibin_psf = 0; ibin_psf < hd2d_psf->GetNbin(); ibin_psf++){
        double xval_psf = hd2d_psf->GetHi2d()->GetBinCenterXFromIbin(ibin_psf);
        double yval_psf = hd2d_psf->GetHi2d()->GetBinCenterYFromIbin(ibin_psf);
        double xval_add = xval_max + xval_psf;
        double yval_add = yval_max + yval_psf;
        double oval_add = hd2d->GetOvalElmAtXY(xval_add, yval_add);
        hd2d_psf->Fill(xval_psf, yval_psf, oval_add);
    }
    // normalize
    double sum_psf = MirMath::GetSum(hd2d_psf->GetOvalArr()->GetNdata(),
                                     hd2d_psf->GetOvalArr()->GetVal());
    for(long ibin_psf = 0; ibin_psf < hd2d_psf->GetNbin(); ibin_psf++){
        long ibinx = hd2d_psf->GetHi2d()->GetIbinX(ibin_psf);
        long ibiny = hd2d_psf->GetHi2d()->GetIbinY(ibin_psf);
        double val = hd2d_psf->GetOvalArr()->GetValElm(ibin_psf);
        val /= sum_psf;
        hd2d_psf->SetOvalElm(ibinx, ibiny, val);
    }

    char file_psf[kLineSize];
    sprintf(file_psf, "%s/%s_psf.dat",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str()); 
    hd2d_psf->Save(file_psf, "x,y,z");


    MifImgInfo* img_info_psf = new MifImgInfo;
    img_info_psf->InitSetImg(1, 1, hi2d_psf->GetNbinX(), hi2d_psf->GetNbinY());
    
    MifFits::OutFitsImageD(argval->GetOutdir(),
                           argval->GetOutfileHead(),
                           "psf",
                           2, 
                           bitpix,
                           img_info_psf->GetNaxesArr(),
                           hd2d_psf->GetOvalArr()->GetVal());
    
    // ------- psf convolve
    double** conv_arr = new double* [ntime];
    for(long itime = 0; itime < ntime; itime ++){
        conv_arr[itime] = new double[img_info_subimg->GetNpixelTotal()];
    }
    for(long itime = 0; itime < ntime; itime ++){
        // push image to 2d hist array
        HistDataNerr2d* hd2d_tmp = new HistDataNerr2d;
        hd2d_tmp->Init(img_info_subimg->GetNaxesArrElm(0), 0, img_info_subimg->GetNaxesArrElm(0),
                       img_info_subimg->GetNaxesArrElm(1), 0, img_info_subimg->GetNaxesArrElm(1));
        hd2d_tmp->SetOvalArr(img_info_subimg->GetNpixelImg(), debias_arr[itime]);

        HistDataNerr2d* hd2d_conv = new HistDataNerr2d;
        hd2d_conv->Init(img_info_subimg->GetNaxesArrElm(0), 0, img_info_subimg->GetNaxesArrElm(0),
                        img_info_subimg->GetNaxesArrElm(1), 0, img_info_subimg->GetNaxesArrElm(1));
        for(long ibin = 0; ibin < hd2d_tmp->GetNbin(); ibin ++){
            long ibinx = hd2d_tmp->GetHi2d()->GetIbinX(ibin);
            long ibiny = hd2d_tmp->GetHi2d()->GetIbinY(ibin);
            double xval = hd2d_tmp->GetHi2d()->GetBinCenterXFromIbin(ibin);
            double yval = hd2d_tmp->GetHi2d()->GetBinCenterYFromIbin(ibin);
            double oval = 0.0;
            for(long ibin_psf = 0; ibin_psf < hd2d_psf->GetNbin(); ibin_psf++){
                double oval_psf = hd2d_psf->GetOvalArr()->GetValElm(ibin_psf);
                double xval_psf = hd2d_psf->GetHi2d()->GetBinCenterXFromIbin(ibin_psf);
                double yval_psf = hd2d_psf->GetHi2d()->GetBinCenterYFromIbin(ibin_psf);
                double xval_add = xval + xval_psf;
                double yval_add = yval + yval_psf;
                if(hd2d_conv->GetHi2d()->GetLoX() <= xval_add &&
                   hd2d_conv->GetHi2d()->GetUpX() >= xval_add &&
                   hd2d_conv->GetHi2d()->GetLoY() <= yval_add &&
                   hd2d_conv->GetHi2d()->GetUpY() >= yval_add){
                    double oval_hd2d = hd2d_tmp->GetOvalElmAtXY(xval_add, yval_add);
                    oval += oval_hd2d * oval_psf;
                }
            }
            hd2d_conv->SetOvalElm(ibinx, ibiny, oval);
        }

        for(long iarr = 0; iarr < img_info_subimg->GetNpixelTotal(); iarr ++){
            conv_arr[itime][iarr] = hd2d_conv->GetOvalArr()->GetValElm(iarr);
        }
        char tag[kLineSize];
        sprintf(tag, "conv_%2.2ld", itime);
        MifFits::OutFitsImageD(argval->GetOutdir(),
                               argval->GetOutfileHead(),
                               tag,
                               2, 
                               bitpix,
                               img_info_subimg->GetNaxesArr(),
                               conv_arr[itime]);
        delete hd2d_tmp;
        delete hd2d_conv;
    }



    
    
    double time_ed = MiTime::GetTimeSec();
    printf("time_ed - time_st = %e\n", time_ed - time_st);
    
    return status_prog;
}
