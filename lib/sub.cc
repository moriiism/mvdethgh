#include "mir_math_util.h"
#include "mir_data1d_ope.h"
#include "mir_hist2d_ope.h"
#include "mir_hist_info.h"
#include "mi_str.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "mifc_gen.h"
#include "sub.h"

void OpenLogfile(string outdir,
                 string outfile_head,
                 string progname,
                 FILE** fp_log_ptr)
{
    char logfile[kLineSize];
    if( MiIolib::TestFileExist(outdir) ){
        char cmd[kLineSize];
        sprintf(cmd, "mkdir -p %s", outdir.c_str());
        system(cmd);
    }
    sprintf(logfile, "%s/%s_%s.log",
            outdir.c_str(),
            outfile_head.c_str(),
            progname.c_str());
    FILE* fp_log = fopen(logfile, "w");
    MiIolib::Printf2(fp_log, "-----------------------------\n");
    *fp_log_ptr = fp_log;
}

void ReadDataList(string data_list,
                  string subimg_dat,
                  long* const ntime_ptr,
                  double** const time_arr_ptr,
                  double*** const data_arr_ptr,
                  MifImgInfo** const img_info_subimg_ptr,
                  int* const bitpix_ptr)
{
    string* line_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(data_list, &line_arr, &nline);
    printf("nline = %ld\n", nline);
    string* fitsfile_arr = new string[nline];
    double* time_arr     = new double[nline];
    for(long iline = 0; iline < nline; iline ++){
        int nsplit = 0;
        string* split_arr = NULL;
        MiStr::GenSplit(line_arr[iline], &nsplit, &split_arr);
        if(2 != nsplit){
            printf("Bad data_list(=%s): nsplit(%d) != 2 @ iline(%ld).\n",
                   data_list.c_str(), nsplit, iline);
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

    printf("--- read 2d images ---\n");
    MifImgInfo* img_info_subimg = new MifImgInfo;
    if("none" == subimg_dat){
        img_info_subimg->InitSetImg(1, 1, npx_amean, npy_amean);
    } else {
        img_info_subimg->Load(subimg_dat);
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

    *time_arr_ptr = time_arr;
    *data_arr_ptr = data_arr;
    *img_info_subimg_ptr = img_info_subimg;
    *ntime_ptr = ntime;
    *bitpix_ptr = bitpix;
}

void MkStdImgArr(long ntime,
                 const double* const* const data_img_arr,
                 const MifImgInfo* const img_info_subimg,
                 int nclip, double significance_clip,
                 double*** std_img_arr_ptr,
                 double*** debias_img_arr_ptr,
                 double** mean_time_arr_ptr,
                 double** stddev_time_arr_ptr)
{
    double** std_img_arr    = new double* [ntime];
    double** debias_img_arr = new double* [ntime];
    double* mean_time_arr   = new double [ntime];
    double* stddev_time_arr = new double [ntime];
    for(long itime = 0; itime < ntime; itime ++){
        std_img_arr[itime] = new double[img_info_subimg->GetNpixelTotal()];
        debias_img_arr[itime] = new double[img_info_subimg->GetNpixelTotal()];
    }
    for(long itime = 0; itime < ntime; itime ++){
        double mean = 0.0;
        double stddev = 0.0;
        MirMathUtil::GetMeanStddevClip(img_info_subimg->GetNpixelTotal(),
                                       data_img_arr[itime],
                                       nclip, significance_clip,
                                       &mean, &stddev);
        for(long iarr = 0; iarr < img_info_subimg->GetNpixelTotal(); iarr ++){
            std_img_arr[itime][iarr] = (data_img_arr[itime][iarr] - mean) / stddev;
            debias_img_arr[itime][iarr] = data_img_arr[itime][iarr] - mean;
        }
        mean_time_arr[itime] = mean;
        stddev_time_arr[itime] = stddev;
    }
    *std_img_arr_ptr = std_img_arr;
    *debias_img_arr_ptr = debias_img_arr;
    *mean_time_arr_ptr = mean_time_arr;
    *stddev_time_arr_ptr = stddev_time_arr;
}


void MkMedianImg(long ntime,
                 const double* const* const debias_img_arr,
                 const MifImgInfo* const img_info_subimg,
                 double** median_img_arr_ptr)
{
    double* median_img_arr = new double [img_info_subimg->GetNpixelImg()];
    for(long iarr = 0; iarr < img_info_subimg->GetNpixelImg(); iarr ++){
        median_img_arr[iarr] = 0.0;
    }
    double* time_tmp_arr = new double [ntime];
    for(long iarr = 0; iarr < img_info_subimg->GetNpixelImg(); iarr ++){
        for(long itime = 0; itime < ntime; itime ++){
            time_tmp_arr[itime] = debias_img_arr[itime][iarr];
        }
        double median = MirMath::GetMedian(ntime, time_tmp_arr);
        median_img_arr[iarr] = median;
    }
    delete [] time_tmp_arr;
    *median_img_arr_ptr = median_img_arr;
}


void GetPeakXY(const MifImgInfo* const img_info_subimg,
               const double* const median_img_arr,
               string func_name, string func_par_name, int nbin_kernel_half,
               double* xval_peak_ptr, double* yval_peak_ptr)
{    
    // push image to 2d hist
    HistDataNerr2d* hd2d = new HistDataNerr2d;
    hd2d->Init(img_info_subimg->GetNaxesArrElm(0), 0, img_info_subimg->GetNaxesArrElm(0),
               img_info_subimg->GetNaxesArrElm(1), 0, img_info_subimg->GetNaxesArrElm(1));
    hd2d->SetOvalArr(img_info_subimg->GetNpixelImg(), median_img_arr);

    MirFunc* func = MifcGen::GenFunc(func_name);
    MirFuncPar* func_par = new MirFuncPar;
    func_par->Load(func_par_name);

    // kernel
    HistInfo2d* hi2d_kernel = new HistInfo2d;
    hi2d_kernel->InitSetByMidPoint(0.0, 1.0, nbin_kernel_half, "ceil",
                                   0.0, 1.0, nbin_kernel_half, "ceil");
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

    long index = MirMath::GetLocMax(hd2d_smooth->GetNbin(), hd2d_smooth->GetOvalArr()->GetVal());
    double xval_peak = hd2d_smooth->GetHi2d()->GetBinCenterXFromIbin(index);
    double yval_peak = hd2d_smooth->GetHi2d()->GetBinCenterYFromIbin(index);

//    MifFits::OutFitsImageD(argval->GetOutdir(),
//                           argval->GetOutfileHead(),
//                           "smooth",
//                           2, 
//                           bitpix,
//                           img_info_subimg->GetNaxesArr(),
//                           hd2d_smooth->GetOvalArr()->GetVal());
    
    delete hd2d;
    delete func;
    delete func_par;
    delete hi2d_kernel;
    delete hd2d_kernel_tmp;
    delete hd2d_kernel;
    delete hd2d_smooth;
    
    *xval_peak_ptr = xval_peak;
    *yval_peak_ptr = yval_peak;
}


void MkPsf(const MifImgInfo* const img_info_subimg,
           const double* const median_img_arr,
           int nbin_psf_half,
           double xval_peak, double yval_peak,
           HistDataNerr2d** hd2d_psf_ptr)
{
    HistInfo2d* hi2d_psf = new HistInfo2d;
    hi2d_psf->InitSetByMidPoint(0.0, 1.0, nbin_psf_half, "ceil",
                                0.0, 1.0, nbin_psf_half, "ceil");
    HistDataNerr2d* hd2d_psf = new HistDataNerr2d;
    hd2d_psf->Init(hi2d_psf);
    hd2d_psf->PrintInfo(stdout);

    // push image to 2d hist
    HistDataNerr2d* hd2d = new HistDataNerr2d;
    hd2d->Init(img_info_subimg->GetNaxesArrElm(0), 0, img_info_subimg->GetNaxesArrElm(0),
               img_info_subimg->GetNaxesArrElm(1), 0, img_info_subimg->GetNaxesArrElm(1));
    hd2d->SetOvalArr(img_info_subimg->GetNpixelImg(), median_img_arr);

    for(long ibin_psf = 0; ibin_psf < hd2d_psf->GetNbin(); ibin_psf++){
        double xval_psf = hd2d_psf->GetHi2d()->GetBinCenterXFromIbin(ibin_psf);
        double yval_psf = hd2d_psf->GetHi2d()->GetBinCenterYFromIbin(ibin_psf);
        double xval_add = xval_peak + xval_psf;
        double yval_add = yval_peak + yval_psf;
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

    delete hi2d_psf;
    delete hd2d;

    *hd2d_psf_ptr = hd2d_psf;
}
