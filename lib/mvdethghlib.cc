#include "mvdethghlib.h"

void OpenLogfile(string outdir,
                 string outfile_head,
                 string progname,
                 FILE** const fp_log_ptr)
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

void LoadData(string data_list,
              string subimg_dat,
              long* const ntime_ptr,
              double** const time_arr_ptr,
              double*** const data_arr_ptr,
              MifImgInfo** const img_info_ptr,
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
    MifImgInfo* img_info = new MifImgInfo;
    if("none" == subimg_dat){
        img_info->InitSetImg(1, 1, npx_amean, npy_amean);
    } else {
        img_info->Load(subimg_dat);
    }
    img_info->PrintInfo();

    long ntime = nline;
    double** data_arr = new double* [ntime];
    int bitpix = 0;
    for(int itime = 0; itime < ntime; itime ++){
        MifFits::InFitsImageD(fitsfile_arr[itime], img_info,
                              &bitpix, &data_arr[itime]);
    }
    printf("bitpix = %d\n", bitpix);
    printf("=== read 2d images ===\n");

    *time_arr_ptr = time_arr;
    *data_arr_ptr = data_arr;
    *img_info_ptr = img_info;
    *ntime_ptr = ntime;
    *bitpix_ptr = bitpix;
}


void LoadDataXYT(string data_list,
                 long* const ntime_ptr,
                 double** const time_arr_ptr,
                 long** const npos_arr_ptr,
                 double*** const xpos_arr_ptr,
                 double*** const ypos_arr_ptr,
                 double* const xpos_lo_ptr,
                 double* const xpos_up_ptr,
                 double* const ypos_lo_ptr,
                 double* const ypos_up_ptr,
                 double* const rho_up_ptr)
{
    string* line_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(data_list, &line_arr, &nline);
    printf("nline = %ld\n", nline);
    string* xytfile_arr = new string[nline];
    for(long iline = 0; iline < nline; iline ++){
        int nsplit = 0;
        string* split_arr = NULL;
        MiStr::GenSplit(line_arr[iline], &nsplit, &split_arr);
        if(1 != nsplit){
            printf("Bad data_list(=%s): nsplit(%d) != 1 @ iline(%ld).\n",
                   data_list.c_str(), nsplit, iline);
            abort();
        }
        xytfile_arr[iline] = split_arr[0];
        MiStr::DelSplit(split_arr);
    }
    delete [] line_arr;
    for(long iline = 0; iline < nline; iline ++){
        printf("%s\n", xytfile_arr[iline].c_str());
    }

    long ntime = nline;
    double* mjd_arr = new double [ntime];
    for(int itime = 0; itime < ntime; itime ++){
        string buf;
        ifstream ifs(xytfile_arr[itime].c_str());
        getline(ifs, buf);
        ifs.close();
        int nsplit = 0;
        string* split_arr = NULL;
        MiStr::GenSplit(buf, &nsplit, &split_arr);
        mjd_arr[itime] = atof(split_arr[0].c_str());
        MiStr::DelSplit(split_arr);
    }

    double mjd_min = MirMath::GetMin(ntime, mjd_arr);
    double mjd_max = MirMath::GetMax(ntime, mjd_arr);
    double mjd_mid = (mjd_max + mjd_min) / 2.0;
    double* time_arr = new double[ntime];
    for(long itime = 0; itime < ntime; itime ++){
        time_arr[itime] = mjd_arr[itime] - mjd_mid;
    }
    double time_min = MirMath::GetMin(ntime, time_arr);
    double time_max = MirMath::GetMax(ntime, time_arr);

    double xpos_min = 1.0e10;
    double xpos_max = -1.0e10;
    double ypos_min = 1.0e10;
    double ypos_max = -1.0e10;
    double rho_min  = 1.0e10;
    double rho_max  = -1.0e10;
    long* npos_arr    = new long[ntime];
    double** xpos_arr = new double* [ntime];
    double** ypos_arr = new double* [ntime];
    double** rho_arr  = new double* [ntime];
    for(int itime = 0; itime < ntime; itime ++){
        string* lines_arr = NULL;
        long nline = 0;
        MiIolib::GenReadFileSkipComment(xytfile_arr[itime],
                                        &lines_arr, &nline);
        npos_arr[itime] = nline;
        xpos_arr[itime] = new double[nline];
        ypos_arr[itime] = new double[nline];
        rho_arr[itime]  = new double[nline];
        for(long iline = 0; iline < nline; iline ++){
            int nsplit = 0;
            string* split_arr = NULL;
            MiStr::GenSplit(lines_arr[iline], &nsplit, &split_arr);
            if(3 != nsplit){
                abort();
            }
            xpos_arr[itime][iline] = atof(split_arr[1].c_str());
            ypos_arr[itime][iline] = atof(split_arr[2].c_str());
            MiStr::DelSplit(split_arr);
            rho_arr[itime][iline] = sqrt(pow(xpos_arr[itime][iline], 2)
                                  + pow(ypos_arr[itime][iline], 2));
            if(xpos_arr[itime][iline] < xpos_min){
                xpos_min = xpos_arr[itime][iline];
            }
            if(xpos_arr[itime][iline] > xpos_max){
                xpos_max = xpos_arr[itime][iline];
            }
            if(ypos_arr[itime][iline] < ypos_min){
                ypos_min = ypos_arr[itime][iline];
            }
            if(ypos_arr[itime][iline] > ypos_max){
                ypos_max = ypos_arr[itime][iline];
            }
            if(rho_arr[itime][iline] < rho_min){
                rho_min = rho_arr[itime][iline];
            }
            if(rho_arr[itime][iline] > rho_max){
                rho_max = rho_arr[itime][iline];
            }
        }
    }

    printf("(mjd_min, mjd_max)   = (%e, %e)\n", mjd_min, mjd_max);
    printf("(time_min, time_max) = (%e, %e)\n", time_min, time_max);
    printf("(xpos_min, xpos_max) = (%e, %e)\n", xpos_min, xpos_max);
    printf("(ypos_min, ypos_max) = (%e, %e)\n", ypos_min, ypos_max);
    printf("(rho_min, rho_max)   = (%e, %e)\n", rho_min, rho_max);

    double time_lo = 0.0;
    double time_up = 0.0;
    double xpos_lo = 0.0;
    double xpos_up = 0.0;
    double ypos_lo = 0.0;
    double ypos_up = 0.0;
    double rho_lo  = 0.0;
    double rho_up  = 0.0;
    double margin_ratio = 0.01;
    MirMathUtil::GetMarginLoUp(time_min, time_max, margin_ratio,
                               &time_lo, &time_up);
    MirMathUtil::GetMarginLoUp(xpos_min, xpos_max, margin_ratio,
                               &xpos_lo, &xpos_up);
    MirMathUtil::GetMarginLoUp(ypos_min, ypos_max, margin_ratio,
                               &ypos_lo, &ypos_up);
    MirMathUtil::GetMarginLoUp(rho_min, rho_max, margin_ratio,
                               &rho_lo, &rho_up);

    *ntime_ptr = ntime;
    *time_arr_ptr = time_arr;
    *npos_arr_ptr = npos_arr;
    *xpos_arr_ptr = xpos_arr;
    *ypos_arr_ptr = ypos_arr;
    *xpos_lo_ptr = xpos_lo;
    *xpos_up_ptr = xpos_up;
    *ypos_lo_ptr = ypos_lo;
    *ypos_up_ptr = ypos_up;
    *rho_up_ptr  = rho_up;
}



void GenStdImgArr(long ntime,
                  const double* const* const data_img_arr,
                  const MifImgInfo* const img_info,
                  int nclip, double significance_clip,
                  double*** const std_img_arr_ptr,
                  double*** const debias_img_arr_ptr,
                  double** const mean_time_arr_ptr,
                  double** const stddev_time_arr_ptr)
{
    double** std_img_arr    = new double* [ntime];
    double** debias_img_arr = new double* [ntime];
    double* mean_time_arr   = new double [ntime];
    double* stddev_time_arr = new double [ntime];
    for(long itime = 0; itime < ntime; itime ++){
        std_img_arr[itime] = new double[img_info->GetNpixelTotal()];
        debias_img_arr[itime] = new double[img_info->GetNpixelTotal()];
    }
    for(long itime = 0; itime < ntime; itime ++){
        double mean = 0.0;
        double stddev = 0.0;
        MirMathUtil::GetMeanStddevClip(img_info->GetNpixelTotal(),
                                       data_img_arr[itime],
                                       nclip, significance_clip,
                                       &mean, &stddev);
        for(long iarr = 0; iarr < img_info->GetNpixelTotal(); iarr ++){
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


void GenMedianImg(long ntime,
                  const double* const* const debias_img_arr,
                  const MifImgInfo* const img_info,
                  double** const median_img_arr_ptr)
{
    double* median_img_arr = new double [img_info->GetNpixelImg()];
    for(long iarr = 0; iarr < img_info->GetNpixelImg(); iarr ++){
        median_img_arr[iarr] = 0.0;
    }
    double* time_tmp_arr = new double [ntime];
    for(long iarr = 0; iarr < img_info->GetNpixelImg(); iarr ++){
        for(long itime = 0; itime < ntime; itime ++){
            time_tmp_arr[itime] = debias_img_arr[itime][iarr];
        }
        double median = MirMath::GetMedian(ntime, time_tmp_arr);
        median_img_arr[iarr] = median;
    }
    delete [] time_tmp_arr;
    *median_img_arr_ptr = median_img_arr;
}

void GenMedianImgByNthElement(long ntime,
                              const double* const* const debias_img_arr,
                              const MifImgInfo* const img_info,
                              double** const median_img_arr_ptr)
{
    double* median_img_arr = new double [img_info->GetNpixelImg()];
    for(long iarr = 0; iarr < img_info->GetNpixelImg(); iarr ++){
        median_img_arr[iarr] = 0.0;
    }
    double* time_tmp_arr = new double [ntime];
    for(long iarr = 0; iarr < img_info->GetNpixelImg(); iarr ++){
        for(long itime = 0; itime < ntime; itime ++){
            time_tmp_arr[itime] = debias_img_arr[itime][iarr];
        }
        double median = GetMedianByNthElement(ntime, time_tmp_arr);
        median_img_arr[iarr] = median;
    }
    delete [] time_tmp_arr;
    *median_img_arr_ptr = median_img_arr;
}

double GetMedianByNthElement(long narr,
                             const double* const val_arr)
{
    long half = (long) (narr / 2.0);
    long* index_arr = new long [narr];
    MiSort::KthElement<double, long>(half, narr, val_arr, index_arr, 1);
    double median = val_arr[index_arr[half]];
    return median;
}

void GenImgAboveThArr(long ntime,
                      const double* const* const data_img_arr,
                      const double* const* const std_img_arr,
                      const MifImgInfo* const img_info,
                      double threshold,
                      double*** const out_img_arr_ptr,
                      long* const ndet_ptr)
{
    double** out_img_arr    = new double* [ntime];
    for(long itime = 0; itime < ntime; itime ++){
        out_img_arr[itime] = new double[img_info->GetNpixelTotal()];
    }
    long ndet = 0;
    for(long itime = 0; itime < ntime; itime ++){
        for(long iarr = 0; iarr < img_info->GetNpixelTotal(); iarr ++){
            if(std_img_arr[itime][iarr] >= threshold){
                out_img_arr[itime][iarr] = data_img_arr[itime][iarr];
                ndet ++;
            } else{
                out_img_arr[itime][iarr] = 0.0;
            }
        }
    }
    *ndet_ptr = ndet;
    *out_img_arr_ptr = out_img_arr;
}


void GetPeakXY(const MifImgInfo* const img_info,
               const double* const median_img_arr,
               int nbin_kernel_half, double val_smooth,
               double* const xval_peak_ptr, double* const yval_peak_ptr)
{    
    // push image to 2d hist
    HistDataNerr2d* hd2d = new HistDataNerr2d;
    hd2d->Init(img_info->GetNaxesArrElm(0), 0, img_info->GetNaxesArrElm(0),
               img_info->GetNaxesArrElm(1), 0, img_info->GetNaxesArrElm(1));
    hd2d->SetOvalArr(img_info->GetNpixelImg(), median_img_arr);

    MirFunc* func = MifcGen::GenFunc("Gauss2dFunc");
    MirFuncPar* func_par = new MirFuncPar;
    func_par->Init(7);
    func_par->SetElm(0, "sigma_xp", val_smooth);
    func_par->SetElm(1, "sigma_yp", val_smooth);
    func_par->SetElm(2, "norm", 1.0);
    func_par->SetElm(3, "rot_angle", 0.0);
    func_par->SetElm(4, "mu_xp", 0.0);
    func_par->SetElm(5, "mu_yp", 0.0);
    func_par->SetElm(6, "shift_z", 0.0);

    // kernel
    HistInfo2d* hi2d_kernel = new HistInfo2d;
    hi2d_kernel->InitSetByMidPoint(0.0, 1.0, nbin_kernel_half, "ceil",
                                   0.0, 1.0, nbin_kernel_half, "ceil");
    HistDataNerr2d* hd2d_kernel_tmp = new HistDataNerr2d;
    hd2d_kernel_tmp->Init(hi2d_kernel);
    hd2d_kernel_tmp->SetByFunc(func, func_par->GetPar());
    double sum = MirMath::GetSum(hd2d_kernel_tmp->GetNbin(),
                                 hd2d_kernel_tmp->GetOvalArr()->GetVal());
    HistDataNerr2d* hd2d_kernel = new HistDataNerr2d;
    HistData2dOpe::GetScale(hd2d_kernel_tmp, 1./sum, 0.0, hd2d_kernel);

    // smoothing
    HistDataNerr2d* hd2d_smooth = new HistDataNerr2d;
    hd2d_smooth->Init(img_info->GetNaxesArrElm(0), 0, img_info->GetNaxesArrElm(0),
                      img_info->GetNaxesArrElm(1), 0, img_info->GetNaxesArrElm(1));
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

    long index = MirMath::GetLocMax(hd2d_smooth->GetNbin(),
                                    hd2d_smooth->GetOvalArr()->GetVal());
    double xval_peak = hd2d_smooth->GetHi2d()->GetBinCenterXFromIbin(index);
    double yval_peak = hd2d_smooth->GetHi2d()->GetBinCenterYFromIbin(index);

//    MifFits::OutFitsImageD(argval->GetOutdir(),
//                           argval->GetOutfileHead(),
//                           "smooth",
//                           2, 
//                           bitpix,
//                           img_info->GetNaxesArr(),
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


void GenPsf(const MifImgInfo* const img_info,
            const double* const median_img_arr,
            int nbin_psf_half,
            double xval_peak, double yval_peak,
            HistDataNerr2d** const hd2d_psf_ptr)
{
    HistInfo2d* hi2d_psf = new HistInfo2d;
    hi2d_psf->InitSetByMidPoint(0.0, 1.0, nbin_psf_half, "ceil",
                                0.0, 1.0, nbin_psf_half, "ceil");
    HistDataNerr2d* hd2d_psf = new HistDataNerr2d;
    hd2d_psf->Init(hi2d_psf);
    hd2d_psf->PrintInfo(stdout);

    // push image to 2d hist
    HistDataNerr2d* hd2d = new HistDataNerr2d;
    hd2d->Init(img_info->GetNaxesArrElm(0), 0, img_info->GetNaxesArrElm(0),
               img_info->GetNaxesArrElm(1), 0, img_info->GetNaxesArrElm(1));
    hd2d->SetOvalArr(img_info->GetNpixelImg(), median_img_arr);

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

void SavePsf(const HistDataNerr2d* const hd2d_psf,
             int bitpix,
             string outdir, string outfile_head)
{
    char file_psf[kLineSize];
    sprintf(file_psf, "%s/%s_psf.dat",
            outdir.c_str(),
            outfile_head.c_str()); 
    hd2d_psf->Save(file_psf, "x,y,z");

    MifImgInfo* img_info_psf = new MifImgInfo;
    img_info_psf->InitSetImg(1, 1,
                             hd2d_psf->GetHi2d()->GetNbinX(),
                             hd2d_psf->GetHi2d()->GetNbinY());
    MifFits::OutFitsImageD(outdir,
                           outfile_head,
                           "psf",
                           2, 
                           bitpix,
                           img_info_psf->GetNaxesArr(),
                           hd2d_psf->GetOvalArr()->GetVal());
    delete img_info_psf;
}


void SaveCube(long ntime, const MifImgInfo* const img_info,
              const double* const* const data_img_arr,
              int bitpix, string outdir, string outfile_head, string tag)
{
    MifImgInfo* img_info_cube = new MifImgInfo;
    img_info_cube->InitSetCube(1, 1, 1,
                               img_info->GetNaxesArrElm(0),
                               img_info->GetNaxesArrElm(1),
                               ntime);
    img_info_cube->PrintInfo();
    double* cube_arr = new double [img_info_cube->GetNpixelTotal()];
    for(long itime = 0; itime < ntime; itime ++){
        for(long iarr = 0; iarr < img_info_cube->GetNpixelImg(); iarr ++){
            long index = iarr + itime * img_info_cube->GetNpixelImg();
            cube_arr[index] = data_img_arr[itime][iarr];
        }
    }
    MifFits::OutFitsCubeD(outdir, outfile_head, tag,
                          3, bitpix,
                          img_info_cube->GetNaxesArr(),
                          cube_arr);
    delete [] cube_arr;
    delete img_info_cube;
}


void GenMvobjImgArr(long ntime,
                    const MifImgInfo* const img_info,
                    const double* const* const debias_img_arr,
                    const double* const median_img_arr,
                    double*** const mvobj_img_arr_ptr)
{
    printf("--- make mvobj ---\n");
    double** mvobj_img_arr = new double* [ntime];
    for(int itime = 0; itime < ntime; itime ++){
        mvobj_img_arr[itime] = new double[img_info->GetNpixelImg()];
        for(long iarr = 0; iarr < img_info->GetNpixelImg(); iarr ++){
            mvobj_img_arr[itime][iarr] = debias_img_arr[itime][iarr] - median_img_arr[iarr];
        }
    }
    printf("=== make mvobj ===\n");

    *mvobj_img_arr_ptr = mvobj_img_arr;
}

void SaveImgArr(long ntime, const MifImgInfo* const img_info,
                const double* const* const data_img_arr,
                int bitpix, string outdir, string outfile_head, string tag)
{
    for(long itime = 0; itime < ntime; itime ++){
        char tag_this[kLineSize];
        sprintf(tag_this, "%s_%2.2ld", tag.c_str(), itime);
        MifFits::OutFitsImageD(outdir, outfile_head, tag_this,
                               2, bitpix,
                               img_info->GetNaxesArr(),
                               data_img_arr[itime]);
    }
}

void GenConvPsf(long ntime,
                const MifImgInfo* const img_info,
                const double* const* const mvobj_img_arr,
                string psf_dat,
                double*** const conv_img_arr_ptr)
{
    // load psf
    HistDataNerr2d* hd2d_psf = new HistDataNerr2d;
    hd2d_psf->Load(psf_dat);
    
    // conv by psf
    double** conv_img_arr = new double* [ntime];
    for(long itime = 0; itime < ntime; itime ++){
        conv_img_arr[itime] = new double[img_info->GetNpixelTotal()];
    }
    for(long itime = 0; itime < ntime; itime ++){
        // push image to 2d hist array
        HistDataNerr2d* hd2d_tmp = new HistDataNerr2d;
        hd2d_tmp->Init(img_info->GetNaxesArrElm(0), 0, img_info->GetNaxesArrElm(0),
                       img_info->GetNaxesArrElm(1), 0, img_info->GetNaxesArrElm(1));
        hd2d_tmp->SetOvalArr(img_info->GetNpixelImg(), mvobj_img_arr[itime]);
        
        HistDataNerr2d* hd2d_conv = new HistDataNerr2d;
        hd2d_conv->Init(img_info->GetNaxesArrElm(0), 0, img_info->GetNaxesArrElm(0),
                        img_info->GetNaxesArrElm(1), 0, img_info->GetNaxesArrElm(1));
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
        for(long iarr = 0; iarr < img_info->GetNpixelTotal(); iarr ++){
            conv_img_arr[itime][iarr] = hd2d_conv->GetOvalArr()->GetValElm(iarr);
        }
        delete hd2d_tmp;
        delete hd2d_conv;
    }
    delete hd2d_psf;
    *conv_img_arr_ptr = conv_img_arr;
}


void LoadHi1dVel(string vel_dat,
                 HistInfo1d* const hi1d_vel)
{
    printf("--- read vel_dat ---\n");
    string* line_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(vel_dat, &line_arr, &nline);
    if(1 != nline){
        printf("nline(%ld) != 1, then abort.", nline);
        abort();
    }
    long nbin_vel = 0;
    double vel_lo = 0;
    double vel_up = 0;
    sscanf(line_arr[0].c_str(), "%ld  %lf  %lf", &nbin_vel, &vel_lo, &vel_up);
    delete [] line_arr;
    printf("=== read vel_dat ===\n");
    
    hi1d_vel->InitSetByNbin(vel_lo, vel_up, nbin_vel);
}


void LoadHi1dPar(string res_dat,
                 const MifImgInfo* const img_info,
                 HistInfo1d* const hi1d_rho,
                 HistInfo1d* const hi1d_phi,
                 HistInfo1d* const hi1d_psi)
{
    printf("--- read res_dat ---\n");
    string* line_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(res_dat, &line_arr, &nline);
    if(1 != nline){
        printf("nline(%ld) != 1, then abort.", nline);
        abort();
    }
    long nbin_rho = 0;
    long nbin_phi = 0;
    long nbin_psi = 0;
    sscanf(line_arr[0].c_str(), "%ld  %ld  %ld",
           &nbin_rho, &nbin_phi, &nbin_psi);
    delete [] line_arr;
    printf("=== read res_dat ===\n");    
    
    double rho_lo = 0.0;
    double rho_up = sqrt( pow(img_info->GetNaxesArrElm(0), 2)
                          + pow(img_info->GetNaxesArrElm(1), 2) );
    double phi_lo = 0.0;
    double phi_up = 2 * M_PI;
    double psi_lo = 0.0;
    double psi_up = 2 * M_PI;

    hi1d_rho->InitSetByNbin(rho_lo, rho_up, nbin_rho);
    hi1d_phi->InitSetByNbin(phi_lo, phi_up, nbin_phi);
    hi1d_psi->InitSetByNbin(psi_lo, psi_up, nbin_psi);
}


void LoadHi1dPar(string res_dat,
                 double rho_up,
                 HistInfo1d* const hi1d_rho,
                 HistInfo1d* const hi1d_phi,
                 HistInfo1d* const hi1d_psi)
{
    printf("--- read res_dat ---\n");
    string* line_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(res_dat, &line_arr, &nline);
    if(1 != nline){
        printf("nline(%ld) != 1, then abort.", nline);
        abort();
    }
    long nbin_rho = 0;
    long nbin_phi = 0;
    long nbin_psi = 0;
    sscanf(line_arr[0].c_str(), "%ld  %ld  %ld",
           &nbin_rho, &nbin_phi, &nbin_psi);
    delete [] line_arr;
    printf("=== read res_dat ===\n");    
    
    double rho_lo = 0.0;
    double phi_lo = 0.0;
    double phi_up = 2 * M_PI;
    double psi_lo = 0.0;
    double psi_up = 2 * M_PI;

    hi1d_rho->InitSetByNbin(rho_lo, rho_up, nbin_rho);
    hi1d_phi->InitSetByNbin(phi_lo, phi_up, nbin_phi);
    hi1d_psi->InitSetByNbin(psi_lo, psi_up, nbin_psi);
}


void LoadHi1dTime(string time_dat,
                  HistInfo1d* const hi1d_time)
{
    printf("--- read time_dat ---\n");
    string* line_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(time_dat, &line_arr, &nline);
    if(1 != nline){
        printf("nline(%ld) != 1, then abort.", nline);
        abort();
    }
    long nbin_time = 0;
    double time_lo = 0;
    double time_up = 0;
    sscanf(line_arr[0].c_str(), "%ld  %lf  %lf", &nbin_time, &time_lo, &time_up);
    delete [] line_arr;    
    printf("=== read time_dat ===\n");

    hi1d_time->InitSetByNbin(time_lo, time_up, nbin_time);
}


void GenDetImg(const HistDataNerr2d* const* const hd2d_arr,
               const double* const time_arr, long ntime,
               int nbin_detimg_half,
               double theta, double rho, double phi, double psi,
               HistDataNerr2d** const hd2d_img_ptr)
{
    long itime0 = 0;
    double zval0  = time_arr[itime0];
    double xval0 = rho * (cos(phi) * cos(psi) - sin(phi) * sin(psi) / cos(theta) )
        + zval0 * tan(theta) * sin(psi);
    double yval0 = rho * (cos(phi) * sin(psi) + sin(phi) * cos(psi) / cos(theta) )
        - zval0 * tan(theta) * cos(psi);

    HistInfo2d* hi2d_img = new HistInfo2d;
    hi2d_img->InitSetByMidPoint(xval0, 1.0, nbin_detimg_half, "ceil",
                                yval0, 1.0, nbin_detimg_half, "ceil");
    HistDataNerr2d* hd2d_img = new HistDataNerr2d;
    hd2d_img->Init(hi2d_img);
        
    for(long itime = 0; itime < ntime; itime ++){
        double zval = time_arr[itime];
        double xval = rho * (cos(phi) * cos(psi) - sin(phi) * sin(psi) / cos(theta) )
            + zval * tan(theta) * sin(psi);
        double yval = rho * (cos(phi) * sin(psi) + sin(phi) * cos(psi) / cos(theta) )
            - zval * tan(theta) * cos(psi);
        double xval_shift = xval - xval0;
        double yval_shift = yval - yval0;
        for(long ibin = 0; ibin < hd2d_img->GetNbin(); ibin ++){
            double xval_bin = hi2d_img->GetBinCenterXFromIbin(ibin);
            double yval_bin = hi2d_img->GetBinCenterYFromIbin(ibin);
            double xval_bin_new = xval_bin + xval_shift;
            double yval_bin_new = yval_bin + yval_shift;
            double oval = 0.0;
            if(xval_bin_new < hd2d_arr[itime]->GetXvalLo() ||
               hd2d_arr[itime]->GetXvalUp() < xval_bin_new ||
               yval_bin_new < hd2d_arr[itime]->GetYvalLo() ||
               hd2d_arr[itime]->GetYvalUp() < yval_bin_new  ){
                oval = 0.0;
            } else {
                oval = hd2d_arr[itime]->GetOvalElmAtXY(xval_bin_new, yval_bin_new);
            }
            hd2d_img->Fill(xval_bin, yval_bin, oval);
        }
    }
    delete hi2d_img;
    *hd2d_img_ptr = hd2d_img;
}

void GenDetImg(const HistDataNerr2d* const* const hd2d_arr,
               const double* const time_arr, long ntime,
               double vel, double psi,
               HistDataNerr2d** const hd2d_img_ptr)
{
    long itime0 = 0;
    double zval0  = time_arr[itime0];
    HistDataNerr2d* hd2d_img = new HistDataNerr2d;
    hd2d_img->Init(hd2d_arr[itime0]->GetHi2d());
    for(long itime = 0; itime < ntime; itime ++){
        double zval = time_arr[itime];
        double delta_zval = zval - zval0;
        double xval_shift = vel * delta_zval * sin(psi);
        double yval_shift = - vel * delta_zval * cos(psi);
        for(long ibin = 0; ibin < hd2d_img->GetNbin(); ibin ++){
            double xval_bin = hd2d_img->GetHi2d()->GetBinCenterXFromIbin(ibin);
            double yval_bin = hd2d_img->GetHi2d()->GetBinCenterYFromIbin(ibin);
            double xval_bin_new = xval_bin + xval_shift;
            double yval_bin_new = yval_bin + yval_shift;
            double oval = 0.0;
            if(xval_bin_new < hd2d_arr[itime]->GetXvalLo() ||
               hd2d_arr[itime]->GetXvalUp() < xval_bin_new ||
               yval_bin_new < hd2d_arr[itime]->GetYvalLo() ||
               hd2d_arr[itime]->GetYvalUp() < yval_bin_new  ){
                oval = 0.0;
            } else {
                oval = hd2d_arr[itime]->GetOvalElmAtXY(xval_bin_new, yval_bin_new);
            }
            hd2d_img->Fill(xval_bin, yval_bin, oval);
        }
    }
    *hd2d_img_ptr = hd2d_img;
}

