#include "mir_math.h"
#include "mir_data1d_ope.h"
#include "mir_hist_info.h"
#include "mi_str.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_mvdethgh.h"
#include "sub_mvdethgh.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValMvdethgh* argval = new ArgValMvdethgh;
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

    printf("--- output subcube ---\n");
    MifImgInfo* img_info_rec = new MifImgInfo;
    img_info_rec->InitSetCube(1, 1, 1, img_info_subimg->GetNaxesArrElm(0),
                              img_info_subimg->GetNaxesArrElm(1), ntime);
    img_info_rec->PrintInfo();
    double* subcube_arr = new double [img_info_rec->GetNpixelTotal()];
    for(long itime = 0; itime < ntime; itime ++){
        for(long iarr = 0; iarr < img_info_rec->GetNpixelImg(); iarr ++){
            long index = iarr + itime * img_info_rec->GetNpixelImg();
            subcube_arr[index] = data_arr[itime][iarr];
        }
    }
    MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), "subcube",
                          3, bitpix,
                          img_info_rec->GetNaxesArr(),
                          subcube_arr);
    delete [] subcube_arr;
    printf("=== output subcube ===\n");

    printf("--- load parameters ---\n");
    HistInfo1d* hi1d_vel   = new HistInfo1d;
    HistInfo1d* hi1d_rho   = new HistInfo1d;
    HistInfo1d* hi1d_phi   = new HistInfo1d;
    HistInfo1d* hi1d_psi   = new HistInfo1d;
    LoadHi1dVel(argval->GetVelDat(), hi1d_vel);
    LoadHi1dPar(argval->GetResDat(), img_info_subimg,
                hi1d_rho, hi1d_phi, hi1d_psi);
    long nbin_vel   = hi1d_vel->GetNbin();
    long nbin_rho   = hi1d_rho->GetNbin();
    long nbin_phi   = hi1d_phi->GetNbin();
    long nbin_psi   = hi1d_psi->GetNbin();
    hi1d_vel->Print(stdout);
    hi1d_rho->Print(stdout);
    hi1d_phi->Print(stdout);
    hi1d_psi->Print(stdout);
    printf("=== load parameters ===\n");
    

    printf("--- standardize images ---\n");
    double** std_arr = new double* [ntime];
    for(long itime = 0; itime < ntime; itime ++){
        std_arr[itime] = new double[img_info_subimg->GetNpixelTotal()];
    }
    int nclip = 30;
    double significance_clip = 5.0;
    for(long itime = 0; itime < ntime; itime ++){
        double mean = 0.0;
        double stddev = 0.0;
        GetMeanStddevClip(img_info_subimg->GetNpixelTotal(), data_arr[itime],
                          nclip, significance_clip,
                          &mean, &stddev);
        for(long iarr = 0; iarr < img_info_subimg->GetNpixelTotal(); iarr ++){
            std_arr[itime][iarr] = (data_arr[itime][iarr] - mean) / stddev;
        }
    }
    printf("=== standardize images ===\n");

    printf("--- output standardized subcube ---\n");
    double* subcube_std_arr = new double [img_info_rec->GetNpixelTotal()];
    for(long itime = 0; itime < ntime; itime ++){
        for(long iarr = 0; iarr < img_info_rec->GetNpixelImg(); iarr ++){
            long index = iarr + itime * img_info_rec->GetNpixelImg();
            subcube_std_arr[index] = std_arr[itime][iarr];
        }
    }
    MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), "subcube_std",
                          3, bitpix,
                          img_info_rec->GetNaxesArr(),
                          subcube_std_arr);
    delete [] subcube_std_arr;
    printf("=== output standardized subcube ===\n");

    
    MifImgInfo* img_info_par_cube = new MifImgInfo;
    img_info_par_cube->InitSetCube(1, 1, 1, nbin_rho, nbin_phi, nbin_psi);
    img_info_par_cube->PrintInfo();
    MifImgInfo* img_info_par_img = new MifImgInfo;
    img_info_par_img->InitSetImg(1, 1, nbin_rho, nbin_phi);
    img_info_par_img->PrintInfo();

    
    HistInfo1d* hi1d_xval = new HistInfo1d;
    HistInfo1d* hi1d_yval = new HistInfo1d;
    HistInfo1d* hi1d_zval = new HistInfo1d;
    hi1d_xval->InitSetByNbin(0.0,
                             img_info_rec->GetNaxesArrElm(0),
                             img_info_rec->GetNaxesArrElm(0));
    hi1d_yval->InitSetByNbin(0.0,
                             img_info_rec->GetNaxesArrElm(1),
                             img_info_rec->GetNaxesArrElm(1));
    LoadHi1dTime(argval->GetTimeDat(), hi1d_zval);
    printf("--- hi1d_xval, yval, zval ---\n");
    hi1d_xval->Print(stdout);
    hi1d_yval->Print(stdout);
    hi1d_zval->Print(stdout);
    printf("=== hi1d_xval, yval, zval ===\n");


    // velocity zero 
    double theta_zero = atan(0.0);
    double* cube_flat_vel_zero_arr = new double [nbin_rho * nbin_phi * nbin_psi];
    double* cube_vel_zero_arr = new double [nbin_rho * nbin_phi * nbin_psi];
    for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
        cube_vel_zero_arr[iarr] = 0.0;
        cube_flat_vel_zero_arr[iarr] = 0.0;
    }
    for(long itime = 0; itime < ntime; itime ++){
        for(long iposx = 0; iposx < hi1d_xval->GetNbin(); iposx ++){
            for(long iposy = 0; iposy < hi1d_yval->GetNbin(); iposy ++){
                double xval = hi1d_xval->GetBinCenter(iposx);
                double yval = hi1d_yval->GetBinCenter(iposy);
                double zval = time_arr[itime];
                for(long ipsi = 0; ipsi < nbin_psi; ipsi ++){
                    double psi = hi1d_psi->GetBinCenter(ipsi);
                    double rho_cos_phi = xval * cos(psi) + yval * sin(psi);
                    double rho_sin_phi = -1 * xval * sin(psi) * cos(theta_zero)
                        + yval * cos(psi) * cos(theta_zero)
                        + zval * sin(theta_zero);
                    double rho = sqrt(pow(rho_cos_phi, 2) + pow(rho_sin_phi, 2));
                    double phi = 0.0;
                    if(rho_sin_phi >= 0.0){
                        phi = acos(rho_cos_phi / rho);
                    } else {
                        phi = -1 * acos(rho_cos_phi / rho) + 2 * M_PI;
                    }
                    long irho = hi1d_rho->GetIbin(rho);
                    long iphi = hi1d_phi->GetIbin(phi);
                    long i_rho_phi_psi = irho + iphi * nbin_rho + ipsi * (nbin_rho * nbin_phi);
                    if(std_arr[itime][iposx + iposy * hi1d_xval->GetNbin()] > argval->GetSig()){
                        cube_vel_zero_arr[i_rho_phi_psi] ++;
                    }
                    cube_flat_vel_zero_arr[i_rho_phi_psi] ++;
                }
            }
        }
    }
    for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
        if(cube_flat_vel_zero_arr[iarr] > 1.0e-10){
            cube_vel_zero_arr[iarr] /= cube_flat_vel_zero_arr[iarr];
        }
    }
    char tag[kLineSize];
    sprintf(tag, "vel_zero");
    MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                          3, bitpix,
                          img_info_par_cube->GetNaxesArr(),
                          cube_vel_zero_arr);

    // velocity zero image
    double* img_vel_zero_arr = new double [nbin_rho * nbin_phi];
    double* img_flat_vel_zero_arr = new double [nbin_rho * nbin_phi];
    for(long iarr = 0; iarr < nbin_rho * nbin_phi; iarr ++){
        img_vel_zero_arr[iarr] = 0.0;
        img_flat_vel_zero_arr[iarr] = 0.0;
    }
    for(long itime = 0; itime < ntime; itime ++){
        for(long iposx = 0; iposx < hi1d_xval->GetNbin(); iposx ++){
            for(long iposy = 0; iposy < hi1d_yval->GetNbin(); iposy ++){
                double xval = hi1d_xval->GetBinCenter(iposx);
                double yval = hi1d_yval->GetBinCenter(iposy);

                double rho_cos_phi = xval;
                double rho_sin_phi = yval;
                double rho = sqrt(pow(rho_cos_phi, 2) + pow(rho_sin_phi, 2));
                double phi = 0.0;
                if(rho_sin_phi >= 0.0){
                    phi = acos(rho_cos_phi / rho);
                } else {
                    phi = -1 * acos(rho_cos_phi / rho) + 2 * M_PI;
                }
                long irho = hi1d_rho->GetIbin(rho);
                long iphi = hi1d_phi->GetIbin(phi);
                long i_rho_phi = irho + iphi * nbin_rho;
                if(std_arr[itime][iposx + iposy * hi1d_xval->GetNbin()] > argval->GetSig()){
                    img_vel_zero_arr[i_rho_phi] ++;
                }
                img_flat_vel_zero_arr[i_rho_phi] ++;
            }
        }
    }
    for(long iarr = 0; iarr < nbin_rho * nbin_phi; iarr ++){
        if(img_flat_vel_zero_arr[iarr] > 1.0e-10){
            img_vel_zero_arr[iarr] /= img_flat_vel_zero_arr[iarr];
        }
    }
    sprintf(tag, "vel_zero_img");
    MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                           2, bitpix,
                           img_info_par_img->GetNaxesArr(),
                           img_vel_zero_arr);
    sprintf(tag, "vel_zero_img_flat");
    MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                           2, bitpix,
                           img_info_par_img->GetNaxesArr(),
                           img_flat_vel_zero_arr);
    ///////
    
    double* par_cube_arr = new double [nbin_rho * nbin_phi * nbin_psi];
    double* par_cube_flat_arr = new double [nbin_rho * nbin_phi * nbin_psi];
    long* par_index_arr = new long [nbin_rho * nbin_phi * nbin_psi];    
    for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
        par_cube_arr[iarr] = 0.0;
        par_cube_flat_arr[iarr] = 0.0;
        par_index_arr[iarr] = 0;
    }
    double* large_arr = new double [img_info_rec->GetNpixelTotal()];
    for(long iarr = 0; iarr < img_info_rec->GetNpixelTotal(); iarr++){
        large_arr[iarr] = 0.0;
    }
    for(long ivel = 0; ivel < nbin_vel; ivel ++){
        double theta =  atan(hi1d_vel->GetBinCenter(ivel));
        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
            par_cube_arr[iarr] = 0.0;
            par_cube_flat_arr[iarr] = 0.0;
            par_index_arr[iarr] = 0;
        }
        for(long itime = 0; itime < ntime; itime ++){
            for(long iposx = 0; iposx < hi1d_xval->GetNbin(); iposx ++){
                for(long iposy = 0; iposy < hi1d_yval->GetNbin(); iposy ++){
                    double xval = hi1d_xval->GetBinCenter(iposx);
                    double yval = hi1d_yval->GetBinCenter(iposy);
                    double zval = time_arr[itime];
                    for(long ipsi = 0; ipsi < nbin_psi; ipsi ++){
                        double psi = hi1d_psi->GetBinCenter(ipsi);
                        double rho_cos_phi = xval * cos(psi) + yval * sin(psi);
                        double rho_sin_phi = -1 * xval * sin(psi) * cos(theta)
                            + yval * cos(psi) * cos(theta)
                            + zval * sin(theta);
                        double rho = sqrt(pow(rho_cos_phi, 2) + pow(rho_sin_phi, 2));
                        double phi = 0.0;
                        if(rho_sin_phi >= 0.0){
                            phi = acos(rho_cos_phi / rho);
                        } else {
                            phi = -1 * acos(rho_cos_phi / rho) + 2 * M_PI;
                        }
                        long irho = hi1d_rho->GetIbin(rho);
                        long iphi = hi1d_phi->GetIbin(phi);
                        long i_rho_phi_psi = irho + iphi * nbin_rho + ipsi * (nbin_rho * nbin_phi);

                        if(std_arr[itime][iposx + iposy * hi1d_xval->GetNbin()] > argval->GetSig()){
                            par_cube_arr[i_rho_phi_psi] ++;
                        }
                        par_cube_flat_arr[i_rho_phi_psi] ++;
                        // cube_arr[i_rho_phi_psi] += std_arr[itime][iposx + iposy * hi1d_xval->GetNbin()];
                    }
                }
            }
        }
        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
            if(par_cube_flat_arr[iarr] > 1.0e-10){
                par_cube_arr[iarr] /= par_cube_flat_arr[iarr];
            }
        }

        sprintf(tag, "vel_%2.2ld", ivel);
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              3, bitpix,
                              img_info_par_cube->GetNaxesArr(),
                              par_cube_arr);

        double* hough_ratio_arr = new double[nbin_rho * nbin_phi * nbin_psi];
        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
            hough_ratio_arr[iarr] = 0.0;
        }
        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
            if(par_cube_arr[iarr] < 1.0e-10){
                continue;
            }
            
            long ipsi = iarr / (nbin_rho * nbin_phi);
            long index2 = iarr % (nbin_rho * nbin_phi);
            long iphi = index2 / nbin_rho;
            long irho = index2 % nbin_rho;
            double rho = hi1d_rho->GetBinCenter(irho);
            double phi = hi1d_phi->GetBinCenter(iphi);
            double psi = hi1d_psi->GetBinCenter(ipsi);

            double zval_md = (hi1d_zval->GetUp() - hi1d_zval->GetLo()) / 2.0;
            double tpar_md = ( zval_md - rho * sin(theta) * sin(phi) ) / cos(theta);
            double xval_md = rho * (cos(psi) * cos(phi)
                                    - sin(psi) * cos(theta) * sin(phi) )
                + tpar_md * sin(psi) * sin(theta);
            double yval_md = rho * (sin(psi) * cos(phi)
                                    + cos(psi) * cos(theta) * sin(phi) )
                - tpar_md * cos(psi) * sin(theta);

            double rho_cos_phi = xval_md;
            double rho_sin_phi = yval_md;
            double rho_md = sqrt(pow(rho_cos_phi, 2) + pow(rho_sin_phi, 2));
            double phi_md = 0.0;
            if(rho_sin_phi >= 0.0){
                phi_md = acos(rho_cos_phi / rho_md);
            } else {
                phi_md = -1 * acos(rho_cos_phi / rho_md) + 2 * M_PI;
            }
            long irho_md = hi1d_rho->GetIbin(rho_md);
            long iphi_md = hi1d_phi->GetIbin(phi_md);
            long i_rho_phi = irho_md + iphi_md * nbin_rho;
            
            hough_ratio_arr[iarr] = par_cube_arr[iarr] - cube_vel_zero_arr[i_rho_phi];
        }

        long nlarge = 10000;
        // MiSort::KthElement<double, long>(nlarge, nbin_rho * nbin_phi * nbin_psi, cube_arr, index_arr, 1);
        //MiSort::Sort<double, long>(nbin_rho * nbin_phi * nbin_psi, par_cube_arr, par_index_arr, 1);
        MiSort::Sort<double, long>(nbin_rho * nbin_phi * nbin_psi, hough_ratio_arr, par_index_arr, 1);
        
        // k-th element
        for(long iarr = 0; iarr < img_info_rec->GetNpixelTotal(); iarr++){
            large_arr[iarr] = 0.0;
        }

        int ndetect = 0;
        for(int ilarge = 0; ilarge < nlarge; ilarge ++){
            long ipsi = par_index_arr[ilarge] / (nbin_rho * nbin_phi);
            long index2 = par_index_arr[ilarge] % (nbin_rho * nbin_phi);
            long iphi = index2 / nbin_rho;
            long irho = index2 % nbin_rho;

            double rho = hi1d_rho->GetBinCenter(irho);
            double phi = hi1d_phi->GetBinCenter(iphi);
            double psi = hi1d_psi->GetBinCenter(ipsi);
            //printf("ivel = %ld, vel = %e, theta = %e, rho = %e, phi = %e, psi = %e \n",
            //       ivel, tan(theta), theta, rho, phi, psi);

            double zval_md = (hi1d_zval->GetUp() - hi1d_zval->GetLo()) / 2.0;
            double tpar_md = ( zval_md - rho * sin(theta) * sin(phi) ) / cos(theta);
            double xval_md = rho * (cos(psi) * cos(phi)
                                    - sin(psi) * cos(theta) * sin(phi) )
                + tpar_md * sin(psi) * sin(theta);
            double yval_md = rho * (sin(psi) * cos(phi)
                                    + cos(psi) * cos(theta) * sin(phi) )
                - tpar_md * cos(psi) * sin(theta);

            double rho_cos_phi = xval_md;
            double rho_sin_phi = yval_md;
            double rho_md = sqrt(pow(rho_cos_phi, 2) + pow(rho_sin_phi, 2));
            double phi_md = 0.0;
            if(rho_sin_phi >= 0.0){
                phi_md = acos(rho_cos_phi / rho_md);
            } else {
                phi_md = -1 * acos(rho_cos_phi / rho_md) + 2 * M_PI;
            }
            long irho_md = hi1d_rho->GetIbin(rho_md);
            long iphi_md = hi1d_phi->GetIbin(phi_md);
            long i_rho_phi = irho_md + iphi_md * nbin_rho;

            printf("par_cube = %e, vel_zero = %e\n",
                   par_cube_arr[par_index_arr[ilarge]], cube_vel_zero_arr[i_rho_phi]);

            double zval_lo = hi1d_zval->GetLo();
            double zval_up = hi1d_zval->GetUp();
            double tpar_lo = ( zval_lo - rho * sin(theta) * sin(phi) ) / cos(theta);
            double tpar_up = ( zval_up - rho * sin(theta) * sin(phi) ) / cos(theta);
            long ntpar = 100;
            double delta_tpar = (tpar_up - tpar_lo) / ntpar;

            if(hough_ratio_arr[par_index_arr[ilarge]] < 0.0){
                continue;
            }
            if(cube_vel_zero_arr[i_rho_phi] < 1.0e-10){
                continue;
            }
            
            for(long itpar = 0; itpar < ntpar; itpar ++){
                double tpar = tpar_lo + itpar * delta_tpar;
                double xval = rho * (cos(psi) * cos(phi)
                                     - sin(psi) * cos(theta) * sin(phi) )
                    + tpar * sin(psi) * sin(theta);
                double yval = rho * (sin(psi) * cos(phi)
                                     + cos(psi) * cos(theta) * sin(phi) )
                    - tpar * cos(psi) * sin(theta);
                double zval = rho * sin(theta) * sin(phi) + tpar * cos(theta);

                if(xval < hi1d_xval->GetLo() || hi1d_xval->GetUp() < xval){
                    continue;
                }
                if(yval < hi1d_yval->GetLo() || hi1d_yval->GetUp() < yval){
                    continue;
                }            
                if(zval < hi1d_zval->GetLo() || hi1d_zval->GetUp() < zval){
                    continue;
                }            
                long iposx = hi1d_xval->GetIbin(xval);
                long iposy = hi1d_yval->GetIbin(yval);
                long iposz = hi1d_zval->GetIbin(zval);
                long index = iposx + iposy * hi1d_xval->GetNbin()
                    + iposz * (hi1d_xval->GetNbin() * hi1d_yval->GetNbin());

                large_arr[index] = 1.0;
            }
            ndetect ++;

            if(ndetect > 10){
                break;
            }
        }
        printf("------\n");
        
        sprintf(tag, "large_%2.2ld", ivel);
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              3, bitpix,
                              img_info_rec->GetNaxesArr(),
                              large_arr);
    }


    double time_ed = MiTime::GetTimeSec();
    printf("time_ed - time_st = %e\n", time_ed - time_st);
    
    return status_prog;
}
