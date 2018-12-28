#include "sub_mvdethgh.h"

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
                 const MifImgInfo* const img_info_subimg,
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
    double rho_up = sqrt( pow(img_info_subimg->GetNaxesArrElm(0), 2)
                          + pow(img_info_subimg->GetNaxesArrElm(1), 2) );
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

// get mean and stddev by clipping
void GetMeanStddevClip(long narr, const double* const val_arr,
                       int nclip, double significance,
                       double* const mean_ptr,
                       double* const stddev_ptr)
{
    double mean   = 0.0;
    double stddev = 0.0;
    int* use_arr = new int[narr];
    for(long iarr = 0; iarr < narr; iarr++){
        use_arr[iarr] = 1;
    }
    
    for(int iclip = 0; iclip < nclip; iclip ++){
        // get mean, stddev
        vector<double> val_vec;
        for(long iarr = 0; iarr < narr; iarr++){
            if(use_arr[iarr] > 0){
                val_vec.push_back(val_arr[iarr]);
            }
        }
        mean   = MirMath::GetAMean(val_vec);
        stddev = MirMath::GetSqrtOfUnbiasedVariance(val_vec);
        printf("iclip, mean, stddev, nsize = %d, %e, %e, %d\n",
               iclip, mean, stddev, (int) val_vec.size());
        long ninc_unuse = 0;
        for(long iarr = 0; iarr < narr; iarr++){
            if( use_arr[iarr] > 0 &&
                fabs(val_arr[iarr] - mean) > significance * stddev){
                use_arr[iarr] = 0;
                ninc_unuse ++;
            }
        }
        if(0 == ninc_unuse){
            printf("clipping end.\n");
            break;
        }
    }
    *mean_ptr = mean;
    *stddev_ptr = stddev;
}
    

//////////////////////////////////////////




//        char out_tmp[kLineSize];
//        sprintf(out_tmp, "tmp_%2.2ld.qdp", ivel);
//        FILE* fp_tmp = fopen(out_tmp, "w");
//        
//        double* lc_mean_arr = new double [img_info_rec->GetNpixelImg()];        
//        double* lc_stddev_arr = new double [img_info_rec->GetNpixelImg()];
//        double* lc_max_arr = new double [img_info_rec->GetNpixelImg()];
//        double* lc_ratio_max_arr = new double [img_info_rec->GetNpixelImg()];
//        for(long iarr = 0; iarr < img_info_rec->GetNpixelImg(); iarr ++){
//            lc_mean_arr[iarr] = 0.0;
//            lc_stddev_arr[iarr] = 0.0;
//            lc_max_arr[iarr] = 0.0;
//            lc_ratio_max_arr[iarr] = 0.0;
//        }
//        for(long iposx = 0; iposx < hi1d_xval->GetNbin(); iposx ++){
//            for(long iposy = 0; iposy < hi1d_yval->GetNbin(); iposy ++){
//                long itime_mid = ntime / 2;
//                double xval = hi1d_xval->GetBinCenter(iposx);
//                double yval = hi1d_yval->GetBinCenter(iposy);
//                double zval = time_arr[itime_mid];
//
//                double* lc_arr = new double[nbin_psi];
//                for(long ipsi = 0; ipsi < nbin_psi; ipsi ++){
//                    double psi = hi1d_psi->GetBinCenter(ipsi);
//                    double rho_cos_phi = xval * cos(psi) + yval * sin(psi);
//                    double rho_sin_phi = -1 * xval * sin(psi) * cos(theta)
//                        + yval * cos(psi) * cos(theta)
//                        + zval * sin(theta);
//                    double rho = sqrt(pow(rho_cos_phi, 2) + pow(rho_sin_phi, 2));
//                    double phi = 0.0;
//                    if(rho_sin_phi >= 0.0){
//                        phi = acos(rho_cos_phi / rho);
//                    } else {
//                        phi = -1 * acos(rho_cos_phi / rho) + 2 * M_PI;
//                    }
//                    long irho = hi1d_rho->GetIbin(rho);
//                    long iphi = hi1d_phi->GetIbin(phi);
//                    long i_rho_phi_psi = irho + iphi * nbin_rho + ipsi * (nbin_rho * nbin_phi);
//                    lc_arr[ipsi] = cube_arr[i_rho_phi_psi];
//                }
//
//                
//                double lc_mean   = MirMath::GetAMean(nbin_psi, lc_arr);                
//                double lc_stddev = MirMath::GetStddev(nbin_psi, lc_arr);
//                double lc_max    = MirMath::GetMax(nbin_psi, lc_arr);
//                long index = iposx + iposy * hi1d_xval->GetNbin();
//                lc_mean_arr[index] = lc_mean;
//                lc_stddev_arr[index] = lc_stddev;
//                lc_max_arr[index] = lc_max;
//                double ratio_max = 0.0;
//                if(lc_mean < 1e-10){
//                    ratio_max = 0.0;
//                } else {
//                    ratio_max = lc_max / lc_mean;
//
//                    if(100 < xval && xval < 200 &&
//                       150 < yval && yval < 250){
//                        for(long ipsi = 0; ipsi < nbin_psi; ipsi ++){
//                            fprintf(fp_tmp, "%ld %e\n", ipsi, lc_arr[ipsi]);
//                        }
//                    }
//                }
//                lc_ratio_max_arr[index] = ratio_max * lc_max;
//                delete [] lc_arr;
//
//                fprintf(fp_tmp, "\n");
//            }
//        }
//
//        fclose(fp_tmp);        
//
//        sprintf(tag, "lc_mean_%2.2ld", ivel);
//        MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
//                             2, bitpix,
//                             img_info_subimg->GetNaxesArr(),
//                             lc_mean_arr);
//        sprintf(tag, "lc_stddev_%2.2ld", ivel);
//        MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
//                             2, bitpix,
//                             img_info_subimg->GetNaxesArr(),
//                             lc_stddev_arr);
//        sprintf(tag, "lc_max_%2.2ld", ivel);
//        MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
//                               2, bitpix,
//                               img_info_subimg->GetNaxesArr(),
//                               lc_max_arr);
//        sprintf(tag, "lc_ratio_max_%2.2ld", ivel);
//        MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
//                               2, bitpix,
//                               img_info_subimg->GetNaxesArr(),
//                               lc_ratio_max_arr);
//        delete [] lc_mean_arr;
//        delete [] lc_stddev_arr;
//        delete [] lc_max_arr;
//        delete [] lc_ratio_max_arr;        

//        // renormalize by initial brightness
//        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
//            long ipsi = iarr  / (nbin_rho * nbin_phi);
//            long index2 = iarr % (nbin_rho * nbin_phi);
//            long iphi = index2 / nbin_rho;
//            long irho = index2 % nbin_rho;
//            double rho = hi1d_rho->GetBinCenter(irho);
//            double phi = hi1d_phi->GetBinCenter(iphi);
//            double psi = hi1d_psi->GetBinCenter(ipsi);
//            double zval_lo = hi1d_zval->GetLo();
//            double tpar_lo = ( zval_lo - rho * sin(theta) * sin(phi) ) / cos(theta);
//           
//            double xval = rho * (cos(psi) * cos(phi)
//                                 - sin(psi) * cos(theta) * sin(phi) )
//                + tpar_lo * sin(psi) * sin(theta);
//            double yval = rho * (sin(psi) * cos(phi)
//                                 + cos(psi) * cos(theta) * sin(phi) )
//                - tpar_lo * cos(psi) * sin(theta);
//            if(xval < hi1d_xval->GetLo() || hi1d_xval->GetUp() < xval){
//                continue;
//            }
//            if(yval < hi1d_yval->GetLo() || hi1d_yval->GetUp() < yval){
//                continue;
//            }            
//            long iposx = hi1d_xval->GetIbin(xval);
//            long iposy = hi1d_yval->GetIbin(yval);
//            if(std_arr[0][iposx + iposy * hi1d_xval->GetNbin()] > argval->GetSig()){
//                cube_arr[iarr] /= std_arr[0][iposx + iposy * hi1d_xval->GetNbin()];
//            }
//        }











//            // calc norm
//            double zval_md = (hi1d_zval->GetUp() - hi1d_zval->GetLo()) / 2.0;
//            double tpar_md = ( zval_md - rho * sin(theta) * sin(phi) ) / cos(theta);
//            double xval_md = rho * (cos(psi) * cos(phi)
//                                    - sin(psi) * cos(theta) * sin(phi) )
//                + tpar_md * sin(psi) * sin(theta);
//            double yval_md = rho * (sin(psi) * cos(phi)
//                                    + cos(psi) * cos(theta) * sin(phi) )
//                - tpar_md * cos(psi) * sin(theta);
//            double* lc_arr = new double[nbin_psi];
//            for(long ipsi = 0; ipsi < nbin_psi; ipsi ++){
//                double psi = hi1d_psi->GetBinCenter(ipsi);
//                double rho_cos_phi = xval_md * cos(psi) + yval_md * sin(psi);
//                double rho_sin_phi = -1 * xval_md * sin(psi) * cos(theta)
//                    + yval_md * cos(psi) * cos(theta)
//                    + zval_md * sin(theta);
//                double rho = sqrt(pow(rho_cos_phi, 2) + pow(rho_sin_phi, 2));
//                double phi = 0.0;
//                if(rho_sin_phi >= 0.0){
//                    phi = acos(rho_cos_phi / rho);
//                } else {
//                    phi = -1 * acos(rho_cos_phi / rho) + 2 * M_PI;
//                }
//                long irho = hi1d_rho->GetIbin(rho);
//                long iphi = hi1d_phi->GetIbin(phi);
//                long i_rho_phi_psi = irho + iphi * nbin_rho + ipsi * (nbin_rho * nbin_phi);
//                lc_arr[ipsi] = cube_arr[i_rho_phi_psi];
//            }
//            
//            double lc_mean   = MirMath::GetAMean(nbin_psi, lc_arr);                
//            double lc_stddev = MirMath::GetStddev(nbin_psi, lc_arr);
//            double lc_max    = MirMath::GetMax(nbin_psi, lc_arr);
//            double ratio_max = 0.0;
//            if(lc_mean < 1e-10){
//                ratio_max = 0.0;
//            } else {
//                ratio_max = lc_max / lc_mean;
//            }
//            delete [] lc_arr;
//
//            double norm = ratio_max;
//            
//            // =====
