#include "mi_time.h"
#include "mvdethghlib.h"
#include "arg_mvdethgh.h"

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

    double** data_img_above_th_arr    = NULL;
    long ndet = 0;
    GenImgAboveThArr(ntime,
                     data_img_arr,
                     std_img_arr,
                     img_info,
                     argval->GetSig(),
                     &data_img_above_th_arr,
                     &ndet);
    printf("ndet above threshold = %ld\n", ndet);
    SaveCube(ntime, img_info, data_img_above_th_arr,
             bitpix,
             argval->GetOutdir(),
             argval->GetOutfileHead(),
             "cube_above_th");

    printf("--- load parameters ---\n");
    HistInfo1d* hi1d_vel   = new HistInfo1d;
    HistInfo1d* hi1d_rho   = new HistInfo1d;
    HistInfo1d* hi1d_phi   = new HistInfo1d;
    HistInfo1d* hi1d_psi   = new HistInfo1d;
    LoadHi1dVel(argval->GetVelDat(), hi1d_vel);
    LoadHi1dPar(argval->GetResDat(), img_info,
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

    MifImgInfo* img_info_par3d = new MifImgInfo;
    img_info_par3d->InitSetCube(1, 1, 1, nbin_rho, nbin_phi, nbin_psi);
    img_info_par3d->PrintInfo();
    
    HistInfo1d* hi1d_xval = new HistInfo1d;
    HistInfo1d* hi1d_yval = new HistInfo1d;
    HistInfo1d* hi1d_zval = new HistInfo1d;
    hi1d_xval->InitSetByNbin(0.0,
                             img_info->GetNaxesArrElm(0),
                             img_info->GetNaxesArrElm(0));
    hi1d_yval->InitSetByNbin(0.0,
                             img_info->GetNaxesArrElm(1),
                             img_info->GetNaxesArrElm(1));
    LoadHi1dTime(argval->GetTimeDat(), hi1d_zval);
    printf("--- hi1d_xval, yval, zval ---\n");
    hi1d_xval->Print(stdout);
    hi1d_yval->Print(stdout);
    hi1d_zval->Print(stdout);
    printf("=== hi1d_xval, yval, zval ===\n");

    long nbin_par = nbin_rho * nbin_phi * nbin_psi * nbin_vel;
    double* par4d_arr = new double [nbin_par];
    double* par4d_weight_arr = new double [nbin_par];
    long* par_index_arr = new long [nbin_par];
    for(long iarr = 0; iarr < nbin_par; iarr ++){
        par4d_arr[iarr] = 0.0;
        par4d_weight_arr[iarr] = 0.0;
        par_index_arr[iarr] = 0;
    }
    MifImgInfo* img_info_rec = new MifImgInfo;
    img_info_rec->InitSetCube(1, 1, 1,
                              img_info->GetNaxesArrElm(0),
                              img_info->GetNaxesArrElm(1),
                              ntime);
    img_info_rec->PrintInfo();

    for(long itime = 0; itime < ntime; itime ++){
        for(long iposx = 0; iposx < hi1d_xval->GetNbin(); iposx ++){
            for(long iposy = 0; iposy < hi1d_yval->GetNbin(); iposy ++){
                double xval = hi1d_xval->GetBinCenter(iposx);
                double yval = hi1d_yval->GetBinCenter(iposy);
                double zval = time_arr[itime];
                if(std_img_arr[itime][iposx + iposy * hi1d_xval->GetNbin()]
                   < argval->GetSig()){
                    continue;
                }
                for(long ivel = 0; ivel < nbin_vel; ivel ++){
                    double theta =  atan(hi1d_vel->GetBinCenter(ivel));
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
                        long i_rho_phi_psi_vel = irho + iphi * nbin_rho
                            + ipsi * (nbin_rho * nbin_phi)
                            + ivel * (nbin_rho * nbin_phi * nbin_psi);
                        par4d_arr[i_rho_phi_psi_vel] ++;
                        par4d_weight_arr[i_rho_phi_psi_vel] +=
                            std_img_arr[itime][iposx + iposy * hi1d_xval->GetNbin()];
                    }
                }
            }
        }
    }


    // load flat
    double* par4d_flat_arr = new double [nbin_par];
    for(long ivel = 0; ivel < nbin_vel; ivel ++){
        double* par3d_flat_arr = NULL;
        int bitpix = 0;
        char flat_file[kLineSize];
        sprintf(flat_file, "%s/%s_vel_%2.2ld.fits", argval->GetOutdir().c_str(), "flat", ivel);
        MifFits::InFitsCubeD(flat_file,
                             img_info_par3d,
                             &bitpix,
                             &par3d_flat_arr);
        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
            long i_rho_phi_psi_vel = iarr + ivel * (nbin_rho * nbin_phi * nbin_psi);
            par4d_flat_arr[i_rho_phi_psi_vel] = par3d_flat_arr[iarr];
        }
        delete [] par3d_flat_arr;
    }
    // divide by flat
    for(long iarr = 0; iarr < nbin_par; iarr ++){
        if(par4d_flat_arr[iarr] > 1.0e-10){
            par4d_arr[iarr] /= par4d_flat_arr[iarr];
            par4d_weight_arr[iarr] /= par4d_flat_arr[iarr];
        }
    }
    delete [] par4d_flat_arr;

    // save par3d_arr
    for(long ivel = 0; ivel < nbin_vel; ivel ++){
        double* par3d_arr = new double [nbin_rho * nbin_phi * nbin_psi];
        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
            long i_rho_phi_psi_vel = iarr + ivel * (nbin_rho * nbin_phi * nbin_psi);
            par3d_arr[iarr] = par4d_weight_arr[i_rho_phi_psi_vel];
        }
        char tag[kLineSize];
        sprintf(tag, "ivel_%2.2ld", ivel);
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              3, bitpix,
                              img_info_par3d->GetNaxesArr(),
                              par3d_arr);
        delete [] par3d_arr;
    }

    // MiSort::Sort<double, long>(nbin_par, par4d_arr, par_index_arr, 1);
    MiSort::Sort<double, long>(nbin_par, par4d_weight_arr, par_index_arr, 1);


    char outfile[kLineSize];
    sprintf(outfile, "%s/%s_par.dat",
            argval->GetOutdir().c_str(), argval->GetOutfileHead().c_str());
    FILE* fp_outfile = fopen(outfile, "w");

    //////////
    double* rec_arr    = new double [img_info_rec->GetNpixelTotal()];
    for(long iarr = 0; iarr < img_info_rec->GetNpixelTotal(); iarr++){
        rec_arr[iarr] = 0.0;
    }
    double margin = argval->GetMargin();
    int ndetect_max = argval->GetNdet();
    int ndetect = 0;
    for(long iarr = 0; iarr < nbin_par; iarr ++){
        long ivel = par_index_arr[iarr] / (nbin_rho * nbin_phi * nbin_psi);
        long index2 = par_index_arr[iarr] % (nbin_rho * nbin_phi * nbin_psi);
        long ipsi = index2 / (nbin_rho * nbin_phi);
        long index3 = index2 % (nbin_rho * nbin_phi);
        long iphi = index3 / nbin_rho;
        long irho = index3 % nbin_rho;
        double rho = hi1d_rho->GetBinCenter(irho);
        double phi = hi1d_phi->GetBinCenter(iphi);
        double psi = hi1d_psi->GetBinCenter(ipsi);
        double vel = hi1d_vel->GetBinCenter(ivel);

        // pixel value along line, which must be constant
        int flag_bad = 0;
        double* line_arr = new double[ntime];
        long* index_line_arr = new long[ntime];
        double sum2_std_img = 0.0;
        for(long itime = 0; itime < hi1d_zval->GetNbin(); itime ++){
            double zval  = hi1d_zval->GetBinCenter(itime);
            double theta =  atan(vel);
            double xval = rho * (cos(phi) * cos(psi) - sin(phi) * sin(psi) / cos(theta) )
                + zval * tan(theta) * sin(psi);
            double yval = rho * (cos(phi) * sin(psi) + sin(phi) * cos(psi) / cos(theta) )
                - zval * tan(theta) * cos(psi);
            if(xval < hi1d_xval->GetLo() + margin || hi1d_xval->GetUp() - margin < xval){
                line_arr[itime] = 0.0;
                flag_bad ++;
                continue;
            }
            if(yval < hi1d_yval->GetLo() + margin || hi1d_yval->GetUp() - margin < yval){
                line_arr[itime] = 0.0;
                flag_bad ++;
                continue;
            }
            long iposx = hi1d_xval->GetIbin(xval);
            long iposy = hi1d_yval->GetIbin(yval);
            index_line_arr[itime] = iposx + iposy * hi1d_xval->GetNbin()
                + itime * hi1d_xval->GetNbin() * hi1d_yval->GetNbin();
            line_arr[itime] = std_img_arr[itime][iposx + iposy * hi1d_xval->GetNbin()];
            sum2_std_img += pow(std_img_arr[itime][iposx + iposy * hi1d_xval->GetNbin()], 2);
        }
        double stddev = MirMath::GetStddev(ntime, line_arr);
        double mean   = MirMath::GetAMean(ntime, line_arr);
        double sigma_std_img = sqrt(sum2_std_img);

        if(mean <= DBL_EPSILON){
            flag_bad ++;
        }
        if(stddev / mean > argval->GetVarRatio()){
            flag_bad ++;
        }
      
        if(0 == flag_bad){
            ndetect ++;
            printf("ndetect = %d, iarr = %ld, mean = %e, stddev / mean = %e, par4d = %e, sigma_line = %e  "
                   "!  %e  %e  %e  %e  ! vel, rho, phi, psi\n",
                   ndetect, iarr, mean, stddev/mean, par4d_weight_arr[par_index_arr[iarr]],
                   sigma_std_img,
                   vel, rho, phi, psi);

            fprintf(fp_outfile, "%e  %e  %e  %e  ! vel, rho, phi, psi\n", vel, rho, phi, psi);
            
            for(long itime = 0; itime < hi1d_zval->GetNbin(); itime ++){
                long index = index_line_arr[itime];
                rec_arr[index] = line_arr[itime];
            }
        }
        delete [] line_arr;
        delete [] index_line_arr;

        if(sigma_std_img < argval->GetSigLine()){
            printf("break by 1\n");
            break;
        }
        if(ndetect >= ndetect_max){
            printf("break by 2\n");
            break;
        }
    }
    fclose(fp_outfile);


    















    
        
    MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), "rec",
                          3, bitpix,
                          img_info_rec->GetNaxesArr(),
                          rec_arr);

    double time_ed = MiTime::GetTimeSec();
    printf("time_ed - time_st = %e\n", time_ed - time_st);
    
    return status_prog;
}
