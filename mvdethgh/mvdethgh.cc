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

    MifImgInfo* img_info_par_cube = new MifImgInfo;
    img_info_par_cube->InitSetCube(1, 1, 1, nbin_rho, nbin_phi, nbin_psi);
    img_info_par_cube->PrintInfo();
    
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

    // velocity zero
    printf("--- velocity zero ---\n");
    double theta_zero = atan(0.0);
    double* cube_vel_zero_arr = new double [nbin_rho * nbin_phi * nbin_psi];
    for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
        cube_vel_zero_arr[iarr] = 0.0;
    }
    for(long itime = 0; itime < ntime; itime ++){
        for(long iposx = 0; iposx < hi1d_xval->GetNbin(); iposx ++){
            for(long iposy = 0; iposy < hi1d_yval->GetNbin(); iposy ++){
                double xval = hi1d_xval->GetBinCenter(iposx);
                double yval = hi1d_yval->GetBinCenter(iposy);
                double zval = time_arr[itime];
                for(long ipsi = 0; ipsi < nbin_psi; ipsi ++){
                    if(data_img_arr[itime][iposx + iposy * hi1d_xval->GetNbin()] < argval->GetSig()){
                        continue;
                    }
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
                    if(data_img_arr[itime][iposx + iposy * hi1d_xval->GetNbin()] > argval->GetSig()){
                        cube_vel_zero_arr[i_rho_phi_psi] ++;
                    }
                }
            }
        }
    }

    // load flat    
    double* cube_flat_vel_zero_arr = NULL;
    char flat_zero_file[kLineSize];
    sprintf(flat_zero_file, "%s/%s_vel_zero.fits", argval->GetOutdir().c_str(), "flat");
    MifFits::InFitsCubeD(flat_zero_file,
                         img_info_par_cube,
                         &bitpix,
                         &cube_flat_vel_zero_arr);
    for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
        if(cube_flat_vel_zero_arr[iarr] > 1.0e-10){
            cube_vel_zero_arr[iarr] /= cube_flat_vel_zero_arr[iarr];
        }
    }
    delete [] cube_flat_vel_zero_arr;
    
    MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), "vel_zero",
                          3, bitpix,
                          img_info_par_cube->GetNaxesArr(),
                          cube_vel_zero_arr);

    
    double* par_cube_arr = new double [nbin_rho * nbin_phi * nbin_psi];
    long* par_index_arr = new long [nbin_rho * nbin_phi * nbin_psi];    
    for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
        par_cube_arr[iarr] = 0.0;
        par_index_arr[iarr] = 0;
    }


    MifImgInfo* img_info_rec = new MifImgInfo;
    img_info_rec->InitSetCube(1, 1, 1,
                              img_info->GetNaxesArrElm(0),
                              img_info->GetNaxesArrElm(1),
                              ntime);
    img_info_rec->PrintInfo();
    
    double* large_arr = new double [img_info_rec->GetNpixelTotal()];
    for(long iarr = 0; iarr < img_info_rec->GetNpixelTotal(); iarr++){
        large_arr[iarr] = 0.0;
    }
    for(long ivel = 0; ivel < nbin_vel; ivel ++){
        double theta =  atan(hi1d_vel->GetBinCenter(ivel));
        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
            par_cube_arr[iarr] = 0.0;
            par_index_arr[iarr] = 0;
        }
        for(long itime = 0; itime < ntime; itime ++){
            for(long iposx = 0; iposx < hi1d_xval->GetNbin(); iposx ++){
                for(long iposy = 0; iposy < hi1d_yval->GetNbin(); iposy ++){
                    double xval = hi1d_xval->GetBinCenter(iposx);
                    double yval = hi1d_yval->GetBinCenter(iposy);
                    double zval = time_arr[itime];
                    for(long ipsi = 0; ipsi < nbin_psi; ipsi ++){

                        if(data_img_arr[itime][iposx + iposy * hi1d_xval->GetNbin()] < argval->GetSig()){
                            continue;
                        }
                        
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
                        if(data_img_arr[itime][iposx + iposy * hi1d_xval->GetNbin()] > argval->GetSig()){
                            par_cube_arr[i_rho_phi_psi] ++;
                        }
                    }
                }
            }
        }

        // load flat
        double* par_cube_flat_arr = NULL;
        int bitpix = 0;
        char flat_file[kLineSize];
        sprintf(flat_file, "%s/%s_vel_%2.2ld.fits", argval->GetOutdir().c_str(), "flat", ivel);
        MifFits::InFitsCubeD(flat_file,
                             img_info_par_cube,
                             &bitpix,
                             &par_cube_flat_arr);
        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
            if(par_cube_flat_arr[iarr] > 1.0e-10){
                par_cube_arr[iarr] /= par_cube_flat_arr[iarr];
            }
        }
        delete [] par_cube_flat_arr;

        char tag[kLineSize];
        sprintf(tag, "vel_%2.2ld", ivel);
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              3, bitpix,
                              img_info_par_cube->GetNaxesArr(),
                              par_cube_arr);

        printf("kkkkkk\n");
        ///////
        
        
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

            printf("ipsi iphi irho = %ld %ld %ld\n", ipsi, iphi, irho);
            
            double rho = hi1d_rho->GetBinCenter(irho);
            double phi = hi1d_phi->GetBinCenter(iphi);
            double psi = hi1d_psi->GetBinCenter(ipsi);

            printf("psi phi rho = %e %e %e\n", psi, phi, rho);
            
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

            printf("rho_md phi_md = %e %e\n", rho_md, phi_md);

            if(rho_md  < hi1d_rho->GetLo() || hi1d_rho->GetUp() < rho_md){
                hough_ratio_arr[iarr] = 0.0;
                continue;
            }
            if(phi_md  < hi1d_phi->GetLo() || hi1d_phi->GetUp() < phi_md){
                hough_ratio_arr[iarr] = 0.0;
                continue;
            }
            long irho_md = hi1d_rho->GetIbin(rho_md);
            long iphi_md = hi1d_phi->GetIbin(phi_md);
            printf("irho_md iphi_md = %ld %ld\n", irho_md, iphi_md);
            
            long i_rho_phi = irho_md + iphi_md * nbin_rho;
            hough_ratio_arr[iarr] = par_cube_arr[iarr] - cube_vel_zero_arr[i_rho_phi];
        }

        printf("jjjjj\n");
        
        MiSort::Sort<double, long>(nbin_rho * nbin_phi * nbin_psi, hough_ratio_arr, par_index_arr, 1);
        
        // k-th element
        for(long iarr = 0; iarr < img_info_rec->GetNpixelTotal(); iarr++){
            large_arr[iarr] = 0.0;
        }

        long nlarge = 10000;
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
