#include "mi_time.h"
#include "mvdethghlib.h"
#include "arg_mkflatxyt.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValMkflatxyt* argval = new ArgValMkflatxyt;
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
    long*   npos_arr = NULL;
    double** xpos_arr = NULL;
    double** ypos_arr = NULL;
    double xpos_lo = 0.0;
    double xpos_up = 0.0;
    double ypos_lo = 0.0;
    double ypos_up = 0.0;    
    double rho_up = 0.0;
    LoadDataXYT(argval->GetDataList(),
                &ntime,
                &time_arr,
                &npos_arr,
                &xpos_arr,
                &ypos_arr,
                &xpos_lo,
                &xpos_up,
                &ypos_lo,
                &ypos_up,
                &rho_up);
    printf("rho_up = %e\n", rho_up);

    printf("--- load parameters ---\n");
    HistInfo1d* hi1d_vel   = new HistInfo1d;
    HistInfo1d* hi1d_rho   = new HistInfo1d;
    HistInfo1d* hi1d_phi   = new HistInfo1d;
    HistInfo1d* hi1d_psi   = new HistInfo1d;
    LoadHi1dVel(argval->GetVelDat(), hi1d_vel);
    LoadHi1dPar(argval->GetResDat(), rho_up,
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

    // scale zpos
    double* zpos_arr = new double[ntime];
    double scale = hi1d_vel->GetUp() / time_arr[ntime - 1];
    printf("scale = %e\n", scale);
    for(long itime = 0; itime < ntime; itime ++){
        zpos_arr[itime] = time_arr[itime] * scale;
        printf("%ld  %e  %e\n", itime, time_arr[itime], zpos_arr[itime]);
    }
    
    MifImgInfo* img_info_par3d = new MifImgInfo;
    img_info_par3d->InitSetCube(1, 1, 1, nbin_rho, nbin_phi, nbin_psi);
    img_info_par3d->PrintInfo();
    
    HistInfo1d* hi1d_xval = new HistInfo1d;
    HistInfo1d* hi1d_yval = new HistInfo1d;
    hi1d_xval->InitSetByNbin(xpos_lo, xpos_up, (long) ceil(xpos_up - xpos_lo));
    hi1d_yval->InitSetByNbin(ypos_lo, ypos_up, (long) ceil(ypos_up - ypos_lo));
    printf("--- hi1d_xval, yval, zval ---\n");
    hi1d_xval->Print(stdout);
    hi1d_yval->Print(stdout);
    printf("=== hi1d_xval, yval, zval ===\n");



    for(long ivel = 0; ivel < nbin_vel; ivel ++){
        double theta =  atan(hi1d_vel->GetBinCenter(ivel) /
                             (time_arr[ntime - 1] * scale) );
        printf("theta = %e (deg)\n", theta / M_PI * 180.0);
        long nbin_par = nbin_rho * nbin_phi * nbin_psi;
        double* par3d_arr = new double [nbin_par];
        for(long iarr = 0; iarr < nbin_par; iarr ++){
            par3d_arr[iarr] = 0.0;
        }
        for(long itime = 0; itime < ntime; itime ++){
            for(long ixpos = 0; ixpos < hi1d_xval->GetNbin(); ixpos ++){
                for(long iypos = 0; iypos < hi1d_yval->GetNbin(); iypos ++){
                    double xval = hi1d_xval->GetBinCenter(ixpos);
                    double yval = hi1d_yval->GetBinCenter(iypos);
                    // double zval = time_arr[itime];
                    double zval = zpos_arr[itime];
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
                        if(rho < hi1d_rho->GetLo() || rho > hi1d_rho->GetUp()){
                            printf("continue: rho = %e, phi = %e\n", rho, phi);
                            continue;
                        }
                        long irho = hi1d_rho->GetIbin(rho);
                        long iphi = hi1d_phi->GetIbin(phi);
                        long i_rho_phi_psi = irho + iphi * nbin_rho
                            + ipsi * (nbin_rho * nbin_phi);
                        par3d_arr[i_rho_phi_psi] ++;
                    }
                }
            }
        }
        char tag[kLineSize];
        int bitpix = -32;
        sprintf(tag, "vel_%2.2ld", ivel);
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              3, bitpix,
                              img_info_par3d->GetNaxesArr(),
                              par3d_arr);
        // theory
        double* par3d_theory_arr = new double [nbin_par];
        for(long iarr = 0; iarr < nbin_par; iarr ++){
            par3d_theory_arr[iarr] = 0.0;
        }
        for(long itime = 0; itime < ntime; itime ++){
            // double zval = time_arr[itime];
            double zval = zpos_arr[itime];
            for(long ipsi = 0; ipsi < nbin_psi; ipsi ++){
                double psi = hi1d_psi->GetBinCenter(ipsi);
                for(long irho = 0; irho < nbin_rho; irho ++){
                    double rho = hi1d_rho->GetBinCenter(irho);
                    for(long iphi = 0; iphi < nbin_phi; iphi ++){
                        double phi = hi1d_phi->GetBinCenter(iphi);
                        double xval = rho * ( cos(phi) * cos(psi) - sin(phi) * sin(psi) / cos(theta) )
                            + zval * tan(theta) * sin(psi);
                        double yval = rho * ( cos(phi) * sin(psi) + sin(phi) * cos(psi) / cos(theta) )
                            - zval * tan(theta) * cos(psi);
                        if(hi1d_xval->GetLo() < xval && xval < hi1d_xval->GetUp() &&
                           hi1d_yval->GetLo() < yval && yval < hi1d_yval->GetUp()){
                            long i_rho_phi_psi = irho + iphi * nbin_rho
                                + ipsi * (nbin_rho * nbin_phi);                        
                            par3d_theory_arr[i_rho_phi_psi] += rho / theta;
                        }
                    }
                }
            }
        }
        sprintf(tag, "vel_th_%2.2ld", ivel);
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              3, bitpix,
                              img_info_par3d->GetNaxesArr(),
                              par3d_theory_arr);

        double* par3d_ratio_arr = new double [nbin_par];
        for(long iarr = 0; iarr < nbin_par; iarr ++){
            par3d_ratio_arr[iarr] = 0.0;
        }

        for(long iarr = 0; iarr < nbin_par; iarr ++){
            if(par3d_theory_arr[iarr] > 1e-10){
                par3d_ratio_arr[iarr] = par3d_arr[iarr] / par3d_theory_arr[iarr];
            }
        }
        sprintf(tag, "vel_ratio_%2.2ld", ivel);
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              3, bitpix,
                              img_info_par3d->GetNaxesArr(),
                              par3d_ratio_arr);

        delete [] par3d_arr;
        delete [] par3d_theory_arr;
        delete [] par3d_ratio_arr;
    }
    

    double time_ed = MiTime::GetTimeSec();
    printf("time_ed - time_st = %e\n", time_ed - time_st);
    
    return status_prog;
}
