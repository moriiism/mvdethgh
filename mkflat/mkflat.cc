#include "mi_time.h"
#include "mvdethghlib.h"
#include "arg_mkflat.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValMkflat* argval = new ArgValMkflat;
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
    double* cube_flat_vel_zero_arr = new double [nbin_rho * nbin_phi * nbin_psi];
    for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
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
                    cube_flat_vel_zero_arr[i_rho_phi_psi] ++;
                }
            }
        }
    }
    char tag[kLineSize];
    sprintf(tag, "vel_zero");
    MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                          3, bitpix,
                          img_info_par_cube->GetNaxesArr(),
                          cube_flat_vel_zero_arr);
    ///////
    
    double* par_cube_flat_arr = new double [nbin_rho * nbin_phi * nbin_psi];
    for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
        par_cube_flat_arr[iarr] = 0.0;
    }
    for(long ivel = 0; ivel < nbin_vel; ivel ++){
        double theta =  atan(hi1d_vel->GetBinCenter(ivel));
        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
            par_cube_flat_arr[iarr] = 0.0;
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
                        par_cube_flat_arr[i_rho_phi_psi] ++;
                    }
                }
            }
        }
        sprintf(tag, "vel_%2.2ld", ivel);
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              3, bitpix,
                              img_info_par_cube->GetNaxesArr(),
                              par_cube_flat_arr);
    }


    double time_ed = MiTime::GetTimeSec();
    printf("time_ed - time_st = %e\n", time_ed - time_st);
    
    return status_prog;
}
