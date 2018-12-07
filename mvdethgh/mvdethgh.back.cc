#include "mir_math.h"
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


    string* line_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(argval->GetDatalist(), &line_arr, &nline);
    string* fitsfile_arr = new string[nline];
    double* time_arr     = new double[nline];
    int ntime = nline;
    for(long iline = 0; iline < nline; iline ++){
        int nsplit = 0;
        string* split_arr = NULL;
        MiStr::GenSplit(line_arr[iline], &nsplit, &split_arr);
        fitsfile_arr[iline] = split_arr[0];
        time_arr[iline] = atof(split_arr[1].c_str());
        MiStr::DelSplit(split_arr);
    }

    // show datalist
    for(int itime = 0; itime < ntime; itime ++){
        printf("%s  %e\n", fitsfile_arr[itime].c_str(), time_arr[itime]);
    }
    
    printf("--- img_info_in ---\n");
    MifImgInfo* img_info_in = new MifImgInfo;
    img_info_in->Load(argval->GetImgInfoDat());
    img_info_in->PrintInfo();
    printf("=== img_info_in ===\n");

    double** data_arr = new double* [ntime];
    int bitpix = 0;
    for(int itime = 0; itime < ntime; itime ++){
        MifFits::InFitsImageD(fitsfile_arr[itime], img_info_in,
                              &bitpix, &data_arr[itime]);
    }

    // theta : 0<= theta <= pi : ntheta
    // phi   : 0<= phi <= 2 pi : nphi
    // psi   : 0<= psi <= pi
    // rho   : 0<= rho <= 50 

    long ntheta = 100;
    long nphi   = 100;
    long nrho   = 100;
    long npsi   = 100;
    double delta_theta = M_PI / ntheta;
    double delta_phi   = 2 * M_PI / nphi;
    double delta_rho   = (double) img_info_in->GetNaxesArrElm(0) / (double) nrho;
    double delta_psi   = M_PI / npsi;

    printf("delta_rho = %e, delta_psi = %e\n", delta_rho, delta_psi);

    MifImgInfo* img_info_rho = new MifImgInfo;
    img_info_rho->InitSetCube(1, 1, 1, ntheta, nphi, nrho);
    img_info_rho->PrintInfo();

    MifImgInfo* img_info_psi = new MifImgInfo;
    img_info_psi->InitSetCube(1, 1, 1, ntheta, nphi, npsi);
    img_info_psi->PrintInfo();
    
    double* rho_arr = new double [ntheta * nphi * nrho];
    double* psi_arr = new double [ntheta * nphi * npsi];

    double posx_lo = 0.0;
    double posy_lo = 0.0;
    double posz_lo = 0.0;
    
    long nposx = img_info_in->GetNaxesArrElm(0);
    long nposy = img_info_in->GetNaxesArrElm(1);

    printf("ntime = %d\n", ntime);
    printf("nposx = %ld\n", nposx);
    printf("nposy = %ld\n", nposy);
    
    for(int itime = 0; itime < ntime; itime ++){
        for(long iposx = 0; iposx < nposx; iposx ++){
            for(long iposy = 0; iposy < nposy; iposy ++){
                double xval = iposx;
                double yval = iposy;
                double zval = itime;

                if(data_arr[itime][iposx + iposy * nposx] < 1e-10){
                    continue;
                }

                for(long itheta = 0; itheta < ntheta; itheta ++){
                    for(long iphi = 0; iphi < nphi; iphi ++){
                        double theta = delta_theta * itheta;
                        double phi   = delta_phi * iphi;
                        double rho_cos_psi
                            = xval * cos(theta) * cos(phi)
                            + yval * cos(theta) * sin(phi)
                            - zval * sin(theta);
                        double rho_sin_psi
                            = -1 * xval * sin(phi) + yval * cos(phi);
                        double rho = sqrt(rho_cos_psi * rho_cos_psi + rho_sin_psi * rho_sin_psi);

                        double psi = 0.0;
                        if(rho_cos_psi >= 0 && rho_sin_psi >= 0){
                            psi = atan( rho_sin_psi / rho_cos_psi);
                        } else if(rho_cos_psi < 0 && rho_sin_psi >= 0) {
                            psi = atan( -1 * rho_sin_psi / rho_cos_psi );
                        }
                        

                        long irho = (long) floor( rho / delta_rho );
                        long i_theta_phi_rho = itheta + iphi * ntheta + irho * (ntheta * nphi);
                        rho_arr[i_theta_phi_rho] += data_arr[itime][iposx + iposy * nposx];

                        long ipsi = (long) floor( ( psi + M_PI/2 ) / delta_psi );
                        long i_theta_phi_psi = itheta + iphi * ntheta + ipsi * (ntheta * nphi);
                        psi_arr[i_theta_phi_psi] += data_arr[itime][iposx + iposy * nposx];
                    }
                }
            }
        }
        printf("itime = %d\n", itime);
    }

    MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), "rho",
                          3, bitpix,
                          img_info_rho->GetNaxesArr(),
                          rho_arr);
    MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), "psi",
                          3, bitpix,
                          img_info_psi->GetNaxesArr(),
                          psi_arr);
    
    for(int itime = 0; itime < ntime; itime ++){
        char tag[kLineSize];
        sprintf(tag, "%2.2d", itime);
        int naxis = 2;
        MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              naxis, bitpix,
                              img_info_in->GetNaxesArr(),
                              data_arr[itime]);
    }
    
    return status_prog;
}
