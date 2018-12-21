#include "sub_mvdethgh.h"

void LoadParDat(string par_dat, const MifImgInfo* const img_info_in,
                HistInfo1d* const hi1d_rho,
                HistInfo1d* const hi1d_phi,
                HistInfo1d* const hi1d_theta,
                HistInfo1d* const hi1d_psi)
{
    printf("--- read par.dat ---\n");
    string* line_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(par_dat, &line_arr, &nline);
    if(4 != nline){
        printf("nline(%ld) != 4, then abort.", nline);
        abort();
    }
    long nbin_rho = 0;
    double rho_lo = 0.0;
    char dummy[kLineSize];
    long nbin_phi = 0;
    double phi_lo = 0.0;
    double phi_up = 0.0;
    long nbin_theta = 0;
    double theta_lo = 0.0;
    double theta_up = 0.0;
    long nbin_psi = 0;
    double psi_lo = 0.0;
    double psi_up = 0.0;
    
    sscanf(line_arr[0].c_str(), "%ld  %lf  %s", &nbin_rho, &rho_lo, dummy);
    sscanf(line_arr[1].c_str(), "%ld  %lf  %lf", &nbin_phi, &phi_lo, &phi_up);
    sscanf(line_arr[2].c_str(), "%ld  %lf  %lf", &nbin_theta, &theta_lo, &theta_up);
    sscanf(line_arr[3].c_str(), "%ld  %lf  %lf", &nbin_psi, &psi_lo, &psi_up);

    double rho_up = sqrt( pow(img_info_in->GetNaxesArrElm(0), 2)
                          + pow(img_info_in->GetNaxesArrElm(1), 2) );
    hi1d_rho->InitSetByNbin(rho_lo, rho_up, nbin_rho);
    hi1d_phi->InitSetByNbin(phi_lo, phi_up, nbin_phi);
    hi1d_theta->InitSetByNbin(theta_lo, theta_up, nbin_theta);
    hi1d_psi->InitSetByNbin(psi_lo, psi_up, nbin_psi);

    delete [] line_arr;    
}
