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
    
