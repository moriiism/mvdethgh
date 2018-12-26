#ifndef MORIIISM_MVDETHGH_SUB_MVDETHGH_H_
#define MORIIISM_MVDETHGH_SUB_MVDETHGH_H_

#include "mi_rand.h"
#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"
#include "mir_hist_info.h"

void LoadHi1dVel(string vel_dat,
                 HistInfo1d* const hi1d_vel);

void LoadHi1dPar(string res_dat,
                 const MifImgInfo* const img_info_subimg,
                 HistInfo1d* const hi1d_rho,
                 HistInfo1d* const hi1d_phi,
                 HistInfo1d* const hi1d_psi);

void LoadHi1dTime(string time_dat,
                  HistInfo1d* const hi1d_time);

void GetMeanStddevClip(long narr, const double* const val_arr,
                       int nclip, double significance,
                       double* const mean_ptr,
                       double* const stddev_ptr);

#endif // MORIIISM_MVDETHGH_SUB_MVDETHGH_H_
