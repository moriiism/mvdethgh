#ifndef MORIIISM_MVDETHGH_SUB_MVDETHGH_H_
#define MORIIISM_MVDETHGH_SUB_MVDETHGH_H_

#include "mi_rand.h"
#include "mib_blas.h"
#include "mi_sort.h"
#include "mi_time.h"
#include "mif_fits.h"
#include "mir_math.h"
#include "mir_hist_info.h"

void LoadParDat(string par_dat, const MifImgInfo* const img_info_in,
                HistInfo1d* const hi1d_rho,
                HistInfo1d* const hi1d_phi,
                HistInfo1d* const hi1d_theta,
                HistInfo1d* const hi1d_psi);

#endif // MORIIISM_MVDETHGH_SUB_MVDETHGH_H_
