#ifndef MORIIISM_MVDETHGH_LIB_MVDETHGHLIB_H_
#define MORIIISM_MVDETHGH_LIB_MVDETHGHLIB_H_

#include "mir_math_util.h"
#include "mir_hist2d_ope.h"
#include "mif_fits.h"
#include "mifc_gen.h"

void OpenLogfile(string outdir,
                 string outfile_head,
                 string progname,
                 FILE** const fp_log_ptr);

void LoadData(string data_list,
              string subimg_dat,
              long* const ntime_ptr,
              double** const time_arr_ptr,
              double*** const data_arr_ptr,
              MifImgInfo** const img_info_ptr,
              int* const bitpix_ptr);

void LoadDataXYT(string data_list,
                 long* const ntime_ptr,
                 double** const time_arr_ptr,
                 long** const npos_arr_ptr,
                 double*** const xpos_arr_ptr,
                 double*** const ypos_arr_ptr,
                 double* const xpos_lo_ptr,
                 double* const xpos_up_ptr,
                 double* const ypos_lo_ptr,
                 double* const ypos_up_ptr,
                 double* const rho_up_ptr);

void GenStdImgArr(long ntime,
                  const double* const* const data_img_arr,
                  const MifImgInfo* const img_info,
                  int nclip, double significance_clip,
                  double*** const std_img_arr_ptr,
                  double*** const debias_img_arr_ptr,
                  double** const mean_time_arr_ptr,
                  double** const stddev_time_arr_ptr);

void GenMedianImg(long ntime,
                  const double* const* const debias_img_arr,
                  const MifImgInfo* const img_info,
                  double** const median_img_arr_ptr);

void GenMedianImgByNthElement(long ntime,
                              const double* const* const debias_img_arr,
                              const MifImgInfo* const img_info,
                              double** const median_img_arr_ptr);

double GetMedianByNthElement(long narr,
                             const double* const val_arr);

void GenImgAboveThArr(long ntime,
                      const double* const* const data_img_arr,
                      const double* const* const std_img_arr,
                      const MifImgInfo* const img_info,
                      double threshold,
                      double*** const out_img_arr_ptr,
                      long* const ndet_ptr);

void GetPeakXY(const MifImgInfo* const img_info,
               const double* const median_img_arr,
               int nbin_kernel_half, double val_smooth,
               double* const xval_peak_ptr, double* const yval_peak_ptr);

void GenPsf(const MifImgInfo* const img_info,
            const double* const median_img_arr,
            int nbin_psf_half,
            double xval_peak, double yval_peak,
            HistDataNerr2d** const hd2d_psf_ptr);

void SavePsf(const HistDataNerr2d* const hd2d_psf,
             int bitpix,
             string outdir, string outfile_head);

void SaveCube(long ntime, const MifImgInfo* const img_info,
              const double* const* const data_img_arr,
              int bitpix, string outdir, string outfile_head, string tag);

void GenMvobjImgArr(long ntime,
                    const MifImgInfo* const img_info,
                    const double* const* const debias_img_arr,
                    const double* const median_img_arr,
                    double*** const mvobj_img_arr_ptr);

void SaveImgArr(long ntime, const MifImgInfo* const img_info,
                const double* const* const data_img_arr,
                int bitpix, string outdir, string outfile_head, string tag);

void GenConvPsf(long ntime,
                const MifImgInfo* const img_info,
                const double* const* const mvobj_img_arr,
                string psf_dat,
                double*** const conv_img_arr_ptr);

void LoadHi1dVel(string vel_dat,
                 HistInfo1d* const hi1d_vel);

void LoadHi1dPar(string res_dat,
                 const MifImgInfo* const img_info,
                 HistInfo1d* const hi1d_rho,
                 HistInfo1d* const hi1d_phi,
                 HistInfo1d* const hi1d_psi);

void LoadHi1dPar(string res_dat,
                 double rho_up,
                 HistInfo1d* const hi1d_rho,
                 HistInfo1d* const hi1d_phi,
                 HistInfo1d* const hi1d_psi);

void LoadHi1dTime(string time_dat,
                  HistInfo1d* const hi1d_time);

void GenDetImg(const HistDataNerr2d* const* const hd2d_arr,
               const double* const time_arr, long ntime,
               int nbin_detimg_half,
               double theta, double rho, double phi, double psi,
               HistDataNerr2d** const hd2d_img_ptr);

void GenDetImg(const HistDataNerr2d* const* const hd2d_arr,
               const double* const time_arr, long ntime,
               double vel, double psi,
               HistDataNerr2d** const hd2d_img_ptr);

#endif // MORIIISM_MVDETHGH_LIB_MVDETHGHLIB_H_
