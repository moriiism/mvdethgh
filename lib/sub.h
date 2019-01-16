#ifndef MORIIISM_MVDETHGH_LIB_SUB_H_
#define MORIIISM_MVDETHGH_LIB_SUB_H_

#include "mi_base.h"

void OpenLogfile(string outdir,
                 string outfile_head,
                 string progname,
                 FILE** fp_log_ptr);

void ReadDataList(string data_list,
                  string subimg_dat,
                  long* const ntime_ptr,
                  double** time_arr_ptr,
                  double*** data_arr_ptr,
                  MifImgInfo** img_info_subimg_ptr,
                  int* const bitpix_ptr);

void MkStdImgArr(long ntime,
                 const double* const* const data_img_arr,
                 const MifImgInfo* const img_info_subimg,
                 int nclip, double significance_clip,
                 double*** std_img_arr_ptr,
                 double*** debias_img_arr_ptr,
                 double** mean_time_arr_ptr,
                 double** stddev_time_arr_ptr);

void MkMedianImg(long ntime,
                 const double* const* const debias_img_arr,
                 const MifImgInfo* const img_info_subimg,
                 double** median_img_arr_ptr);

void GetPeakXY(const MifImgInfo* const img_info_subimg,
               const double* const median_img_arr,
               string func_name, string func_par_name, int nbin_kernel_half,
               double* xval_peak_ptr, double* yval_peak_ptr);

void MkPsf(const MifImgInfo* const img_info_subimg,
           const double* const median_img_arr,
           int nbin_psf_half,
           double xval_peak, double yval_peak,
           HistDataNerr2d** hd2d_psf_ptr);

#endif // MORIIISM_MVDETHGH_LIB_SUB_H_
