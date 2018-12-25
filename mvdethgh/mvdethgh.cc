#include "mir_math.h"
#include "mir_data1d_ope.h"
#include "mir_hist_info.h"
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

    printf("--- read data list ---\n");
    string* line_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(argval->GetDatalist(), &line_arr, &nline);
    string* fitsfile_arr = new string[nline];
    double* time_arr     = new double[nline];
    for(long iline = 0; iline < nline; iline ++){
        int nsplit = 0;
        string* split_arr = NULL;
        MiStr::GenSplit(line_arr[iline], &nsplit, &split_arr);
        fitsfile_arr[iline] = split_arr[0];
        time_arr[iline] = atof(split_arr[1].c_str());
        MiStr::DelSplit(split_arr);
    }
    delete [] line_arr;
    for(int iline = 0; iline < nline; iline ++){
        printf("%s  %e\n", fitsfile_arr[iline].c_str(), time_arr[iline]);
    }
    printf("=== read data list ===\n");
   
    printf("--- img_info_in ---\n");
    MifImgInfo* img_info_in = new MifImgInfo;
    img_info_in->Load(argval->GetImgInfoDat());
    img_info_in->PrintInfo();
    printf("=== img_info_in ===\n");

    printf("--- read 2d images ---\n");
    long ntime = nline;
    double** data_arr = new double* [ntime];
    int bitpix = 0;
    for(int itime = 0; itime < ntime; itime ++){
        MifFits::InFitsImageD(fitsfile_arr[itime], img_info_in,
                              &bitpix, &data_arr[itime]);
    }
    printf("bitpix = %d\n", bitpix);
    printf("=== read 2d images ===\n");

    printf("--- Init output cube ---\n");
    MifImgInfo* img_info_rec = new MifImgInfo;
    img_info_rec->InitSetCube(1, 1, 1, img_info_in->GetNaxesArrElm(0),
                              img_info_in->GetNaxesArrElm(1), ntime);
    img_info_rec->PrintInfo();
    printf("=== Init output cube ===\n");

    printf("--- Output subcube ---\n");
    double* data_subcube_arr = new double [img_info_rec->GetNpixelTotal()];
    for(long itime = 0; itime < ntime; itime ++){
        for(long iarr = 0; iarr < img_info_rec->GetNpixelImg(); iarr ++){
            long index = iarr + itime * img_info_rec->GetNpixelImg();
            if(data_arr[itime][iarr] < 1000){
                continue;
            }
            data_subcube_arr[index] = data_arr[itime][iarr];
        }
    }
    MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), "subcube",
                          3, bitpix,
                          img_info_rec->GetNaxesArr(),
                          data_subcube_arr);
    printf("=== Output subcube ===\n");


    printf("--- statistic of input data_arr ---\n");
    DataArray1d** da1d_arr = new DataArray1d* [ntime];
    for(long itime = 0; itime < ntime; itime ++){
        da1d_arr[itime] = new DataArrayNerr1d;
        da1d_arr[itime]->Init(img_info_rec->GetNpixelImg());
        da1d_arr[itime]->SetVal(img_info_rec->GetNpixelImg(), data_arr[itime]);
    }
    DataArrayNerr1d* da1d_amean = new DataArrayNerr1d;
    DataArray1dOpe::GetAMean(da1d_arr, ntime, da1d_amean);
    DataArrayNerr1d* da1d_stddev = new DataArrayNerr1d;
    DataArray1dOpe::GetSqrtOfUnbiasedVariance(da1d_arr, ntime, da1d_stddev);
    
    printf("=== statistic of input data_arr ===\n");





    
    HistInfo1d* hi1d_rho   = new HistInfo1d;
    HistInfo1d* hi1d_phi   = new HistInfo1d;
    HistInfo1d* hi1d_theta = new HistInfo1d;
    HistInfo1d* hi1d_psi   = new HistInfo1d;
    LoadParDat(argval->GetParDat(), img_info_in,
               hi1d_rho, hi1d_phi, hi1d_theta, hi1d_psi);
    double rho_lo      = hi1d_rho->GetLo();
    double rho_up      = hi1d_rho->GetUp();
    long   nbin_rho    = hi1d_rho->GetNbin();
    double delta_rho   = hi1d_rho->GetBinWidth();
    
    double phi_lo      = hi1d_phi->GetLo();
    double phi_up      = hi1d_phi->GetUp();
    long   nbin_phi    = hi1d_phi->GetNbin();
    double delta_phi   = hi1d_phi->GetBinWidth();
    
    double theta_lo    = hi1d_theta->GetLo();
    double theta_up    = hi1d_theta->GetUp();
    long   nbin_theta  = hi1d_theta->GetNbin();
    double delta_theta = hi1d_theta->GetBinWidth();
    
    double psi_lo      = hi1d_psi->GetLo();
    double psi_up      = hi1d_psi->GetUp();
    long   nbin_psi    = hi1d_psi->GetNbin();
    double delta_psi   = hi1d_psi->GetBinWidth();

    printf("rho_lo, rho_up, nbin_rho, delta_rho = %e %e %ld %e\n",
           rho_lo, rho_up, nbin_rho, delta_rho);
    printf("phi_lo, phi_up, nbin_phi, delta_phi = %e %e %ld %e\n",
           phi_lo, phi_up, nbin_phi, delta_phi);
    printf("theta_lo, theta_up, nbin_theta, delta_theta = %e %e %ld %e\n",
           theta_lo, theta_up, nbin_theta, delta_theta);
    printf("psi_lo, psi_up, nbin_psi, delta_psi = %e %e %ld %e\n",
           psi_lo, psi_up, nbin_psi, delta_psi);        
    
    MifImgInfo* img_info_cube = new MifImgInfo;
    img_info_cube->InitSetCube(1, 1, 1, nbin_rho, nbin_phi, nbin_psi);
    img_info_cube->PrintInfo();

    double posx_lo = 0.0;
    double posy_lo = 0.0;
    double posz_lo = 0.0;
    long nposx = img_info_rec->GetNaxesArrElm(0);
    long nposy = img_info_rec->GetNaxesArrElm(1);
    long nposz = img_info_rec->GetNaxesArrElm(2);
    printf("nposx = %ld\n", nposx);
    printf("nposy = %ld\n", nposy);
    printf("nposz = %ld\n", nposz);

    double* hough_max_theta_arr = new double [nbin_theta];
    long*   hough_max_index_arr = new long   [nbin_theta];
    for(long itheta = 0; itheta < nbin_theta; itheta ++){
        hough_max_theta_arr[itheta] = 0.0;
        hough_max_index_arr[itheta] = 0;
    }
    double* cube_arr = new double [nbin_rho * nbin_phi * nbin_psi];
    for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
        cube_arr[iarr] = 0.0;
    }
    double* flat_arr = new double [nbin_rho * nbin_phi * nbin_psi];
    for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
        flat_arr[iarr] = 0.0;
    }
    double* rec_arr = new double [img_info_rec->GetNpixelTotal()];
    for(long iarr = 0; iarr < img_info_rec->GetNpixelTotal(); iarr++){
        rec_arr[iarr] = 0.0;
    }
    for(long itheta = 0; itheta < nbin_theta; itheta ++){
        double theta_rad = (theta_lo + delta_theta * (itheta + 0.5) ) / 180.0 * M_PI ;
        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
            cube_arr[iarr] = 0.0;
            flat_arr[iarr] = 0.0;
        }
        for(long iarr = 0; iarr < img_info_rec->GetNpixelTotal(); iarr++){
            rec_arr[iarr] = 0.0;
        }
        for(long itime = 0; itime < ntime; itime ++){
            for(long iposx = 0; iposx < nposx; iposx ++){
                for(long iposy = 0; iposy < nposy; iposy ++){
                    double xval = 1.0 * (iposx + 0.5);
                    double yval = 1.0 * (iposy + 0.5);
                    double zval = time_arr[itime];

                    if(data_arr[itime][iposx + iposy * nposx] < 1000){
                        continue;
                    }

                    for(long ipsi = 0; ipsi < nbin_psi; ipsi ++){
                        double psi_rad = (psi_lo + delta_psi * (ipsi + 0.5) ) / 180. * M_PI;
                        double rho_cos_phi = xval * cos(psi_rad) + yval * sin(psi_rad);
                        double rho_sin_phi = -1 * xval * sin(psi_rad) * cos(theta_rad)
                            + yval * cos(psi_rad) * cos(theta_rad)
                            + zval * sin(theta_rad);
                        double rho = sqrt(pow(rho_cos_phi, 2) + pow(rho_sin_phi, 2));
                        double phi_rad = 0.0;
                        if(rho_sin_phi >= 0.0){
                            phi_rad = acos(rho_cos_phi / rho);
                        } else {
                            phi_rad = -1 * acos(rho_cos_phi / rho) + 2 * M_PI;
                        }
                        long irho = (long) floor( rho / delta_rho );
                        long iphi = (long) floor( phi_rad / M_PI * 180.0 / delta_phi );
                        if(irho < 0 || nbin_rho <= irho){
                            printf("!!!! irho = %ld\n", irho);
                        }
                        if(iphi < 0 || nbin_phi <= iphi){
                            printf("!!!! iphi = %ld phi(deg) = %e\n", iphi, phi_rad/M_PI * 180.);
                        }
                        long i_rho_phi_psi = irho + iphi * nbin_rho + ipsi * (nbin_rho * nbin_phi);
                        if(data_arr[itime][iposx + iposy * nposx] >= 1000){
                            cube_arr[i_rho_phi_psi] += 1.0;
                        }
                        // cube_arr[i_rho_phi_psi] += data_arr[itime][iposx + iposy * nposx];
                    }
                }
            }
        }

        char tag[kLineSize];
        sprintf(tag, "theta_%2.2ld", itheta);
        int naxis = 2;
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              3, bitpix,
                              img_info_cube->GetNaxesArr(),
                              cube_arr);
        sprintf(tag, "theta_flat_%2.2ld", itheta);
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              3, bitpix,
                              img_info_cube->GetNaxesArr(),
                              flat_arr);
        
        double hough_max = 0.0;
        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
            if(cube_arr[iarr] > hough_max){
                hough_max = cube_arr[iarr];
                hough_max_theta_arr[itheta] = cube_arr[iarr];
                hough_max_index_arr[itheta] = iarr;
            }
        }

        long ipsi = hough_max_index_arr[itheta] / (nbin_rho * nbin_phi);
        long index2 = hough_max_index_arr[itheta] % (nbin_rho * nbin_phi);
        long iphi = index2 / nbin_rho;
        long irho = index2 % nbin_rho;

        double rho     = rho_lo + (irho + 0.5) * delta_rho;
        double phi_deg = phi_lo + (iphi + 0.5) * delta_phi;
        double psi_deg = psi_lo + (ipsi + 0.5) * delta_psi;
        double phi_rad = phi_deg / 180.0 * M_PI;
        double psi_rad = psi_deg / 180.0 * M_PI;

        double theta_deg = theta_rad / M_PI * 180.0;
        printf("itheta = %ld, theta_deg = %e, rho = %e, phi_deg = %e, psi_deg = %e \n",
               itheta, theta_deg, rho, phi_deg, psi_deg);

        int fill_flag = 0;
        double zval_lo = time_arr[0];
        double zval_up = time_arr[ntime - 1];
        double tpar_lo = ( zval_lo - rho * sin(theta_rad) * sin(phi_rad) ) / cos(theta_rad );
        double tpar_up = ( zval_up - rho * sin(theta_rad) * sin(phi_rad) ) / cos(theta_rad );
        long ntpar = 100;
        double delta_tpar = (tpar_up - tpar_lo) / ntpar;
        for(long itpar = 0; itpar < ntpar; itpar ++){
            double tpar = tpar_lo + itpar * delta_tpar;
            double xval = rho * (cos(psi_rad) * cos(phi_rad)
                                 - sin(psi_rad) * cos(theta_rad) * sin(phi_rad) )
                + tpar * sin(psi_rad) * sin(theta_rad);
            double yval = rho * (sin(psi_rad) * cos(phi_rad)
                                 + cos(psi_rad) * cos(theta_rad) * sin(phi_rad) )
                - tpar * cos(psi_rad) * sin(theta_rad);
            double zval = rho * sin(theta_rad) * sin(phi_rad) + tpar * cos(theta_rad);

            long iposx = (long) xval;
            long iposy = (long) yval;
            long iposz = (long) zval;
            long index = iposx + iposy * nposx + iposz * (nposx * nposy);

            double val11 = cos(psi_rad) * cos(phi_rad);
            double val12 = - sin(psi_rad) * cos(theta_rad) * sin(phi_rad);
            double val2 = sin(psi_rad) * sin(theta_rad);

            
            printf("val11, val12, val2 %f, %f, %f, tpar = %f, xval, yval, zval = %f, %f, %f\n",
                   val11, val12, val2, tpar, xval, yval, zval);
            
            if(iposx < 0 || nposx - 1 < iposx){
                continue;
            }
            if(iposy < 0 || nposy - 1 < iposy){
                continue;
            }
            if(iposz < 0 || nposz - 1 < iposz){
                continue;
            }
            rec_arr[index] = 1.0;
            fill_flag = 1;
        }
        
        sprintf(tag, "rec_%2.2ld_%d", itheta, fill_flag);
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              3, bitpix,
                              img_info_rec->GetNaxesArr(),
                              rec_arr);
    }

    FILE* fp_out = fopen("out.qdp", "w");
    for(long itheta = 0; itheta < nbin_theta; itheta++){
        fprintf(fp_out, "%ld %e\n", itheta, hough_max_theta_arr[itheta]);
    }
    fclose(fp_out);


//    {
//        // fake
//
//        // theta = 4.600000e+01, rho = 8.485281e+01, phi_deg = 2.916000e+02, psi_deg = 1.620000e+02 
//
//        double rho_fake = 100;
//        double phi_deg_fake = 275.0;
//        double theta_deg_fake = 45.0;
//        double psi_deg_fake = 180.0;
//        double phi_rad_fake = phi_deg_fake / 180. * M_PI;
//        double theta_rad_fake = theta_deg_fake / 180. * M_PI;
//        double psi_rad_fake = psi_deg_fake / 180. * M_PI;
//    
//        double* fake_arr = new double [img_info_rec->GetNpixelTotal()];
//        for(long iarr = 0; iarr < img_info_rec->GetNpixelTotal(); iarr++){
//            fake_arr[iarr] = 0.0;
//        }
//        double zval_lo = time_arr[0];
//        double zval_up = time_arr[ntime - 1];
//        double tpar_lo = ( zval_lo - rho_fake * sin(theta_rad_fake) * sin(phi_rad_fake) )
//            / cos(theta_rad_fake );
//        double tpar_up = ( zval_up - rho_fake * sin(theta_rad_fake) * sin(phi_rad_fake) )
//            / cos(theta_rad_fake );
//        long ntpar = 100;
//        double delta_tpar = (tpar_up - tpar_lo) / ntpar;
//        for(long itpar = 0; itpar < ntpar; itpar ++){
//            double tpar = tpar_lo + itpar * delta_tpar;
//            double xval = rho_fake * (cos(psi_rad_fake) * cos(phi_rad_fake)
//                                      - sin(psi_rad_fake) * cos(theta_rad_fake) * sin(phi_rad_fake) )
//                + tpar * sin(psi_rad_fake) * sin(theta_rad_fake);
//            double yval = rho_fake * (sin(psi_rad_fake) * cos(phi_rad_fake)
//                                      + cos(psi_rad_fake) * cos(theta_rad_fake) * sin(phi_rad_fake) )
//                - tpar * cos(psi_rad_fake) * sin(theta_rad_fake);
//            double zval = rho_fake * sin(theta_rad_fake) * sin(phi_rad_fake)
//                + tpar * cos(theta_rad_fake);
//            long iposx = (long) xval;
//            long iposy = (long) yval;
//            long iposz = (long) zval;
//            long index = iposx + iposy * nposx + iposz * (nposx * nposy);
//            if(iposx < 0 || nposx - 1 < iposx){
//                continue;
//            }
//            if(iposy < 0 || nposy - 1 < iposy){
//                continue;
//            }
//            if(iposz < 0 || nposz - 1 < iposz){
//                continue;
//            }
//            fake_arr[index] = 1.0;
//        }
//        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), "fake",
//                              3, bitpix,
//                              img_info_rec->GetNaxesArr(),
//                              fake_arr);
//    }    
//
//






    
    
    return status_prog;
}
