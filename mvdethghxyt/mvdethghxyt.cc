#include "mi_time.h"
#include "mvdethghlib.h"
#include "arg_mvdethghxyt.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValMvdethghxyt* argval = new ArgValMvdethghxyt;
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

    // for this CCD chip
    xpos_lo = -1024.5;
    xpos_up = +1024.5;
    ypos_lo = -2088.5;
    ypos_up = +2088.5;

   
    printf("--- load parameters ---\n");
    HistInfo1d* hi1d_vel   = new HistInfo1d;
    HistInfo1d* hi1d_rho   = new HistInfo1d;
    HistInfo1d* hi1d_phi   = new HistInfo1d;
    HistInfo1d* hi1d_psi   = new HistInfo1d;
    LoadHi1dVel(argval->GetVelDat(), hi1d_vel);
    
    // scale zpos
    double* zpos_arr = new double[ntime];
    double scale = hi1d_vel->GetUp() / time_arr[ntime - 1] * 10.0;
    printf("scale = %e\n", scale);
    for(long itime = 0; itime < ntime; itime ++){
        zpos_arr[itime] = time_arr[itime] * scale;
        printf("%ld  %e  %e\n", itime, time_arr[itime], zpos_arr[itime]);
    }


    rho_up = 0.0;
    for(long ipsi = 0; ipsi < 100; ipsi ++){
        double theta =  atan(hi1d_vel->GetBinCenter(0) /
                             (time_arr[ntime - 1] * scale) );
        double psi = 2.0 * M_PI / 100.0 * ipsi;
        double rho_cos_phi = xpos_up * cos(psi) + ypos_up * sin(psi);
        double rho_sin_phi = -1 * xpos_up * sin(psi) * cos(theta)
            + ypos_up * cos(psi) * cos(theta)
            + zpos_arr[ntime - 1] * sin(theta);
        double rho_up_this = sqrt(pow(rho_cos_phi, 2) + pow(rho_sin_phi, 2));
        if(rho_up_this > rho_up){
            rho_up = rho_up_this;
        }
    }
    printf("rho_up = %e\n", rho_up);
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


    // plot input sources
    MifImgInfo* img_info_3d = new MifImgInfo;
    img_info_3d->InitSetCube(1, 1, 1, 
                             hi1d_xval->GetNbin(),
                             hi1d_yval->GetNbin(),
                             ntime);
    img_info_3d->PrintInfo();
    double* cube_arr = new double[img_info_3d->GetNpixelTotal()];
    for(long iarr = 0; iarr < img_info_3d->GetNpixelTotal(); iarr++){
        cube_arr[iarr] = 0.0;
    }
    HistDataNerr2d** hd2d_arr = new HistDataNerr2d* [ntime];
    for(long itime = 0; itime < ntime; itime ++){
        // plot input sources
        MifImgInfo* img_info_2d = new MifImgInfo;
        img_info_2d->InitSetImg(1, 1,
                                hi1d_xval->GetNbin(),
                                hi1d_yval->GetNbin());
        img_info_2d->PrintInfo();

        // fill 2d hist
        hd2d_arr[itime] = new HistDataNerr2d;
        hd2d_arr[itime]->Init(hi1d_xval->GetNbin(), hi1d_xval->GetLo(), hi1d_xval->GetUp(), 
                              hi1d_yval->GetNbin(), hi1d_yval->GetLo(), hi1d_yval->GetUp());
        for(long ipos = 0; ipos < npos_arr[itime]; ipos ++){
            hd2d_arr[itime]->Fill(xpos_arr[itime][ipos], ypos_arr[itime][ipos]);
        }

        for(long ibin = 0; ibin < hd2d_arr[itime]->GetNbin(); ibin ++){
            if(hd2d_arr[itime]->GetOvalArr()->GetValElm(ibin) > 1e-10){
                long ibinx = hd2d_arr[itime]->GetHi2d()->GetIbinX(ibin);
                long ibiny = hd2d_arr[itime]->GetHi2d()->GetIbinY(ibin);
                hd2d_arr[itime]->SetOvalElm(ibinx, ibiny, 1.0);
            }
        }
        char tag[kLineSize];
        sprintf(tag, "img_%2.2ld", itime);
        int bitpix = -32;
        MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                               2, bitpix,
                               img_info_2d->GetNaxesArr(),
                               hd2d_arr[itime]->GetOvalArr()->GetVal());

        for(long iarr = 0; iarr < hd2d_arr[itime]->GetNbin(); iarr ++){
            long icube = iarr + itime * hd2d_arr[itime]->GetNbin();
            cube_arr[icube] = hd2d_arr[itime]->GetOvalArr()->GetValElm(iarr);
        }
        delete img_info_2d;
    }

    {
        int bitpix = -32;
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), "cube",
                              3, bitpix,
                              img_info_3d->GetNaxesArr(),
                              cube_arr);
    }
    delete [] cube_arr;
    
    long nbin_par = nbin_rho * nbin_phi * nbin_psi * nbin_vel;
    double* par4d_arr = new double [nbin_par];
    long* par_index_arr = new long [nbin_par];
    for(long iarr = 0; iarr < nbin_par; iarr ++){
        par4d_arr[iarr] = 0.0;
        par_index_arr[iarr] = 0;
    }

    printf("par4d_arr init\n");

   

    for(long ivel = 0; ivel < nbin_vel; ivel ++){
        double theta =  atan(hi1d_vel->GetBinCenter(ivel) /
                             (time_arr[ntime - 1] * scale) );
        printf("theta = %e (deg)\n", theta / M_PI * 180.0);


//        // fake
//        double theta_fake = atan(hi1d_vel->GetBinCenter(0) /
//                                 (time_arr[ntime - 1] * scale) );
//        double rho_fake   = hi1d_rho->GetBinCenter(100);
//        double psi_fake   = hi1d_psi->GetBinCenter(10);
//        double phi_fake   = hi1d_phi->GetBinCenter(250);
//
//        printf("theta_fake, rho_fake, phi_fake, psi_fake = %e %e %e %e\n",
//               theta_fake, rho_fake, phi_fake, psi_fake);
//
//        
//        double* xpos_fake_arr = new double[ntime];
//        double* ypos_fake_arr = new double[ntime];
//        for(long itime = 0; itime < ntime; itime ++){
//            xpos_fake_arr[itime] = 0.0;
//            ypos_fake_arr[itime] = 0.0;
//        }
//        for(long itime = 0; itime < ntime; itime ++){
//            double zval  = zpos_arr[itime];
//            double xval = rho_fake * (cos(phi_fake) * cos(psi_fake)
//                                      - sin(phi_fake) * sin(psi_fake) / cos(theta_fake) )
//                + zval * tan(theta_fake) * sin(psi_fake);
//            double yval = rho_fake * (cos(phi_fake) * sin(psi_fake)
//                                      + sin(phi_fake) * cos(psi_fake) / cos(theta_fake) )
//                - zval * tan(theta_fake) * cos(psi_fake);
//            if(hi1d_xval->GetLo() < xval && xval < hi1d_xval->GetUp() &&
//               hi1d_yval->GetLo() < yval && yval < hi1d_yval->GetUp()){            
//                xpos_fake_arr[itime] = xval;
//                ypos_fake_arr[itime] = yval;
//            }
//            printf("xval_fake + 1024, yval_fake + 2088= %e, %e\n",
//                   xpos_fake_arr[itime] + 1024, ypos_fake_arr[itime] + 2088);
//            printf("circle(%e, %e, 20)\n",
//                   xpos_fake_arr[itime] + 1024, ypos_fake_arr[itime] + 2088);            
//        }
//        
//        // fill for fake
//        for(long ipsi = 0; ipsi < nbin_psi; ipsi ++){
//            printf("ipsi (nbin_psi) = %ld (%ld)\n", ipsi, nbin_psi);
//            
//            double psi = hi1d_psi->GetBinCenter(ipsi);
//            for(long itime = 0; itime < ntime; itime ++){
//                double xval = xpos_fake_arr[itime];
//                double yval = ypos_fake_arr[itime];
//                double zval = zpos_arr[itime];
//
//                // cluster
//                int nextx = 5;
//                int nexty = 5;
//                for(int iextx = 0; iextx < nextx; iextx ++){
//                    double extx = -2.0 + iextx * 1.0;
//                    for(int iexty = 0; iexty < nexty; iexty ++){
//                        double exty = -2.0 + iexty * 1.0;
//                        double xval_this = xval + extx;
//                        double yval_this = yval + exty;
//                        
//                        double rho_cos_phi = xval_this * cos(psi) + yval_this * sin(psi);
//                        double rho_sin_phi = -1 * xval_this * sin(psi) * cos(theta)
//                            + yval_this * cos(psi) * cos(theta)
//                            + zval * sin(theta);
//                        double rho = sqrt(pow(rho_cos_phi, 2) + pow(rho_sin_phi, 2));
//                        double phi = 0.0;
//                        if(rho_sin_phi >= 0.0){
//                            phi = acos(rho_cos_phi / rho);
//                        } else {
//                            phi = -1 * acos(rho_cos_phi / rho) + 2 * M_PI;
//                        }
//                        if(rho <= hi1d_rho->GetLo() || hi1d_rho->GetUp() <= rho){
//                            continue;
//                        }
//                        //if(phi <= hi1d_phi->GetLo() || hi1d_phi->GetUp() <= phi){
//                        //    continue;
//                        //}
//                        long irho = hi1d_rho->GetIbin(rho);
//                        long iphi = hi1d_phi->GetIbin(phi);
//                        long i_rho_phi_psi_vel = irho + iphi * nbin_rho
//                            + ipsi * (nbin_rho * nbin_phi)
//                            + ivel * (nbin_rho * nbin_phi * nbin_psi);
//                        par4d_arr[i_rho_phi_psi_vel] += 1.0;
//                    }
//                }
//                // cluster
//            }
//        }
//
//        // === end of fake
        
        for(long ipsi = 0; ipsi < nbin_psi; ipsi ++){
            printf("ipsi (nbin_psi) = %ld (%ld)\n", ipsi, nbin_psi);
            
            double psi = hi1d_psi->GetBinCenter(ipsi);
            for(long itime = 0; itime < ntime; itime ++){
                for(long ipos = 0; ipos < npos_arr[itime]; ipos ++){
                    double xval = xpos_arr[itime][ipos];
                    double yval = ypos_arr[itime][ipos];
                    // double zval = time_arr[itime];
                    double zval = zpos_arr[itime];

//                    // cluster
//                    int nextx = 5;
//                    int nexty = 5;
//                    for(int iextx = 0; iextx < nextx; iextx ++){
//                        double extx = -2.0 + iextx * 1.0;
//                        for(int iexty = 0; iexty < nexty; iexty ++){
//                            double exty = -2.0 + iexty * 1.0;
//                            double xval_this = xval + extx;
//                            double yval_this = yval + exty;
//                        
//                            double rho_cos_phi = xval_this * cos(psi) + yval_this * sin(psi);
//                            double rho_sin_phi = -1 * xval_this * sin(psi) * cos(theta)
//                                + yval_this * cos(psi) * cos(theta)
//                                + zval * sin(theta);
//                            double rho = sqrt(pow(rho_cos_phi, 2) + pow(rho_sin_phi, 2));
//                            double phi = 0.0;
//                            if(rho_sin_phi >= 0.0){
//                                phi = acos(rho_cos_phi / rho);
//                            } else {
//                                phi = -1 * acos(rho_cos_phi / rho) + 2 * M_PI;
//                            }
//                            if(rho <= hi1d_rho->GetLo() || hi1d_rho->GetUp() <= rho){
//                                continue;
//                            }
//                            //if(phi <= hi1d_phi->GetLo() || hi1d_phi->GetUp() <= phi){
//                            //    continue;
//                            //}
//                            long irho = hi1d_rho->GetIbin(rho);
//                            long iphi = hi1d_phi->GetIbin(phi);
//                            long i_rho_phi_psi_vel = irho + iphi * nbin_rho
//                                + ipsi * (nbin_rho * nbin_phi)
//                                + ivel * (nbin_rho * nbin_phi * nbin_psi);
//                            par4d_arr[i_rho_phi_psi_vel] += 1.0;
//                        }
//                    }
//                    // cluster

                    double xval_this = xval;
                    double yval_this = yval;
                    double rho_cos_phi = xval_this * cos(psi) + yval_this * sin(psi);
                    double rho_sin_phi = -1 * xval_this * sin(psi) * cos(theta)
                        + yval_this * cos(psi) * cos(theta)
                        + zval * sin(theta);
                    double rho = sqrt(pow(rho_cos_phi, 2) + pow(rho_sin_phi, 2));
                    double phi = 0.0;
                    if(rho_sin_phi >= 0.0){
                        phi = acos(rho_cos_phi / rho);
                    } else {
                        phi = -1 * acos(rho_cos_phi / rho) + 2 * M_PI;
                    }
                    if(rho <= hi1d_rho->GetLo() || hi1d_rho->GetUp() <= rho){
                        continue;
                    }
                    //if(phi <= hi1d_phi->GetLo() || hi1d_phi->GetUp() <= phi){
                    //    continue;
                    //}
                    long irho = hi1d_rho->GetIbin(rho);
                    long iphi = hi1d_phi->GetIbin(phi);
                    long i_rho_phi_psi_vel = irho + iphi * nbin_rho
                        + ipsi * (nbin_rho * nbin_phi)
                        + ivel * (nbin_rho * nbin_phi * nbin_psi);
                    par4d_arr[i_rho_phi_psi_vel] += 1.0;

                    
                }
            }
        }


// calc flat
        
        double* flat_arr = new double[nbin_rho * nbin_phi * nbin_psi];
        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr++){
            flat_arr[iarr] = 0.0;
        }

        for(long ipsi = 0; ipsi < nbin_psi; ipsi ++){

            printf("for flat: ipsi (nbin_psi) = %ld (%ld)\n", ipsi, nbin_psi);
            
            double psi = hi1d_psi->GetBinCenter(ipsi);
            for(long itime = 0; itime < ntime; itime ++){
                //                double zval = time_arr[itime];
                double zval = zpos_arr[itime];
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
                            flat_arr[i_rho_phi_psi] += rho / theta;
                        }
                    }
                }
            }
        }

//        // load flat data
//        double* flat_arr = NULL;
//        int bitpix = 0;
//        char flat_file[kLineSize];
//        sprintf(flat_file, "flat/flat_vel_%2.2ld.fits", ivel);
//        MifFits::InFitsCubeD(flat_file,
//                             img_info_par3d,
//                             &bitpix,
//                             &flat_arr);
        
        // div by flat
        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr++){
            long i_rho_phi_psi = iarr;
            long i_rho_phi_psi_vel = iarr
                + ivel * (nbin_rho * nbin_phi * nbin_psi);
            if(flat_arr[i_rho_phi_psi] > 1e-10){
                par4d_arr[i_rho_phi_psi_vel] /= flat_arr[i_rho_phi_psi];
            }
        }
        delete [] flat_arr;
    }
    printf("par4d_arr filled\n");

    // save par3d_arr
    for(long ivel = 0; ivel < nbin_vel; ivel ++){
        double* par3d_arr = new double [nbin_rho * nbin_phi * nbin_psi];
        for(long iarr = 0; iarr < nbin_rho * nbin_phi * nbin_psi; iarr ++){
            long i_rho_phi_psi_vel = iarr + ivel * (nbin_rho * nbin_phi * nbin_psi);
            par3d_arr[iarr] = par4d_arr[i_rho_phi_psi_vel];
        }
        char tag[kLineSize];
        sprintf(tag, "ivel_%2.2ld", ivel);
        int bitpix = -32;
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              3, bitpix,
                              img_info_par3d->GetNaxesArr(),
                              par3d_arr);
        delete [] par3d_arr;
    }

    MiSort::Sort<double, long>(nbin_par, par4d_arr, par_index_arr, 1);

    printf("sort done.\n");

    char outfile[kLineSize];
    sprintf(outfile, "%s/%s_par.dat",
            argval->GetOutdir().c_str(), argval->GetOutfileHead().c_str());
    FILE* fp_outfile = fopen(outfile, "w");

    char qdpfile[kLineSize];
    sprintf(qdpfile, "%s/%s_plot.qdp",
            argval->GetOutdir().c_str(), argval->GetOutfileHead().c_str());
    FILE* fp_qdp = fopen(qdpfile, "w");

    char regfile[kLineSize];
    sprintf(regfile, "%s/%s_plot.reg",
            argval->GetOutdir().c_str(), argval->GetOutfileHead().c_str());
    FILE* fp_reg = fopen(regfile, "w");
    fprintf(fp_reg, "# Region file format: DS9 version 4.1\n");
    fprintf(fp_reg, "global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 "
            "highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n");
    fprintf(fp_reg, "image\n");
    
    // reconstructed movie
    MifImgInfo* img_info_rec_3d = new MifImgInfo;
    img_info_rec_3d->InitSetCube(1, 1, 1,
                                 hi1d_xval->GetNbin(),
                                 hi1d_yval->GetNbin(),
                                 ntime);
    img_info_rec_3d->PrintInfo();

    
    double* rec_arr    = new double [img_info_rec_3d->GetNpixelTotal()];
    for(long iarr = 0; iarr < img_info_rec_3d->GetNpixelTotal(); iarr++){
        rec_arr[iarr] = 0.0;
    }
    double margin = 5.0;
    int ndetect_max = 100;
    int ndetect = 0;
    for(long iarr = 0; iarr < nbin_par; iarr ++){
        long ivel = par_index_arr[iarr] / (nbin_rho * nbin_phi * nbin_psi);
        long index2 = par_index_arr[iarr] % (nbin_rho * nbin_phi * nbin_psi);
        long ipsi = index2 / (nbin_rho * nbin_phi);
        long index3 = index2 % (nbin_rho * nbin_phi);
        long iphi = index3 / nbin_rho;
        long irho = index3 % nbin_rho;
        double rho = hi1d_rho->GetBinCenter(irho);
        double phi = hi1d_phi->GetBinCenter(iphi);
        double psi = hi1d_psi->GetBinCenter(ipsi);
        double vel = hi1d_vel->GetBinCenter(ivel);

        // pixel value along line, which must be constant
        int flag_bad = 0;
        double* line_arr = new double[ntime];
        long* index_line_arr = new long[ntime];
        double* xval_arr = new double[ntime];
        double* yval_arr = new double[ntime];
        for(long itime = 0; itime < ntime; itime ++){
            line_arr[itime] = 0.0;
            index_line_arr[itime] = 0;
            xval_arr[itime] = 0.0;
            yval_arr[itime] = 0.0;
        }
        double sum2_std_img = 0.0;
        int ninfov = 0;
        for(long itime = 0; itime < ntime; itime ++){
            // double zval  = time_arr[itime];
            double zval  = zpos_arr[itime];
            double theta =  atan(vel / (time_arr[ntime - 1] * scale) );
            double xval = rho * (cos(phi) * cos(psi) - sin(phi) * sin(psi) / cos(theta) )
                + zval * tan(theta) * sin(psi);
            double yval = rho * (cos(phi) * sin(psi) + sin(phi) * cos(psi) / cos(theta) )
                - zval * tan(theta) * cos(psi);

//            if(hi1d_rho->GetIbin(rho) < 10){
//                line_arr[itime] = 0.0;
//                flag_bad ++;
//                continue;
//            }
//            if(-2 < xval && xval < 2){
//                line_arr[itime] = 0.0;
//                flag_bad ++;
//                continue;
//            }
//            if(-2 < yval && yval < 2){
//                line_arr[itime] = 0.0;
//                flag_bad ++;
//                continue;
//            }

            xval_arr[itime] = xval;
            yval_arr[itime] = yval;
            
            if(xval < hi1d_xval->GetLo() + margin || hi1d_xval->GetUp() - margin < xval){
                line_arr[itime] = 0.0;
                // flag_bad ++;
                continue;
            }
            if(yval < hi1d_yval->GetLo() + margin || hi1d_yval->GetUp() - margin < yval){
                line_arr[itime] = 0.0;
                // flag_bad ++;
                continue;
            }
            ninfov ++;
            
            // get hist array data
            
            double oval = hd2d_arr[itime]->GetOvalElmAtXY(xval, yval);
            long iposx = hi1d_xval->GetIbin(xval);
            long iposy = hi1d_yval->GetIbin(yval);
            index_line_arr[itime] = iposx + iposy * hi1d_xval->GetNbin()
                + itime * hi1d_xval->GetNbin() * hi1d_yval->GetNbin();
            line_arr[itime] = oval;
            sum2_std_img += 1.0;
        }
        double sum = MirMath::GetSum(ntime, line_arr);
        //double stddev = MirMath::GetStddev(ntime, line_arr);
        //double mean   = MirMath::GetAMean(ntime, line_arr);
        //double sigma_std_img = sqrt(sum2_std_img);
        
        if(sum < 3){
            flag_bad ++;
        }
        if( fabs(ninfov - sum) > 1e-10){
            flag_bad ++;
        }
        // fprintf(fp_qdp, "\n");
        //fprintf(fp_qdp, "no\n");
        //fprintf(fp_qdp, "\n");
        // fprintf(fp_reg, "\n");
        
        if(0 == flag_bad){
            ndetect ++;
            printf("ndetect = %d, iarr = %ld, par4d = %e, sum = %e, ninfov = %d, "
                   "!  %e  %e  %e  %e ! %ld, %ld, %ld "
                   "! vel, rho, phi, psi ! irho, iphi, ipsi\n",
                   ndetect, iarr, par4d_arr[par_index_arr[iarr]], sum, ninfov,
                   vel, rho, phi, psi,
                   hi1d_rho->GetIbin(rho),
                   hi1d_phi->GetIbin(phi),
                   hi1d_psi->GetIbin(psi));
            fprintf(fp_outfile, "%e  %e  %e  %e  %e ! vel, rho, phi, psi, par4d \n",
                    vel, rho, phi, psi,
                    par4d_arr[par_index_arr[iarr]]);
            for(long itime = 0; itime < ntime; itime ++){
                long index = index_line_arr[itime];
                rec_arr[index] = line_arr[itime];

                fprintf(fp_reg, "circle(%e, %e, 10)  # text={%d - %ld} \n",
                        xval_arr[itime] + hi1d_xval->GetFullWidth()/2.0,
                        yval_arr[itime] + hi1d_yval->GetFullWidth()/2.0,
                        ndetect, itime + 1);
            }
        }
        
        delete [] line_arr;
        delete [] index_line_arr;
        delete [] xval_arr;
        delete [] yval_arr;

        if(ndetect >= ndetect_max){
            printf("iarr = %ld, iarr/nbin_par = %e, par4d = %e\n",
                   iarr, (double) iarr / (double) nbin_par, par4d_arr[par_index_arr[iarr]]);
            printf("break by 2\n");
            break;
        }
///        if(iarr > nbin_par * 0.5){
///            printf("iarr = %ld, iarr/nbin_par = %e, par4d = %e\n",
///                   iarr, (double) iarr / (double) nbin_par, par4d_arr[par_index_arr[iarr]]);
///            printf("break by 3\n");
///            break;
///        }
        if(par4d_arr[par_index_arr[iarr]] < 1.0e-10){
            printf("iarr = %ld, iarr/nbin_par = %e, par4d = %e\n",
                   iarr, (double) iarr / (double) nbin_par, par4d_arr[par_index_arr[iarr]]);
            printf("break by 4\n");
            break;
        }
    }
    fclose(fp_outfile);
    fclose(fp_qdp);
    fclose(fp_reg);
    
    int bitpix = -32;
    MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), "rec",
                          3, bitpix,
                          img_info_rec_3d->GetNaxesArr(),
                          rec_arr);

    double time_ed = MiTime::GetTimeSec();
    printf("time_ed - time_st = %e\n", time_ed - time_st);
    
    return status_prog;
}
