#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "load_dat.hpp"
#include "Euler.hpp"
using namespace std;

float random_float(float mid_val, float variance) {
    float val = mid_val + (0.5 - ((float)rand() / (float)(RAND_MAX))) * (2. * variance);
    return val;
}

int main(int argc, char *argv[]) {

    //Variable declarations
    FILE *in, *error;
    char input_file[128], buffer[256], keyword[64], spec_file[128];
    int l, sites = 0, nspec = 0, steps = 200, max_steps = 100000, tries = 8, simplex=1000;
    float rate = 0.5, S=1.5;
    vector<int> TD;
    vector<int> MAS;
    vector<int> sat;
    vector<float> vs;
    vector<vector<float> > bin(nspec, vector<float>(1, 0.));
    vector<vector<float> > ppm(nspec, vector<float>(1, 0.));

    //Terminate the program if to input file is provided.
    if (argc < 2) {
        error = fopen("error.txt", "w");
        fprintf(error, "\nERROR: Missing input file declaration\n");
        fclose(error);
        exit(1);
    }

    sprintf(input_file, "%s", argv[1]);
    in = fopen(input_file, "r");

    //Terminate the program if the input file is missing.
    if (in == NULL) {
        error = fopen("error.txt", "w");
        fprintf(error, "\nERROR: Input file '%s' not found\n", input_file);
        fclose(error);
        exit(1);
    }

    //Read the input file
    while (fgets(buffer, sizeof(buffer), in) != NULL) {
        sscanf(buffer, "%s", keyword);

        //Load a spectrum to memory
        if (strcmp(keyword, "spectrum") == 0) {
            TD.push_back(1);
            MAS.push_back(0);
            sat.push_back(0);
            vs.push_back(0.);

            sscanf(buffer, "%s %s %d %s", keyword, keyword, &TD[nspec], spec_file);
            bin.push_back(vector<float>(TD[nspec], 0.));
            ppm.push_back(vector<float>(TD[nspec], 0.));

            vs[nspec] = load_dat(spec_file, TD[nspec], bin[nspec], ppm[nspec]);

            if (strcmp(keyword, "MAS") == 0) {
                MAS[nspec] = 1;
            }
            if (strcmp(keyword, "satellites") == 0) {
                sat[nspec] = 1;
            }

            nspec++;
            sprintf(keyword, "void");
        } else if (strcmp(keyword, "sites") == 0) {
            sscanf(buffer, "%s %d", keyword, &sites);
            sprintf(keyword, "void");
        } else if (strcmp(keyword, "spin") == 0) {
            sscanf(buffer, "%s %f", keyword, &S);
            sprintf(keyword, "void");
        } else if (strcmp(keyword, "steps") == 0) {
            sscanf(buffer, "%s %d", keyword, &steps);
            sprintf(keyword, "void");
        } else if (strcmp(keyword, "max_steps") == 0) {
            sscanf(buffer, "%s %d", keyword, &max_steps);
            sprintf(keyword, "void");
        } else if (strcmp(keyword, "tries") == 0) {
            sscanf(buffer, "%s %d", keyword, &tries);
            sprintf(keyword, "void");
        } else if (strcmp(keyword, "simplex") == 0) {
            sscanf(buffer, "%s %d", keyword, &simplex);
            sprintf(keyword, "void");
        } else if (strcmp(keyword, "rate") == 0) {
            sscanf(buffer, "%s %f", keyword, &rate);
            sprintf(keyword, "void");
        }
    }

    fclose(in);

    // Descriptions of the experimental data
    float width[nspec];
    for (l = 0; l < nspec; l++) {
        width[l] = ppm[l][0] - ppm[l][TD[l] - 1];
    }

    //create the best-fits comma-delimited output file.
    FILE *fp2;
    fp2 = fopen("best-fits.csv", "w");
    fprintf(fp2, "site,RMSD,intensity,diso/ppm,CQ/MHz,eta,span/ppm,skew,alpha,beta,gamma,dani,etaCS,phi,chi,psi\n");
    fclose(fp2);

    //The number of tries is set to a minimum of the number of threads
    #pragma omp parallel
    { l = omp_get_num_threads(); }
    if (tries < l) tries = l;

    //parallel loop over the various tries.
    //A while loop is used that terminates when the right number of tries have been done.
    //This ensures maximal hardware utilization.


    #pragma omp parallel
    {
        while (true) {
            //thread-specific parameter declarations
            float diso_range[sites][3], span_range[sites][3], skew_range[sites][3], chi_range[sites][3], eta_range[sites][3], GB_range[nspec][sites][3], LB_range[nspec][sites][3], alpha_range[sites][3], beta_range[sites][3], gamma_range[sites][3], intensity_range[nspec][sites][3];
            int i, j, k;
			FILE *in2, *fp3;

            // initializing parameters
            for (i = 0; i < sites; i++) {
                diso_range[i][0] = diso_range[i][1] = 0.;
                span_range[i][0] = span_range[i][1] = 0.;
                skew_range[i][0] = skew_range[i][1] = 0.;
                chi_range[i][0] = chi_range[i][1] = 0.;
                eta_range[i][0] = eta_range[i][1] = 0.;
                alpha_range[i][0] = alpha_range[i][1] = 0.;
                beta_range[i][0] = beta_range[i][1] = 0.;
                gamma_range[i][0] = gamma_range[i][1] = 0.;
                for (j = 0; j < nspec; j++) {
                    GB_range[j][i][0] = GB_range[j][i][1] = 0.;
                    LB_range[j][i][0] = LB_range[j][i][1] = 0.;
                    intensity_range[j][i][0] = intensity_range[j][i][1] = 0.;
                }
            }

            //Each thread reads the input file again to extract the mean parameter values and their ranges
            if (sites > 0) {
                char buffer2[256], keyword2[64];
                in2 = fopen(input_file, "r");
                while (fgets(buffer2, sizeof(buffer2), in2) != NULL) {
                    sscanf(buffer2, "%s", keyword2);
                    if (strcmp(keyword2, "diso") == 0) {
                        sscanf(buffer2, "%s %d", keyword2, &i);
                        sscanf(buffer2, "%s %d %f %f", keyword2, &i, &diso_range[i - 1][0], &diso_range[i - 1][1]);
                        sprintf(keyword2, "void");
                    } else if (strcmp(keyword2, "CQ") == 0) {
                        sscanf(buffer2, "%s %d", keyword2, &i);
                        sscanf(buffer2, "%s %d %f %f", keyword2, &i, &chi_range[i - 1][0], &chi_range[i - 1][1]);
                        chi_range[i - 1][0] *= 1.e6;
                        chi_range[i - 1][1] *= 1.e6;
                        sprintf(keyword2, "void");
                    }

                    else if (strcmp(keyword2, "eta") == 0) {
                        sscanf(buffer2, "%s %d", keyword2, &i);
                        sscanf(buffer2, "%s %d %f %f", keyword2, &i, &eta_range[i - 1][0], &eta_range[i - 1][1]);
                        sprintf(keyword2, "void");
                    }

                    else if (strcmp(keyword2, "span") == 0) {
                        sscanf(buffer2, "%s %d", keyword2, &i);
                        sscanf(buffer2, "%s %d %f %f", keyword2, &i, &span_range[i - 1][0], &span_range[i - 1][1]);
                        sprintf(keyword2, "void");
                    }

                    else if (strcmp(keyword2, "skew") == 0) {
                        sscanf(buffer2, "%s %d", keyword2, &i);
                        sscanf(buffer2, "%s %d %f %f", keyword2, &i, &skew_range[i - 1][0], &skew_range[i - 1][1]);
                        sprintf(keyword2, "void");
                    }

                    else if (strcmp(keyword2, "alpha") == 0) {
                        sscanf(buffer2, "%s %d", keyword2, &i);
                        sscanf(buffer2, "%s %d %f %f", keyword2, &i, &alpha_range[i - 1][0], &alpha_range[i - 1][1]);
                        alpha_range[i - 1][0] *= 0.017453293;
                        alpha_range[i - 1][1] *= 0.017453293;
                        sprintf(keyword2, "void");
                    }

                    else if (strcmp(keyword2, "beta") == 0) {
                        sscanf(buffer2, "%s %d", keyword2, &i);
                        sscanf(buffer2, "%s %d %f %f", keyword2, &i, &beta_range[i - 1][0], &beta_range[i - 1][1]);
                        beta_range[i - 1][0] *= 0.017453293;
                        beta_range[i - 1][1] *= 0.017453293;
                        sprintf(keyword2, "void");
                    }

                    else if (strcmp(keyword2, "gamma") == 0) {
                        sscanf(buffer2, "%s %d", keyword2, &i);
                        sscanf(buffer2, "%s %d %f %f", keyword2, &i, &gamma_range[i - 1][0], &gamma_range[i - 1][1]);
                        gamma_range[i - 1][0] *= 0.017453293;
                        gamma_range[i - 1][1] *= 0.017453293;
                        sprintf(keyword2, "void");
                    }

                    //Line broadening is specified in ppm, so here they are reduced to correspond to the same Hz values for all spectra
                    else if (strcmp(keyword2, "GB") == 0) {
                        sscanf(buffer2, "%s %d", keyword2, &i);
                        sscanf(buffer2, "%s %d %f %f", keyword2, &i, &GB_range[0][i - 1][0], &GB_range[0][i - 1][1]);
                        for (j = 0; j < nspec; j++) {
                            GB_range[j][i - 1][0] = 1.e6 * GB_range[0][i - 1][0] / vs[j];
                            GB_range[j][i - 1][1] = 1.e6 * GB_range[0][i - 1][1] / vs[j];
                        }
                        sprintf(keyword2, "void");
                    } else if (strcmp(keyword2, "LB") == 0) {
                        sscanf(buffer2, "%s %d", keyword2, &i);
                        sscanf(buffer2, "%s %d %f %f", keyword2, &i, &LB_range[0][i - 1][0], &LB_range[0][i - 1][1]);

                        for (j = 0; j < nspec; j++) {
                            LB_range[j][i - 1][0] = 1.e6 * LB_range[0][i - 1][0] / vs[j];
                            LB_range[j][i - 1][1] = 1.e6 * LB_range[0][i - 1][1] / vs[j];
                        }
                        sprintf(keyword2, "void");
                    }
                }
                fclose(in2);
            }
            for (i = 0; i < sites; i++) {
                for (j = 0; j < nspec; j++) {
                    intensity_range[j][i][0] = 1.0 / sites;
                    intensity_range[j][i][1] = .5 / sites;
                }
            }

            srand(time(NULL) * (1 + omp_get_thread_num()));

            float diso[sites], span[sites], skew[sites], chi[sites], eta[sites];//, GB[nspec][sites], LB[nspec][sites], intensity[nspec][sites];
            float alpha[sites], beta[sites], gamma[sites];

            vector<vector<float> > GB(nspec,vector<float>(sites, 0.));
            vector<vector<float> > LB(nspec,vector<float>(sites, 0.));
            vector<vector<float> > intensity(nspec,vector<float>(sites, 0.));

            gsl_vector *ss, *var;
            var = gsl_vector_alloc (8*sites+sites*3*nspec);
            ss = gsl_vector_alloc (8*sites+sites*3*nspec);
            struct par parameters;
            create_params(&parameters, nspec, TD, MAS, sat, 1, sites, bin, ppm, S, vs, width);

            // RMSD calculation
            float RMSD, RMSD_min;
            RMSD_min = 1.e10;
            float temp = 1.;
            float RMSD_log[3];
            RMSD_log[0] = RMSD_min;
            RMSD_log[2] = RMSD_min;
            i = 0;
            int fast=1;

            do {
                i++;
                for (j = 0; j < sites; j++) {
                    for (k = 0; k < nspec; k++) {
                        GB[k][j] = fabs(random_float(GB_range[k][j][0], temp * GB_range[k][j][1]));
                        LB[k][j] = fabs(random_float(LB_range[k][j][0], temp * LB_range[k][j][1]));
                        intensity[k][j] = fabs(random_float(intensity_range[k][j][0], temp * intensity_range[k][j][1]));
                    }

                    if(fast){
                        //In the first part of the fitting, the same Hz broadening is applied to all spectra, with is being 1/3 lower for the MAS data
                        for (k = 0; k < nspec; k++) {
                            GB[k][j] = GB[0][j]*vs[0]/vs[k];
                            LB[k][j] = LB[0][j]*vs[0]/vs[k];
                            if(!MAS[k]){
                                GB[k][j]*=3.;
                                LB[k][j]*=3.;
                            }
                        }
                    }

                    diso[j] = random_float(diso_range[j][0], temp * diso_range[j][1]);
                    chi[j] = random_float(chi_range[j][0], temp * chi_range[j][1]);
                    eta[j] = fabs(random_float(eta_range[j][0], temp * eta_range[j][1]));
                    span[j] = fabs(random_float(span_range[j][0], temp * span_range[j][1]));
                    skew[j] = random_float(skew_range[j][0], temp * skew_range[j][1]);
                    alpha[j] = random_float(alpha_range[j][0], temp * alpha_range[j][1]);
                    beta[j] = fabs(random_float(beta_range[j][0], temp * beta_range[j][1]));
                    gamma[j] = random_float(gamma_range[j][0], temp * gamma_range[j][1]);

                    //Correcting the parameters to ensure proper values
                    if (eta[j] > 1.0) eta[j] = 1.0;
                    if (skew[j] > 1.0)
                        skew[j] = 1.0;
                    else if (skew[j] < -1.0)
                        skew[j] = -1.0;
                    if (beta[j] > Pi) {
                        beta[j] -= Pi;
                        alpha[j] = -alpha[j];
                    }
                    if (alpha[j] < 0) alpha[j] += Pi;
                    if (gamma[j] < 0) gamma[j] += Pi;
                    if (alpha[j] > Pi) alpha[j] -= Pi;
                    if (gamma[j] > Pi) gamma[j] -= Pi;
                }

                arrange_variables(var, nspec, sites, diso, span, skew, chi, eta, alpha, beta, gamma, LB, GB, intensity);

                //Calculation of the total RMSD
                RMSD=total_RMSD(var, &parameters);

                //Commands that are done when a better fit is found, namely all mean parameters are updated
                if (RMSD < RMSD_min) {
                    RMSD_min = RMSD;
                    if(fast)
                        printf("%d %f (ASSR fast interpolation) (%d/%d) \n",i, RMSD, l, tries);
                    else
                        printf("%d %f (ASSR slow interpolation) (%d/%d) \n",i, RMSD, l, tries);
                    for (j = 0; j < sites; j++) {
                        for (k = 0; k < nspec; k++) {
                            GB_range[k][j][0] = GB[k][j];
                            LB_range[k][j][0] = LB[k][j];
                            intensity_range[k][j][0] = intensity[k][j];
                        }
                        diso_range[j][0] = diso[j];
                        chi_range[j][0] = chi[j];
                        eta_range[j][0] = eta[j];
                        span_range[j][0] = span[j];
                        skew_range[j][0] = skew[j];
                        alpha_range[j][0] = alpha[j];
                        beta_range[j][0] = beta[j];
                        gamma_range[j][0] = gamma[j];
                    }
                }

                //Checking so see if improvement has continued and whether the step size should be reduced.
                if (i % steps == 0) {
                    RMSD_log[1] = RMSD_min;
                    if (RMSD_log[1] == RMSD_log[0]) {
                        temp = temp * rate;
                        if(temp<0.05){
                            fast=0;
                            parameters.fast=0;
                            RMSD=total_RMSD(var, &parameters);
                            RMSD_log[1] = RMSD_min = RMSD;
                        }
                        if (RMSD_log[2] == RMSD_log[0])
                            break;
                        else
                            RMSD_log[2] = RMSD_log[0];
                    }
                    RMSD_log[0] = RMSD_log[1];
                }

            } while (i < max_steps);



            //selecting the minimizer, here is a O(N) simplex optimizer
            parameters.fast=0;
            const gsl_multimin_fminimizer_type *T =  gsl_multimin_fminimizer_nmsimplex2;
            gsl_multimin_fminimizer *s = NULL;
            gsl_multimin_function minex_func;
            minex_func.n = 8*sites+sites*3*nspec;
            minex_func.f = total_RMSD; //function assignment
            minex_func.params = &parameters; //parameter assignment
            set_step_size(var,ss,sites,nspec);
            s = gsl_multimin_fminimizer_alloc (T, minex_func.n); //giving the minimizer (s) the algorithm (T) and dimensionality (2)
            gsl_multimin_fminimizer_set (s, &minex_func, var, ss); //giving the minimizer (s) the function details (minex_func), variables (x), and step sizes (ss).

            size_t iter = 0;
            int status;
            double size, init_size=gsl_multimin_fminimizer_size (s);

            do{
                iter++;
                status = gsl_multimin_fminimizer_iterate(s);

                if (status)
                    break;

                size = gsl_multimin_fminimizer_size (s);
                status = gsl_multimin_test_size (size/init_size, 1e-4);

                if (status == GSL_SUCCESS){
                    break;
                }
                RMSD=s->fval;

                if(RMSD<RMSD_min){
                    printf ("%d %f (simplex of size %lf) (%d/%d)\n",i+iter,s->fval,size,l,tries);
                    RMSD_min=RMSD;
                }

            }  while (status == GSL_CONTINUE && iter < simplex);

            disentangle_variables(s->x,nspec,sites,diso,span,skew,chi,eta,alpha,beta,gamma,LB,GB,intensity);
            RMSD=RMSD_min=total_RMSD(s->x,&parameters);
            printf ("simplex optimized to %f (%d/%d)\n",RMSD,l,tries);


            //One fit has completed, the resulting parameters are written to the best-fits file.
            for (j = 0; j < sites; j++) {
                //Euler angles range is corrected for consistency
                Fix_Angles(alpha[j],beta[j],gamma[j]);
                //Conversions for those that want to use DMFIT
                float phi,chii,psi,dani,etaCS;
                float d22 = skew[j] * span[j] / 3. + diso[j];
                float d33 = (3. * diso[j] - d22 - span[j]) / 2.;
                float d11 = 3. * diso[j] - d22 - d33;
                if(skew[j]>0){
                    dani=d33-diso[j];
                    etaCS=(d22-d11)/dani;
                }
                else{
                    dani=d11-diso[j];
                    etaCS=(d22-d33)/dani;
                }
                Calculate_DMFIT_angles(alpha[j],beta[j],gamma[j],skew[j],phi,chii,psi);
				fp3 = fopen("best-fits.csv", "a");
                fprintf(fp3, "%d,%f,%.1f,%.2f,%.2f,%.2f,%.2f,%.2f,%.0f,%.0f,%.0f,%.2f,%.2f,%.0f,%.0f,%.0f\n", j+1, RMSD_min, intensity[0][j] * 100., diso[j], chi[j] / 1e6, eta[j], span[j], skew[j], alpha[j] * 180. / Pi, beta[j] * 180. / Pi, gamma[j] * 180. / Pi,dani,etaCS,phi * 180. / Pi,chii * 180. / Pi,psi * 180. / Pi);
                fclose(fp3);
                //simulations are written in a separate CSV file.
                write_fits(RMSD_min,nspec,TD,MAS,sat,sites,bin,ppm,diso,span,skew,chi,eta,S,vs,alpha,beta,gamma,LB,GB,intensity,width);
            }
            l++;
            if (l > tries) {
                l--;
                break;
            }
        }  // tries
    }      // OMP
    return 0;
}
