// Copyright 2018 @ Chun Shen

#include "eos_PCE.h"
#include "util.h"

#include <sstream>
#include <fstream>
#include <cmath>

using std::stringstream;
using std::string;

EOS_PCE::EOS_PCE(const InitData &DATA_in) : DATA(DATA_in) {
    set_EOS_id(18);
    set_number_of_tables(0);
    set_eps_max(1e5);
    set_flag_muB(false);
    set_flag_muS(false);
    set_flag_muC(false);
}


void EOS_PCE::initialize_eos() {
    // read the lattice EOS pressure, temperature, and 
    music_message.info("reading PCE EOS ...");
    
    auto envPath = get_hydro_env_path();
    stringstream slocalpath;
    slocalpath << envPath << "/EOS/PCE";

    string path = slocalpath.str();
    music_message << "from path " << path;
    music_message.flush("info");
   
    set_number_of_tables(2);
    resize_table_info_arrays();

    int ntables = get_number_of_tables();

    pressure_tb    = new double** [ntables];
    temperature_tb = new double** [ntables]; 
    cs2_tb         = new double** [ntables];
    energy_tb      = new double** [ntables];
    entropy_tb     = new double** [ntables];

    std::ifstream qcd_eos(path + "/QGP_eos_Vovchenko_e0p25spacing.dat", std::ios::binary);
    std::ifstream gluon_eos(path + "/gluon_eos_Vovchenko.dat", std::ios::binary);
    
    for (int itable = 0; itable < ntables; itable++) { // 0 = QCD, 1 = gluon
        std::ifstream* eos_file;
	if (itable == 0) {
	  eos_file = &qcd_eos;
	}
	else {
	  eos_file = &gluon_eos;
	}

        if (!*eos_file) {
            music_message.error("Can not find the EoS file.");
            exit(1);
        }
	
        e_length[itable]  = 100000;
        nb_length[itable] = 1;
        // allocate memory for pressure arrays
        pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
	cs2_tb[itable] = Util::mtx_malloc(nb_length[itable], e_length[itable]);
        energy_tb[itable] = Util::mtx_malloc(nb_length[itable], e_length[itable]);
	entropy_tb[itable] = Util::mtx_malloc(nb_length[itable], e_length[itable]);

        double temp;
        for (int ii = 0; ii < e_length[itable]; ii++) {
	    eos_file->read((char*)&temp, sizeof(double));  // e^(1/4) (1/fm)
            //temp /= Util::hbarc;      // 1/fm^4
            if (ii == 0) e_bounds[itable] = temp;
            if (ii == 1) e_spacing[itable] = temp - e_bounds[itable];
            if (ii == e_length[itable] - 1) set_eps_max(std::max(get_eps_max(),std::pow(temp,4)));
	    energy_tb[itable][0][ii] = temp;

            eos_file->read((char*)&temp, sizeof(double));  // P
            pressure_tb[itable][0][ii] = temp/Util::hbarc;      // 1/fm^4

            eos_file->read((char*)&temp, sizeof(double));  // s
	    entropy_tb[itable][0][ii] = temp;                   // 1/fm^3

            eos_file->read((char*)&temp, sizeof(double));  // T
            temperature_tb[itable][0][ii] = temp/Util::hbarc;   // 1/fm

	    eos_file->read((char*)&temp, sizeof(double));  // cs2
	    cs2_tb[itable][0][ii] = temp;
        }

	// if (itable == 1) {
        //     acc = gsl_interp_accel_alloc();
	    
        //     pressure_spline = gsl_spline_alloc(gsl_interp_linear, e_length[itable]);
	//     temperature_spline = gsl_spline_alloc(gsl_interp_linear, e_length[itable]);
	//     cs2_spline = gsl_spline_alloc(gsl_interp_linear, e_length[itable]);
	//     entropy_spline = gsl_spline_alloc(gsl_interp_linear, e_length[itable]);

	//     gsl_spline_init(pressure_spline, energy_tb[itable][0], pressure_tb[itable][0], e_length[itable]);
	//     gsl_spline_init(temperature_spline, energy_tb[itable][0], temperature_tb[itable][0], e_length[itable]);
	//     gsl_spline_init(cs2_spline, energy_tb[itable][0], cs2_tb[itable][0], e_length[itable]);
	//     gsl_spline_init(entropy_spline, energy_tb[itable][0], entropy_tb[itable][0], e_length[itable]);
	// }
    }
    music_message.info("Done reading EOS.");
}


// double EOS_PCE::interpolate1D(double e, int table_idx, double ***table, gsl_spline *spline) const {
// // This is a generic linear interpolation routine for EOS at zero mu_B
// // it assumes the class has already read in
// //        P(e), T(e), s(e)
// // as one-dimensional arrays on an equally spacing lattice grid
// // units: e is in 1/fm^4
//     //double local_ed = e*hbarc;  // [GeV/fm^3]
//     double local_ed = e;

//     const double e0       = e_bounds[table_idx];
//     const double delta_e  = e_spacing[table_idx];
//     const int N_e         = e_length[table_idx];

//     // compute the indices
//     int idx_e  = static_cast<int>((local_ed - e0)/delta_e);

//     // treatment for overflow, use the last two points to do extrapolation
//     idx_e  = std::min(N_e - 2, idx_e);

//     // check underflow
//     idx_e  = std::max(0, idx_e);

//     const double frac_e = (local_ed - (idx_e*delta_e + e0))/delta_e;

//     double result;
    
//     double temp1 = table[table_idx][0][idx_e];
//     double temp2 = table[table_idx][0][idx_e + 1];
//     result = temp1*(1. - frac_e) + temp2*frac_e;

//     if (idx_e > 1 and table_idx == 1) {result = gsl_spline_eval(spline, e, acc);}
//     //    std::cout << "Spacing: " << table_idx << "  " << e0 << "  " << delta_e << "  " << N_e << std::endl;
//     //   std::cout << "Interpolate1D: " << table_idx << "  " << e << "  " << result << "  " << get_eps_max() << std::endl; 
//     return(result);
// }


double EOS_PCE::p_e_func(double e, double rhob, double proper_tau) const {
    return(get_dpOverde3(e, rhob, proper_tau));
}

double EOS_PCE::get_dpOverde3(double e, double rhob, double proper_tau) const {
    return(get_cs2(e, rhob, proper_tau));
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_PCE::get_temperature(double e, double rhob, double proper_tau) const {
    double T_QCD = interpolate1D(std::pow(e,0.25), 0, temperature_tb);  // 1/fm
    double T_gluon = interpolate1D(std::pow(e,0.25), 1, temperature_tb);
    double fugacity = get_fugacity(proper_tau);

    //T_QCD = std::pow(e * 30./std::pow(M_PI,2) * 1/(2*(3*3.-1) + 3.5*3*3), 0.25);
    double T = fugacity * T_QCD + (1 - fugacity) * T_gluon;
    return(std::max(1e-15, T));
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS_PCE::get_pressure(double e, double rhob, double proper_tau) const {
    double f_QCD = interpolate1D(std::pow(e,0.25), 0, pressure_tb);  // 1/fm^4
    double f_gluon = interpolate1D(std::pow(e,0.25), 1, pressure_tb);

    //    f_QCD = e/3.;
    double fugacity = get_fugacity(proper_tau);
    double f = fugacity * f_QCD + (1 - fugacity) * f_gluon;

    return(std::max(1e-15, f));
}


double EOS_PCE::get_s2e(double s, double rhob, double proper_tau) const {
    double e = get_s2e_finite_rhob(s, 0.0, proper_tau);
    return(e);
}

double EOS_PCE::get_T2e(double T, double rhob, double proper_tau) const {
    double e = get_T2e_finite_rhob(T, 0.0, proper_tau);
    return(e);
}

double EOS_PCE::get_fugacity(double proper_tau) const {
    double fugacity;
    double tau0 = DATA.tau0;
    double tau_eq = DATA.tau_eq;
    if (tau0 <= 0 or tau_eq <= 0) {
        fugacity = 1;
    }
    else {
        fugacity = 1 - exp(2*(tau0 - proper_tau)/tau_eq);

	if (fugacity < 0) return(0);
    }
    return(fugacity);
}

double EOS_PCE::get_cs2(double e, double rhob, double proper_tau) const {
    double cs2_gluon = interpolate1D(std::pow(e,0.25), 1, cs2_tb);
    double cs2_QCD = interpolate1D(std::pow(e,0.25), 0, cs2_tb);
    double fugacity = get_fugacity(proper_tau);
    //cs2_QCD = 1/3.;
    double cs2 = fugacity * cs2_QCD + (1 - fugacity) * cs2_gluon;
    return(std::max(1e-15, cs2));
}

// double EOS_PCE::get_entropy(double e, double rhob, double proper_tau) const {
//     double T = get_temperature(e, rhob, proper_tau);
//     double P_QCD = get_pressure(e, rhob, 0);
//     double P_gluon = get_pressure(e, rhob, 1e10);
//     double entropy_QCD = interpolate1D(std::pow(e,0.25), 0, entropy_tb);
//     double entropy_gluon = interpolate1D(std::pow(e,0.25), 1, entropy_tb);
//     double fugacity = get_fugacity(proper_tau);
//     //entropy_QCD = (e + get_pressure(e, rhob, proper_tau))/(get_temperature(e, rhob, proper_tau));
//     double entropy = fugacity * entropy_QCD + (1 - fugacity) * entropy_gluon - fugacity*log(fugacity)/T * (P_QCD - P_gluon);
//     return(std::max(1e-15,entropy));
// }
  
