// Copyright 2018 @ Chun Shen

#include "eos_PCE.h"
#include "util.h"

#include <sstream>
#include <fstream>
#include <cmath>

using std::stringstream;
using std::string;

EOS_PCE::EOS_PCE() {
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
    energy_tb      = new double** [ntables]; 
    cs2_tb         = new double** [ntables];

    std::ifstream qcd_eos(path + "/hrg_hotqcd_eos_binary.dat", std::ios::binary);
    std::ifstream gluon_eos(path + "/gluon_eos_binary.dat", std::ios::binary);
    
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
        energy_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
	cs2_tb[itable] = Util::mtx_malloc(nb_length[itable], e_length[itable]);
        double temp;
        for (int ii = 0; ii < e_length[itable]; ii++) {
	    eos_file->read((char*)&temp, sizeof(double));  // e
            energy_tb[itable][0][ii] = temp/Util::hbarc;      // 1/fm^4
	    //    if (ii == 0) e_bounds[itable] = temp;
            //    if (ii == 1) e_spacing[itable] = temp - e_bounds[itable];
            if (ii == e_length[itable] - 1) set_eps_max(temp/Util::hbarc);
	    

            eos_file->read((char*)&temp, sizeof(double));  // P
            pressure_tb[itable][0][ii] = temp/Util::hbarc;      // 1/fm^4

            eos_file->read((char*)&temp, sizeof(double));  // s

            eos_file->read((char*)&temp, sizeof(double));  // T
            if (ii == 0) e_bounds[itable] = temp/Util::hbarc;   // 1/fm
	    if (ii == 1) e_spacing[itable] = temp/Util::hbarc - e_bounds[itable];
	    if (ii == e_length[itable] -1) set_T_max(temp/Util::hbarc);

	    eos_file->read((char*)&temp, sizeof(double));  // cs2
	    cs2_tb[itable][0][ii] = temp;
        }
    }
    std::cout << e_bounds[0] << " " << e_spacing[0] << " " << e_bounds[1] << " " << e_spacing[1] << std::endl;
    music_message.info("Done reading EOS.");
}


double EOS_PCE::p_e_func(double e, double rhob, double proper_tau) const {
    return(get_dpOverde3(e, rhob, proper_tau));
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
// Uses binary search as in get_T2e_finite_rhob - Andrew
double EOS_PCE::get_temperature(double e, double rhob, double proper_tau) const {
    double e_goal = e/Util::hbarc;         // convert to 1/fm^4
    double T_lower = 1e-15;  
    double T_upper = T_max;
    double T_mid  = (T_upper + T_lower)/2.;
    double e_lower   = get_T2e(T_lower, rhob, proper_tau);
    double e_upper   = 10*get_T2e(T_upper, rhob, proper_tau);
    int ntol         = 1000;
    
    if (e_goal < 0.0 || e_goal > e_upper) {
       std::cout << "get_temperature:: e is out of bounds, "
                 << "e = " << e << ", e_upper = " << e_upper*Util::hbarc
                 << ", e_lower = " << e_lower*Util::hbarc << " " << get_eps_max()*Util::hbarc << " " << T_max*Util::hbarc << std::endl;
       exit(1);
    }
    
    if (e_goal < e_lower) return(T_lower);

    double rel_accuracy = 1e-8;
    double abs_accuracy = 1e-15;
    double e_mid;
    int iter = 0;
    while (((T_upper - T_lower)/T_mid > rel_accuracy
            && (T_upper - T_lower) > abs_accuracy) && iter < ntol) {
        e_mid = get_T2e(T_mid, rhob, proper_tau);
        if (e_goal < e_mid)
            T_upper = T_mid;
        else 
            T_lower = T_mid;
        T_mid = (T_upper + T_lower)/2.;
        iter++;
    }
    if (iter == ntol) {
        std::cout << "get_temperature:: max iteration reached, "
	     	  << "e = " << e << ", rhob = " << rhob << std::endl;;
	std::cout << "e_upper = " << e_upper*Util::hbarc
		  << " , e_lower = " << e_lower*Util::hbarc << std::endl;
	std::cout << "T_upper = " << T_upper
                  << " , T_lower = " << T_lower
		  << ", diff = " << (T_upper - T_lower) << std::endl;
        exit(1);
    }
    return (T_mid);
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS_PCE::get_pressure(double e, double rhob, double proper_tau) const {
    double T = get_temperature(e, rhob, proper_tau);
    double f_QCD = interpolate1D(T, 0, pressure_tb);  // 1/fm^4
    double f_gluon = interpolate1D(T, 1, pressure_tb);
    double fugacity = get_fugacity(proper_tau);
    double f = fugacity * f_QCD + (1 - fugacity) * f_gluon;
    return(std::max(1e-15, f));
}


double EOS_PCE::get_s2e(double s, double rhob, double proper_tau) const {
    double e = get_s2e_finite_rhob(s, 0.0, proper_tau);
    return(e);
}

double EOS_PCE::get_T2e(double T, double rhob, double proper_tau) const {
    double e_QCD = interpolate1D(T, 0, energy_tb); // 1/fm^4
    double e_gluon = interpolate1D(T, 1, energy_tb);
    double fugacity = get_fugacity(proper_tau);
    double e = fugacity * e_QCD + (1 - fugacity) * e_gluon;
    return(e);
}

double EOS_PCE::get_fugacity(double proper_tau) const {
    double fugacity;
    double tau0 = 0;
    double tau_eq = 5.0;
    if (tau0 <= 0) {
        fugacity = 1;
    }
    else {
        fugacity = 1 - exp((tau0 - proper_tau)/tau_eq);
    }
    return(fugacity);
}

double EOS_PCE::get_cs2(double e, double rhob, double proper_tau) const {
    double T = get_temperature(e, rhob, proper_tau);
    double cs2_QCD = interpolate1D(T, 0, cs2_tb);  // 1/fm
    double cs2_gluon = interpolate1D(T, 1, cs2_tb);
    double fugacity = get_fugacity(proper_tau);
    double cs2 = fugacity * cs2_QCD + (1 - fugacity) * cs2_gluon;
    return(std::max(1e-15, cs2));
}
