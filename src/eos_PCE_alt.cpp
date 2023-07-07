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
   
    set_number_of_tables(4);
    resize_table_info_arrays();

    int ntables = get_number_of_tables();

    pressure_tb    = new double** [ntables];
    energy_tb      = new double** [ntables]; 
    entropy_tb     = new double** [ntables];
    temperature_tb = new double** [ntables];
    
    std::ifstream qcd_eos(path + "/QGP_eos_Vovchenko_Tspacing.dat", std::ios::binary);
    std::ifstream gluon_eos(path + "/gluon_eos_binary_Tspacing.dat", std::ios::binary);
    std::ifstream qcd_eos_espacing(path + "/QGP_eos_Vovchenko_e0p25spacing.dat", std::ios::binary); // For fast T(e) lookup
    std::ifstream gluon_eos_espacing(path + "/gluon_eos_binary_e0p25spacing.dat", std::ios::binary); 
    
    for (int itable = 0; itable < ntables; itable++) { // 0 = QCD, 1 = gluon
        std::ifstream* eos_file;
        switch(itable) {
	    case 0:
	         eos_file = &qcd_eos;
		 break;
	    case 1:
	         eos_file = &gluon_eos;
		 break;
            case 2:
	         eos_file = &qcd_eos_espacing;
                 break;
	    case 3:
	         eos_file = &gluon_eos_espacing;
                 break;
	}

        if (!*eos_file) {
            music_message.error("Can not find the EoS file.");
            exit(1);
        }

        e_length[itable]  = 100000;
        nb_length[itable] = 1;
        // allocate memory for pressure arrays
	if (itable == 0 or itable == 1) {
            pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
            energy_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
	    entropy_tb[itable] = Util::mtx_malloc(nb_length[itable], e_length[itable]);
	
            double temp;
            for (int ii = 0; ii < e_length[itable]; ii++) {
	        eos_file->read((char*)&temp, sizeof(double));  // e
                energy_tb[itable][0][ii] = temp/Util::hbarc;      // 1/fm^4
	        //    if (ii == 0) e_bounds[itable] = temp;
                //    if (ii == 1) e_spacing[itable] = temp - e_bounds[itable];
                if (ii == e_length[itable] - 1) set_eps_max(std::max(get_eps_max(),temp/Util::hbarc));

                eos_file->read((char*)&temp, sizeof(double));  // P
                pressure_tb[itable][0][ii] = temp/Util::hbarc;      // 1/fm^4

                eos_file->read((char*)&temp, sizeof(double));  // s
	        entropy_tb[itable][0][ii] = temp;

                eos_file->read((char*)&temp, sizeof(double));  // T
                if (ii == 0) e_bounds[itable] = temp/Util::hbarc;   // 1/fm
	        if (ii == 1) e_spacing[itable] = temp/Util::hbarc - e_bounds[itable];
	        if (ii == e_length[itable] -1) set_T_max(std::max(T_max,temp/Util::hbarc));
	    }
        }
	
	else {
            pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
            temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
	    entropy_tb[itable] = Util::mtx_malloc(nb_length[itable], e_length[itable]);

            double temp;
            for (int ii = 0; ii < e_length[itable]; ii++) {
	        eos_file->read((char*)&temp, sizeof(double));  // e^(1/4) (1/fm)
                //temp /= Util::hbarc;      // 1/fm^4
                if (ii == 0) e_bounds[itable] = temp;
                if (ii == 1) e_spacing[itable] = temp - e_bounds[itable];
                if (ii == e_length[itable] - 1) set_eps_max(std::max(get_eps_max(),std::pow(temp,4)));

                eos_file->read((char*)&temp, sizeof(double));  // P
                pressure_tb[itable][0][ii] = temp/Util::hbarc;      // 1/fm^4

                eos_file->read((char*)&temp, sizeof(double));  // s
	        entropy_tb[itable][0][ii] = temp;                   // 1/fm^3

                eos_file->read((char*)&temp, sizeof(double));  // T
                temperature_tb[itable][0][ii] = temp/Util::hbarc;   // 1/fm

	        eos_file->read((char*)&temp, sizeof(double));  // cs2
	    }
	}
    }
    //std::cout << e_bounds[0] << " " << e_spacing[0] << " " << e_bounds[1] << " " << e_spacing[1] << std::endl;
    music_message.info("Done reading EOS.");
}


double EOS_PCE::p_e_func(double e, double rhob, double proper_tau) const {
    return(get_dpOverde3(e, rhob, proper_tau));
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
// Uses binary search as in get_T2e_finite_rhob - Andrew
double EOS_PCE::get_temperature(double e, double rhob, double proper_tau) const {
    double e_goal = e;         // convert to 1/fm^4
    double T_QCD = interpolate1D(std::pow(e,0.25), 2, temperature_tb);  // 1/fm
    double T_gluon = interpolate1D(std::pow(e,0.25), 3, temperature_tb);
    double T_lower = std::min(T_QCD,T_gluon);//1e-15;  
    double T_upper = std::max(T_QCD,T_gluon);//T_max;
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

    double rel_accuracy = 1e-15;
    double abs_accuracy = 1e-8;
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
    double tau0 = DATA.tau0;
    double tau_eq = DATA.tau_eq;
    if (tau0 <= 0 or tau_eq <= 0) {
        fugacity = 1;
    }
    else {
        fugacity = 1 - exp((tau0 - proper_tau)/tau_eq);
	if (fugacity < 0) return(0);
    }
    return(fugacity);
}

// double EOS_PCE::get_cs2(double e, double rhob, double proper_tau) const {
//     double T = get_temperature(e, rhob, proper_tau);
//     double cs2_QCD = interpolate1D(T, 0, cs2_tb);  // 1/fm
//     double cs2_gluon = interpolate1D(T, 1, cs2_tb);
//     double fugacity = get_fugacity(proper_tau);
//     double cs2 = fugacity * cs2_QCD + (1 - fugacity) * cs2_gluon;
//     return(std::max(1e-15, cs2));
// }

double EOS_PCE::get_entropy(double e, double rhob, double proper_tau) const {
    double T = get_temperature(e, rhob, proper_tau);
    double P_QCD = interpolate1D(T, 0, pressure_tb);
    double P_gluon = interpolate1D(T, 1, pressure_tb);
    double entropy_QCD = interpolate1D(T, 0, entropy_tb);
    double entropy_gluon = interpolate1D(T, 1, entropy_tb);
    double fugacity = get_fugacity(proper_tau);
    //entropy_QCD = (e + get_pressure(e, rhob, proper_tau))/(get_temperature(e, rhob, proper_tau));
    double entropy = fugacity * entropy_QCD + (1 - fugacity) * entropy_gluon - fugacity*log(fugacity)/T * (P_QCD - P_gluon);
    return(std::max(1e-15,entropy));
}
