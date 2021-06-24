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
    temperature_tb = new double** [ntables]; 
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
        temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
	cs2_tb[itable] = Util::mtx_malloc(nb_length[itable], e_length[itable]);
        double temp;
        for (int ii = 0; ii < e_length[itable]; ii++) {
	    eos_file->read((char*)&temp, sizeof(double));  // e^(1/4)
            //temp /= Util::hbarc;      // 1/fm
            if (ii == 0) e_bounds[itable] = temp;
            if (ii == 1) e_spacing[itable] = temp - e_bounds[itable];
            if (ii == e_length[itable] - 1) set_eps_max(temp);

            eos_file->read((char*)&temp, sizeof(double));  // P
            pressure_tb[itable][0][ii] = temp/Util::hbarc;      // 1/fm^4

            eos_file->read((char*)&temp, sizeof(double));  // s

            eos_file->read((char*)&temp, sizeof(double));  // T
            temperature_tb[itable][0][ii] = temp/Util::hbarc;   // 1/fm

	    eos_file->read((char*)&temp, sizeof(double));  // cs2
	    cs2_tb[itable][0][ii] = temp;
        }
    }
    music_message.info("Done reading EOS.");
}


double EOS_PCE::p_e_func(double e, double rhob, double proper_tau) const {
    return(get_dpOverde3(e, rhob, proper_tau));
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_PCE::get_temperature(double e, double rhob, double proper_tau) const {
    double T_QCD;
    //   if (e < 0.1) {  // This corresponds to T < ~0.1 GeV
    //    T_QCD = (0.18175818963347465 * std::pow(e*Util::hbarc, 0.1433156325035367))/Util::hbarc;
    //}
    //else T_QCD = interpolate1D(e, 0, temperature_tb);  // 1/fm
    T_QCD = interpolate1D(std::pow(e,0.25), 0, temperature_tb);  // 1/fm
    double T_gluon;
    //if (e < 0.1) {  // This corresponds to T < ~0.2 GeV
    //   T_gluon = (0.545827055067574 * std::pow(e*Util::hbarc, 0.2237373605050102))/Util::hbarc;
    //}
    //else T_gluon = interpolate1D(e, 1, temperature_tb);
    T_gluon = interpolate1D(std::pow(e,0.25), 1, temperature_tb);
    double fugacity = get_fugacity(proper_tau);
    double T = fugacity * T_QCD + (1 - fugacity) * T_gluon;

    //T = (.15274608127810668 * std::pow(e*Util::hbarc, .2478330589996049))/Util::hbarc;
    return(std::max(1e-15, T));
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS_PCE::get_pressure(double e, double rhob, double proper_tau) const {
  // double T = get_temperature(e, rhob, proper_tau);
    double f_QCD = interpolate1D(std::pow(e,0.25), 0, pressure_tb);  // 1/fm^4
    double f_gluon = interpolate1D(std::pow(e,0.25), 1, pressure_tb);
    double fugacity = get_fugacity(proper_tau);
    double f = fugacity * f_QCD + (1 - fugacity) * f_gluon;

    //f = (-.0000000775447 * std::pow(e*Util::hbarc,0.247833) + 0.329492 * e*Util::hbarc)/Util::hbarc;
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
    double tau0 = 0.0;
    double tau_eq = 5.0;
    if (tau0 <= 0 or tau_eq <= 0) {
        fugacity = 0;
    }
    else {
        fugacity = 1 - exp((tau0 - proper_tau)/tau_eq);
    }
    return(fugacity);
}

double EOS_PCE::get_cs2(double e, double rhob, double proper_tau) const {
  // double T = get_temperature(e, rhob, proper_tau);
    double cs2_QCD;
    //   if (e < 0.08445) {  // This corresponds to T < 0.1 GeV
    //    double T = interpolate1D(e, 0, temperature_tb)*Util::hbarc;
    ///    cs2_QCD = std::max(285.383289 * std::pow(T,3) - 91.6171247 * std::pow(T,2) + 8.44067621 * T - 0.00689668358, 1e-10); // Cubic fit for low temp cs2
	//  	std::cout << "cs2 debug: " << e << " " << T << " " << cs2_QCD << std::endl;
    //}
    //else
    cs2_QCD = interpolate1D(std::pow(e,0.25), 0, cs2_tb);  // 1/fm
    double cs2_gluon = interpolate1D(std::pow(e,0.25), 1, cs2_tb);
    double fugacity = get_fugacity(proper_tau);
    double cs2 = fugacity * cs2_QCD + (1 - fugacity) * cs2_gluon;
    //double T = interpolate1D(e,0,temperature_tb)*Util::hbarc;
    //if (T < 0.18410179) cs2 = 246.0490286078087 * std::pow(T,3) - 83.01937867454492 * std::pow(T,2) + 7.936821067345355 * T - 0.00025794218162239854;
    //else cs2 = 1/3. - 0.008183471854138312/(T - 0.12997412301945538);
    return(std::max(1e-15, cs2));
}
