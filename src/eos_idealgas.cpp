// Copyright 2018 @ Chun Shen
#include "eos_idealgas.h"

#include <cmath>

EOS_idealgas::EOS_idealgas() {
    set_EOS_id(0);
    set_number_of_tables(0);
    set_eps_max(1e5);
    Nc = 3.;
    Nf = 2.5;
    set_flag_muB(false);
    set_flag_muS(false);
    set_flag_muC(false);
}

void EOS_idealgas::initialize_eos() {
    music_message.info("initialze EOS ideal gas ...");
}

double EOS_idealgas::get_temperature(double eps, double rhob, double proper_tau) const {
    return pow(90.0/M_PI/M_PI*(eps/3.0)/(2*(Nc*Nc-1)+7./2*Nc*Nf), .25);
}

double EOS_idealgas::get_s2e(double s, double rhob, double proper_tau) const {
    return(3./4.*s*pow(3.*s/4./(M_PI*M_PI*3.0*(2*(Nc*Nc-1)+7./2*Nc*Nf)/90.0), 1./3.));  // in 1/fm^4
}


double EOS_idealgas::get_T2e(double T, double rhob, double proper_tau) const {
    return 3*T*T*T*T*M_PI*M_PI/90*(2*(Nc*Nc-1)+7./2*Nc*Nf);
}
