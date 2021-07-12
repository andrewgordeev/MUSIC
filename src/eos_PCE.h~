// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_PCE_H_
#define SRC_EOS_PCE_H_

#include "eos_base.h"

class EOS_PCE : public EOS_base {
 private:
   
 public:
    EOS_PCE();
    
    void initialize_eos();
    double p_e_func       (double e, double rhob, double proper_tau) const;
    double get_temperature(double e, double rhob, double proper_tau) const;
    double get_pressure   (double e, double rhob, double proper_tau) const;
    double get_s2e        (double s, double rhob, double proper_tau) const;
    double get_T2e        (double T, double rhob, double proper_tau) const;
    double get_fugacity  (double proper_tau) const;
    double get_cs2        (double e, double rhob, double proper_tau) const;

    void check_eos() const {check_eos_no_muB();}
};

#endif  // SRC_EOS_PCE_H_
