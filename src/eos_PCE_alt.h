// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_PCE_H_
#define SRC_EOS_PCE_H_

#include "eos_base.h"
#include "data.h"

class EOS_PCE : public EOS_base {
 private:
    double T_max;
    const InitData &DATA;
   
 public:
    EOS_PCE(const InitData &DATA_in);
    
    void initialize_eos();
    double p_e_func       (double e, double rhob, double proper_tau) const;
    double get_temperature(double e, double rhob, double proper_tau) const;
    double get_pressure   (double e, double rhob, double proper_tau) const;
    double get_s2e        (double s, double rhob, double proper_tau) const;
    double get_T2e        (double T, double rhob, double proper_tau) const;
    double get_fugacity  (double proper_tau) const;
  //double get_cs2        (double e, double rhob, double proper_tau) const;
    double get_entropy   (double e, double rhob, double proper_tau) const;
    // returns maximum temperature of the EoS table
    // in the unit of [1/fm]
    void   set_T_max(double T_max_in) {T_max = T_max_in;}
    double get_T_max() const {return(T_max);}

    void check_eos() const {check_eos_no_muB();}
};

#endif  // SRC_EOS_PCE_H_
