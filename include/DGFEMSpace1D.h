/**
 * @file DGFEMSpace1D.h
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2018-04-26
 */

#ifndef DGFEMSPACE1D_H
#define DGFEMSPACE1D_H

#include "Eigen/Dense"
#include "para.h"
#include "BasFun.h"
#include "Quadrature.h"
#include "interval_crd_trs.h"
#include "engine.h"

class DGFEMSpace1D {
  private:
    u_int Nx;
    double xl, xr;
    double h;
    VEC<double> mesh;
    TemplateQuadrature TemQuad_Gauss;
    VEC<Quadrature> QUADINFO;
    VEC<double> mu;
    VEC<double> wgt;
    double sum_wm;
    EMAT Kt, Kxl, Fxl, Kxr, Fxr, Ft, MM;//universal matrices
    SOL I, I1, In, I_last;//dimension: Nx*M*K
    VEC<VEC<double>> I_av;
    VEC<ST_ele> W_list, W_list_last;
    BM bml, bmr;
    double BD_L, BD_R;
    EVEC FLUX;
    EMAT A;
    EVEC rhs;
    Eigen::ColPivHouseholderQR<EMAT> solver;

  public:
    DGFEMSpace1D(u_int Nx, double xl, double xr);
    void BuildQuad(u_int np);
    void st2coe(const u_int q, u_int& qx, u_int& qt);
    void BuildUniMAT();
    void Projection(u_int cell, func I0, double t, bU&);
    double Composition(const SOL&, u_int cell, u_int m, double x);
    void init(func I0);
    double cal_dt(const SOL&, const double);
    void update_sol(const VEC<ST_ele>& W, SOL& I);
    double RAD_BE_unsteady(const SOL& In, SOL& I,
        func, const double, const double, func, func);
    void Boundary(func BL, func BR, double& BD_L, double& BD_R, const double mu, const double t, const double dt);
    EVEC STDG_BD(const u_int, const u_int m, const double t, const double dt, func BL);
    void solve_leqn(const EMAT&, const EVEC&, EVEC&);
    void STDG(const SOL&,const double t, const double dt, const u_int, const u_int, func, func, ST_ele& W);
    void STDG(const SOL&,const double t, const double dt, const u_int, const u_int, const ST_ele&, func, ST_ele& W);
    EVEC source(func q, u_int cell, u_int m,
        const double t, const double dt);
    void DG2av(const SOL& I, VEC<VEC<double>>& I_av);
    void find_trouble_cell(const SOL& I, const u_int m, VEC<u_int>& TC);
    void Reconstruct(const SOL& I0, SOL& I);
    void Reconstruct_TC(const SOL& I0, SOL& I, std::ostream&);
    void run_unsteady(func, func, func, double t_end);
    double cal_res(const SOL& s1, const SOL& s2);
    double cal_res(const VEC<ST_ele>& w1, const VEC<ST_ele>& w2);
    void print_solution_integral(const SOL&, std::ostream&);
    void print_solution_average(std::ostream&);
    void print_sol_ex_integral(std::ostream& os, func, const double t);
    int plot_I(engine*, const VEC<double>& mesh, const SOL& In);
    void scaling_limiter(const u_int, EVEC& u);
};

#endif //DGFEMSPACE1D_H

