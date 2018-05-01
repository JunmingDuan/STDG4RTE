#include "DGFEMSpace1D.h"
#include "conWENO3.h"
#include "conWENO5.h"
#include "conWENO7.h"

double minmod(const double a1, const double a2, const double a3) {
  if(a1 > 0 && a2 > 0 && a3 > 0) return std::min(std::min(a1,a2),a3);
  else if(a1 < 0 && a2 < 0 && a3 < 0) return std::max(std::max(a1,a2),a3);
  else return 0;
}

double TVB(const double a1, const double a2, const double a3, const double M) {
  if(fabs(a1) < M) return a1;
  else return minmod(a1, a2, a3);
}

void DGFEMSpace1D::find_trouble_cell(const SOL& I, const u_int m, VEC<u_int>& TC) {
  DG2av(I, I_av);
  const double M = 50;
  double ul, ur, ut1, ut2, ut1_mod, ut2_mod, FD, BD;
  for(u_int i = 0; i < Nx; ++i) {
    ul = Composition(I, i, m, mesh[i]);
    ur = Composition(I, i, m, mesh[i+1]);
    if(i == 0) {
      FD = I_av[i+1][m]-I_av[i][m];
      BD = FD;
    }
    else if(i == Nx-1) {
      BD = I_av[i][m]-I_av[i-1][m];
      FD = BD;
    }
    else {
      FD = I_av[i+1][m]-I_av[i][m];
      BD = I_av[i][m]-I_av[i-1][m];
    }
    ut1 = ur-I_av[i][m];
    ut2 = I_av[i][m]-ul;
    ut1_mod = TVB(ut1, FD, BD, M*h*h);
    ut2_mod = TVB(ut2, FD, BD, M*h*h);
    if(fabs(ut1-ut1_mod) + fabs(ut2-ut2_mod) > 1e-14) TC.push_back(i);
  }
}

void DGFEMSpace1D::Reconstruct_TC(const SOL& I0, SOL& I, std::ostream& os) {
  I = I0;
  if(K == 1) {
  }
  else if(K == 2) {//WENO3, period
    u_int n1, n2, n3;
    DG2av(I0, I_av);
    for(u_int m = 0; m < M; ++m) {
      VEC<u_int> TC;
      find_trouble_cell(I, m, TC);
      std::cout << "Number of the trouble cells: " << TC.size() << std::endl;
      for(u_int i = 0; i < TC.size(); ++i) {
        if(I_av[i][m] < 0) {
          std::cout << "Negative cell average: " << i << " " << I_av[i][m] << std::endl;
          abort();
        }
        u_int cell = TC[i];
        n1 = (cell-1+Nx)%Nx;
        n2 = (cell+Nx)%Nx;
        n3 = (cell+1+Nx)%Nx;
        I[cell][m][0] = WENONI21(I_av[n1][m], I_av[n2][m], I_av[n3][m]);
        I[cell][m][1] = WENONI22(I_av[n1][m], I_av[n2][m], I_av[n3][m]);
      }
    }
  }
  else if(K == 3) {//WENO5, period
    u_int n1, n2, n3, n4, n5;
    DG2av(I0, I_av);
    for(u_int m = 0; m < M; ++m) {
      VEC<u_int> TC;
      find_trouble_cell(I, m, TC);
      std::cout << "Number of the trouble cells: " << TC.size() << std::endl;
      for(u_int i = 0; i < TC.size(); ++i) {
        if(I_av[i][m] < 0) {
          std::cout << "Negative cell average: " << i << " " << I_av[i][m] << std::endl;
          abort();
        }
        u_int cell = TC[i];
        n1 = (cell-2+Nx)%Nx;
        n2 = (cell-1+Nx)%Nx;
        n3 = (cell+Nx)%Nx;
        n4 = (cell+1+Nx)%Nx;
        n5 = (cell+2+Nx)%Nx;
        I[cell][m][0] = WENONI31(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m]);
        I[cell][m][1] = WENONI32(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m]);
        I[cell][m][2] = WENONI33(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m]);
      }
    }
  }
  else if(K == 4) {//WENO7, period
    u_int n1, n2, n3, n4, n5, n6, n7;
    DG2av(I0, I_av);
    for(u_int m = 0; m < M; ++m) {
      VEC<u_int> TC;
      find_trouble_cell(I, m, TC);
      std::cout << "Number of the trouble cells: " << TC.size() << std::endl;
      os << TC << std::endl;
      for(u_int i = 0; i < TC.size(); ++i) {
        if(I_av[i][m] < 0) {
          std::cout << "Negative cell average: " << i << " " << I_av[i][m] << std::endl;
          abort();
        }
        u_int cell = TC[i];
        n1 = (cell-3+Nx)%Nx;
        n2 = (cell-2+Nx)%Nx;
        n3 = (cell-1+Nx)%Nx;
        n4 = (cell+Nx)%Nx;
        n5 = (cell+1+Nx)%Nx;
        n6 = (cell+2+Nx)%Nx;
        n7 = (cell+3+Nx)%Nx;
        I[cell][m][0] = WENONI41(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m], I_av[n6][m], I_av[n7][m]);
        I[cell][m][1] = WENONI42(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m], I_av[n6][m], I_av[n7][m]);
        I[cell][m][2] = WENONI43(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m], I_av[n6][m], I_av[n7][m]);
        I[cell][m][3] = WENONI44(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m], I_av[n6][m], I_av[n7][m]);
      }
    }
  }
}


void DGFEMSpace1D::Reconstruct(const SOL& I0, SOL& I) {
  if(K == 1) {
    for(u_int m = 0; m < M; ++m) {
      for(u_int i = 0; i < Nx; ++i) {
        I[i][m][0] = I0[i][m][0];
      }
    }
  }
  else if(K == 2) {//WENO3, period
    u_int n1, n2, n3;
    DG2av(I0, I_av);
    for(u_int m = 0; m < M; ++m) {
      for(u_int i = 0; i < Nx; ++i) {
        n1 = (i-1+Nx)%Nx;
        n2 = (i+Nx)%Nx;
        n3 = (i+1+Nx)%Nx;
        I[i][m][0] = WENONI21(I_av[n1][m], I_av[n2][m], I_av[n3][m]);
        I[i][m][1] = WENONI22(I_av[n1][m], I_av[n2][m], I_av[n3][m]);
      }
    }
  }
  else if(K == 3) {//WENO5, period
    u_int n1, n2, n3, n4, n5;
    DG2av(I0, I_av);
    for(u_int m = 0; m < M; ++m) {
      for(u_int i = 0; i < Nx; ++i) {
        n1 = (i-2+Nx)%Nx;
        n2 = (i-1+Nx)%Nx;
        n3 = (i+Nx)%Nx;
        n4 = (i+1+Nx)%Nx;
        n5 = (i+2+Nx)%Nx;
        I[i][m][0] = WENONI31(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m]);
        I[i][m][1] = WENONI32(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m]);
        I[i][m][2] = WENONI33(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m]);
      }
    }
  }
  else if(K == 4) {//WENO7, period
    u_int n1, n2, n3, n4, n5, n6, n7;
    DG2av(I0, I_av);
    for(u_int m = 0; m < M; ++m) {
      for(u_int i = 0; i < Nx; ++i) {
        n1 = (i-3+Nx)%Nx;
        n2 = (i-2+Nx)%Nx;
        n3 = (i-1+Nx)%Nx;
        n4 = (i+Nx)%Nx;
        n5 = (i+1+Nx)%Nx;
        n6 = (i+2+Nx)%Nx;
        n7 = (i+3+Nx)%Nx;
        I[i][m][0] = WENONI41(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m], I_av[n6][m], I_av[n7][m]);
        I[i][m][1] = WENONI42(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m], I_av[n6][m], I_av[n7][m]);
        I[i][m][2] = WENONI43(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m], I_av[n6][m], I_av[n7][m]);
        I[i][m][3] = WENONI44(I_av[n1][m], I_av[n2][m], I_av[n3][m], I_av[n4][m], I_av[n5][m], I_av[n6][m], I_av[n7][m]);
      }
    }
  }
}

