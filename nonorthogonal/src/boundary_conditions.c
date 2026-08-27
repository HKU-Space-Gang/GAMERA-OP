#include <omp.h>
#include "common.h"
#include "config.h"
#include "setup_mpi.h"
#include "utils.h"
#include "math.h"
#include "problem.h"

void gas_bc_symmetric_i_low() {
  if (nbr_low[0] == MPI_PROC_NULL) {
    // symmetric boundary condition for i lowest boundary
    #pragma omp parallel for collapse(5) schedule(static)
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < NG; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < config.NK; k++) {
              gas[s][f][i][j][k] =
                  gas[s][f][2 * NG - 1 - i + gas_onface_i[f]][j][k];
            }
          }
        }
      }
    }
  }
}
void gas_bc_symmetric_i_high() {
  if (nbr_high[0] == MPI_PROC_NULL) {
    // symmetric boundary condition for i highest boundary
    #pragma omp parallel for collapse(5) schedule(static)
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < NG; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < config.NK; k++) {
              gas[s][f][config.NI - 1 - NG + i + gas_onface_i[f]][j][k] =
                  gas[s][f][config.NI - 1 - NG - 1 - i][j][k];
            }
          }
        }
      }
    }
  }
}
void gas_bc_symmetric_j_low() {
  if (nbr_low[1] == MPI_PROC_NULL) {
    // symmetric boundary condition for j lowest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < NG; j++) {
            for (int k = 0; k < config.NK; k++) {
              gas[s][f][i][j][k] =
                  gas[s][f][i][2 * NG - 1 - j + gas_onface_j[f]][k];
            }
          }
        }
      }
    }
  }
}
void gas_bc_symmetric_j_high() {
  if (nbr_high[1] == MPI_PROC_NULL) {
    // symmetric boundary condition for j highest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < NG; j++) {
            for (int k = 0; k < config.NK; k++) {
              gas[s][f][i][config.NJ - 1 - NG + j + gas_onface_j[f]][k] =
                  gas[s][f][i][config.NJ - 1 - NG - 1 - j][k];
            }
          }
        }
      }
    }
  }
}
void gas_bc_symmetric_k_low() {
  if (nbr_low[2] == MPI_PROC_NULL) {
    // symmetric boundary condition for k lowest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < NG; k++) {
              gas[s][f][i][j][k] =
                  gas[s][f][i][j][2 * NG - 1 - k + gas_onface_k[f]];
            }
          }
        }
      }
    }
  }
}
void gas_bc_symmetric_k_high() {
  if (nbr_high[2] == MPI_PROC_NULL) {
    // symmetric boundary condition for k highest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < NG; k++) {
              gas[s][f][i][j][config.NK - 1 - NG + k + gas_onface_k[f]] =
                  gas[s][f][i][j][config.NK - 1 - NG - 1 - k];
            }
          }
        }
      }
    }
  }
}

void gas_bc_reflective_i_low() {
  if (nbr_low[0] == MPI_PROC_NULL) {
    // reflective boundary condition for i lowest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < NG; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < config.NK; k++) {
              if (f == gas_v1) {  // for Vx
                gas[s][f][i][j][k] =
                    -gas[s][f][2 * NG - 1 - i + gas_onface_i[f]][j][k];
              } else {
                gas[s][f][i][j][k] =
                    gas[s][f][2 * NG - 1 - i + gas_onface_i[f]][j][k];
              }
            }
          }
        }
      }
    }
  }
}
void gas_bc_reflective_i_high() {
  if (nbr_high[0] == MPI_PROC_NULL) {
    // reflective boundary condition for i highest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < NG; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < config.NK; k++) {
              if (f == gas_v1) {
                gas[s][f][config.NI - 1 - NG + i + gas_onface_i[f]][j][k] =
                    -gas[s][f][config.NI - 1 - NG - 1 - i][j][k];
              } else {
                gas[s][f][config.NI - 1 - NG + i + gas_onface_i[f]][j][k] =
                    gas[s][f][config.NI - 1 - NG - 1 - i][j][k];
              }
            }
          }
        }
      }
    }
  }
}
void gas_bc_reflective_j_low() {
  if (nbr_low[1] == MPI_PROC_NULL) {
    // reflective boundary condition for j lowest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < NG; j++) {
            for (int k = 0; k < config.NK; k++) {
              if (f == gas_v2) {
                gas[s][f][i][j][k] =
                    -gas[s][f][i][2 * NG - 1 - j + gas_onface_j[f]][k];
              } else {
                gas[s][f][i][j][k] =
                    gas[s][f][i][2 * NG - 1 - j + gas_onface_j[f]][k];
              }
            }
          }
        }
      }
    }
  }
}
void gas_bc_reflective_j_high() {
  if (nbr_high[1] == MPI_PROC_NULL) {
    // reflective boundary condition for j highest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < NG; j++) {
            for (int k = 0; k < config.NK; k++) {
              if (f == gas_v2) {
                gas[s][f][i][config.NJ - 1 - NG + j + gas_onface_j[f]][k] =
                    -gas[s][f][i][config.NJ - 1 - NG - 1 - j][k];
              } else {
                gas[s][f][i][config.NJ - 1 - NG + j + gas_onface_j[f]][k] =
                    gas[s][f][i][config.NJ - 1 - NG - 1 - j][k];
              }
            }
          }
        }
      }
    }
  }
}
void gas_bc_reflective_k_low() {
  if (nbr_low[2] == MPI_PROC_NULL) {
    // reflective boundary condition for k lowest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < NG; k++) {
              if (f == gas_v3) {
                gas[s][f][i][j][k] =
                    -gas[s][f][i][j][2 * NG - 1 - k + gas_onface_k[f]];
              } else {
                gas[s][f][i][j][k] =
                    gas[s][f][i][j][2 * NG - 1 - k + gas_onface_k[f]];
              }
            }
          }
        }
      }
    }
  }
}
void gas_bc_reflective_k_high() {
  if (nbr_high[2] == MPI_PROC_NULL) {
    // symmetric boundary condition for k highest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < NG; k++) {
              if (f == gas_v3) {
                gas[s][f][i][j][config.NK - 1 - NG + k + gas_onface_k[f]] =
                    -gas[s][f][i][j][config.NK - 1 - NG - 1 - k];
              } else {
                gas[s][f][i][j][config.NK - 1 - NG + k + gas_onface_k[f]] =
                    gas[s][f][i][j][config.NK - 1 - NG - 1 - k];
              }
            }
          }
        }
      }
    }
  }
}

void gas_bc_fixed_i_low() {
  if (nbr_low[0] == MPI_PROC_NULL) {
    // constant fixed boundary condition for i lowest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < NG; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < config.NK; k++) {
              if (f == gas_rho) {  // for rho
                gas[s][f][i][j][k] = 1.0;
              } else if (f == gas_v1) {  // for Vx
                gas[s][f][i][j][k] = 0.0;
              } else if (f == gas_v2) {  // for Vy
                gas[s][f][i][j][k] = 0.0;
              } else if (f == gas_v3) {  // for Vz
                gas[s][f][i][j][k] = 0.0;
              } else if (f == gas_p) {  // for P
                gas[s][f][i][j][k] = 1.0;
              }
            }
          }
        }
      }
    }
  }
}
void gas_bc_fixed_i_high() {
  if (nbr_high[0] == MPI_PROC_NULL) {
    // constant fixed boundary condition for i highest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < NG; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < config.NK; k++) {
              if (f == gas_rho) {  // for rho
                gas[s][f][config.NI - 1 - NG + i][j][k] = 1.0;
              } else if (f == gas_v1) {  // for Vx
                gas[s][f][config.NI - 1 - NG + i][j][k] = 0.0;
              } else if (f == gas_v2) {  // for Vy
                gas[s][f][config.NI - 1 - NG + i][j][k] = 0.0;
              } else if (f == gas_v3) {  // for Vz
                gas[s][f][config.NI - 1 - NG + i][j][k] = 0.0;
              } else if (f == gas_p) {  // for P
                gas[s][f][config.NI - 1 - NG + i][j][k] = 1.0;
              }
            }
          }
        }
      }
    }
  }
}
void gas_bc_fixed_j_low() {
  if (nbr_low[1] == MPI_PROC_NULL) {
    // constant fixed boundary condition for j lowest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < NG; j++) {
            for (int k = 0; k < config.NK; k++) {
              if (f == gas_rho) {  // for rho
                gas[s][f][i][j][k] = 2.0;
              } else if (f == gas_v1) {  // for Vx
                gas[s][f][i][j][k] = 0.0;
              } else if (f == gas_v2) {  // for Vy
                gas[s][f][i][j][k] = 0.0;
              } else if (f == gas_v3) {  // for Vz
                gas[s][f][i][j][k] = 0.0;
              } else if (f == gas_p) {  // for P
                gas[s][f][i][j][k] = 1.0;
              }
            }
          }
        }
      }
    }
  }
}
void gas_bc_fixed_j_high() {
  if (nbr_high[1] == MPI_PROC_NULL) {
    // constant fixed boundary condition for j highest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < NG; j++) {
            for (int k = 0; k < config.NK; k++) {
              if (f == gas_rho) {  // for rho
                gas[s][f][i][config.NJ - 1 - NG + j][k] = 1.0;
              } else if (f == gas_v1) {  // for Vx
                gas[s][f][i][config.NJ - 1 - NG + j][k] = 0.0;
              } else if (f == gas_v2) {  // for Vy
                gas[s][f][i][config.NJ - 1 - NG + j][k] = 0.0;
              } else if (f == gas_v3) {  // for Vz
                gas[s][f][i][config.NJ - 1 - NG + j][k] = 0.0;
              } else if (f == gas_p) {  // for P
                gas[s][f][i][config.NJ - 1 - NG + j][k] = 2.5;
              }
            }
          }
        }
      }
    }
  }
}
void gas_bc_fixed_k_low() {
  if (nbr_low[2] == MPI_PROC_NULL) {
    // constant fixed boundary condition for k lowest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < NG; k++) {
              if (f == gas_rho) {  // for rho
                gas[s][f][i][j][k] = 1.0;
              } else if (f == gas_v1) {  // for Vx
                gas[s][f][i][j][k] = 0.0;
              } else if (f == gas_v2) {  // for Vy
                gas[s][f][i][j][k] = 0.0;
              } else if (f == gas_v3) {  // for Vz
                gas[s][f][i][j][k] = 0.0;
              } else if (f == gas_p) {  // for P
                gas[s][f][i][j][k] = 1.0;
              }
            }
          }
        }
      }
    }
  }
}
void gas_bc_fixed_k_high() {
  if (nbr_high[2] == MPI_PROC_NULL) {
    // constant fixed boundary condition for k highest boundary
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < NG; k++) {
              if (f == gas_rho) {  // for rho
                gas[s][f][i][j][config.NK - 1 - NG + k] = 1.0;
              } else if (f == gas_v1) {  // for Vx
                gas[s][f][i][j][config.NK - 1 - NG + k] = 0.0;
              } else if (f == gas_v2) {  // for Vy
                gas[s][f][i][j][config.NK - 1 - NG + k] = 0.0;
              } else if (f == gas_v3) {  // for Vz
                gas[s][f][i][j][config.NK - 1 - NG + k] = 0.0;
              } else if (f == gas_p) {  // for P
                gas[s][f][i][j][config.NK - 1 - NG + k] = 1.0;
              }
            }
          }
        }
      }
    }
  }
}

void gas_bc_extrapolated_i_low() {
  if (nbr_low[0] == MPI_PROC_NULL) {
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < NG; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < config.NK; k++) {
              gas[s][f][i][j][k] = gas[s][f][NG][j][k];
            }
          }
        }
      }
    }
  }
}
void gas_bc_extrapolated_i_high() {
  if (nbr_high[0] == MPI_PROC_NULL) {
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < NG; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < config.NK; k++) {
              gas[s][f][config.NI - 1 - NG + i + gas_onface_i[f]][j][k] =
                  gas[s][f][config.NI - 1 - NG - 1 + gas_onface_i[f]][j][k];
            }
          }
        }
      }
    }
  }
}
void gas_bc_extrapolated_j_low() {
  if (nbr_low[1] == MPI_PROC_NULL) {
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < NG; j++) {
            for (int k = 0; k < config.NK; k++) {
              gas[s][f][i][j][k] = gas[s][f][i][NG][k];
            }
          }
        }
      }
    }
  }
}
void gas_bc_extrapolated_j_high() {
  if (nbr_high[1] == MPI_PROC_NULL) {
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < NG; j++) {
            for (int k = 0; k < config.NK; k++) {
              gas[s][f][i][config.NJ - 1 - NG + j + gas_onface_j[f]][k] =
                  gas[s][f][i][config.NJ - 1 - NG - 1 + gas_onface_j[f]][k];
            }
          }
        }
      }
    }
  }
}
void gas_bc_extrapolated_k_low() {
  if (nbr_low[2] == MPI_PROC_NULL) {
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < NG; k++) {
              gas[s][f][i][j][k] = gas[s][f][i][j][NG];
            }
          }
        }
      }
    }
  }
}
void gas_bc_extrapolated_k_high() {
  if (nbr_high[2] == MPI_PROC_NULL) {
    for (int s = 0; s < NS1; s++) {
      for (int f = 0; f < NF_gas_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < NG; k++) {
              gas[s][f][i][j][config.NK - 1 - NG + k + gas_onface_k[f]] =
                  gas[s][f][i][j][config.NK - 1 - NG - 1 + gas_onface_k[f]];
            }
          }
        }
      }
    }
  }
}

void gas_bc_none_i_low(){}
void gas_bc_none_i_high(){}
void gas_bc_none_j_low(){}
void gas_bc_none_j_high(){}
void gas_bc_none_k_low(){}
void gas_bc_none_k_high(){}

void gem_bc_symmetric_i_low() {
  if (nbr_low[0] == MPI_PROC_NULL) {
    // symmetric boundary condition for i lowest boundary
    #pragma omp parallel for collapse(4) schedule(static)
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < NG; i++) {
        for (int j = 0; j < config.NJ; j++) {
          for (int k = 0; k < config.NK; k++) {
            gem[f][i][j][k] = gem[f][2 * NG - 1 - i + gem_onface_i[f]][j][k];
          }
        }
      }
    }
  }
}
void gem_bc_symmetric_i_high() {
  if (nbr_high[0] == MPI_PROC_NULL) {
    // symmetric boundary condition for i highest boundary
    #pragma omp parallel for collapse(4) schedule(static)
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < NG; i++) {
        for (int j = 0; j < config.NJ; j++) {
          for (int k = 0; k < config.NK; k++) {
            gem[f][config.NI - 1 - NG + i + gem_onface_i[f]][j][k] =
                gem[f][config.NI - 1 - NG - 1 - i][j][k];
          }
        }
      }
    }
  }
}
void gem_bc_symmetric_j_low() {
  if (nbr_low[1] == MPI_PROC_NULL) {
    // symmetric boundary condition for j lowest boundary
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < config.NI; i++) {
        for (int j = 0; j < NG; j++) {
          for (int k = 0; k < config.NK; k++) {
            gem[f][i][j][k] = gem[f][i][2 * NG - 1 - j + gem_onface_j[f]][k];
          }
        }
      }
    }
  }
}
void gem_bc_symmetric_j_high() {
  if (nbr_high[1] == MPI_PROC_NULL) {
    // symmetric boundary condition for j highest boundary
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < config.NI; i++) {
        for (int j = 0; j < NG; j++) {
          for (int k = 0; k < config.NK; k++) {
            gem[f][i][config.NJ - 1 - NG + j + gem_onface_j[f]][k] =
                gem[f][i][config.NJ - 1 - NG - 1 - j][k];
          }
        }
      }
    }
  }
}
void gem_bc_symmetric_k_low() {
  if (nbr_low[2] == MPI_PROC_NULL) {
    // symmetric boundary condition for k lowest boundary
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < config.NI; i++) {
        for (int j = 0; j < config.NJ; j++) {
          for (int k = 0; k < NG; k++) {
            gem[f][i][j][k] = gem[f][i][j][2 * NG - 1 - k + gem_onface_k[f]];
          }
        }
      }
    }
  }
}
void gem_bc_symmetric_k_high() {
  if (nbr_high[2] == MPI_PROC_NULL) {
    // symmetric boundary condition for k highest boundary
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < config.NI; i++) {
        for (int j = 0; j < config.NJ; j++) {
          for (int k = 0; k < NG; k++) {
            gem[f][i][j][config.NK - 1 - NG + k + gem_onface_k[f]] =
                gem[f][i][j][config.NK - 1 - NG - 1 - k];
          }
        }
      }
    }
  }
}

void gem_bc_reflective_i_low() {}
void gem_bc_reflective_i_high() {}
void gem_bc_reflective_j_low() {}
void gem_bc_reflective_j_high() {}
void gem_bc_reflective_k_low() {}
void gem_bc_reflective_k_high() {}

void gem_bc_fixed_i_low() {
  if (nbr_low[0] == MPI_PROC_NULL) {
    // constant fixed boundary condition for i lowest boundary
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < NG; i++) {
        for (int j = 0; j < config.NJ; j++) {
          for (int k = 0; k < config.NK; k++) {
            gem[f][i][j][k] = 0.0;
          }
        }
      }
    }
  }
}
void gem_bc_fixed_i_high() {
  if (nbr_high[0] == MPI_PROC_NULL) {
    // constant fixed boundary condition for i highest boundary
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < NG; i++) {
        for (int j = 0; j < config.NJ; j++) {
          for (int k = 0; k < config.NK; k++) {
            gem[f][config.NI - 1 - NG + i + gem_onface_i[f]][j][k] = 0.0;
          }
        }
      }
    }
  }
}
void gem_bc_fixed_j_low() {
  if (nbr_low[1] == MPI_PROC_NULL) {
    // constant fixed boundary condition for j lowest boundary
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < config.NI; i++) {
        for (int j = 0; j < NG; j++) {
          for (int k = 0; k < config.NK; k++) {
            gem[f][i][j][k] = 0.0;
          }
        }
      }
    }
  }
}
void gem_bc_fixed_j_high() {
  if (nbr_high[1] == MPI_PROC_NULL) {
    // constant fixed boundary condition for j highest boundary
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < config.NI; i++) {
        for (int j = 0; j < NG; j++) {
          for (int k = 0; k < config.NK; k++) {
            gem[f][i][config.NJ - 1 - NG + j + gem_onface_j[f]][k] = 0.0;
          }
        }
      }
    }
  }
}
void gem_bc_fixed_k_low() {
  if (nbr_low[2] == MPI_PROC_NULL) {
    // constant fixed boundary condition for k lowest boundary
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < config.NI; i++) {
        for (int j = 0; j < config.NJ; j++) {
          for (int k = 0; k < NG; k++) {
            gem[f][i][j][k] = 0.0;
          }
        }
      }
    }
  }
}
void gem_bc_fixed_k_high() {
  if (nbr_high[2] == MPI_PROC_NULL) {
    // constant fixed boundary condition for k highest boundary
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < config.NI; i++) {
        for (int j = 0; j < config.NJ; j++) {
          for (int k = 0; k < NG; k++) {
            gem[f][i][j][config.NK - 1 - NG + k + gem_onface_k[f]] = 0.0;
          }
        }
      }
    }
  }
}

void gem_bc_extrapolated_i_low() {
  if (nbr_low[0] == MPI_PROC_NULL) {
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < NG; i++) {
        for (int j = 0; j < config.NJ; j++) {
          for (int k = 0; k < config.NK; k++) {
            gem[f][i][j][k] = gem[f][NG][j][k];
          }
        }
      }
    }
  }
}
void gem_bc_extrapolated_i_high() {
  if (nbr_high[0] == MPI_PROC_NULL) {
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < NG; i++) {
        for (int j = 0; j < config.NJ; j++) {
          for (int k = 0; k < config.NK; k++) {
            gem[f][config.NI - 1 - NG + i + gem_onface_i[f]][j][k] =
                gem[f][config.NI - 1 - NG - 1 + gem_onface_i[f]][j][k];
          }
        }
      }
    }
  }
}
void gem_bc_extrapolated_j_low() {
  if (nbr_low[1] == MPI_PROC_NULL) {
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < config.NI; i++) {
        for (int j = 0; j < NG; j++) {
          for (int k = 0; k < config.NK; k++) {
            gem[f][i][j][k] = gem[f][i][NG][k];
          }
        }
      }
    }
  }
}
void gem_bc_extrapolated_j_high() {
  if (nbr_high[1] == MPI_PROC_NULL) {
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < config.NI; i++) {
        for (int j = 0; j < NG; j++) {
          for (int k = 0; k < config.NK; k++) {
            gem[f][i][config.NJ - 1 - NG + j + gem_onface_j[f]][k] =
                gem[f][i][config.NJ - 1 - NG - 1 + gem_onface_j[f]][k];
          }
        }
      }
    }
  }
}
void gem_bc_extrapolated_k_low() {
  if (nbr_low[2] == MPI_PROC_NULL) {
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < config.NI; i++) {
        for (int j = 0; j < config.NJ; j++) {
          for (int k = 0; k < NG; k++) {
            gem[f][i][j][k] = gem[f][i][j][NG];
          }
        }
      }
    }
  }
}
void gem_bc_extrapolated_k_high() {
  if (nbr_high[2] == MPI_PROC_NULL) {
    for (int f = 0; f < NF_gem_prim; f++) {
      for (int i = 0; i < config.NI; i++) {
        for (int j = 0; j < config.NJ; j++) {
          for (int k = 0; k < NG; k++) {
            gem[f][i][j][config.NK - 1 - NG + k + gem_onface_k[f]] =
                gem[f][i][j][config.NK - 1 - NG - 1 + gem_onface_k[f]];
          }
        }
      }
    }
  }
}

void gem_bc_none_i_low(){}
void gem_bc_none_i_high(){}
void gem_bc_none_j_low(){}
void gem_bc_none_j_high(){}
void gem_bc_none_k_low(){}
void gem_bc_none_k_high(){}