#include <stdio.h>
#include <mpi.h>
#include "ftcs.h"

const int T_MAX = 13;

void init(Formura_Navigator &navi) {
  for(int x = navi.lower_x; x < navi.upper_x/2; ++x) {
    U[x] = 1.0;
  }
}

int main(int argc, char **argv) {
  Formura_Navigator navi;
  MPI_Init(&argc, &argv);
  Formura_Init(&navi, MPI_COMM_WORLD);

  FILE *fp = fopen("data/lax-0.5.dat", "w");

  init(navi);
  while(navi.time_step < T_MAX) {
    printf("t = %d\n", navi.time_step);
    for(int x = navi.lower_x; x < navi.upper_x; ++x) {
      int x1 = x + navi.offset_x;
      fprintf(fp, "%d %f\n", x1, U[x]);
    }
    fprintf(fp, "\n\n");
    Formura_Forward(&navi);
  }

  fclose(fp);

  MPI_Finalize();
}


