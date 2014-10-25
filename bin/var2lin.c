
#include "varian.h"

int main (int argc, char **argv) {
  hx_array x;

  if (argc < 2) {
    printf("error: need input file\n");
    return 1;
  }

  if (!varian_read(argv[1], &x)) {
    printf("error: failed to read data from '%s'\n", argv[1]);
    return 1;
  }

  int *sz = hx_array_index_alloc(x.k);
  memcpy(sz, x.sz, x.k * sizeof(int));
  sz[0] = 1;
  while (sz[0] < x.sz[0]) sz[0] *= 2;
  hx_array_resize(&x, x.d, x.k, sz);

  if (!hx_array_fft(&x, 0, 0)) {
    printf("error: failed to fourier transform data\n");
    return 1;
  }

  if (!hx_array_print(&x, NULL)) {
    printf("error: failed to print array data\n");
    return 1;
  }

  hx_array_free(&x);
  return 0;
}

