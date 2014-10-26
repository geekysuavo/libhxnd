
#include <hxnd/nmr-bruker.h>

int main (int argc, char **argv) {
  unsigned int end, nblk, szblk;
  hx_array x;

  if (argc < 2) {
    printf("error: need input file\n");
    return 1;
  }

  end = BYTES_ENDIAN_LITTLE;
  if (argc >= 3 && strcmp(argv[2], "big") == 0)
    end = BYTES_ENDIAN_BIG;

  if (argc >= 5) {
    nblk = atol(argv[3]);
    szblk = atol(argv[4]);
  }
  else {
    nblk = 1;
    szblk = bytes_size(argv[1]);
  }

  if (!bruker_read(argv[1], end, nblk, szblk, &x)) {
    printf("error: failed to read data from '%s'\n", argv[1]);
    return 1;
  }

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

