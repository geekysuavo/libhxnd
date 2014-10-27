
#include <hxnd/nmr-bruker.h>

int main (int argc, char **argv) {
  datum D;

  if (argc < 2) {
    printf("error: need input directory\n");
    return 1;
  }

  if (!bruker_datum(argv[1], &D)) {
    printf("error: failed to read data from '%s'\n", argv[1]);
    return 1;
  }

  if (!datum_print(&D, NULL)) {
    printf("error: failed to print datum header\n");
    return 1;
  }

  return 0;
}

