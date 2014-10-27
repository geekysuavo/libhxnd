
#include <hxnd/nmr-varian.h>

int main (int argc, char **argv) {
  datum D;

  if (argc < 2) {
    printf("usage: %s [data-directory]\n", argv[0]);
    return 1;
  }

  datum_init(&D);
  if (!varian_fill_datum(argv[1], &D)) {
    traceback_print();
    return 1;
  }

  if (!datum_read_array(&D)) {
    traceback_print();
    return 1;
  }

  if (!datum_print(&D, NULL)) {
    traceback_print();
    return 1;
  }

  return 0;
}

