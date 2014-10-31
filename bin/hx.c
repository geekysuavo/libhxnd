
/* include the datum header. */
#include <hxnd/nmr-datum.h>

/* include the portable option parsing header. */
#include <hxnd/opts.h>

/* correction: structure that defines parameter corrections that need to
 * be made to datum metadata fields.
 */
struct correction {
  /* @d: dimension index.
   * @key: parameter key string.
   * @val: parameter value string.
   */
  unsigned int d;
  char key[32];
  char val[32];
};

/* main(): application entry point.
 * @argc: argument count.
 * @argv: argument array.
 */
int main (int argc, char **argv) {
  /* declare a datum structure for loading data into. */
  datum D;

  /* declare variables used for argument parsing:
   * @argi: current argument array index in opts_get()
   * @c: returned option character from opts_get()
   */
  int argi = 0;
  int c;

  /* declare variables for file input/output:
   * @fname_out: output filename.
   * @fname_in: input filename.
   * @fh: output file handle.
   */
  char *fname_out = NULL;
  char *fname_in = NULL;
  FILE *fh;

  /* declare variables for behavior determination. */
  int pretend = 0;

  /* declare variables for parameter correction:
   * @kbuf, vbuf, ibuf: temporary parsing locations.
   * @corrs: array of correction values.
   * @n_corrs: array length.
   */
  struct correction *corrs = NULL;
  unsigned int n_corrs = 0;
  char kbuf[32], vbuf[32];
  unsigned int i, ibuf;

  /* loop until the arguments are exhausted. */
  while ((c = opts_get(argc, argv, "i:o:pf:v:", &argi)) != -1) {
    /* determine which option was specified. */
    switch ((char) c) {
      /* i: input filename. */
      case 'i':
        fname_in = argv[argi - 1];
        break;

      /* o: output filename. */
      case 'o':
        fname_out = argv[argi - 1];
        break;

      /* p: pretend mode. */
      case 'p':
        pretend = 1;
        break;

      /* f: processing function. */
      case 'f':
        /* FIXME: add processing functions. */
        break;

      /* v: dimension parameter correction. */
      case 'v':
        /* initialize the parsing locations. */
        strncpy(kbuf, "", 32);
        strncpy(vbuf, "", 32);
        ibuf = 0;

        /* attempt to parse the parameter values. */
        if (sscanf(argv[argi - 1], " %31[^[] [ %u ] = %31s ",
                   kbuf, &ibuf, vbuf) != 3) {
          /* raise an error and end execution. */
          raise("failed to parse correction '%s'", argv[argi - 1]);
          traceback_print();
          return 1;
        }

        /* ensure string termination. */
        kbuf[31] = vbuf[31] = '\0';

        /* reallocate the array of corrections. */
        corrs = (struct correction*)
          realloc(corrs, ++n_corrs * sizeof(struct correction));

        /* check that the reallocation succeeded. */
        if (!corrs) {
          /* raise an error and end execution. */
          raise("failed to reallocate correction array");
          traceback_print();
          return 1;
        }

        /* store the correction information. */
        corrs[n_corrs - 1].d = ibuf;
        strcpy(corrs[n_corrs - 1].key, kbuf);
        strcpy(corrs[n_corrs - 1].val, vbuf);
        break;

      /* unknown option or argument. */
      default:
        /* raise an error and quit parsing. */
        raise("failed to parse arguments");
        traceback_print();
        return 1;
    }
  }

  /* initialize the datum structure. */
  datum_init(&D);

  /* determine how to read data in. */
  if (fname_in) {
    /* determine the type of input data. */
    D.type = datum_guess_type(fname_in);

    /* ensure that the file format is supported. */
    if (D.type == DATUM_TYPE_UNDEFINED) {
      /* raise an error and end execution. */
      raise("unsupported data type in '%s'", fname_in);
      traceback_print();
      return 1;
    }

    /* fill the datum structure */
    if (!datum_fill(&D, fname_in)) {
      /* raise an error and end execution. */
      raise("failed to generate hx-format metadata from '%s'", fname_in);
      traceback_print();
      return 1;
    }

    /* load the array data, if required. */
    if (!pretend && !datum_read_array(&D)) {
      /* raise an error and end execution. */
      raise("failed to generate hx-format data from '%s'", fname_in);
      traceback_print();
      return 1;
    }
  }
  else {
    /* standard input mode only supports hx-formatted binary files. */
    D.type = DATUM_TYPE_HXND;

    /* fill the datum structure from standard input. */
    if (!datum_fread(&D, stdin, !pretend)) {
      /* raise an error and end execution. */
      raise("failed to read hx-format data from stdin");
      traceback_print();
      return 1;
    }
  }

  /* apply parameter corrections at this point. */
  for (i = 0; i < n_corrs; i++) {
    /* apply the currently indexed correction. */
    if (!datum_set_dim_parameter(&D, corrs[i].key, corrs[i].d - 1,
                                 corrs[i].val)) {
      /* raise an exception. */
      raise("failed to correct %s[%u] to '%s'",
            corrs[i].key, corrs[i].d, corrs[i].val);

      /* end execution. */
      traceback_print();
      return 1;
    }
  }

  /* check if we're only pretending to read the data. */
  if (pretend) {
    /* make sure the array has been dropped from the structure. */
    datum_free_array(&D);

    /* print the datum metadata to the terminal. */
    if (!datum_print(&D, fname_out)) {
      /* raise an error. */
      raise("failed to print hx-format metadata to %s%s%s",
            fname_out ? "'" : "",
            fname_out ? fname_out : "stdout",
            fname_out ? "'" : "");

      /* end execution. */
      traceback_print();
      return 1;
    }

    /* return successfully. */
    return 0;
  }

  /* FIXME: apply processing functions at this point. */

  /* determine how to write data out. */
  if (fname_out)
    fh = fopen(fname_out, "wb");
  else
    fh = stdout;

  /* write the data out. */
  if (!datum_fwrite(&D, fh)) {
    /* raise an error. */
    raise("failed to write hx-format data to %s%s%s",
          fname_out ? "'" : "",
          fname_out ? fname_out : "stdout",
          fname_out ? "'" : "");

    /* end execution. */
    traceback_print();
    return 1;
  }

  /* return successfully. */
  return 0;
}

