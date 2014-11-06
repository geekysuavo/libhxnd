
/* hxnd: A framework for n-dimensional hypercomplex calculations for NMR.
 * Copyright (C) 2014  Bradley Worley  <geekysuavo@gmail.com>.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to:
 *
 *   Free Software Foundation, Inc.
 *   51 Franklin Street, Fifth Floor
 *   Boston, MA  02110-1301, USA.
 */

/* include the datum and function headers. */
#include <hxnd/nmr-datum.h>
#include <hxnd/fn.h>

/* include the portable option parsing header. */
#include <hxnd/opts.h>

/* define a short help message string to display if the help flag is thrown.
 */
#define HX_HELPSTRING "\
 hx: A command-line multi-tool for handling NMR data\n\
 Copyright (C) 2014 Bradley Worley. Released under the GNU GPL 2.0.\n\
\n\
 Usage:\n\
   hx [OPTIONS]\n\
\n\
 Options:\n\
   -h, --help             Display this help message\n\
   -i, --input FIN        Specify an input filename [stdin]\n\
   -o, --output FOUT      Specify an output filename [stdout]\n\
   -p, --pretend          Perform no actual processing\n\
   -f, --function FNDEF   Apply a processing function (optional)\n\
   -v, --value VALDEF     Change a parameter value (optional)\n\
\n\
 The hx tool performs all functions required to convert and process NMR\n\
 time-domain and spectral data, based on the libhxnd framework for using\n\
 multidimensional arrays of hypercomplex numbers.\n\
\n\
 For more information on available processing functions and their syntax,\n\
 see the manual page for hx(1).\n\
\n\
"

/* processor: structure that defines functions that need to be applied to
 * datum content.
 */
struct processor {
  /* @d: dimension index.
   * @name: function name string.
   * @args: function arguments string.
   */
  int d;
  char *name;
  char *args;
};

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

  /* declare the option definition array used by opts_get() for parsing.
   */
  const opts_def long_options[] = {
    { "help",     0, 'h' },
    { "input",    1, 'i' },
    { "output",   1, 'o' },
    { "pretend",  0, 'p' },
    { "function", 1, 'f' },
    { "value",    1, 'v' },
    { NULL, 0, '\0' }
  };

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

  /* declare variables for function execution:
   * @fnname: function name temporary string.
   * @fnargs: function args temporary string.
   * @procs: array of processing functions.
   * @n_procs: array length.
   */
  struct processor *procs = NULL;
  unsigned int n, fnc, n_procs = 0;
  char *fnname, *fnargs;
  char **fnv;
  int dbuf;

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
  while ((c = opts_get(argc, argv, long_options, &argi)) != -1) {
    /* determine which option was specified. */
    switch ((char) c) {
      /* h: help mode. */
      case 'h':
        fprintf(stdout, HX_HELPSTRING);
        return 0;

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
        /* allocate the buffer strings for storing function name/args. */
        n = strlen(argv[argi - 1]);
        fnname = (char*) malloc(n * sizeof(char));
        fnargs = (char*) malloc(n * sizeof(char));

        /* initialize the buffer strings. */
        strcpy(fnname, "");
        strcpy(fnargs, "");

        /* initialize the dimension index. this value will result in
         * an 'invalid dimension' error for any function that requires
         * a valid dimension index, because the call to fn_execute()
         * will subtract 1.
         */
        dbuf = 0;

        /* check that allocation succeeded. */
        if (!fnname || !fnargs) {
          /* raise an error and end execution. */
          raise("failed to allocate temporary buffers");
          traceback_print();
          return 1;
        }

        /* try to split the function and arguments. */
        fnv = strsplit(argv[argi - 1], ":", &fnc);

        /* check that the split succeeded. */
        if (!fnv || fnc < 1) {
          /* raise an error and end execution. */
          raise("failed to split function string");
          traceback_print();
          return 1;
        }

        /* attempt to parse the function name and dimension. */
        if (sscanf(fnv[0], " %[^ [] [ %d ] ", fnname, &dbuf) != 2 &&
            sscanf(fnv[0], " %s ", fnname) != 1) {
          /* raise an error and end execution. */
          raise("failed to parse function '%s'", fnv[0]);
          traceback_print();
          return 1;
        }

        /* attempt to parse the function arguments. */
        if (fnc >= 2)
          strcpy(fnargs, fnv[1]);

        /* free the function string array. */
        strvfree(fnv, fnc);

        /* reallocate the array of functions. */
        procs = (struct processor*)
          realloc(procs, ++n_procs * sizeof(struct processor));

        /* check that the reallocation succeeded. */
        if (!procs) {
          /* raise an error and end execution. */
          raise("failed to reallocate processing function array");
          traceback_print();
          return 1;
        }

        /* store the processing information. */
        procs[n_procs - 1].name = fnname;
        procs[n_procs - 1].args = fnargs;
        procs[n_procs - 1].d = dbuf;
        break;

      /* v: dimension parameter correction. */
      case 'v':
        /* attempt to parse the parameter values. */
        if (sscanf(argv[argi - 1], " %31[^ [] [ %u ] = %31s ",
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

    /* free the parameter corrections array. */
    if (n_corrs) {
      /* only if corrections have been allocated. */
      free(corrs);
      n_corrs = 0;
    }

    /* load the array data. */
    if (!datum_read_array(&D)) {
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
    if (!datum_fread(&D, stdin, 1)) {
      /* raise an error and end execution. */
      raise("failed to read hx-format data from stdin");
      traceback_print();
      return 1;
    }
  }

  /* check if we're only pretending to read the data. */
  if (pretend) {
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

  /* apply processing functions at this point. */
  for (i = 0; i < n_procs; i++) {
    /* apply the currently indexed processing function. */
    if (!fn_execute(&D, procs[i].name, procs[i].d - 1, procs[i].args)) {
      /* raise an exception. */
      raise("failed to apply function '%s' (#%u)", procs[i].name, i);

      /* end execution. */
      traceback_print();
      return 1;
    }

    /* free the allocated processing function strings. */
    free(procs[i].name);
    free(procs[i].args);
  }

  /* free the processing functions array. */
  if (n_procs) {
    /* but only if it's been allocated. */
    free(procs);
    n_procs = 0;
  }

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

