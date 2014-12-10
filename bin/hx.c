
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
   -n, --new              Create a new datum, do not read input\n\
   -i, --input FIN        Specify an input filename [stdin]\n\
   -o, --output FOUT      Specify an output filename [stdout]\n\
   -F, --format FMT       Specify an output format [hx]\n\
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

/* parsed_arg: structure that defines...
 *  (1) parameter corrections that need to be made to datum metadata fields.
 *  (2) processing functions that need to be applied to datum content.
 */
struct parsed_arg {
  /* @d: dimension index.
   * @lstr: argument left-value string.
   * @rstr: argument right-valuestring.
   */
  int d;
  char *lstr;
  char *rstr;
};

/* function declarations: */

int main_apply_procs (datum *D, struct parsed_arg *lst, unsigned int n);

int main_apply_corrs (datum *D, struct parsed_arg *lst, unsigned int n);

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
    { "new",      1, 'n' },
    { "input",    1, 'i' },
    { "output",   1, 'o' },
    { "format",   1, 'F' },
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
   * @fmt_out: output file format.
   * @fname_out: output filename.
   * @fname_in: input filename.
   */
  enum datum_type fmt_out = DATUM_TYPE_HXND;
  char *fname_out = NULL;
  char *fname_in = NULL;

  /* declare variables for behavior determination:
   * @ndnew: number of dimensions to create if @mknew is raised.
   * @pretend: whether to continue writing to @fh or print header info.
   * @mknew: whether to create a new datum or read from the input file.
   */
  unsigned int ndnew = 0;
  int pretend = 0;
  int mknew = 0;

  /* declare variables for function execution:
   * @procs: array of processing functions.
   * @n_procs: array length.
   */
  struct parsed_arg *procs = NULL;
  unsigned int n_procs = 0;

  /* declare variables for parameter correction:
   * @corrs: array of correction values.
   * @n_corrs: array length.
   */
  struct parsed_arg *corrs = NULL;
  unsigned int n_corrs = 0;

  /* declare variables for general option argument string parsing:
   * @lval: left-hand value string.
   * @rval: right-hand value string.
   * @dval: dimension index.
   */
  char *lval, *rval;
  int dval;

  /* loop until the arguments are exhausted. */
  while ((c = opts_get(argc, argv, long_options, &argi)) != -1) {
    /* determine which option was specified. */
    switch ((char) c) {
      /* h: help mode. */
      case 'h':
        fprintf(stdout, HX_HELPSTRING);
        return 0;

      /* n: new mode. */
      case 'n':
        ndnew = atoi(argv[argi - 1]);
        mknew = 1;
        break;

      /* i: input filename. */
      case 'i':
        fname_in = argv[argi - 1];
        break;

      /* o: output filename. */
      case 'o':
        fname_out = argv[argi - 1];
        break;

      /* F: output format. */
      case 'F':
        /* determine which output type was specified. */
        fmt_out = datum_type_lookup(argv[argi - 1]);

        /* check that the datum type is supported. */
        if (fmt_out == DATUM_TYPE_UNDEFINED)
          trace("unsupported output format '%s'", argv[argi - 1]);

        /* break the switch. */
        break;

      /* p: pretend mode. */
      case 'p':
        pretend = 1;
        break;

      /* f: processing function. */
      case 'f':
        /* parse the processing function argument string. */
        if (!opts_parse_arg(argv[argi - 1], ":", &lval, &rval, &dval))
          trace("failed to parse function argument '%s'", argv[argi - 1]);

        /* reallocate the array of functions. */
        procs = (struct parsed_arg*)
          realloc(procs, ++n_procs * sizeof(struct parsed_arg));

        /* check that the reallocation succeeded. */
        if (!procs)
          trace("failed to reallocate processing function array");

        /* store the processing information. */
        procs[n_procs - 1].lstr = lval;
        procs[n_procs - 1].rstr = rval;
        procs[n_procs - 1].d = dval;
        break;

      /* v: dimension parameter correction. */
      case 'v':
        /* parse the modifier argument string. */
        if (!opts_parse_arg(argv[argi - 1], "=", &lval, &rval, &dval))
          trace("failed to parse value argument '%s'", argv[argi - 1]);

        /* check for the required right-value string. */
        if (strcmp(rval, "") == 0)
          trace("value argument '%s' lacks required right-hand value", lval);

        /* reallocate the array of corrections. */
        corrs = (struct parsed_arg*)
          realloc(corrs, ++n_corrs * sizeof(struct parsed_arg));

        /* check that the reallocation succeeded. */
        if (!corrs)
          trace("failed to reallocate correction array");

        /* store the correction information. */
        corrs[n_corrs - 1].lstr = lval;
        corrs[n_corrs - 1].rstr = rval;
        corrs[n_corrs - 1].d = dval;
        break;

      /* unknown option or argument. */
      default:
        trace("failed to parse arguments");
    }
  }

  /* initialize the datum structure. */
  datum_init(&D);

  /* determine whether and how to read data in. */
  if (mknew) {
    /* initialize the datum type field. */
    D.type = DATUM_TYPE_HXND;

    /* allocate the required datum dimensions. */
    if (!datum_dims_realloc(&D, ndnew))
      throw("failed to allocate dimension array");

    /* apply parameter corrections. */
    if (!main_apply_corrs(&D, corrs, n_corrs))
      trace("failed to apply parameter corrections");

    /* allocate the array to match the datum parameters. */
    if (!datum_array_alloc(&D))
      trace("failed to allocate new datum array");
  }
  else {
    /* determine the type of input data. */
    if (fname_in)
      D.type = datum_type_guess(fname_in);
    else
      D.type = DATUM_TYPE_HXND;

    /* ensure that the file format is supported. */
    if (D.type == DATUM_TYPE_UNDEFINED)
      trace("unsupported data type in '%s'", fname_in);

    /* decode the file parameters into the datum structure. */
    if (!datum_type_decode(&D, fname_in, D.type))
      trace("failed to read %s-format data from '%s'",
            datum_type_name(D.type), fname_in);

    /* apply parameter corrections at this point. */
    if (!main_apply_corrs(&D, corrs, n_corrs))
      trace("failed to apply parameter corrections");

    /* load the array data. */
    if (!datum_array_read(&D))
      trace("failed to generate hx-format data from '%s'", fname_in);
  }

  /* check if we're only pretending to read the data. */
  if (pretend) {
    /* print the datum metadata to the terminal. */
    if (!datum_print(&D, fname_out))
      trace("failed to print hx-format metadata to %s%s%s",
            fname_out ? "'" : "",
            fname_out ? fname_out : "standard output",
            fname_out ? "'" : "");

    /* free the datum structure. */
    datum_free(&D);

    /* return successfully. */
    return 0;
  }

  /* apply processing functions at this point. */
  if (!main_apply_procs(&D, procs, n_procs))
    trace("failed to apply processing functions");

  /* encode the data into the output file. datum formats that support writing
   * to standard output must accept @fname as NULL.
   */
  if (!datum_type_encode(&D, fname_out, fmt_out))
    trace("failed to write '%s'-format data to %s%s%s",
          datum_type_name(fmt_out),
          fname_out ? "'" : "",
          fname_out ? fname_out : "standard output",
          fname_out ? "'" : "");

  /* free the datum structure. */
  datum_free(&D);

  /* return successfully. */
  return 0;
}

/* main_apply_procs(): apply processing functions to a datum structure.
 * @D: pointer to the datum to manipulate.
 * @lst: list of processing function arguments.
 * @n: number of processing functions.
 */
int main_apply_procs (datum *D, struct parsed_arg *lst, unsigned int n) {
  /* declare a required loop counter. */
  unsigned int i;

  /* check that the list is non-empty. */
  if (!n)
    return 1;

  /* loop over the array of function arguments. */
  for (i = 0; i < n; i++) {
    /* apply the currently indexed processing function. */
    if (!fn_execute(D, lst[i].lstr, lst[i].d - 1, lst[i].rstr))
      throw("failed to apply function '%s' (#%u)", lst[i].lstr, i);

    /* free the allocated processing function strings. */
    free(lst[i].lstr);
    free(lst[i].rstr);
  }

  /* free the processing functions array. */
  free(lst);

  /* return success. */
  return 1;
}

/* main_apply_corrs(): apply parameter corrections to a datum structure.
 * @D: pointer to the datum to manipulate.
 * @lst: list of correction arguments.
 * @n: number of corrections.
 */
int main_apply_corrs (datum *D, struct parsed_arg *lst, unsigned int n) {
  /* declare a required loop counter. */
  unsigned int i;

  /* check that the list is non-empty. */
  if (!n)
    return 1;

  /* loop over the array of corrections. */
  for (i = 0; i < n; i++) {
    /* apply the currently indexed correction. */
    if (!datum_dims_setparm(D, lst[i].lstr,
                            lst[i].d - 1,
                            lst[i].rstr)) {
      /* raise an exception. */
      throw("failed to correct %s[%u] to '%s'",
            lst[i].lstr, lst[i].d, lst[i].rstr);
    }

    /* free the allocated correction strings. */
    free(lst[i].lstr);
    free(lst[i].rstr);
  }

  /* free the corrections array. */
  free(lst);

  /* return success. */
  return 1;
}

