
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

/* include the nmr data header. */
#include <hxnd/nmr-datum.h>

/* datum_type_def: datum type definition structure for holding all available
 * datum type names, values and encode/decode/guess functions.
 */
struct datum_type_def {
  /* identifying fields:
   * @type: the value of the datum type.
   * @name: the string name of the datum type.
   * @desc: short description of the datum type.
   */
  enum datum_type type;
  const char *name;
  const char *desc;

  /* function pointers for format determination and conversion:
   * @guessfn: pointer to the function for format checking.
   * @encfn: pointer to the function for formatted writing.
   * @decfn: pointer to the function for formatted parameter reading.
   * @arrayfn: pointer to the function for formatted array reading.
   * @postfn: pointer to the function for post-array read fixes.
   */
  int (*guessfn) (const char *fname);
  int (*encfn) (datum *D, const char *fname);
  int (*decfn) (datum *D, const char *fname);
  int (*arrayfn) (datum *D);
  int (*postfn) (datum *D);
};

/* datum_types: local table of all available datum type names, values
 * and function handles.
 */
static const struct datum_type_def datum_types[] = {
  { DATUM_TYPE_HXND, "hx",
    "Native hypercomplex",
    &hxnd_guess,
    &hxnd_encode,
    &hxnd_decode,
    &hxnd_array,
    NULL },
  { DATUM_TYPE_TEXT, "text",
    "Plain text",
    NULL,
    &text_encode,
    NULL,
    NULL,
    NULL },
  { DATUM_TYPE_BRUKER, "bruker",
    "Bruker unprocessed",
    &bruker_guess,
    NULL,
    &bruker_decode,
    &bruker_array,
    &bruker_post },
  { DATUM_TYPE_VARIAN, "varian",
    "Varian/Agilent unprocessed",
    &varian_guess,
    NULL,
    &varian_decode,
    &varian_array,
    NULL },
  { DATUM_TYPE_PIPE, "pipe",
    "NMRPipe",
    &pipe_guess,
    &pipe_encode,
    &pipe_decode,
    &pipe_array,
    NULL },
  { DATUM_TYPE_UCSF, "ucsf",
    "UCSF/Sparky",
    &ucsf_guess,
    &ucsf_encode,
    &ucsf_decode,
    &ucsf_array,
    NULL },
  { DATUM_TYPE_NV, "nv",
    "NMRView/NMRViewJ",
    &nv_guess,
    &nv_encode,
    &nv_decode,
    &nv_array,
    &nv_post },
  { DATUM_TYPE_RNMRTK, "rnmrtk",
    "Rowland NMR Toolkit",
    &rnmrtk_guess,
    &rnmrtk_encode,
    &rnmrtk_decode,
    &rnmrtk_array,
    NULL },
  { DATUM_TYPE_UNDEFINED, NULL, NULL, NULL, NULL, NULL, NULL, NULL }
};

/* datum_type_name(): return the string short name of an enumerated datum
 * type. do not attempt to free the returned string.
 * @typ: the datum enumerated type.
 */
const char *datum_type_name (enum datum_type typ) {
  /* declare a required variable. */
  unsigned int i;

  /* loop over all supported datum types. */
  for (i = 0; datum_types[i].name; i++) {
    /* break if the datum type matches. */
    if (typ == datum_types[i].type)
      break;
  }

  /* return the identified type name. */
  return (datum_types[i].name ? datum_types[i].name : "unknown");
}

/* datum_type_desc(): return the string long name of an enumerated datum
 * type. do not attempt to free the returned string.
 * @typ: the datum enumerated type.
 */
const char *datum_type_desc (enum datum_type typ) {
  /* declare a required variable. */
  unsigned int i;

  /* loop over all supported datum types. */
  for (i = 0; datum_types[i].name; i++) {
    /* break if the datum type matches. */
    if (typ == datum_types[i].type)
      break;
  }

  /* return the identified type description. */
  return (datum_types[i].desc ? datum_types[i].desc : "Unknown");
}

/* datum_type_lookup(): return an enumerated datum type based on its string
 * representation.
 * @name: the datum type string.
 */
enum datum_type datum_type_lookup (const char *name) {
  /* declare a required variable. */
  unsigned int i;

  /* loop over all supported datum types. */
  for (i = 0; datum_types[i].name; i++) {
    /* break if the datum name matches. */
    if (strcmp(name, datum_types[i].name) == 0)
      break;
  }

  /* return the identified type value. */
  return datum_types[i].type;
}

/* datum_type_guess(): determine the type of datum format contained in a
 * specified file or directory.
 * @fname: the file or directory name to check.
 */
enum datum_type datum_type_guess (const char *fname) {
  /* declare a required variable. */
  unsigned int i;

  /* loop over all supported datum types. */
  for (i = 0; datum_types[i].name; i++) {
    /* skip formats with no available guess function. */
    if (!datum_types[i].guessfn)
      continue;

    /* break if the format guess function returns true. */
    if (datum_types[i].guessfn(fname))
      break;
  }

  /* clear any errors that may have come up during type guessing. */
  traceback_clear();

  /* return the identified type value. */
  return datum_types[i].type;
}

/* datum_type_encode(): write an NMR datum structure to a file in a specified
 * format.
 * @D: pointer to the source datum structure.
 * @fname: output filename string.
 * @typ: the output format type.
 */
int datum_type_encode (datum *D, const char *fname, enum datum_type typ) {
  /* declare a required variable. */
  unsigned int i;

  /* loop over all supported datum types. */
  for (i = 0; datum_types[i].name; i++) {
    /* skip formats with no available encode function. */
    if (!datum_types[i].encfn)
      continue;

    /* run the encode function if the type matches. */
    if (typ == datum_types[i].type)
      return datum_types[i].encfn(D, fname);
  }

  /* return failure. */
  throw("datum format '%s' does not support file encoding",
        datum_type_name(typ));
}

/* datum_type_decode(): populate an NMR datum structure from a file that
 * follows a specified format.
 * @D: pointer to the destination datum structure.
 * @fname: input filename string.
 * @typ: the input format type.
 */
int datum_type_decode (datum *D, const char *fname, enum datum_type typ) {
  /* declare a required variable. */
  unsigned int i;

  /* loop over all supported datum types. */
  for (i = 0; datum_types[i].name; i++) {
    /* skip formats with no available decode function. */
    if (!datum_types[i].decfn)
      continue;

    /* run the decode function if the type matches. */
    if (typ == datum_types[i].type)
      return datum_types[i].decfn(D, fname);
  }

  /* return failure. */
  throw("datum format '%s' does not support file decoding",
        datum_type_name(typ));
}

/* datum_type_array(): read array data into an NMR datum structure from a file
 * that follows a specified format.
 * @D: pointer to the destination datum structure.
 * @typ: the input format type.
 */
int datum_type_array (datum *D, enum datum_type typ) {
  /* declare a required variable. */
  unsigned int i;

  /* loop over all supported datum types. */
  for (i = 0; datum_types[i].name; i++) {
    /* skip formats with no available array function. */
    if (!datum_types[i].arrayfn)
      continue;

    /* run the array function if the type matches. */
    if (typ == datum_types[i].type)
      return datum_types[i].arrayfn(D);
  }

  /* return failure. */
  throw("datum format '%s' does not support array reading",
        datum_type_name(typ));
}

/* datum_type_post(): post-process array data in an NMR datum structure that
 * was loaded from a file that follows a specified format.
 * @D: pointer to the destination datum structure.
 * @typ: the input format type.
 */
int datum_type_post (datum *D, enum datum_type typ) {
  /* declare a required variable. */
  unsigned int i;

  /* loop over all supported datum types. */
  for (i = 0; datum_types[i].name; i++) {
    /* skip formats with no available post-process function. */
    if (!datum_types[i].postfn)
      continue;

    /* run the array function if the type matches. */
    if (typ == datum_types[i].type)
      return datum_types[i].postfn(D);
  }

  /* return success. not all formats require post-processing. */
  return 1;
}

