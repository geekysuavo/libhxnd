
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

/* include the traceback header. */
#include <hxnd/trace.h>

/* traceback: data structure for holding stack trace information.
 */
struct traceback {
  /* @line: the line number where the error was emitted.
   * @file: the source file that emitted the error.
   * @num: the emitted error code.
   */
  unsigned int line;
  char *file;
  int num;

  /* @msg: the custom error string, or NULL if none. */
  char *msg;
};

/* tb: a private variable in the trace_* namespace that holds the current
 * stack trace of any thrown errors.
 */
struct traceback *tb;
int n_tb = -1;

/* traceback_init(): initializes the traceback array for use.
 */
void traceback_init (void) {
  /* check if the traceback array has been initialized. */
  if (n_tb >= 0)
    return;

  /* initialize the traceback array and count. */
  tb = NULL;
  n_tb = 0;
}

/* traceback_print(): prints the contents of the stack trace to standard out.
 */
void traceback_print (void) {
  /* declare a loop index. */
  int i;

  /* ensure the traceback array is initialized. */
  traceback_init();

  /* loop over the traceback elements. */
  for (i = 0; i < n_tb; i++) {
    /* print the traceback line. */
    fprintf(stderr, "%s:%u: (%s) %s\n",
      tb[i].file, tb[i].line,
      strerror(tb[i].num),
      tb[i].msg);
  }
}

/* traceback_clear(): clears the current stack trace in preparation for the
 * next error, ugh...
 */
void traceback_clear (void) {
  /* declare a loop index. */
  int i;

  /* ensure the traceback array is initialized. */
  traceback_init();

  /* loop over the traceback elements. */
  for (i = 0; i < n_tb; i++) {
    /* free the allocated strings. */
    if (tb[i].file) free(tb[i].file);
    if (tb[i].msg) free(tb[i].msg);
  }

  /* free the traceback array. */
  free(tb);

  /* prepare for the next error. */
  tb = NULL;
  n_tb = 0;
}

/* traceback_throw(): main function that appends another frame to the
 * stack trace array.
 * @f: the emitter filename.
 * @l: the emitter line number.
 * @format: printf-style format string for custom messages.
 * @...: arguments that go with the format string.
 */
int traceback_throw (const char *f, const unsigned int l,
                     const char *format, ...) {
  /* declare a few required variables. */
  unsigned int n_msg;
  va_list vl;
  int i;

  /* ensure the traceback array is initialized. */
  traceback_init();

  /* increment the traceback array size. */
  i = n_tb;
  n_tb++;

  /* reallocate the traceback array. */
  tb = (struct traceback*) realloc(tb, n_tb * sizeof(struct traceback));
  if (!tb)
    return 0;

  /* add the filename string to the traceback. */
  n_msg = 2 * strlen(f);
  tb[i].file = (char*) malloc(n_msg * sizeof(char));
  if (tb[i].file)
    strcpy(tb[i].file, f);

  /* add the information into the traceback. */
  tb[i].num = errno;
  tb[i].line = l;

  /* build the custom message string. */
  if (format) {
    /* allocate memory for the message. */
    n_msg = 2 * strlen(format);
    tb[i].msg = (char*) malloc(n_msg * sizeof(char));
    if (!tb[i].msg)
      return 0;

    /* write the formatted message string. */
    va_start(vl, format);
    vsnprintf(tb[i].msg, n_msg, format, vl);
    va_end(vl);
  }

  /* always return failure. this allows the throw() macro to be used
   * easily in tail-calls.
   */
  return 0;
}

