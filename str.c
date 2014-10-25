
/* ndmath: A framework for n-dimensional hypercomplex calculations for NMR
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

/* include the strings header. */
#include "str.h"

/* strsplit(): splits a string @s1 into a string array by @s2 tokens.
 * @s1: the string to split.
 * @s2: the token string.
 * @ntok: pointer to the output token count.
 */
char **strsplit (const char *s1, const char *s2, unsigned int *ntok) {
  /* declare a few required variables:
   * @count: number of tokens parsed.
   * @n: number of chars to buffer.
   * @pa: token start pointer.
   * @pb: token end pointer.
   * @strv: output string array.
   * @done: completion flag.
   */
  unsigned int count, n;
  char *pa, *pb;
  char **strv;
  char *buf;
  int done;

  /* allocate memory for the temporary token buffer. */
  buf = (char*) malloc((strlen(s1) + 1) * sizeof(char));

  /* check that the buffer was allocated. */
  if (!buf)
    return NULL;

  /* initialize the start pointer. */
  pa = (char*) s1;

  /* initialize the string array. */
  strv = NULL;
  count = 0;

  /* initialize the completion flag. */
  done = 0;

  /* loop until completion is reached. */
  do {
    /* find the next token end pointer. */
    pb = strstr(pa, s2);

    /* check if a new pointer was found. */
    if (!pb) {
      /* no. set up for a copy of the remaining string contents. */
      done = 1;
      pb = pa + strlen(pa);
    }

    /* copy the current token into the buffer. */
    n = pb - pa;
    strncpy(buf, pa, n);
    buf[n] = '\0';

    /* reallocate memory for the next token. */
    count++;
    strv = (char**) realloc(strv, count * sizeof(char*));

    /* check that the allocation was successful. */
    if (!strv)
      return NULL;

    /* allocate memory for the token string. */
    strv[count - 1] = (char*) malloc((n + 1) * sizeof(char));

    /* check that the allocation was successful. */
    if (!strv[count - 1])
      return NULL;

    /* copy the token and move past the delimiter. */
    strcpy(strv[count - 1], buf);
    pa = pb + strlen(s2);
  }
  while (!done);

  /* free the temporary token buffer. */
  free(buf);

  /* set the number of read tokens and return the string array. */
  *ntok = count;
  return strv;
}

/* strvfree(): free a string array.
 * @strv: the string array to free.
 * @n: the number of strings in the array.
 */
void strvfree (char **strv, unsigned int n) {
  /* declare a loop counter. */
  unsigned int i;

  /* free the individual strings in the array. */
  for (i = 0; i < n; i++)
    free(strv[i]);

  /* free the array memory. */
  free(strv);
}

