
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

/* include the strings header. */
#include <hxnd/str.h>

/* strtrim(): trim spaces from the beginning and end of a string.
 * @s: the strig to trim.
 */
void strtrim (char *s) {
  /* declare required variables. */
  int i, j, n;
  char *stmp;

  /* check that the string is non-null. */
  if (!s)
    return;

  /* get the string length. */
  n = strlen(s);

  /* loop from the beginning until non-whitespace is found. */
  for (i = 0; i < n && s[i] == ' '; i++) {}

  /* loop from the end until non-whitespace is found. */
  for (j = n - 1; j >= 0 && s[j] == ' '; j--) {}

  /* if the whole string is whitespace, null it out. */
  if (i > j)
    strcpy(s, "");

  /* allocate a temporary string. */
  stmp = (char*) malloc((j - i + 2) * sizeof(char));
  if (!stmp)
    return;

  /* build the trimmed substring. */
  strncpy(stmp, s + i, j - i + 1);
  stmp[j - i + 1] = '\0';

  /* copy the trimmed substring and free the temporary string. */
  strcpy(s, stmp);
  free(stmp);
}

/* strltrim(): trim spaces from the beginning of a string.
 * @s: the string to trim.
 */
void strltrim (char *s) {
  /* declare required variables. */
  int i, n;
  char *stmp;

  /* check that the string is non-null. */
  if (!s)
    return;

  /* get the string length. */
  n = strlen(s);

  /* loop from the beginning until non-whitespace is found. */
  for (i = 0; i < n && s[i] == ' '; i++) {}

  /* if the whole string is whitespace, null it out. */
  if (i >= n)
    strcpy(s, "");

  /* allocate a temporary string. */
  stmp = (char*) malloc((n - i + 1) * sizeof(char));
  if (!stmp)
    return;

  /* build the trimmed substring. */
  strncpy(stmp, s + i, n - i + 1);
  stmp[n - i + 1] = '\0';

  /* copy the trimmed substring and free the temporary string. */
  strcpy(s, stmp);
  free(stmp);
}

/* strrtrim(): trim spaces from the end of a string.
 * @s: the string to trim.
 */
void strrtrim (char *s) {
  /* declare a variable to hold the current string length. */
  int n;

  /* check that the string is non-null. */
  if (!s)
    return;

  /* loop until the end of the string contains no undesireables.
   */
  while ((n = strlen(s)) && s[n - 1] == ' ')
    s[n - 1] = '\0';
}

/* strnltrim(): trim newlines and carriage returns from the end of a string.
 * @s: the string to trim.
 */
void strnltrim (char *s) {
  /* declare a variable to hold the current string length. */
  int n;

  /* check that the string is non-null. */
  if (!s)
    return;

  /* loop until the end of the string contains no undesireables.
   */
  while ((n = strlen(s)) &&
         (s[n - 1] == '\n' ||
          s[n - 1] == '\r'))
    s[n - 1] = '\0';
}

/* strtolower(): convert a string to all lowercase characters.
 * @s: the string to case-fold.
 */
void strtolower (char *s) {
  /* declare a few required variables. */
  int i;

  /* check that the string is non-null. */
  if (!s)
    return;

  /* loop over the characters of the string. */
  for (i = 0; i < strlen(s); i++)
    s[i] = tolower(s[i]);
}

/* strtoupper(): convert a string to all uppercase characters.
 * @s: the string to case-fold.
 */
void strtoupper (char *s) {
  /* declare a few required variables. */
  int i;

  /* check that the string is non-null. */
  if (!s)
    return;

  /* loop over the characters of the string. */
  for (i = 0; i < strlen(s); i++)
    s[i] = toupper(s[i]);
}

/* strbool(): read any rational string name for a boolean value and return
 * an integer representation.
 * @s: the string to parse.
 */
int strbool (const char *s) {
  /* check that the string is allocated. */
  if (!s || strlen(s) < 1)
    return 0;

  /* check the acceptable values for truth. */
  if (strcmp(s, "true") == 0 ||
      strcmp(s, "True") == 0 ||
      strcmp(s, "TRUE") == 0 ||
      strcmp(s, "yes") == 0 ||
      strcmp(s, "Yes") == 0 ||
      strcmp(s, "YES") == 0 ||
      strcmp(s, "t") == 0 ||
      strcmp(s, "T") == 0 ||
      strcmp(s, "y") == 0 ||
      strcmp(s, "Y") == 0 ||
      strcmp(s, "1") == 0) {
    /* an acceptable value was identified. return 'true'. */
    return 1;
  }

  /* return 'false' if no 'true' value was identified. */
  return 0;
}

/* strsplit(): split a string @s1 into a string array by @s2 tokens.
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

/* strvtrim(): trim leading and trailing spaces from each string in an array.
 * @strv: the string array to trim.
 * @n: the number of strings in the array.
 */
void strvtrim (char **strv, unsigned int n) {
  /* declare a loop counter. */
  unsigned int i;

  /* trim the individual strings in the array. */
  for (i = 0; i < n; i++)
    strtrim(strv[i]);
}

/* strvcompact(): compact a string array by removing all empty strings.
 * @strv: the string array to compact.
 * @pn: pointer to the array length.
 */
void strvcompact (char **strv, unsigned int *pn) {
  /* declare a few required variables. */
  unsigned int i, n;

  /* loop over the array until we reach a non-empty member. */
  for (i = 0, n = 0; i < *pn; i++) {
    /* check if the current member is non-empty. */
    if (strlen(strv[i])) {
      /* check if it may be shifted left. */
      if (i > n) {
        /* yes. resize the destination string. */
        strv[n] = (char*)
          realloc(strv[n], (strlen(strv[i]) + 1) * sizeof(char));

        /* move the string. */
        strcpy(strv[n], strv[i]);
        strcpy(strv[i], "");
      }

      /* increment the adjusted index. */
      n++;
    }
  }

  /* check if fields were moved. */
  if (n < *pn) {
    /* free the trailing empty strings. */
    for (i = n; i < *pn; i++)
      free(strv[i]);

    /* reallocate the string array. */
    strv = (char**) realloc(strv, n * sizeof(char*));
    *pn = n;
  }
}

/* strvtolower(): convert each string in an array to lowercase.
 * @strv: the string array to convert.
 * @n: the number of strings in the array.
 */
void strvtolower (char **strv, unsigned int n) {
  /* declare a loop counter. */
  unsigned int i;

  /* convert the individual strings in the array. */
  for (i = 0; i < n; i++)
    strtolower(strv[i]);
}

/* strvtoupper(): convert each string in an array to uppercase.
 * @strv: the string array to convert.
 * @n: the number of strings in the array.
 */
void strvtoupper (char **strv, unsigned int n) {
  /* declare a loop counter. */
  unsigned int i;

  /* convert the individual strings in the array. */
  for (i = 0; i < n; i++)
    strtoupper(strv[i]);
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

