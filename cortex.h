/*
 cortex.h
 project: Cortex Library
 author: Isaac Turner <turner.isaac@gmail.com>

 Copyright (c) 2012, Isaac Turner
 All rights reserved.

 see: README

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CORTEX_H_SEEN
#define CORTEX_H_SEEN

#include "string_buffer.h"

enum CORTEX_FILE_TYPE {UNKNOWN,BUBBLE_FILE,ALIGNMENT_FILE};

typedef struct CORTEX_FILE CORTEX_FILE;
typedef struct CORTEX_ALIGNMENT CORTEX_ALIGNMENT;
typedef struct CORTEX_BUBBLE CORTEX_BUBBLE;
typedef struct CORTEX_BUBBLE_PATH CORTEX_BUBBLE_PATH;
typedef struct COLOUR_COVG COLOUR_COVG;

struct CORTEX_FILE
{
  // For reading the file
  char *path;
  gzFile *file;
  STRING_BUFFER *buffer;
  unsigned long line_number; // line currently in buffer (starting at 1)

  // Syntax of the file
  enum CORTEX_FILE_TYPE filetype;
  char has_likelihoods;
  unsigned long num_of_colours;
};

struct COLOUR_COVG
{
  unsigned long length, capacity;
  unsigned long *colour_covgs;
};

struct CORTEX_ALIGNMENT
{
  STRING_BUFFER *name, *seq;
  COLOUR_COVG *colour_covgs;
};

struct CORTEX_BUBBLE_PATH
{
  char *seq;
  size_t seq_length;
  float mean_covg;
  unsigned long min_covg, max_covg, fst_covg, lst_covg;
  char *fst_kmer, *fst_r, *fst_f, *lst_kmer, *lst_r, *lst_f;
};

struct CORTEX_BUBBLE
{
  unsigned long var_num;
  CORTEX_BUBBLE_PATH flank_5p, flank_3p, branches[2];

  // Array of pointers to arrays of covgs
  // colour_covgs[colour][kmer_index] = coverage
  COLOUR_COVG **branches_colour_covgs[2];
};

CORTEX_FILE* cortex_open(const char *path);
void cortex_close(CORTEX_FILE *cortex);

CORTEX_BUBBLE* cortex_bubble_create(const CORTEX_FILE *c_file);
void cortex_bubble_free(CORTEX_BUBBLE* bubble, const CORTEX_FILE *c_file);

char cortex_read_alignment(CORTEX_FILE* file, CORTEX_ALIGNMENT* alignment);
void cortex_print_alignment(const CORTEX_ALIGNMENT* alignment);

char cortex_read_bubble(CORTEX_FILE* file, CORTEX_BUBBLE* bubble);
void cortex_print_bubble(const CORTEX_BUBBLE* bubble);

#endif
