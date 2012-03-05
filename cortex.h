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

enum CORTEX_FILE_TYPE {UNKNOWN_FILE,BUBBLE_FILE,ALIGNMENT_FILE};
enum HETEROGENEITY {UNKNOWN_HET,HOM1,HOM2,HET};

typedef enum HETEROGENEITY HETEROGENEITY;
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
  unsigned char has_likelihoods, kmer_size,
                fails_classifier_line, discovery_phase_line, is_diploid;

  unsigned long num_of_colours;
  unsigned long *colour_arr;
};

struct COLOUR_COVG
{
  unsigned long length, capacity;
  unsigned long *colour_covgs;
};

struct CORTEX_ALIGNMENT
{
  STRING_BUFFER *name, *seq;
  COLOUR_COVG **colour_covgs;
};

struct CORTEX_BUBBLE_PATH
{
  STRING_BUFFER *seq;
  size_t seq_length;
  float mean_covg;
  unsigned long min_covg, max_covg, fst_covg, lst_covg;
  char fst_kmer[65], fst_r[65], fst_f[65], lst_kmer[65], lst_r[65], lst_f[65];
};

struct CORTEX_BUBBLE
{
  unsigned long var_num;
  CORTEX_BUBBLE_PATH flank_5p, flank_3p;
  CORTEX_BUBBLE_PATH branches[2];

  // Arrays of colours
  HETEROGENEITY *calls;
  float *llk_hom_br1, *llk_het, *llk_hom_br2;

  // Array of pointers to arrays of covgs
  // colour_covgs[colour][kmer_index] = coverage
  COLOUR_COVG **branches_colour_covgs[2];
};

// path can point to a .colour_covgs or .colour_covgs.gzip file
CORTEX_FILE* cortex_open(const char *path);
// cortex_close frees CORTEX_FILE
void cortex_close(CORTEX_FILE *cortex);

char* cortex_colour_list_str(const CORTEX_FILE* c_file);

long cortex_file_get_colour_index(unsigned long colour,
                                  const CORTEX_FILE* c_file);

//
// Reading bubbles
//

// Before reading any bubbles allocate memory for the result with:
CORTEX_BUBBLE* cortex_bubble_create(const CORTEX_FILE *c_file);
// Reset values
void cortex_bubble_reset(CORTEX_BUBBLE* bubble, const CORTEX_FILE *c_file);
// Once you're done reading, free memory
void cortex_bubble_free(CORTEX_BUBBLE* bubble, const CORTEX_FILE *c_file);

// Read a bubble from a file
char cortex_read_bubble(CORTEX_BUBBLE* bubble, CORTEX_FILE* file);
// Print a bubble that came from a given file
void cortex_print_bubble(const CORTEX_BUBBLE* bubble, const CORTEX_FILE *c_file);

//
// Reading alignments
//

// Before reading any alignments allocate memory for the result with:
CORTEX_ALIGNMENT* cortex_alignment_create(const CORTEX_FILE *c_file);
// Reset values
void cortex_alignment_reset(CORTEX_ALIGNMENT* alignment,
                            const CORTEX_FILE *c_file);
// Once you're done reading, free memory
void cortex_alignment_free(CORTEX_ALIGNMENT* alignment,
                           const CORTEX_FILE *c_file);

// Read an alignment to a (possibly multicoloured) graph from a file
char cortex_read_alignment(CORTEX_ALIGNMENT* alignment, CORTEX_FILE* file);
// Print an alignment that came from a given file
void cortex_print_alignment(const CORTEX_ALIGNMENT* alignment,
                            const CORTEX_FILE* file);

#endif
