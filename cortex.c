/*
 cortex.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "cortex.h"

t_buf_pos _cortex_read_line(CORTEX_FILE* file)
{
  t_buf_pos chars_read = string_buff_reset_gzreadline(file->buffer, file->file);
  string_buff_chomp(file->buffer);
  file->line_number++;
  return chars_read;
}

// If cannot open file or file is empty will print error and return NULL
CORTEX_FILE* cortex_open(const char* path)
{
  CORTEX_FILE* cortex = (CORTEX_FILE*) malloc(sizeof(CORTEX_FILE));

  // Give initial values
  cortex->buffer = string_buff_init(500); // Create read in buffer
  cortex->line_number = 0;
  cortex->filetype = UNKNOWN;
  cortex->has_likelihoods = 0;
  cortex->num_of_colours = 0;

  // Set path
  size_t path_len = strlen(path);
  cortex->path = (char*) malloc(path_len+1);
  strcpy(cortex->path, path);
  cortex->path[path_len] = '\0';

  // Open file
  cortex->file = gzopen(path, "r");

  if(cortex->file == NULL)
  {
    fprintf(stderr, "cortex.c: couldn't open file (%s)\n", path);
    cortex_close(cortex);
    return NULL;
  }

  // Whilst still reading but lines empty (_cortex_read_line does chomp)
  t_buf_pos chars_read;
  
  while((chars_read = _cortex_read_line(cortex)) > 0 &&
        string_buff_strlen(cortex->buffer) == 0) {}

  if(chars_read == 0)
  {
    fprintf(stderr, "cortex.c: file is empty (%s)\n", path);
    cortex_close(cortex);
    return NULL;
  }

  if(strncasecmp(cortex->buffer->buff, "Colour", 6) == 0)
  {
    cortex->filetype = BUBBLE_FILE;
    cortex->has_likelihoods = 1;
  
    // Count lines begining [0-9] to count colours
    while(_cortex_read_line(cortex) > 0)
    {
      // Ignore blank lines
      if(string_buff_strlen(cortex->buffer) > 0)
      {
        char first_char = cortex->buffer->buff[0];

        if(first_char == '>')
        {
          break;
        }
        else if(first_char >= '0' && first_char <= '9')
        {
          cortex->num_of_colours++;
        }
      }
    }

    // Reset file
    gzrewind(cortex->file);
    return cortex;
  }

  if(string_buff_get_char(cortex->buffer, 0) != '>')
  {
    fprintf(stderr, "cortex.c: unrecognised line (%s:%lu)\n",
            path, cortex->line_number);

    cortex_close(cortex);
    return NULL;
  }

  _cortex_read_line(cortex);

  _cortex_read_line(cortex);
  char* third_line = string_buff_as_str(cortex->buffer);

  _cortex_read_line(cortex);
  char* fourth_line = string_buff_as_str(cortex->buffer);

  // else if third_line ends "_colour_0_kmer_coverages" AND
  //         fourth_line starts with a [0-9]
  // then ALIGNMENT FILE
  char *expected_line = "_colour_0_kmer_coverages";
  char *third_line_search = strstr(third_line, expected_line);


  if(third_line_search != NULL &&
     third_line_search - third_line + strlen(expected_line)
       == strlen(third_line) &&
     fourth_line[0] >= '0' && fourth_line[0] <= '9')
  {
    cortex->filetype = ALIGNMENT_FILE;

    // read in pairs on lines
    //  and as long as the second one begins [0-9] then num_of_colours++

    while(_cortex_read_line(cortex) > 0 &&
          _cortex_read_line(cortex) > 0 &&
          cortex->buffer->buff[0] >= '0' && cortex->buffer->buff[0] <= '9')
    {
      cortex->num_of_colours++;
    }

    // Reset file
    gzrewind(cortex->file);
    return cortex;
  }
  
  // else if third_line matches bubble path AND
  //         fourth_line start [ACGT]
  // then BUBBLE_FILE
  // 

  unsigned long var_num, length_5p;
  float mean_covg_5p;
  unsigned long min_covg_5p, max_covg_5p, fst_covg_5p, lst_covg_5p;
  char *fst_kmer_5p, *fst_r_5p, *fst_f_5p,
       *lst_kmer_5p, *lst_r_5p, *lst_f_5p;

  int bubble_items_read
    = sscanf(third_line,
             ">var_%lu_5p_flank length:%lu average_coverage: %f "
             "min_coverage:%lu max_coverage:%lu "
             "fst_coverage:%lu fst_kmer:%s fst_r:%s fst_f:%s "
             "lst_coverage:%lu lst_kmer:%s lst_r:%s lst_f:%s",
             &var_num, &length_5p, &mean_covg_5p,
             &min_covg_5p, &max_covg_5p,
             &fst_covg_5p, fst_kmer_5p, fst_r_5p, fst_f_5p,
             &lst_covg_5p, lst_kmer_5p, lst_r_5p, lst_f_5p);
  
  char first_char = tolower(fourth_line[0]);
    
  if(bubble_items_read == 13 &&
     (first_char == 'a' || first_char == 'c' ||
      first_char == 'g' || first_char == 't'))
  {
    cortex->filetype = BUBBLE_FILE;

    // read 6 lines, then start reading in pairs
    int i;
    for(i = 0; i < 6 && _cortex_read_line(cortex) > 0; i++);

    //  as long as the first line begins 'Covg ..' and
    //             the second line begins [0-9]
    //  then num_of_colours++
    while(_cortex_read_line(cortex) > 0 &&
          strncasecmp(cortex->buffer->buff, "Covg ", 5) == 0 &&
          _cortex_read_line(cortex) > 0 &&
          cortex->buffer->buff[0] >= '0' && cortex->buffer->buff[0] <= '9')
    {
      cortex->num_of_colours++;
    }
    
    // Reset file
    gzrewind(cortex->file);
    return cortex;
  }

  fprintf(stderr, "cortex.c: Couldn't determine file type (%s:%lu)",
          cortex->path, cortex->line_number);

  cortex_close(cortex);
  return NULL;
}

void cortex_close(CORTEX_FILE *cortex)
{
  if(cortex->buffer != NULL)
  {
    string_buff_free(cortex->buffer);
  }

  if(cortex->file != NULL)
  {
    gzclose(cortex->file);
  }

  free(cortex->path);
  free(cortex);
}


CORTEX_BUBBLE* cortex_bubble_create(const CORTEX_FILE *c_file)
{
  CORTEX_BUBBLE* bubble = (CORTEX_BUBBLE*) malloc(sizeof(bubble));

  unsigned long col;
  int branch;

  for(branch = 0; branch < 2; branch++)
  {
    bubble->branches_colour_covgs[branch]
      = (COLOUR_COVG**) malloc(c_file->num_of_colours * sizeof(COLOUR_COVG*));
  
    for(col = 0; col < c_file->num_of_colours; col++)
    {
      COLOUR_COVG** covgs = bubble->branches_colour_covgs[branch] + col;

      *covgs = (COLOUR_COVG*) malloc(sizeof(COLOUR_COVG));
      (*covgs)->colour_covgs = (unsigned long*) malloc(200 * sizeof(unsigned long));
      (*covgs)->capacity = 200;
      (*covgs)->length = 0;
    }
  }

  return bubble;
}

void cortex_bubble_free(CORTEX_BUBBLE* bubble, const CORTEX_FILE *c_file)
{
  unsigned long col;
  int branch;

  for(branch = 0; branch < 2; branch++)
  {
    for(col = 0; col < c_file->num_of_colours; col++)
    {
      COLOUR_COVG** covgs = bubble->branches_colour_covgs[branch] + col;

      free((*covgs)->colour_covgs);
      free(*covgs);
    }

    free(bubble->branches_colour_covgs[branch]);
  }

  free(bubble);
}

// line of numbers should have bean already read into c_file->buffer
char _read_covg(CORTEX_FILE *c_file, COLOUR_COVG *covgs, size_t required_size)
{
  // Ensure size
  if(covgs->capacity < required_size)
  {
    // Enlarge for all colours
    covgs->colour_covgs = realloc(covgs->colour_covgs,
                                  required_size * sizeof(unsigned long));

    covgs->capacity = required_size;
  }

  covgs->length = 0;

  char* pos = c_file->buffer->buff;
  char* new_pos;
  long value;

  while((value = strtol(pos, &new_pos, 10)) && pos != new_pos)
  {
    if(covgs->length == covgs->capacity)
    {
      fprintf(stderr, "cortex.c: more numbers than expected (%s:%lu)\n",
              c_file->path, c_file->line_number);
    }

    covgs->colour_covgs[covgs->length++] = value;
  }

  if(!string_is_all_whitespace(new_pos))
  {
    fprintf(stderr, "cortex.c: unexpected content on the end of line (%s:%lu)\n",
            c_file->path, c_file->line_number);
  }

  return 1;
}

//
// Alignments
//

char cortex_read_alignment(CORTEX_FILE* c_file, CORTEX_ALIGNMENT* alignment)
{
  if(c_file->filetype != ALIGNMENT_FILE)
  {
    fprintf(stderr, "cortex.c: cortex_read_alignment cannot read from "
                    "bubble file (%s:%lu)\n", c_file->path, c_file->line_number);

    return 0;
  }

  if(string_buff_get_char(c_file->buffer, 0) != '>')
  {
    fprintf(stderr, "cortex.c: cortex_read_alignment cannot read from "
                    "bubble file (%s:%lu)\n", c_file->path, c_file->line_number);

    return 0;
  }

  string_buff_copy(alignment->name, 0, c_file->buffer, 1,
                   string_buff_strlen(c_file->buffer));
  _cortex_read_line(c_file);

  string_buff_copy(alignment->seq, 0, c_file->buffer, 0,
                   string_buff_strlen(c_file->buffer));
  _cortex_read_line(c_file);

  unsigned long col;
  for(col = 0; col < c_file->num_of_colours; col++)
  {
    _read_covg(c_file, alignment->colour_covgs+col,
               string_buff_strlen(alignment->seq));
  }

  return 1;
}

void cortex_print_alignment(const CORTEX_ALIGNMENT* alignment)
{
  // DEV:
}

//
// Bubbles
//

// Returns 1 (success) or 0 (failure).  Path argument is where to store result
char _read_bubble_path(CORTEX_FILE *c_file, CORTEX_BUBBLE_PATH *path)
{
  unsigned long var_num;

  int bubble_items_read
    = sscanf(c_file->buffer->buff,
             ">var_%lu_5p_flank length:%lu average_coverage: %f "
             "min_coverage:%lu max_coverage:%lu "
             "fst_coverage:%lu fst_kmer:%s fst_r:%s fst_f:%s "
             "lst_coverage:%lu lst_kmer:%s lst_r:%s lst_f:%s",
             &var_num, &path->seq_length, &path->mean_covg,
             &path->min_covg, &path->max_covg,
             &path->fst_covg, path->fst_kmer, path->fst_r, path->lst_r,
             &path->lst_covg, path->lst_kmer, path->lst_r, path->lst_f);

  if(bubble_items_read != 13)
  {
    fprintf(stderr, "cortex.c: cortex_read_bubble couldn't parse line "
                    "(%s:%lu)\n", c_file->path, c_file->line_number);
    return 0;
  }

  // Read sequence line
  _cortex_read_line(c_file);

  if(gzeof(c_file->file))
  {
    return 0;
  }

  path->seq = string_buff_as_str(c_file->buffer);

  return 1;
}

char cortex_read_bubble(CORTEX_FILE* c_file, CORTEX_BUBBLE* bubble)
{
  if(c_file->filetype != BUBBLE_FILE)
  {
    fprintf(stderr, "cortex.c: cortex_read_bubble cannot read from "
                    "alignment file (%s:%lu)\n",
            c_file->path, c_file->line_number);

    return 0;
  }

  if(c_file->has_likelihoods)
  {
    // DEV: read in likelihoods
  }

  unsigned long var_num1 = _read_bubble_path(c_file, &bubble->flank_5p);
  unsigned long var_num2 = _read_bubble_path(c_file, &bubble->branches[0]);
  unsigned long var_num3 = _read_bubble_path(c_file, &bubble->branches[1]);
  unsigned long var_num4 = _read_bubble_path(c_file, &bubble->flank_3p);

  if(var_num1 != var_num2 || var_num2 != var_num3 || var_num3 != var_num4)
  {
    fprintf(stderr, "cortex.c: cortex_read_bubble() lines have different var "
                    "numbers (%s:%lu)\n", c_file->path, c_file->line_number);
    return 0;
  }

  bubble->var_num = var_num1;

  // Skip whitelines (_cortex_read_line does chomp)
  t_buf_pos chars_read;

  while((chars_read = _cortex_read_line(c_file)) > 0 &&
        string_buff_strlen(c_file->buffer) == 0);

  if(chars_read == 0)
  {
    fprintf(stderr, "cortex.c: file ended prematurely (%s:%lu)\n",
            c_file->path, c_file->line_number);

    return 0;
  }

  // Ensure we have enough space for coverage numbers
  unsigned long col;
  int branch;

  for(branch = 0; branch < 2; branch++)
  {
    unsigned long branch_length = bubble->branches[branch].seq_length;

    for(col = 0; col < c_file->num_of_colours; col++)
    {
      // Skip a line and then read the covg line of numbers
      if(_cortex_read_line(c_file) == 0 || _cortex_read_line(c_file) == 0)
      {
        fprintf(stderr, "cortex.c: missing content from end of file (%s:%lu)\n",
                c_file->path, c_file->line_number);
      }

      // Get coverage of branch 'branch' on colour 'col'
      COLOUR_COVG* colour_covgs = *bubble->branches_colour_covgs[branch] + col;
      _read_covg(c_file, colour_covgs, branch_length);
    }

    _cortex_read_line(c_file) == 0);
  }

  // Read until not whiteline
  while((chars_read = _cortex_read_line(c_file)) > 0 &&
        string_buff_strlen(c_file->buffer) == 0);

  return 1;
}

void cortex_print_bubble(const CORTEX_BUBBLE* bubble)
{
  // DEV:
}
