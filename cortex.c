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

enum PATH_TYPE {FLANK_5P,FLANK_3P,BRANCH1,BRANCH2};

char *fail_msg = "FAILS CLASSIFIER:";
char *discovery_msg = "DISCOVERY PHASE:";

t_buf_pos _cortex_read_line(CORTEX_FILE* c_file)
{
  t_buf_pos chars_read = string_buff_reset_gzreadline(c_file->buffer,
                                                      c_file->file);
  string_buff_chomp(c_file->buffer);
  c_file->line_number++;
  
  //printf("Read: %s\n", c_file->buffer->buff);
  
  return chars_read;
}

t_buf_pos _cortex_read_reset(CORTEX_FILE* c_file)
{
  gzrewind(c_file->file);
  c_file->line_number = 0;
  return _cortex_read_line(c_file);
}

char* _parse_bubble_meta(const char* buff, const char* search,
                         char* result, CORTEX_FILE *c_file)
{
  char *hit;

  if((hit = strstr(buff, search)) == NULL)
  {
    fprintf(stderr, "cortex.c: cortex_read_bubble couldn't parse line "
                    "[no '%s'] (%s:%lu)\n",
            search, c_file->path, c_file->line_number);
    return 0;
  }

  hit = hit + strlen(search);
  
  if(*hit == 'A' || *hit == 'C' || *hit == 'G' || *hit == 'T')
  {
    sscanf(hit, "%64s ", result);
    hit += strlen(result);
  }
  else
  {
    result[0] = '\0';
  }

  while(isspace(*hit))
  {
    hit++;
  }

  return hit;
}

inline void _set_kmer_size(char *buffer, CORTEX_FILE *c_file)
{
  char *first_kmer = (char*)malloc(65*sizeof(char));
  _parse_bubble_meta(buffer, "fst_kmer:", first_kmer, c_file);
  c_file->kmer_size = (unsigned char)strlen(first_kmer);
  free(first_kmer);
}

// If cannot open file or file is empty will print error and return NULL
CORTEX_FILE* cortex_open(const char* path)
{
  CORTEX_FILE* c_file = (CORTEX_FILE*) malloc(sizeof(CORTEX_FILE));

  // Give initial values
  c_file->buffer = string_buff_init(500); // Create read in buffer
  c_file->line_number = 0;
  c_file->filetype = UNKNOWN_FILE;
  c_file->kmer_size = 0;
  c_file->num_of_colours = 0;
  // Likelihoods
  c_file->has_likelihoods = 0;
  c_file->is_diploid = 0;
  c_file->fails_classifier_line = 0;
  c_file->discovery_phase_line = 0;

  // Set path
  size_t path_len = strlen(path);
  c_file->path = (char*) malloc(path_len+1);
  strcpy(c_file->path, path);
  c_file->path[path_len] = '\0';

  // Open file
  c_file->file = gzopen(path, "r");

  if(c_file->file == NULL)
  {
    fprintf(stderr, "cortex.c: couldn't open file (%s)\n", path);
    cortex_close(c_file);
    return NULL;
  }

  // Whilst still reading but lines empty (_cortex_read_line does chomp)
  t_buf_pos chars_read;
  
  while((chars_read = _cortex_read_line(c_file)) > 0 &&
        string_buff_strlen(c_file->buffer) == 0);

  if(chars_read == 0)
  {
    fprintf(stderr, "cortex.c: file is empty (%s)\n", path);
    cortex_close(c_file);
    return NULL;
  }

  // Skip FAILS CLASSIFIER line
  if(strncasecmp(c_file->buffer->buff, fail_msg, strlen(fail_msg)) == 0)
  {
    c_file->fails_classifier_line = 1;

    if(_cortex_read_line(c_file) == 0)
    {
      fprintf(stderr, "cortex.c: file is mostly empty (%s)\n", path);
      cortex_close(c_file);
      return NULL;
    }
  }

  // Skip DISCOVERY PHASE line
  if(strncasecmp(c_file->buffer->buff, discovery_msg, strlen(discovery_msg)) == 0)
  {
    c_file->discovery_phase_line = 1;

    if(_cortex_read_line(c_file) == 0)
    {
      fprintf(stderr, "cortex.c: file is mostly empty (%s)\n", path);
      cortex_close(c_file);
      return NULL;
    }
  }

  if(strncasecmp(c_file->buffer->buff, "Colour", strlen("Colour")) == 0)
  {
    c_file->filetype = BUBBLE_FILE;
    c_file->has_likelihoods = 1;
  
    // Count lines begining [0-9] to count colours
    while(_cortex_read_line(c_file) > 0)
    {
      // Ignore blank lines
      if(string_buff_strlen(c_file->buffer) > 0)
      {
        char first_char = c_file->buffer->buff[0];

        if(first_char == '>')
        {
          break;
        }
        else if(isdigit(first_char))
        {
          c_file->num_of_colours++;
        }
      }
    }

    // Get kmer size
    _set_kmer_size(c_file->buffer->buff, c_file);

    // Reset file
    _cortex_read_reset(c_file);
    return c_file;
  }

  if(string_buff_get_char(c_file->buffer, 0) != '>')
  {
    fprintf(stderr, "cortex.c: unrecognised line (%s:%lu)\n",
            path, c_file->line_number);

    cortex_close(c_file);
    return NULL;
  }

  _cortex_read_line(c_file);
  unsigned long second_line_len = string_buff_strlen(c_file->buffer);

  _cortex_read_line(c_file);
  char* third_line = string_buff_as_str(c_file->buffer);

  _cortex_read_line(c_file);
  char* fourth_line = string_buff_as_str(c_file->buffer);

  // else if third_line ends "_colour_0_kmer_coverages" AND
  //         fourth_line starts with a [0-9]
  // then ALIGNMENT FILE
  char *expected_line = "_colour_0_kmer_coverages";
  char *third_line_search = strstr(third_line, expected_line);

  if(third_line_search != NULL &&
     third_line_search - third_line + strlen(expected_line)
       == strlen(third_line) &&
     isdigit(fourth_line[0]))
  {
    c_file->filetype = ALIGNMENT_FILE;
    c_file->num_of_colours = 1;

    // read in pairs on lines
    //  and as long as the second one begins [0-9] then num_of_colours++

    while(_cortex_read_line(c_file) > 0 && _cortex_read_line(c_file) > 0 &&
          isdigit(c_file->buffer->buff[0]))
    {
      c_file->num_of_colours++;
    }

    // Get kmer size
    char* digit_start = fourth_line;
    unsigned long num_of_kmers = 0;

    while(isdigit(*digit_start))
    {
      num_of_kmers++;
      while(isdigit(*digit_start)) {
        digit_start++;
      }
      while(!isspace(*digit_start)) {
        digit_start++;
      }
    }

    c_file->kmer_size = (unsigned char)(second_line_len - num_of_kmers + 1);

    // Reset file
    free(third_line);
    free(fourth_line);
    _cortex_read_reset(c_file);
    return c_file;
  }
  
  // else if third_line matches bubble path AND
  //         fourth_line start [ACGT]
  // then BUBBLE_FILE
  // 

  unsigned long var_num, length_5p;
  float mean_covg_5p;
  unsigned long min_covg_5p, max_covg_5p;

  int items_read;

  items_read  = sscanf(third_line,
                       ">branch_%lu_1 length:%lu average_coverage: %f "
                       "min_coverage:%lu max_coverage:%lu ",
                       &var_num, &length_5p, &mean_covg_5p,
                       &min_covg_5p, &max_covg_5p);

  char first_char = tolower(fourth_line[0]);
    
  if(items_read == 5 &&
     (first_char == 'a' || first_char == 'c' ||
      first_char == 'g' || first_char == 't'))
  {
    // Get kmer size
    _set_kmer_size(third_line, c_file);

    c_file->filetype = BUBBLE_FILE;

    // read 6 lines, then start reading in pairs
    int i;
    for(i = 0; i < 6 && _cortex_read_line(c_file) > 0; i++);

    //  as long as the first line begins 'Covg ..' and
    //             the second line begins [0-9]
    //  then num_of_colours++
    while(_cortex_read_line(c_file) > 0 &&
          strncasecmp(c_file->buffer->buff, "Covg ", 5) == 0 &&
          _cortex_read_line(c_file) > 0 && isdigit(c_file->buffer->buff[0]))
    {
      c_file->num_of_colours++;
    }

    // Reset file
    free(third_line);
    free(fourth_line);
    _cortex_read_reset(c_file);
    return c_file;
  }

  fprintf(stderr, "cortex.c: Couldn't determine file type (%s:%lu)\n",
          c_file->path, c_file->line_number);

  free(third_line);
  free(fourth_line);
  cortex_close(c_file);
  return NULL;
}

void cortex_close(CORTEX_FILE *c_file)
{
  if(c_file->buffer != NULL)
  {
    string_buff_free(c_file->buffer);
  }

  if(c_file->file != NULL)
  {
    gzclose(c_file->file);
  }

  free(c_file->path);
  free(c_file);
}

COLOUR_COVG* _create_colour_covgs()
{
  COLOUR_COVG* covgs = (COLOUR_COVG*) malloc(sizeof(COLOUR_COVG));
  covgs->colour_covgs = (unsigned long*) malloc(200 * sizeof(unsigned long));
  covgs->capacity = 200;
  covgs->length = 0;

  return covgs;
}

void _free_colour_covgs(COLOUR_COVG* covgs)
{
  free(covgs->colour_covgs);
  free(covgs);
}

CORTEX_BUBBLE* cortex_bubble_create(const CORTEX_FILE *c_file)
{
  CORTEX_BUBBLE* bubble = (CORTEX_BUBBLE*) malloc(sizeof(CORTEX_BUBBLE));

  unsigned long col;
  int branch;

  for(branch = 0; branch < 2; branch++)
  {
    bubble->branches_colour_covgs[branch]
      = (COLOUR_COVG**) malloc(c_file->num_of_colours * sizeof(COLOUR_COVG*));
  
    for(col = 0; col < c_file->num_of_colours; col++)
    {
      bubble->branches_colour_covgs[branch][col] = _create_colour_covgs();
    }
  }

  bubble->flank_5p.seq = string_buff_init(200);
  bubble->flank_3p.seq = string_buff_init(200);
  bubble->branches[0].seq = string_buff_init(200);
  bubble->branches[1].seq = string_buff_init(200);

  if(bubble->branches[0].seq == NULL)
  {
    fprintf(stderr, "cortex.c: mem fail\n");
    exit(EXIT_FAILURE);
  }

  bubble->calls
    = (HETEROGENEITY*) malloc(c_file->num_of_colours * sizeof(HETEROGENEITY));

  bubble->llk_hom_br1 = (float*) malloc(c_file->num_of_colours * sizeof(float));
  bubble->llk_het     = (float*) malloc(c_file->num_of_colours * sizeof(float));
  bubble->llk_hom_br2 = (float*) malloc(c_file->num_of_colours * sizeof(float));

  if(bubble->llk_hom_br1 == NULL || bubble->llk_het == NULL ||
     bubble->llk_hom_br2 == NULL)
  {
    fprintf(stderr, "cortex.c: Couldn't allocate enough memory\n");
    exit(EXIT_FAILURE);
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
      _free_colour_covgs(bubble->branches_colour_covgs[branch][col]);
    }

    free(bubble->branches_colour_covgs[branch]);
  }

  string_buff_free(bubble->flank_5p.seq);
  string_buff_free(bubble->flank_3p.seq);
  string_buff_free(bubble->branches[0].seq);
  string_buff_free(bubble->branches[1].seq);

  free(bubble->calls);
  free(bubble->llk_hom_br1);
  free(bubble->llk_het);
  free(bubble->llk_hom_br2);

  free(bubble);
}

CORTEX_ALIGNMENT* cortex_alignment_create(const CORTEX_FILE *c_file)
{
  CORTEX_ALIGNMENT* alignment
    = (CORTEX_ALIGNMENT*) malloc(sizeof(CORTEX_ALIGNMENT));

  alignment->name = string_buff_init(200);
  alignment->seq = string_buff_init(200);

  alignment->colour_covgs
    = (COLOUR_COVG**) malloc(c_file->num_of_colours * sizeof(COLOUR_COVG*));

  unsigned long col;
  for(col = 0; col < c_file->num_of_colours; col++)
  {
    alignment->colour_covgs[col] = _create_colour_covgs();
  }

  return alignment;
}

void cortex_alignment_free(CORTEX_ALIGNMENT* alignment, const CORTEX_FILE *c_file)
{
  string_buff_free(alignment->name);
  string_buff_free(alignment->seq);

  unsigned long col;
  for(col = 0; col < c_file->num_of_colours; col++)
  {
    _free_colour_covgs(alignment->colour_covgs[col]);
  }

  free(alignment->colour_covgs);
  free(alignment);
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
  unsigned long value;

  value = strtoul(pos, &new_pos, 10);

  while(pos != new_pos)
  {
    if(covgs->length == covgs->capacity)
    {
      fprintf(stderr, "cortex.c: more numbers than expected [%lu] (%s:%lu)\n",
              required_size, c_file->path, c_file->line_number);
    }

    covgs->colour_covgs[covgs->length++] = value;

    pos = new_pos;
    value = strtoul(pos, &new_pos, 10);
  }

  if(!string_is_all_whitespace(new_pos))
  {
    fprintf(stderr, "cortex.c: unexpected content on the end of line ['%s'] "
                    "(%s:%lu)\n",
            new_pos, c_file->path, c_file->line_number);
  }

  return 1;
}

//
// Alignments
//

char cortex_read_alignment(CORTEX_ALIGNMENT* alignment, CORTEX_FILE* c_file)
{
  if(c_file->filetype != ALIGNMENT_FILE)
  {
    fprintf(stderr, "cortex.c: cortex_read_alignment cannot read from "
                    "alignment file (%s:%lu)\n",
            c_file->path, c_file->line_number);

    return 0;
  }

  // Read until not whiteline
  while(string_buff_strlen(c_file->buffer) == 0 &&
        _cortex_read_line(c_file) > 0);

  if(string_buff_strlen(c_file->buffer) == 0)
  {
    // EOF
    return 0;
  }

  if(string_buff_get_char(c_file->buffer, 0) != '>')
  {
    fprintf(stderr, "cortex.c: cortex_read_alignment line doesn't start '>' "
                    "(%s:%lu)\n",
            c_file->path, c_file->line_number);

    return 0;
  }

  string_buff_copy(alignment->name, 0, c_file->buffer, 1,
                   string_buff_strlen(c_file->buffer));

  if(_cortex_read_line(c_file) == 0)
  {
    fprintf(stderr, "cortex.c: alignment ended early (%s:%lu)\n",
            c_file->path, c_file->line_number);
    return 0;
  }

  string_buff_copy(alignment->seq, 0, c_file->buffer, 0,
                   string_buff_strlen(c_file->buffer));

  unsigned long col;
  for(col = 0; col < c_file->num_of_colours; col++)
  {
    if(_cortex_read_line(c_file) == 0 || _cortex_read_line(c_file) == 0)
    {
      fprintf(stderr, "cortex.c: alignment ended early (%s:%lu)\n",
              c_file->path, c_file->line_number);
      return 0;
    }

    _read_covg(c_file, alignment->colour_covgs[col],
               string_buff_strlen(alignment->seq));
  }

  _cortex_read_line(c_file);

  return 1;
}

void cortex_print_alignment(const CORTEX_ALIGNMENT* alignment,
                            const CORTEX_FILE* c_file)
{
  printf(">%s\n", alignment->name->buff);
  printf("%s\n", alignment->seq->buff);

  unsigned long col, covgs_i;

  for(col = 0; col < c_file->num_of_colours; col++)
  {
    COLOUR_COVG* covgs = alignment->colour_covgs[col];

    printf(">%s_colour_%lu_kmer_coverages\n", alignment->name->buff, col);
    printf("%lu", covgs->colour_covgs[0]);

    for(covgs_i = 0; covgs_i < covgs->length; covgs_i++)
    {
      printf(" %lu", covgs->colour_covgs[covgs_i]);
    }

    printf("\n");
  }
}

//
// Bubbles
//

void _print_bubble_path(const unsigned long var_num,
                        const CORTEX_BUBBLE_PATH *bp,
                        enum PATH_TYPE path_type)
{
  switch (path_type)
  {
    case FLANK_5P:
      printf(">var_%lu_5p_flank ", var_num);
      break;
    case BRANCH1:
      printf(">branch_%lu_1 ", var_num);
      break;
    case BRANCH2:
      printf(">branch_%lu_2 ", var_num);
      break;
    case FLANK_3P:
      printf(">var_%lu_3p_flank ", var_num);
      break;
    default:
      break;
  }

  printf("length:%lu average_coverage: %f min_coverage:%lu max_coverage:%lu "
         "fst_coverage:%lu fst_kmer:%s fst_r:%s fst_f:%s "
         "lst_coverage:%lu lst_kmer:%s lst_r:%s lst_f:%s\n",
         bp->seq_length, bp->mean_covg, bp->min_covg, bp->max_covg,
         bp->fst_covg, bp->fst_kmer, bp->fst_r, bp->fst_f,
         bp->lst_covg, bp->lst_kmer, bp->lst_r, bp->lst_f);

  printf("%s\n", bp->seq->buff);
}

// Returns 1 (success) or 0 (failure).  Path argument is where to store result
char _read_bubble_path(CORTEX_FILE *c_file, CORTEX_BUBBLE_PATH *path)
{
  // Line looks like:
  // >var_1_5p_flank length:50 average_coverage: 2.00 min_coverage:2 
  // max_coverage:2 fst_coverage:2 fst_kmer:GACCATAGCAAGGACAC fst_r: fst_f:C 
  // lst_coverage:2 lst_kmer:ACGTTCAACGCCAAGGG lst_r:C lst_f:AT 

  char var_name[100];

  int items_read
    = sscanf(c_file->buffer->buff,
             ">%50s length:%lu average_coverage: %f "
             "min_coverage:%lu max_coverage:%lu "
             "fst_coverage:%lu fst_kmer:%s ",
             var_name, &path->seq_length, &path->mean_covg,
             &path->min_covg, &path->max_covg, &path->fst_covg, path->fst_kmer);

  if(items_read != 7)
  {
    fprintf(stderr, "cortex.c: cortex_read_bubble couldn't parse line "
                    "[%i items] (%s:%lu)\n",
            items_read, c_file->path, c_file->line_number);
    return 0;
  }

  // Parse var_name
  unsigned long var_num;
  
  if(!sscanf(var_name, "var_%lu_5p_flank", &var_num) && 
     !sscanf(var_name, "branch_%lu_1", &var_num) && 
     !sscanf(var_name, "branch_%lu_2", &var_num) && 
     !sscanf(var_name, "var_%lu_3p_flank", &var_num))
  {
    fprintf(stderr, "cortex.c: cortex_read_bubble couldn't parse name "
                    "['%s'] (%s:%lu)\n",
            var_name, c_file->path, c_file->line_number);
    return 0;
  }

  // Read the rest of the line
  char* end;
  
  end = _parse_bubble_meta(c_file->buffer->buff, "fst_r:", path->fst_r, c_file);
  end = _parse_bubble_meta(end, "fst_f:", path->fst_f, c_file);

  items_read = sscanf(end, "lst_coverage:%lu lst_kmer:%s ",
                      &path->lst_covg, path->lst_kmer);

  if(items_read != 2)
  {
    fprintf(stderr, "cortex.c: cortex_read_bubble couldn't parse line "
                    "[%i items] (%s:%lu)\n",
            items_read, c_file->path, c_file->line_number);
    return 0;
  }

  end = _parse_bubble_meta(end, "lst_r:", path->lst_r, c_file);
  end = _parse_bubble_meta(end, "lst_f:", path->lst_f, c_file);

  // Read sequence line
  if(_cortex_read_line(c_file) == 0)
  {
    return 0;
  }

  string_buff_copy(path->seq, 0, c_file->buffer, 0,
                   string_buff_strlen(c_file->buffer));

  // Load up next line
  _cortex_read_line(c_file);

  return 1;
}

char cortex_read_bubble(CORTEX_BUBBLE* bubble, CORTEX_FILE* c_file)
{
  if(c_file->filetype != BUBBLE_FILE)
  {
    fprintf(stderr, "cortex.c: cortex_read_bubble cannot read from "
                    "alignment file (%s:%lu)\n",
            c_file->path, c_file->line_number);

    return 0;
  }

  // Read until not whiteline
  while(string_buff_strlen(c_file->buffer) == 0 &&
        _cortex_read_line(c_file) > 0);

  if(string_buff_strlen(c_file->buffer) == 0)
  {
    // EOF
    return 0;
  }

  if((c_file->fails_classifier_line && (_cortex_read_line(c_file) == 0)) ||
     (c_file->discovery_phase_line && (_cortex_read_line(c_file) == 0)))
  {
    fprintf(stderr, "cortex.c: premature end of call start (%s:%lu)\n",
            c_file->path, c_file->line_number);
    return 0;
  }

  if(c_file->has_likelihoods)
  {
    // Read in likelihoods
    if(strncasecmp(c_file->buffer->buff, "Colour", strlen("Colour")) != 0)
    {
      fprintf(stderr, "cortex.c: premature end of file likelihoods (%s:%lu)\n",
                c_file->path, c_file->line_number);
      return 0;
    }

    // Check if diploid (ie. has het. option)
    if(strstr(c_file->buffer->buff, "llk_het") != NULL)
    {
      c_file->is_diploid = 1;
    }

    unsigned long col;
    for(col = 0; col < c_file->num_of_colours; col++)
    {
      if(_cortex_read_line(c_file) == 0 || !isdigit(c_file->buffer->buff[0]))
      {
        fprintf(stderr, "cortex.c: premature end of call likelihoods (%s:%lu)\n",
                c_file->path, c_file->line_number);
        return 0;
      }

      unsigned long col2;
      char str[6];
      float col_llk_hom_br1, col_llk_het, col_llk_hom_br2;

      int items_read;
      
      if(c_file->is_diploid)
      {
        items_read = sscanf(c_file->buffer->buff, "%lu %4s %f %f %f",
                            &col2, str, &col_llk_hom_br1, &col_llk_het,
                            &col_llk_hom_br2);
      }
      else
      {
        // haploid - can't be het
        items_read = sscanf(c_file->buffer->buff, "%lu %4s %f %f",
                            &col2, str, &col_llk_hom_br1, &col_llk_hom_br2);
      }

      if((c_file->is_diploid && items_read != 5) ||
         (!c_file->is_diploid && items_read != 4))
      {
        fprintf(stderr, "cortex.c: invalid likelihood line ['%s'] (%s:%lu)\n",
                c_file->buffer->buff, c_file->path, c_file->line_number);
        return 0;
      }

      HETEROGENEITY call = UNKNOWN_HET;

      if(strncasecmp(str, "HOM1", strlen("HOM1")) == 0)
      {
        call = HOM1;
      }
      else if(strncasecmp(str, "HET", strlen("HET")) == 0)
      {
        call = HET;
      }
      else if(strncasecmp(str, "HOM2", strlen("HOM2")) == 0)
      {
        call = HOM2;
      }
      else
      {
        fprintf(stderr, "cortex.c: unexpected likelihood line ['%s'] (%s:%lu)\n",
                str, c_file->path, c_file->line_number);
        return 0;
      }

      bubble->calls[col] = call;
      bubble->llk_hom_br1[col] = col_llk_hom_br1;

      if(c_file->is_diploid)
      {
        bubble->llk_het[col] = col_llk_het;
      }

      bubble->llk_hom_br2[col] = col_llk_hom_br2;
    }
  
    // Read first line of the bubble
    _cortex_read_line(c_file);
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
        return 0;
      }

      // Get coverage of branch 'branch' on colour 'col'
      _read_covg(c_file, bubble->branches_colour_covgs[branch][col], branch_length);
    }

    _cortex_read_line(c_file);
  }

  // Read until not whiteline
  while((chars_read = _cortex_read_line(c_file)) > 0 &&
        string_buff_strlen(c_file->buffer) == 0);

  return 1;
}

void cortex_print_bubble(const CORTEX_BUBBLE* bubble, const CORTEX_FILE *c_file)
{
  if(c_file->fails_classifier_line)
  {
    printf("FAILS CLASSIFIER: fits repeat model better than variation model\n");
  }

  if(c_file->discovery_phase_line)
  {
    printf("DISCOVERY PHASE:  VARIANT vs REPEAT MODEL "
           "LOG_LIKELIHOODS:	llk_var:nan	llk_rep:-inf\n");
  }

  if(c_file->has_likelihoods)
  {
    if(c_file->is_diploid)
    {
      printf("Colour/sample	GT_call	llk_hom_br1	llk_het	llk_hom_br2\n");
    }
    else
    {
      printf("Colour/sample	GT_call	llk_hom_br1	llk_hom_br2\n");
    }

    unsigned long col;
    for(col = 0; col < c_file->num_of_colours; col++)
    {
      char *call = NULL;

      switch(bubble->calls[col])
      {
        case HOM1:
          call = "HOM1";
          break;
        case HET:
          call = "HET";
          break;
        case HOM2:
          call = "HOM2";
          break;
        default:
          fprintf(stderr, "cortex.c: unknown call heterogeneity [%i] (%s:%lu)\n",
                  bubble->calls[col], c_file->path, c_file->line_number);
          break;
      }

      if(c_file->is_diploid)
      {
        // Has het
        printf("%lu	%s	%.2f	%.2f	%.2f\n", col, call,
               bubble->llk_hom_br1[col], bubble->llk_het[col],
               bubble->llk_hom_br2[col]);
      }
      else
      {
        printf("%lu	%s	%.2f	%.2f\n", col, call,
               bubble->llk_hom_br1[col],  bubble->llk_hom_br2[col]);
      }
    }
  }

  _print_bubble_path(bubble->var_num, &bubble->flank_5p, FLANK_5P);
  _print_bubble_path(bubble->var_num, &bubble->branches[0], BRANCH1);
  _print_bubble_path(bubble->var_num, &bubble->branches[1], BRANCH2);
  _print_bubble_path(bubble->var_num, &bubble->flank_3p, FLANK_3P);

  printf("\n\n");

  // Print branches
  int branch;
  unsigned long col, covgs_i;

  for(branch = 0; branch < 2; branch++)
  {
    printf("branch%i coverages\n", branch);
    
    for(col = 0; col < c_file->num_of_colours; col++)
    {
      COLOUR_COVG* covgs = bubble->branches_colour_covgs[branch][col];

      printf("Covg in Col %lu:\n", col);
      printf("%lu", covgs->colour_covgs[0]);

      for(covgs_i = 0; covgs_i < covgs->length; covgs_i++)
      {
        printf(" %lu", covgs->colour_covgs[covgs_i]);
      }

      printf("\n");
    }
  }

  printf("\n\n");
}
