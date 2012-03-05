/*
 cortex_test.c
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cortex.h"

int main(int argc, char* argv[])
{
  if(argc != 2)
  {
    printf("Usage: %s <in.colour_covgs>\n", argv[0]);
    return EXIT_FAILURE;
  }

  CORTEX_FILE *cortex_file = cortex_open(argv[1]);

  if(cortex_file == NULL)
  {
    return EXIT_FAILURE;
  }

  fprintf(stderr, "Colours [%lu]: %s\n",
          cortex_file->num_of_colours, cortex_colour_list_str(cortex_file));

  if(cortex_file->filetype == BUBBLE_FILE)
  {
    CORTEX_BUBBLE* bubble = cortex_bubble_create(cortex_file);

    while(cortex_read_bubble(bubble, cortex_file))
    {
      cortex_print_bubble(bubble, cortex_file);
    }
  
    cortex_bubble_free(bubble, cortex_file);
  }
  else if(cortex_file->filetype == ALIGNMENT_FILE)
  {
    CORTEX_ALIGNMENT* alignment = cortex_alignment_create(cortex_file);
    
    while(cortex_read_alignment(alignment, cortex_file))
    {
      cortex_print_alignment(alignment, cortex_file);
    }
    
    cortex_alignment_free(alignment, cortex_file);
  }

  cortex_close(cortex_file);

  return EXIT_SUCCESS;
}
