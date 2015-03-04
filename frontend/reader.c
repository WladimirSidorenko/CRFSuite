/*
 *        Data reader.
 *
 * Copyright (c) 2007-2010, Naoaki Okazaki
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the names of the authors nor the names of its contributors
 *       may be used to endorse or promote products derived from this
 *       software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* $Id$ */

#include <os.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <crfsuite.h>
#include "iwa.h"

static int progress(FILE *fpo, int prev, int current)
{
  while (prev < current) {
    ++prev;
    if (prev % 2 == 0) {
      if (prev % 10 == 0) {
	fprintf(fpo, "%d", prev / 10);
	fflush(fpo);
      } else {
	fprintf(fpo, ".");
	fflush(fpo);
      }
    }
  }
  return prev;
}

/**
 * Function for reading training data.
 *
 * @param fpi - input file
 * @param fpo - output file
 * @param data - pointer to data instance
 * @param group - pointer to data instance
 * @param ftype - type of trained model (affects the way in which input data
 * are interpreted)
 *
 * @return number of instances read
 */
int read_data(FILE *fpi, FILE *fpo, crfsuite_data_t* data, int group, \
	      crfsuite_trainer_t *trainer)
{
  crfsuite_dictionary_t *attrs = data->attrs;
  crfsuite_dictionary_t *labels = data->labels;
  crfsuite_dictionary_t *node_labels = data->node_labels;

  int n = 0;
  int lid = -1;
  int ftype = trainer->ftype;
  unsigned attr_cnt = 0;
  crfsuite_instance_t inst;
  crfsuite_item_t item;
  crfsuite_attribute_t cont;
  iwa_t* iwa = NULL;
  const iwa_token_t *token = NULL;
  long filesize = 0, begin = 0, offset = 0;
  int prev = 0, current = 0, ret = 0;

  /* Initialize instance.*/
  crfsuite_instance_init(&inst);
  inst.group = group;

  /* Obtain file size. */
  begin = ftell(fpi);
  fseek(fpi, 0, SEEK_END);
  filesize = ftell(fpi) - begin;
  fseek(fpi, begin, SEEK_SET);

  /* */
  fprintf(fpo, "0");
  fflush(fpo);
  prev = 0;

  iwa = iwa_reader(fpi);
  while ((token = iwa_read(iwa))) {
    /* Progress report. */
    offset = ftell(fpi);
    current = (int)((offset - begin) * 100.0 / (double)filesize);
    prev = progress(fpo, prev, current);

    switch (token->type) {
    case IWA_BOI:
      /* Initialize an item. */
      lid = -1;
      attr_cnt = 0;
      crfsuite_item_init(&item);
      break;

    case IWA_EOI:
      /* check that node id is specified for tree CRF's */
      if (ftype == FTYPE_CRF1TREE && attr_cnt < 2) {
	fprintf(stderr, "ERROR: Incorrect number of attributes specified for tree (%d instead of %d)",
		attr_cnt, 2);
	n = 5;
	goto clear_exit;
      }

      if (0 <= lid)
	crfsuite_instance_append(&inst, &item, lid);

      crfsuite_item_finish(&item);
      break;

    case IWA_ITEM:
      ++attr_cnt;
      if (lid == -1) {
	lid = labels->get(labels, token->attr);
      } else {
	if (ftype == FTYPE_CRF1TREE) {
	  if (attr_cnt == 2) {
	    // check that same id is not used twice for different nodes within
	    // an instance
	    item.id = node_labels->get(node_labels, token->attr);
	    // remember string label of this node
	    item.node_label = (char *) malloc(sizeof(char) * (strlen(token->attr) + 1));
	    if (item.node_label) {
	      // be sure to delete node_label at the end
	      strcpy(item.node_label, token->attr);
	    } else {
	      fprintf(stderr, "ERROR: Could not allocate memory for storing node label '%s'.\n", token->attr);
	      goto clear_exit;
	    }
	    break;
	  } else if (attr_cnt == 3) {
	    if (strcmp(token->attr, "_") == 0)
	      item.prnt = -1;
	    else
	      item.prnt = node_labels->get(node_labels, token->attr);
	    break;
	  }
	}
	crfsuite_attribute_init(&cont);
	cont.aid = attrs->get(attrs, token->attr);
	if (token->value && *token->value)
	  cont.value = atof(token->value);
	else
	  cont.value = 1.0;
	crfsuite_item_append_attribute(&item, &cont);
      }
      break;

    case IWA_NONE:
    case IWA_EOF:
      /* perform some sanity check and create a tree from nodes */
      if (ftype == FTYPE_CRF1TREE && (ret = crfsuite_tree_init(&inst)) != 0) {
	fprintf(stderr, "ERROR: Could not create tree for training instance '%d'.\n", n);
	n = -1;
	goto clear_exit;
      }
      /* Add training instance to data. */
      crfsuite_data_append(data, &inst);
      crfsuite_instance_finish(&inst);

      inst.group = group;
      ++n;

      /* clear dictionary of node labels so that new instances will
	 have dense representation of node ids again */
      if (ftype == FTYPE_CRF1TREE)
	node_labels->reset(node_labels);

      break;
    }
  }
  progress(fpo, prev, 100);
  fprintf(fpo, "\n");

 clear_exit:
  if (ftype == FTYPE_CRF1TREE)
    node_labels->reset(node_labels);

  return n;
}
