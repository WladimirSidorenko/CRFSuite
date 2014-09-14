/*
 *      Data structures and auxiliary functions for semi-markov CRFs.
 *
 * Copyright (c) 2014-2015, Uladzimir Sidarenka
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

/* Libraries */
#include "ring.h"
#include "semimarkov.h"

/* Macros */
#define CLEAR(a_item)					\
  if ((a_item)) {					\
  free(a_item);						\
  a_item = NULL;					\
  }

/* Implementation */
static void crf1de_semimarkov_finish(crf1de_semimarkov_t *sm)
{
  CLEAR(max_seg_len);
  CLEAR(fs_llabels);
  CLEAR(forward_trans1);
  CLEAR(forward_trans2);

  CLEAR(backward_trans);

  CLEAR(pattern_trans1);
  CLEAR(pattern_trans2);

  if (sm->forward_states) {
    rumavl_destroy(crf1de->forward_states);
    sm->forward_states = NULL;
    sm->num_fs = 0;
  }

  if (crf1de->backward_states) {
    rumavl_destroy(crf1de->backward_states);
    crf1de->backward_states = NULL;
    crf1de->num_bs = 0;
  }
}

crf1de_semimarkov_t *crf1de_create_semimarkov(void);
 crf1de_semimarkov_t *sm = calloc(1, sizeof(crf1de_semimarkov_t));
 if (sm == NULL)
   return sm;

 sm->finish = crf1de_semimarkov_finish;
 return sm;
}

#i

f

#define MYMACRO 1

#ifndef MYMACRO
/* Initialize state dictionaries for crf1de instance.

   @param a_crf1de - crf1de instance for which new affixes should be initialized
   @param L - number of distinct tags
   @param a_wb - workbench for constructing affixes (should provide
                 capacity for at least 2 labels, unless a_max_order <= 0)
   @param a_max_order - maximum order of transition features
   @param a_seg_len_constr - boolean indicating whether maximum
                             segment length is constrained or not

   @return \c 0 if function completed successfully, non-\c 0 otherwise
 */
static int crf1de_sm_init_states(crf1de_t *a_crf1de, const int L,	\
				     int *a_wb, const int a_max_order, \
				    const int a_seg_len_constr)
{
  /* create an AVL for storing prefixes and suffixes (a prefix can
     have a maximum of `max_order` labels; suffixes have length 1
     greater than prefix) */
  a_crf1de->forward_states = rumavl_new(a_max_order * sizeof(int), \
					crf1de_cmp_lseq, NULL, NULL);
  if (a_crf1de->forward_states == NULL)
    return CRFSUITEERR_OUTOFMEMORY;

  a_crf1de->backward_states = rumavl_new((a_max_order + 1) * sizeof(int), \
					 crf1de_cmp_lseq, NULL, NULL);
  if (a_crf1de->backward_states == NULL)
    return CRFSUITEERR_OUTOFMEMORY;

  a_crf1de->num_fs = a_crf1de->num_bs = 0;
  /* if max_order is greater than zero, then unconditionally remember
     all existing tags as prefixes */
  if (a_max_order > 0) {
    /* remember the empty prefix */
    a_wb[0] = -1;
    rumavl_insert(a_crf1de->forward_states, a_wb);
    /* remember all tags */
    int i, j;
    for (i = 0; i < L; ++i) {
      /* insert tag `i` as prefix and suffix */
      a_wb[0] = i;
      rumavl_insert(a_crf1de->forward_states, a_wb);
      rumavl_insert(a_crf1de->backward_states, a_wb);
      /* insert tag sequence `ij` as suffix */
      for (j = 0; j < L; ++j) {
	a_wb[1] = j;
	if (i == j && ! a_seg_len_constr)
	  continue;
	rumavl_insert(a_crf1de->backward_states, a_wb);
      }
      /* reset last element of tag sequence to -1 */
      a_wb[1] = -1;
    }
    a_crf1de->num_fs = L + 1;
    a_crf1de->num_bs = L * (a_seg_len_constr? L + 1: L);
  }
  return 0;
}

/* Generate and add new states to states dictionaries.

   @param a_crf1de - crf1de instance to which new affixes should be added
   @param a_labels - ring of collected labels used to generate new affixes
   @param a_wb - workbench for constructing affixes (its capacity
                 should be one more than current number of items in `a_labels`)
   @param a_wb_size - size of workbench
   @param L - number of distinct tags
   @param a_sl_constr - boolean indicating whether constraint on
                        segment length is imposed or not

   @return \c void
 */
static void crf1de_sm_add_states(crf1de_t *const a_crf1de,			\
			      const crfsuite_ring_t *const a_labels,	\
			       int *const a_wb, const size_t a_wb_size, \
			      const int L, const int a_sl_constr)
{
  memset(a_wb, -1, a_wb_size);
  RUMAVL *forward_states = a_crf1de->forward_states;
  RUMAVL *backward_states = a_crf1de->backward_states;
  crfsuite_chain_link_t *clink = a_labels->head;
  /* produce all prefixes from the ring and generate suffixes for them */
  int i, j, crnt;
  for (i = 0; i < a_labels->num_items; ++i) {
    a_wb[i] = crnt = clink->data;
    /* insert prefix in dictionary */
    if (rumavl_find(forward_states, a_wb) == NULL) {
      rumavl_insert(forward_states, a_wb);
      ++a_crf1de->num_fs;
      /* construct and insert suffixes */
      for (j = 0; j < L; ++j) {
	if (crnt == j && ! a_sl_constr)
	  continue;

	a_wb[i+1] = j;
	rumavl_insert(backward_states, a_wb);
      }
      a_wb[i+1] = -1;
      a_crf1de->num_bs += a_sl_constr? L: L - 1;
    }
    /* proceed to next element */
    clink = clink->next;
  }
}

/* Initialize forward transition tables for semi-Markov CRFs.

   @param crf1de - pointer to crf1de instance for which transition should be
                   generated
   @param L - number of labels

   @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static int crf1de_sm_init_fw_trans(crf1de_t *a_crf1de, const int L) {
  ;
}

/* Initialize backward transition tables for semi-Markov CRFs.

   @param crf1de - pointer to crf1de instance for which transition should be
                   generated
   @param L - number of labels

   @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static int crf1de_sm_init_bw_trans(crf1de_t *a_crf1de, const int L) {
  ;
}

/* Initialize backward transition tables for semi-Markov CRFs.

   @param crf1de - pointer to crf1de instance for which transition should be
                   generated
   @param L - number of labels

   @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static int crf1de_sm_init_ptrn_trans(crf1de_t *a_crf1de, const int L) {
  ;
}

/* Initialize transition tables for semi-Markov CRFs.

   @param crf1de - pointer to crf1de instance for which transition should be
                   generated
   @param L - number of labels

   @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static int crf1de_sm_init_trans(crf1de_t *a_crf1de, const int L)
{
  int ret = (crf1de_sm_init_fw_trans(a_crf1de, L) ||	\
	     crf1de_sm_init_bw_trans(a_crf1de, L) ||	\
	     crf1de_sm_init_ptrn_trans(a_crf1de, L));
  return ret;
}

/* Set data specific to semi-Markov CRF.

   @param crf1de - interface to graphical model
   @param ds - pointer to dataset
   @param L - number of labels
   @param N - number of instances
   @param T - pointer to int storing maximum number of items for given instance

   @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static int crf1de_set_semimarkov(crf1de_t *crf1de, dataset_t *ds, \
				 const int L, const int N, int *T)
{
  int i, j, t, prev, crnt, nitems, seg_len, ret = 0;
  const crfsuite_instance_t *inst = NULL;
  crfsuite_dictionary_t *labels = ds->data->labels;


  /* create and populate initial set of affixes */
  memset(wb, -1, wb_size);
  crf1de_sm_init_states(crf1de, L, wb, max_order, rest_seg_len);

  for (j = 0; j < nitems; ++j) {
    crnt = inst->labels[j];
    if (prev < 0) {		/* maybe, subject to change because -1 is not pushed on ring */
      seg_len = 1;
    } else if (prev != crnt) {
      /* update maximum segment length, if necessary */
      if (seg_len > crf1de->max_seg_len[prev])
	crf1de->max_seg_len[prev] = seg_len;
      /* add new label to the ring of prefixes and remeber the new prefix and suffix */
      prev_labels->push(prev_labels, prev);
      crf1de_sm_add_states(crf1de, prev_labels, wb, wb_size, L, restr_seg_len);
      seg_len = 1;
    } else {
      ++seg_len;
    }
    prev = crnt;
  }
  crf1de_sm_add_states(crf1de, prev_labels, wb, wb_size, L, restr_seg_len);
  if (seg_len > crf1de->max_seg_len[prev])
    crf1de->max_seg_len[prev] = seg_len;

  /* store maximum number of items per instance */
  *T = t;
  /* if user specified positive maximum segment length, store
     specified value as maximum segment length for all labels */
  if (restr_seg_len) {
    for (i = 0; i < L; ++i) {
      crf1de->max_seg_len[i] = crf1de->opt.feature_max_seg_len;
    }
  }
  /* TODO: check prefixes; */
  RUMAVL_NODE *node = NULL;
  while ((node = rumavl_node_next(crf1de->forward_states, node, 1, (void**)&wb)) != NULL) {
    printf("prefix = '");
    for (i = 0; i < max_order && wb[i] != -1; ++i) {
      printf("%d", wb[i]);
    }
      printf("'\n");
  }
  printf("crfde->num_prefixes = %d\n", crf1de->num_fs);

  /* TODO: check suffixes; */
  node = NULL;
  while ((node = rumavl_node_next(crf1de->backward_states, node, 1, (void**)&wb)) != NULL) {
    printf("suffix = '");
    for (i = 0; i <= max_order && wb[i] != -1; ++i) {
      printf("%d", wb[i]);
    }
      printf("'\n");
  }
  /* TODO: check transitions; */
  printf("crfde->num_suffixes = %d\n", crf1de->num_bs);
  /* TODO: generate last labels and transitions; */
  crf1de_sm_init_trans(crf1de, L);
  exit(66);

 final_steps:
  prev_labels->free(prev_labels);
  free(wb);
  return ret;
}
#endif
