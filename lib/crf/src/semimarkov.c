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
#include "semimarkov.h"

#include <assert.h>		/* for assert() */
#include <stdlib.h>		/* for calloc() */
#include <string.h>		/* for memset() */

/* Macros */
#define CLEAR(a_item)					\
  if ((a_item)) {					\
    free(a_item);					\
    a_item = NULL;					\
  }

/* Implementation */


/** Function for comparing label sequences.
 *
 * @param a_lseq1 - first label sequence to compare
 * @param a_lseq2 - second label sequence to compare
 * @param a_size - maximum sequence size
 * @param a_udata - unused parameter (needed for compliance)
 *
 * @return > \c 1 if `a_lseq1` is greater than `a_lseq2`, \c 0 if both
 * sequences are the same, and < 0 if `a_lseq1` is smaller than
 * `a_lseq2`
*/
static int crf1de_cmp_lseq(const void *a_lseq1, const void *a_lseq2,	\
			   size_t a_size, void *a_udata)
{
  int ret = 0;
  size_t n = a_size / sizeof(int);
  const int *el1 = (const int*) a_lseq1;
  const int *el2 = (const int*) a_lseq2;

  for (size_t i = 0; i < n && ret == 0; ++i) {
    ret = el1[i] - el2[i];
    /* -1 terminates tagging sequence */
    if (el1[i] < 0 || el2[i] < 0)
      break;
  }
  return ret;
}

int crf1de_semimarkov_find_lng_sfx(crf1de_semimarkov_t *sm, const int *a_seqstart, \
				   const int *a_seqend)
{

}

/* Initialize forward transition tables.

   @param sm - pointer to semi-markov model data

   @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static int crf1de_semimarkov_build_frw_transitions(crf1de_semimarkov_t *sm)
{
  sm->m_forward_trans1 = calloc(sm->num_fs * sm->max_order, sizeof(int));
  if (sm->m_forward_trans1 == NULL)
    return -1;

  sm->m_forward_trans2 = calloc(sm->num_fs * (sm->max_order + 1), sizeof(int));
  if (sm->m_forward_trans2 == NULL) {
    free(sm->m_forward_trans1);
    return -2;
  }

  int i = 0, j = 0, idx = 0, max_seq_len = sm->m_max_order;
  int *seq_start = NULL, *seq_end = NULL;
  RUMAVL_NODE *node = NULL, *;
  memset(sm->m_wrkbench, -1, sizeof(sm->m_wrkbench));
  while ((node = rumavl_node_next(m_forward_states, node, 1, (void**)&sm->m_wrkbench)) != NULL) {
    seq_start = sm->m_wrkbench;
    fprintf(stderr, "sm->forwards_states[%d] = %d", i, *seq_start);
    while(*seq_start != -1) {
      fprintf(stderr, "%d", *seq_start++);
    }
    seq_start = sm->m_wrkbench;
    seq_end = memchr(seq_start, -1, max_seq_len);
    if (! seq_end)
      seg_end = seq_start + max_seg_len + 1;

    for (j = 0; j < sm->L; ++j) {
      *seg_end = j;
      idx = crf1de_semimarkov_find_lng_sfx(sm, seq_start, seq_end);
      rumavl_find();
    }

    memset(sm->m_wrkbench, -1, sizeof(sm->m_wrkbench));
  }
  exit(66);
  return 0;
}

/* Initialize forward transition tables.

   @param sm - pointer to semi-markov model data

   @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static int crf1de_semimarkov_build_bck_transitions(crf1de_semimarkov_t *sm)
{
  return 0;
}

static int crf1de_semimarkov_initialize(crf1de_semimarkov_t *sm, const int a_max_order, \
					const int a_seg_len_lim, const int L)
{
  int j;

  assert(L > 0);			/* can't allow zero length arrays */
  sm->m_max_order = a_max_order;
  sm->m_seg_len_lim = a_seg_len_lim;
  sm->L = L;
  /* initialize dictionaries of forward and backward states */
  sm->m_num_fs = 0;
  sm->m_forward_states = rumavl_new(a_max_order * sizeof(int), crf1de_cmp_lseq, \
			      NULL, NULL);
  sm->m_num_bs = 0;
  sm->m_backward_states = rumavl_new((1 + a_max_order) * sizeof(int), crf1de_cmp_lseq, \
			       NULL, NULL);

  /* allocate memory for storing maximum segment lengths */
  sm->m_max_seg_len = calloc(L, sizeof(int));
  if (sm->m_max_seg_len == NULL)
    return -1;
  /* create a workbench space */
  sm->m_wrkbench = calloc(a_max_order + 2, sizeof(int));
  if (sm->m_wrkbench == NULL) {
    free(sm->m_max_seg_len);
    return -2;
  }
  /* initialize ring for storing labels */
  if (crfsuite_ring_create_instance(&sm->m_ring, a_max_order)) {
    free(sm->m_max_seg_len);
    free(sm->m_wrkbench);
    return -3;
  }
  /* add labels to forward and backward states maps */
  memset(sm->m_wrkbench, -1, sizeof(sm->m_wrkbench));

  for (int i = 0; i < L; ++i) {
    sm->m_wrkbench[0] = i;

    rumavl_insert(sm->m_forward_states, sm->m_wrkbench);
    ++sm->m_num_fs;
    rumavl_insert(sm->m_backward_states, sm->m_wrkbench);
    ++sm->m_num_bs;

    for (j = 0; j < L; ++j) {
      sm->m_wrkbench[1] = j;
      rumavl_insert(sm->m_backward_states, sm->m_wrkbench);
      ++sm->m_num_bs;
    }
    sm->m_wrkbench[1] = -1;
  }
  return 0;
}

static void crf1de_semimarkov_update(crf1de_semimarkov_t *sm, int a_lbl, int a_seg_len)
{
  int j, k;
  /* update maximum segment length if necessary */
  if (sm->m_max_seg_len[a_lbl] < a_seg_len)
    sm->m_max_seg_len[a_lbl] = a_seg_len;

  /* append new label to the label ring */
  sm->m_ring->push(sm->m_ring, a_lbl);
  /* clear workbench */
  memset(sm->m_wrkbench, -1, sizeof(sm->m_wrkbench));
  /* add forward and backward states (a.k.a prefixes and suffixes) */
  crfsuite_chain_link_t *chlink = sm->m_ring->head;
  sm->m_wrkbench[0] = chlink->data;
  chlink = chlink->next;
  /* append prefixes to forwards states */
  for (int i = 1; i < sm->m_ring->num_items; ++i) {
    sm->m_wrkbench[i] = chlink->data;
    if (rumavl_find(sm->m_forward_states, sm->m_wrkbench) == NULL) {
      rumavl_insert(sm->m_forward_states, sm->m_wrkbench);
      ++sm->m_num_fs;
      j = i + 1;
      for (k = 0; k < sm->L; ++k) {
	/* skip equal tags for semi-markov CRFs */
	if (chlink->data == k && sm->m_seg_len_lim < 0)
	  continue;

	sm->m_wrkbench[j] = k;
	rumavl_insert(sm->m_backward_states, sm->m_wrkbench);
	++sm->m_num_bs;
      }
    }
    chlink = chlink->next;
  }
}

static int crf1de_semimarkov_finalize(crf1de_semimarkov_t *sm)
{
  /* generate forward transitions */
  if (crf1de_semimarkov_build_frw_transitions(sm))
    return -1;

  /* generate backward transitions */
  if (crf1de_semimarkov_build_bck_transitions(sm)) {
    CLEAR(sm->m_forward_trans1);
    CLEAR(sm->m_forward_trans2);
    return -2;
  }
  return 0;
}

static void crf1de_semimarkov_clear(crf1de_semimarkov_t *sm)
{
  CLEAR(sm->m_max_seg_len);
  CLEAR(sm->m_wrkbench);
  CLEAR(sm->m_fs_llabels);
  CLEAR(sm->m_forward_trans1);
  CLEAR(sm->m_forward_trans2);
  CLEAR(sm->m_backward_trans);
  CLEAR(sm->m_pattern_trans1);
  CLEAR(sm->m_pattern_trans2);

  if (sm->m_ring) {
    sm->m_ring->free(sm->m_ring);
    free(sm->m_ring);
    sm->m_ring = NULL;
  }

  if (sm->m_forward_states) {
    rumavl_destroy(sm->m_forward_states);
    sm->m_forward_states = NULL;
    sm->m_num_fs = 0;
  }

  if (sm->m_backward_states) {
    rumavl_destroy(sm->m_backward_states);
    sm->m_backward_states = NULL;
    sm->m_num_bs = 0;
  }
}

crf1de_semimarkov_t *crf1de_create_semimarkov(void) {
  crf1de_semimarkov_t *sm = calloc(1, sizeof(crf1de_semimarkov_t));
  if (sm == NULL)
    return sm;

  sm->initialize = crf1de_semimarkov_initialize;
  sm->update = crf1de_semimarkov_update;
  sm->finalize = crf1de_semimarkov_finalize;
  sm->clear = crf1de_semimarkov_clear;

  return sm;
}
