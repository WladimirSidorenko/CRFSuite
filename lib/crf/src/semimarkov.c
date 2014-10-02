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
#include <stdio.h>		/* for fprintf() */
#include <stdlib.h>		/* for calloc() */
#include <string.h>		/* for memset() */

/* Macros */
#define CLEAR(a_item)					\
  if ((a_item)) {					\
    free(a_item);					\
    a_item = NULL;					\
  }

/* Declaration */

/**
 * Synonym for crf1de_wrkbench struct.
 */
typedef struct crf1de_wrkbench crf1de_wrkbench_t;

/**
 * Auxiliary container for constructing label prefixes.
 */
struct crf1de_wrkbench {
  int *m_lbl;
  int m_id;
  size_t m_size;

  void (*free)(crf1de_wrkbench_t *a_wb);
  void (*reset)(crf1de_wrkbench_t *a_wb);
};


/* Implementation */

/* crf1de_wrkbench */
static void crf1de_wrkbench_free(crf1de_wrkbench_t *a_wb)
{
  if (a_wb->m_lbl) {
    free(a_wb->m_lbl);

    a_wb->m_id = -1;
    a_wb->m_size = 0;
    a_wb->m_lbl = NULL;
    a_wb->free = NULL;
    a_wb->reset = NULL;
  }
}

static void crf1de_wrkbench_reset(crf1de_wrkbench_t *a_wb)
{
  memset(a_wb->m_lbl, -1, a_wb->m_size);
  a_wb->m_id = -1;
}

static void *crf1de_init_wrkbench(const size_t a_size)
{
  crf1de_wrkbench_t *iwrkbench = (crf1de_wrkbench_t *) calloc(1, sizeof(crf1de_wrkbench_t));
  if (! iwrkbench)
    return iwrkbench;

  iwrkbench->m_lbl = (int *) malloc(a_size);
  if (! iwrkbench->m_lbl) {
    free(iwrkbench);
    return NULL;
  }

  iwrkbench->m_id = -1;
  iwrkbench->m_size = a_size;

  iwrkbench->free = crf1de_wrkbench_free;
  iwrkbench->reset = crf1de_wrkbench_reset;

  return (void *) iwrkbench;
}

/* crf1de_semimarkov */

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
static int crf1de_cmp_lseq(const void *a_entry1, const void *a_entry2,	\
			   size_t a_size, void *a_udata)
{
  int ret = 0;
  size_t n = a_size / sizeof(int);
  fprintf(stderr, "crf1de_cmp_lseq: a_size = %d; n = %d\n", a_size, n);
  const int *el1 = (const int*) ((const crf1de_wrkbench_t *) a_entry1)->m_lbl;
  const int *el2 = (const int*) ((const crf1de_wrkbench_t *) a_entry2)->m_lbl;

  for (size_t i = 0; i < n; ++i) {
    ret = el1[i] - el2[i];
    fprintf(stderr, "crf1de_cmp_lseq: el1[%d] = %d\n", i, el1[i]);
    fprintf(stderr, "crf1de_cmp_lseq: el2[%d] = %d\n", i, el2[i]);
    /* -1 terminates tagging sequence */
    if (ret || el1[i] < 0)
      break;
  }
  fprintf(stderr, "crf1de_cmp_lseq: ret = %d\n", ret);
  return ret;
}

int crf1de_semimarkov_find_lng_sfx(crf1de_semimarkov_t *sm, const int *a_seqstart, \
				   const int *a_seqend)
{
  return 0;
}

/* Initialize forward transition tables.

   @param sm - pointer to semi-markov model data

   @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static int crf1de_semimarkov_build_frw_transitions(crf1de_semimarkov_t *sm)
{
  sm->m_forward_trans1 = calloc(sm->m_num_fs * sm->m_max_order, sizeof(int));
  if (sm->m_forward_trans1 == NULL)
    return -1;

  sm->m_forward_trans2 = calloc(sm->m_num_fs * (sm->m_max_order + 1), sizeof(int));
  if (sm->m_forward_trans2 == NULL) {
    free(sm->m_forward_trans1);
    return -2;
  }

  int i = 0, j = 0, idx = 0, max_seq_len = sm->m_max_order;
  int *seq_start = NULL, *seq_end = NULL;
  RUMAVL_NODE *node = NULL;
  memset(sm->m_wrkbench, -1, (sm->m_max_order + 2) * sizeof(int));
  while ((node = rumavl_node_next(sm->m_forward_states, node, 1, (void**)&sm->m_wrkbench)) != NULL) {
    seq_start = sm->m_wrkbench;
    fprintf(stderr, "\nsm->forwards_states[%d] = %d|", i, *seq_start);
    while(*seq_start >= 0) {
      fprintf(stderr, "%d|", *seq_start++);
    }
    seq_start = sm->m_wrkbench;
    seq_end = memchr(seq_start, -1, max_seq_len);
    if (! seq_end)
      seq_end = seq_start + sm->m_max_order + 1;

    for (j = 0; j < sm->L; ++j) {
      *seq_end = j;
      idx = crf1de_semimarkov_find_lng_sfx(sm, seq_start, seq_end);
    }

    memset(sm->m_wrkbench, -1, (sm->m_max_order + 2) * sizeof(int));
    ++i;
  }
  fprintf(stderr, "\n");
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
  if (L <= 0)
    return 0;

  int j;
  size_t max_afx_size = (a_max_order + 2) * sizeof(int);

  sm->L = L;
  sm->m_max_order = a_max_order;
  sm->m_seg_len_lim = a_seg_len_lim;
  sm->m_max_fs_size = sizeof(int) + sm->m_max_order;
  sm->m_max_bs_size = sm->m_max_fs_size + sizeof(int);

  /* allocate memory for storing maximum segment lengths */
  sm->m_max_seg_len = calloc(L, sizeof(int));
  if (sm->m_max_seg_len == NULL)
    return -1;

  /* initialize ring for storing labels */
  if (crfsuite_ring_create_instance(&sm->m_ring, sm->m_max_order)) {
    free(sm->m_max_seg_len);
    return -2;
  }

  /* create workbench space */
  sm->m_wrkbench = crf1de_init_wrkbench(max_afx_size);
  if (sm->m_wrkbench == NULL) {
    free(sm->m_max_seg_len);
    sm->m_ring->free(sm->m_ring);
    free(sm->m_ring);
    return -3;
  }
  crf1de_wrkbench_t *iwrkbench = (crf1de_wrkbench_t *) sm->m_wrkbench;
  iwrkbench->reset(iwrkbench);

  /* initialize dictionaries of forward and backward states */
  sm->m_num_fs = 0;
  sm->m_forward_states = rumavl_new(sizeof(*iwrkbench), crf1de_cmp_lseq, NULL, NULL);
  sm->m_num_bs = 0;
  sm->m_backward_states = rumavl_new(sizeof(*iwrkbench), crf1de_cmp_lseq, NULL, NULL);

  /* add labels to forward and backward states maps */
  for (int i = 0; i < L; ++i) {
    iwrkbench->m_lbl[0] = i;
    fprintf(stderr, "crf1de_semimarkov_initialize: sm->m_wrkbench->m_lbl[0] = %d\n", i);

    iwrkbench->m_id = sm->m_num_fs++;
    rumavl_insert(sm->m_forward_states, iwrkbench);
    iwrkbench->m_id = sm->m_num_bs++;
    rumavl_insert(sm->m_backward_states, iwrkbench);

    for (j = 0; j < L; ++j) {
      if (j == i && sm->m_seg_len_lim <= 0)
	continue;

      iwrkbench->m_id = sm->m_num_bs++;
      iwrkbench->m_lbl[1] = j;
      fprintf(stderr, "crf1de_semimarkov_initialize: sm->m_wrkbench[1] = %d\n", j);
      rumavl_insert(sm->m_backward_states, iwrkbench);
    }
    iwrkbench->reset(iwrkbench);
  }
  return 0;
}

static void crf1de_semimarkov_update(crf1de_semimarkov_t *sm, int a_lbl, int a_seg_len)
{
  int j, k;
  crf1de_wrkbench_t *iwrkbench = (crf1de_wrkbench_t *) sm->m_wrkbench;
  /* update maximum segment length if necessary */
  if (sm->m_max_seg_len[a_lbl] < a_seg_len)
    sm->m_max_seg_len[a_lbl] = a_seg_len;

  fprintf(stderr, "crf1de_semimarkov_update: new label = %d\n", a_lbl);
  /* append new label to the label ring */
  sm->m_ring->push(sm->m_ring, a_lbl);
  fprintf(stderr, "crf1de_semimarkov_update: sm->m_ring->num_items = %d\n", sm->m_ring->num_items);
  /* clear workbench */
  iwrkbench->reset(iwrkbench);
  /* add forward and backward states (a.k.a prefixes and suffixes) */
  crfsuite_chain_link_t *chlink = sm->m_ring->head;
  iwrkbench->m_lbl[0] = chlink->data;
  chlink = chlink->next;
  /* append prefixes to forwards states */
  for (int i = 1; i < sm->m_ring->num_items; ++i) {
    iwrkbench->m_lbl[i] = chlink->data;
    fprintf(stderr, "crf1de_semimarkov_update: sm->m_wrkbench[%d] = %d\n", i, iwrkbench[i]);
    if (rumavl_find(sm->m_forward_states, iwrkbench) == NULL) {
      iwrkbench->m_id = sm->m_num_fs++;
      rumavl_insert(sm->m_forward_states, iwrkbench);
      j = i + 1;
      for (k = 0; k < sm->L; ++k) {
	/* skip equal tags for semi-markov CRFs */
	if (chlink->data == k && sm->m_seg_len_lim <= 0)
	  continue;

	iwrkbench->m_id = sm->m_num_bs++;
	iwrkbench->m_lbl[j] = k;
	rumavl_insert(sm->m_backward_states, iwrkbench);
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

  /* free workbench and ring */
  if (sm->m_ring) {
    sm->m_ring->free(sm->m_ring);
    free(sm->m_ring);
    sm->m_ring = NULL;
  }
  crf1de_wrkbench_t *iwrkbench = (crf1de_wrkbench_t *) sm->m_wrkbench;
  if (iwrkbench) {
    iwrkbench->free(iwrkbench);
    free(iwrkbench);
    sm->m_wrkbench = NULL;
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

  crfsuite_ring_t *iring = (crfsuite_ring_t *) sm->m_ring;
  if (sm->m_ring) {
    sm->m_ring->free(sm->m_ring);
    free(sm->m_ring);
    sm->m_ring = NULL;
  }

  crf1de_wrkbench_t *iwrkbench = (crf1de_wrkbench_t *) sm->m_wrkbench;
  if (iwrkbench) {
    iwrkbench->free(iwrkbench);
    free(iwrkbench);
    sm->m_wrkbench = NULL;
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
