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

/// index at which the id of label sequence is stored
#define F_ID 0
/// index at which the length of label sequence is stored
#define F_LEN 1
/// index of array at which new label sequence starts
#define F_START 2

/// function for freeing an item if it is not null
#define CLEAR(a_item)					\
  if ((a_item)) {					\
    free(a_item);					\
    a_item = NULL;					\
  }

/* Implementation */
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
  /* pointers to whole workbenches */
  const int *entry1 = (const int *) a_entry1, *entry2 = (const int *) a_entry2;
  /* pointers to label sequences */
  const int *lbl_seq1 = &entry1[F_START], *lbl_seq2 = &entry2[F_START];
  /* number of labels in both sequences */
  const int N = entry1[F_LEN];

  if (N > entry2[F_LEN])
    return RUMAVL_ASC;
  else if (N < entry2[F_LEN])
    return RUMAVL_DESC;

  /* consecutively check every tag */
  for (int i = 0; i < N; ++i) {
    if (ret = lbl_seq1[i] - lbl_seq2[i])
      break;
  }
  /* normalize return value */
  if (ret > 0)
    ret = RUMAVL_ASC;
  else if (ret < 0)
    ret = RUMAVL_DESC;

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
  fprintf(stderr, "crf1de_semimarkov_build_frw_transitions started\n");
  sm->m_forward_trans1 = calloc(sm->m_num_fs * sm->m_max_order, sizeof(int));
  if (sm->m_forward_trans1 == NULL)
    return -1;

  sm->m_forward_trans2 = calloc(sm->m_num_fs * (sm->m_max_order + 1), sizeof(int));
  if (sm->m_forward_trans2 == NULL) {
    free(sm->m_forward_trans1);
    return -2;
  }

  /* int i = 0, j = 0, idx = 0, max_seq_len = sm->m_max_order; */
  int i = 0, j = 0, seq_id, seq_len, new_seg_len;
  int *seq_start, *entry = NULL;
  RUMAVL_NODE *node = NULL;
  while ((node = rumavl_node_next(sm->m_forward_states, node, 1, (void**) &entry)) != NULL) {
    seq_id = entry[F_ID];
    seq_len = entry[F_LEN];
    seq_start = &entry[F_START];
    fprintf(stderr, "\nsm->forwards_states[%d] = id = %d; len = %d; seq = %d", i++, seq_id, seq_len, *seq_start);
    for (j = 1; j < seq_len; ++j) {
      fprintf(stderr, "|%d", entry[F_START + j]);
    }
  }

  node = NULL;
  while ((node = rumavl_node_next(sm->m_backward_states, node, 1, (void**) &entry)) != NULL) {
    seq_id = entry[F_ID];
    seq_len = entry[F_LEN];
    seq_start = &entry[F_START];
    fprintf(stderr, "\nsm->backward_states[%d] = id = %d; len = %d; seq = %d", i++, seq_id, seq_len, *seq_start);
    for (j = 1; j < seq_len; ++j) {
      fprintf(stderr, "|%d", entry[F_START + j]);
    }
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

  int j, bs_lbl_idx;

  sm->L = L;
  sm->m_max_order = a_max_order;
  sm->m_seg_len_lim = a_seg_len_lim;
  sm->m_max_fs_size = sizeof(int) * (sm->m_max_order + 2);
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

  /* create workbenches for prefixes and suffixes */
  sm->m_fs_wrkbench = (int *) malloc(sm->m_max_fs_size);
  if (sm->m_fs_wrkbench == NULL) {
    free(sm->m_max_seg_len);
    sm->m_ring->free(sm->m_ring);
    free(sm->m_ring);
    return -3;
  }
  memset(sm->m_fs_wrkbench, -1, sm->m_max_fs_size);

  sm->m_bs_wrkbench = (int *) malloc(sm->m_max_bs_size);
  if (sm->m_bs_wrkbench == NULL) {
    free(sm->m_max_seg_len);
    sm->m_ring->free(sm->m_ring);
    free(sm->m_ring);
    free(sm->m_fs_wrkbench);
    return -4;
  }
  memset(sm->m_bs_wrkbench, -1, sm->m_max_bs_size);

  /* initialize dictionaries of forward and backward states */
  sm->m_num_fs = 0;
  sm->m_forward_states = rumavl_new(sm->m_max_fs_size, crf1de_cmp_lseq, NULL, NULL);

  sm->m_num_bs = 0;
  sm->m_backward_states = rumavl_new(sm->m_max_bs_size, crf1de_cmp_lseq, NULL, NULL);

  /* add labels to forward and backward states maps */
  bs_lbl_idx = F_START + 1;		/* insertion point for the second tag in the sequence */
  sm->m_fs_wrkbench[F_LEN] = 1;	/* 1-st cell holds the length of the prefix */
  for (int i = 0; i < L; ++i) {
    sm->m_fs_wrkbench[F_ID] = sm->m_num_fs++; /* 0-th cell holds the running id of the prefix */
    sm->m_bs_wrkbench[F_ID] = sm->m_num_bs++;
    sm->m_bs_wrkbench[F_LEN] = 1;	/* 1-st cell holds the length of the prefix */
    sm->m_fs_wrkbench[F_START] = sm->m_bs_wrkbench[F_START] = i; /* 2-nd and subsequent cells hold the actual affixes */
    rumavl_insert(sm->m_forward_states, (const void *) sm->m_fs_wrkbench);
    rumavl_insert(sm->m_backward_states, sm->m_bs_wrkbench);

    sm->m_bs_wrkbench[F_LEN] = 2;
    sm->m_bs_wrkbench[bs_lbl_idx] = sm->m_bs_wrkbench[F_START];
    for (j = 0; j < L; ++j) {
      if (j == i && sm->m_seg_len_lim <= 0)
    	continue;

      sm->m_bs_wrkbench[F_ID] = sm->m_num_bs++;
      sm->m_bs_wrkbench[F_START] = j;
      rumavl_insert(sm->m_backward_states, sm->m_bs_wrkbench);
    }
    sm->m_bs_wrkbench[F_START] = -1;
  }
  return 0;
}

static void crf1de_semimarkov_update(crf1de_semimarkov_t *sm, int a_lbl, int a_seg_len)
{
  int j, fs_lbl_idx, bs_lbl_idx, last_label, seg_len;
  /* update maximum segment length if necessary */
  if (sm->m_max_seg_len[a_lbl] < a_seg_len)
    sm->m_max_seg_len[a_lbl] = a_seg_len;

  sm->m_ring->push(sm->m_ring, a_lbl);
  /* reset both workbenches */
  memset(sm->m_fs_wrkbench, -1, sm->m_max_fs_size);
  memset(sm->m_bs_wrkbench, -1, sm->m_max_bs_size);
  /* add forward and backward states (a.k.a prefixes and suffixes) */
  seg_len = 1;
  crfsuite_chain_link_t *chlink = sm->m_ring->tail->prev; /* tail points to the one after last element */
  last_label = chlink->data;
  sm->m_fs_wrkbench[F_START] = sm->m_bs_wrkbench[F_START + 1] = last_label;
  /* append prefixes to forwards states */
  for (int i = 1; i < sm->m_ring->num_items; ++i) {
    ++seg_len;
    fs_lbl_idx = F_START + i;
    bs_lbl_idx = fs_lbl_idx + 1;
    chlink = chlink->prev;
    sm->m_fs_wrkbench[F_LEN] = seg_len;
    sm->m_fs_wrkbench[fs_lbl_idx] = sm->m_bs_wrkbench[bs_lbl_idx] = chlink->data;
    if (rumavl_find(sm->m_forward_states, sm->m_fs_wrkbench) == NULL) {
      sm->m_fs_wrkbench[F_ID] = sm->m_num_fs++;
      rumavl_insert(sm->m_forward_states, sm->m_fs_wrkbench);
      /* increase segment length for backward states */
      sm->m_bs_wrkbench[F_LEN] = seg_len + 1;
      for (j = 0; j < sm->L; ++j) {
      	/* skip equal tags for semi-markov CRFs */
      	if (last_label == j && sm->m_seg_len_lim <= 0)
      	  continue;

      	sm->m_bs_wrkbench[F_ID] = sm->m_num_bs++;
      	sm->m_bs_wrkbench[F_START] = j;
      	rumavl_insert(sm->m_backward_states, sm->m_bs_wrkbench);
      }
    }
  }
}

static int crf1de_semimarkov_finalize(crf1de_semimarkov_t *sm)
{
  /* exit(66); */
  /* generate forward transitions */
  if (crf1de_semimarkov_build_frw_transitions(sm))
    return -1;

  /* generate backward transitions */
  if (crf1de_semimarkov_build_bck_transitions(sm)) {
    CLEAR(sm->m_forward_trans1);
    CLEAR(sm->m_forward_trans2);
    return -2;
  }

  /* clear workbenches */
  CLEAR(sm->m_fs_wrkbench);
  CLEAR(sm->m_bs_wrkbench);

  /* clear ring */
  if (sm->m_ring) {
    sm->m_ring->free(sm->m_ring);
    free(sm->m_ring);
    sm->m_ring = NULL;
  }

  return 0;
}

static void crf1de_semimarkov_clear(crf1de_semimarkov_t *sm)
{
  CLEAR(sm->m_max_seg_len);
  CLEAR(sm->m_fs_wrkbench);
  CLEAR(sm->m_bs_wrkbench);
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
