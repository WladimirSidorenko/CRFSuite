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
#include "crf1d.h"
#include "semimarkov.h"

#include <assert.h>		/* for assert() */
#include <stdio.h>		/* for fprintf() */
#include <stdlib.h>		/* for calloc() */
#include <string.h>		/* for memset() */

/* Macros */

/// macro for accessing entries of 1-st forward transition table
#define FORWARD_TRANS1(sm, y, x)				\
  (&MATRIX(sm->m_forward_trans1, (1 + sm->L), x, y))

/// macro for accessing entries of 2-nd forward transition table
#define FORWARD_TRANS2(sm, y, x)					\
  (&MATRIX(sm->m_forward_trans2, (1 + sm->L), x, y))

/// macro for accessing entries of backward transition table
#define BACKWARD_TRANS(sm, y, x)		\
  (&MATRIX(sm->m_backward_trans, sm->L, x, y))

/// index at which the number of affixes is stored in transitions
#define F_PRFX_N 0
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

/// function for destroying RUMAVL dictionary
#define RUMAVL_CLEAR(a_item)			\
  if (a_item) {					\
    rumavl_destroy(a_item);			\
    a_item = NULL;				\
  }

/* Implementation */

/**
 * Output single state.
 *
 * @param a_fstream - file stream for outputting information
 * @param a_name - symbolic name of state's container
 * @param a_entry - pointer to rumavl entry holding the state
 *
 * @return \c void
 */
static void semimarkov_output_state(FILE *a_fstream, const char *a_name, int *a_entry)
{
  int pk_id = a_entry[F_ID];
  int pk_len = a_entry[F_LEN];
  int *pk_start = &a_entry[F_START];

  if (a_name)
    fprintf(a_fstream, "%s[%d] (%d) = %d", a_name, pk_id, pk_len, *pk_start);
  else
    fprintf(a_fstream, "%d", *pk_start);

  for (int i = 1; i < pk_len; ++i) {
    fprintf(a_fstream, "|%d", *++pk_start);
  }
  if (a_name)
    fprintf(a_fstream, "\n");
}

/**
 * Output states.
 *
 * @param sm - pointer to semi-markov model data
 *
 * @return \c void
 */
static void semimarkov_debug_states(const crf1de_semimarkov_t * const sm)
{
  int i, j, pk_id, pk_len;
  int *pk_entry = NULL, *pk_start;
  RUMAVL_NODE *node = NULL;

  /* output patterns */
  while ((node = rumavl_node_next(sm->m_patterns, node, 1, (void**) &pk_entry)) != NULL) {
    semimarkov_output_state(stderr, "sm->m_patterns", pk_entry);
  }

  /* output forward states */
  node = NULL;
  while ((node = rumavl_node_next(sm->m_forward_states, node, 1, (void**) &pk_entry)) != NULL) {
    semimarkov_output_state(stderr, "sm->m_forward_states", pk_entry);
  }

  /* output backward states */
  node = NULL;
  while ((node = rumavl_node_next(sm->m_backward_states, node, 1, (void**) &pk_entry)) != NULL) {
    semimarkov_output_state(stderr, "sm->m_backward_states", pk_entry);
  }
}

/**
 * Output transition tables.
 *
 * @param sm - pointer to semi-markov model data
 *
 * @return \c void
 */
static void semimarkov_debug_transitions(const crf1de_semimarkov_t * const sm)
{
  int i, j, i_max;
  int pk_id, pky_id;
  const int *pk_entry, *pky_entry;

  /* debug forward transition matrices */
  for (j = 0; j < sm->m_num_fs; ++j) {
    /* obtain number of states for which j is maximum suffix */
    i_max = F_PRFX_N + 1 + *FORWARD_TRANS1(sm, j, F_PRFX_N);
    fprintf(stderr, "sm->forward_trans1[%d (", j);
    semimarkov_output_state(stderr, NULL, sm->m_fsid2fs[j]);
    fprintf(stderr, ")]: ");

    for (i = F_PRFX_N + 1; i < i_max; ++i) {
      pk_id = *FORWARD_TRANS1(sm, j, i);
      fprintf(stderr, "(pk_id = %d) ", pk_id);
      semimarkov_output_state(stderr, NULL, sm->m_fsid2fs[pk_id]);

      if (i_max - i > 1)
	fprintf(stderr, "; ");
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "sm->forward_trans2[%d (", j);
    semimarkov_output_state(stderr, NULL, sm->m_fsid2fs[j]);
    fprintf(stderr, ")]");
    for (i = F_PRFX_N + 1; i < i_max; ++i) {
      pky_id = *FORWARD_TRANS2(sm, j, i);
      fprintf(stderr, "(pky_id = %d) ", pky_id);
      semimarkov_output_state(stderr, NULL, sm->m_bsid2bs[pky_id]);

      if (i_max - i > 1)
  	fprintf(stderr, "; ");
    }
    fprintf(stderr, "\n");
  }

  /* debug backward transition matrices */
  for (j = 0; j < sm->m_num_bs; ++j) {
    for (i = 0; i < sm->L; ++i) {
      fprintf(stderr, "sm->backward_trans1[");
      semimarkov_output_state(stderr, NULL, sm->m_bsid2bs[j]);
      fprintf(stderr, "]");
      fprintf(stderr, "[%d] =", i);

      pk_id = *BACKWARD_TRANS(sm, j, i);
      if (pk_id >= 0)
	semimarkov_output_state(stderr, NULL, sm->m_bsid2bs[pk_id]);
      else
	fprintf(stderr, "-1");

      fprintf(stderr, "\n");
    }
  }
}

/**
 * Function for comparing two label sequences.
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

/**
 * Find maximum known suffix for a given tag sequence.
 *
 * @param a_dic - reference dictionary in which the suffix should be checked
 * @param a_seq - tag sequence for which we should look for suffixes
 *
 * @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static int semimarkov_find_max_sfx(RUMAVL *a_dic, int *a_seq)
{
  int *seq_start = &a_seq[F_START];
  int orig_len = a_seq[F_LEN], ret = -1;
  void *s_id = NULL;

  for (int ilen = orig_len; ilen > 0; --ilen) {
    a_seq[F_LEN] = ilen;
    if (s_id = rumavl_find(a_dic, a_seq))
      break;
  }

  if (s_id)
    ret = ((int *) s_id)[F_ID];

  /* restore original length and return */
  a_seq[F_LEN] = orig_len;
  return ret;
}

/**
 * Initialize forward transition tables.
 *
 * @param sm - pointer to semi-markov model data
 *
 * @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static int semimarkov_build_frw_transitions(crf1de_semimarkov_t *sm)
{
  /* allocate memory for storing last labels and prefix indices */
  sm->m_fs_llabels = calloc(sm->m_num_fs, sizeof(int));
  if (sm->m_fs_llabels == NULL)
    return -1;

  sm->m_forward_trans1 = calloc(sm->m_num_fs * (1 + sm->L), sizeof(int));
  if (sm->m_forward_trans1 == NULL) {
    CLEAR(sm->m_fs_llabels);
    return -2;
  }

  sm->m_forward_trans2 = calloc(sm->m_num_fs * (1 + sm->L), sizeof(int));
  if (sm->m_forward_trans2 == NULL) {
    CLEAR(sm->m_fs_llabels);
    CLEAR(sm->m_forward_trans1);
    return -3;
  }

  sm->m_fsid2fs = calloc(sm->m_num_fs, sizeof(int *));
  if (sm->m_fsid2fs == NULL) {
    CLEAR(sm->m_fs_llabels);
    CLEAR(sm->m_forward_trans1);
    CLEAR(sm->m_forward_trans2);
    return -4;
  }
  /* iterate over each prefix, append a label to it and find the longest
     possible suffix in both forward and backward transitions */
  RUMAVL_NODE *node = NULL;
  int *pk_start, *pk_entry = NULL, *pky_entry = NULL, pk_idx, pky_idx;
  int i = 0, j = 0, idx, pk_id, pky_id, pk_len, pky_len, cnt, last_label;
  memset(sm->m_bs_wrkbench, -1, sm->m_max_bs_size);

  while ((node = rumavl_node_next(sm->m_forward_states, node, 1, (void**) &pk_entry)) != NULL) {
    pk_id = pk_entry[F_ID];
    pk_len = pk_entry[F_LEN];
    pk_start = &pk_entry[F_START];
    sm->m_fsid2fs[pk_id] = pk_entry;

    /* remember last tag of the sequence */
    sm->m_fs_llabels[pk_id] = last_label = *pk_start;
    sm->m_bs_wrkbench[F_LEN] = pk_len + 1;
    memcpy(&sm->m_bs_wrkbench[F_START + 1], pk_start, pk_len * sizeof(int));
    /* append every possible tag and find longest known suffix */
    for (j = 1; j < sm->L; ++j) {
      if (j == last_label && sm->m_seg_len_lim <= 0)
	continue;

      sm->m_bs_wrkbench[F_START] = j;
      pky_entry = (int *) rumavl_find(sm->m_backward_states, sm->m_bs_wrkbench);
      pky_id = pky_entry[F_ID];
      idx = semimarkov_find_max_sfx(sm->m_forward_states, sm->m_bs_wrkbench);

      /* increase the counter of prefixes for which `pk_id` and `pky_id` were
	 the longest suffixes and insert the id's of those suffixes */
      cnt = ++(*FORWARD_TRANS1(sm, idx, F_PRFX_N));
      /* check that we don't overflow */
      assert(cnt <= sm->L);
      cnt += F_PRFX_N;
      *FORWARD_TRANS1(sm, idx, cnt) = pk_id;
      *FORWARD_TRANS2(sm, idx, cnt) = pky_id;
    }
  }
  return 0;
}

/**
 * Initialize backward transition tables.
 *
 * @param sm - pointer to semi-markov model data
 *
 * @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static int semimarkov_build_bkw_transitions(crf1de_semimarkov_t *sm)
{
  fprintf(stderr, "generating backward states\n");
  sm->m_bsid2bs = calloc(sm->m_num_bs, sm->L *sizeof(int *));
  if (sm->m_bsid2bs == NULL)
    return -1;

  sm->m_backward_trans = calloc(sm->m_num_bs, sm->L * sizeof(int));
  if (sm->m_backward_trans == NULL) {
    CLEAR(sm->m_bsid2bs);
    return -2;
  }

  RUMAVL_NODE *node = NULL;
  int pk_id, pk_len, lst_lbl, *pk_start, *pk_entry = NULL;
  int y, sfx_id, pky_len, *pky_start, *pky_entry = sm->m_bs_wrkbench;

  fprintf(stderr, "entering while loop\n");
  while ((node = rumavl_node_next(sm->m_backward_states, node, 1, (void**) &pk_entry)) != NULL) {
    pk_id = pk_entry[F_ID];
    pk_len = pk_entry[F_LEN];
    pk_start = &pk_entry[F_START];

    lst_lbl = *pk_start;
    sm->m_bsid2bs[pk_id] = pk_entry;

    pky_entry[F_LEN] = pk_len + 1;
    memcpy(&pky_entry[F_START + 1], pk_start, pk_len * sizeof(int));

    for (y = 0; y < sm->L; ++y) {
      if (y == lst_lbl && sm->m_seg_len_lim <= 0) {
	*BACKWARD_TRANS(sm, pk_id, y) = -1;
      } else {
	pky_entry[F_START] = y;
	sfx_id = semimarkov_find_max_sfx(sm->m_backward_states, pky_entry);
	*BACKWARD_TRANS(sm, pk_id, y) = sfx_id;
      }
    }
  }
  return 0;
}

/**
 * Generate backward states for given prefix.
 *
 * @param sm - pointer to semi-markov model
 * @param a_wrkbench - pointer to workbench with prefix
 *
 * @return \c void
 */
static void semimarkov_add_bkw_states(crf1de_semimarkov_t *sm, int *a_wrkbench)
{
  semimarkov_output_state(stderr, "semimarkov_add_bkw_states: input state", a_wrkbench);

  int last_label = a_wrkbench[F_START + a_wrkbench[F_LEN] - 1];
  size_t len = a_wrkbench[F_LEN]++;

  /* append all possible tags to given prefix */
  for (int i = 0; i < sm->L; ++i) {
    if (i == last_label && sm->m_seg_len_lim < 0)
      continue;

    fprintf(stderr, "i = %d\n", i);
    /* we assume that backward state is not known */
    a_wrkbench[F_ID] = sm->m_num_bs++;
    a_wrkbench[F_START + len] = i;
    rumavl_insert(sm->m_backward_states, a_wrkbench);
  }
  /* restore original prefix */
  a_wrkbench[F_START + len] = -1;
  --a_wrkbench[F_LEN];
}

/**
 * Generate all possible patterns and add them to the semimarkov model.
 *
 * @param sm - pointer to semi-markov model
 * @param a_wrkbench - pointer to workbench with pattern
 *
 * @return \c void
 */
static void semimarkov_add_patterns(crf1de_semimarkov_t *sm, int *a_wrkbench)
{
  while (a_wrkbench[F_LEN] > 1 && rumavl_find(sm->m_patterns, a_wrkbench) == NULL) {
    a_wrkbench[F_ID] = sm->m_num_patterns++;
    rumavl_insert(sm->m_patterns, a_wrkbench);

    /* reduce pattern for the next loop */
    --a_wrkbench[F_LEN];
    memmove(a_wrkbench[F_START], a_wrkbench[F_START + 1], a_wrkbench[F_LEN] * sizeof(int));
  }
}

/**
 * Generate forward and backward states and store them in the
 * semimarkov model.
 *
 * @param sm - pointer to semi-markov model
 * @param a_wrkbench - pointer to workbench with pattern
 *
 * @return \c void
 */
static void semimarkov_add_states(crf1de_semimarkov_t *sm, int *a_wrkbench)
{

  size_t len = a_wrkbench[F_LEN];

  while (--len > 0) {
    /* erase last element */
    a_wrkbench[F_LEN] = len;
    a_wrkbench[F_START + len] = -1;

    /* ASSUMPTION: if we already know forward state, then we already
       know all backward states and patterns */
    if (rumavl_find(sm->m_forward_states, a_wrkbench) != NULL)
      break;

    sm->m_bs_wrkbench[F_ID] = sm->m_num_fs++;
    rumavl_insert(sm->m_forward_states, a_wrkbench);
    semimarkov_add_bkw_states(sm, a_wrkbench);
  }
}

/**
 * Add known labels to the sets of patterns, forward and backward transitions.
 *
 * @param sm - pointer to semi-markov model data
 *
 * @return \c void
 */
static void semimarkov_initialize_states(crf1de_semimarkov_t *sm)
{
  int j;
  const int lbl_idx = F_START + 1;
  memset(sm->m_bs_wrkbench, -1, sm->m_max_bs_size);

  for (int i = 0; i < sm->L; ++i) {
    sm->m_bs_wrkbench[F_LEN] = 1;
    sm->m_bs_wrkbench[F_START] = i;

    /* add label to the set of prefixes (forward states) */
    sm->m_bs_wrkbench[F_ID] = sm->m_num_fs++;
    rumavl_insert(sm->m_forward_states, sm->m_bs_wrkbench);

    /* add label to the set of backward states and patterns */
    sm->m_bs_wrkbench[F_ID] = sm->m_num_patterns++;
    rumavl_insert(sm->m_patterns, sm->m_bs_wrkbench);

    /* generate first-order transitions */
    sm->m_bs_wrkbench[F_ID] = sm->m_num_bs++;
    rumavl_insert(sm->m_backward_states, sm->m_bs_wrkbench);
    semimarkov_add_bkw_states(sm, sm->m_bs_wrkbench);
  }
}

/**
 * Allocate space and set initial values for semi-markov container.
 *
 * @param sm - pointer to semi-markov model data
 * @param a_max_order - maximum order of transition features
 * @param a_seg_len_lim - imposed limit on the maximum length of a contiguous
 *                        tag segment
 * @param L - number of tags
 *
 * @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static int semimarkov_initialize(crf1de_semimarkov_t *sm, const int a_max_order, \
					const int a_seg_len_lim, const int L)
{
  if (L <= 0)
    return -1;

  sm->L = L;
  sm->m_max_order = a_max_order;
  sm->m_seg_len_lim = a_seg_len_lim;
  /* a prefix (aka forward state can have a maximum length of max_order, plus
     we need two additional ints for storing the id of that prefix and its
     length) */
  sm->m_max_fs_size = sizeof(int) * (sm->m_max_order + 2);
  /* backward states can have one more tag than forward states */
  sm->m_max_bs_size = sizeof(int) + sm->m_max_fs_size;

  /* allocate memory for storing maximum segment lengths */
  sm->m_max_seg_len = calloc(L, sizeof(int));
  if (sm->m_max_seg_len == NULL)
    return -1;

  /* initialize ring for storing labels */
  if (crfsuite_ring_create_instance(&sm->m_ring, sm->m_max_order + 1)) {
    free(sm->m_max_seg_len);
    return -2;
  }

  /* create workbench for states and transitions */
  sm->m_bs_wrkbench = (int *) malloc(sm->m_max_bs_size + sizeof(int));
  if (sm->m_bs_wrkbench == NULL) {
    free(sm->m_max_seg_len);
    sm->m_ring->free(sm->m_ring);
    free(sm->m_ring);
    return -4;
  }
  memset(sm->m_bs_wrkbench, -1, sm->m_max_bs_size);

  /* initialize dictionaries of patterns, forward and backward states */
  sm->m_num_patterns = 0;
  sm->m_patterns = rumavl_new(sm->m_max_bs_size, crf1de_cmp_lseq, NULL, NULL);

  sm->m_num_fs = 0;
  sm->m_forward_states = rumavl_new(sm->m_max_fs_size, crf1de_cmp_lseq, NULL, NULL);

  sm->m_num_bs = 0;
  sm->m_backward_states = rumavl_new(sm->m_max_bs_size, crf1de_cmp_lseq, NULL, NULL);

  /* unconditionally add all labels to the sets of patterns, forward and
     backward states */
  semimarkov_initialize_states(sm);

  return 0;
}


/**
 * Add new patterns to the semimarkov model if needed.
 *
 * @param sm - pointer to semi-markov model data
 * @param a_lbl - new label to be added
 * @param a_seg_len - number of words continuously tagged with `a_lbl`
 *
 * @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static void semimarkov_update(crf1de_semimarkov_t *sm, int a_lbl, int a_seg_len)
{
  /* update maximum segment length if necessary */
  sm->m_ring->push(sm->m_ring, a_lbl);

  if (sm->m_max_seg_len[a_lbl] < a_seg_len)
    sm->m_max_seg_len[a_lbl] = a_seg_len;

  /* since we have already added all labels during initialization, we
     can skip segments of length one */
  if (sm->m_ring->num_items == 1)
    return;

  /* add pattern to the set */

  /* transfer tag sequence from ring to workbench */
  crfsuite_chain_link_t *chlink = sm->m_ring->head;
  memset(sm->m_bs_wrkbench, -1, sm->m_max_bs_size);
  sm->m_bs_wrkbench[F_LEN] = 0;

  for (size_t i = 0; i < sm->m_ring->num_items; ++i) {
    ++sm->m_bs_wrkbench[F_LEN];
    sm->m_bs_wrkbench[F_START + i] = chlink->data;
    chlink = chlink->next;
  }
  /* generate all possible prefixes and multiply them by L */
  if (rumavl_find(sm->m_patterns, sm->m_bs_wrkbench) == NULL) {
    semimarkov_add_states(sm, sm->m_bs_wrkbench);
    semimarkov_add_patterns(sm, sm->m_bs_wrkbench);
  }
}

/**
 * Add new affixes to the semimarkov container if needed.
 *
 * @param sm - pointer to semi-markov model data
 *
 * @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static int semimarkov_finalize(crf1de_semimarkov_t *sm)
{
  /* debug forward and backward transitions */
  semimarkov_debug_states(sm);
  exit(66);

  /* generate forward transitions */
  if (semimarkov_build_frw_transitions(sm))
    return -1;

  /* generate backward transitions */
  if (semimarkov_build_bkw_transitions(sm)) {
    CLEAR(sm->m_forward_trans1);
    CLEAR(sm->m_forward_trans2);
    CLEAR(sm->m_fsid2fs);
    return -2;
  }

  /* generate pattern transitions */
  /* semimarkov_build_ptrn_transitions(sm); */

  semimarkov_debug_transitions(sm);
  exit(66);
  /* clear workbenches */
  CLEAR(sm->m_bs_wrkbench);

  /* clear ring */
  if (sm->m_ring) {
    sm->m_ring->free(sm->m_ring);
    free(sm->m_ring);
    sm->m_ring = NULL;
  }
  return 0;
}

/**
 * Deallocate space and reset values.
 *
 * @param sm - pointer to semi-markov model data
 *
 * @return void
 */
static void semimarkov_clear(crf1de_semimarkov_t *sm)
{
  CLEAR(sm->m_max_seg_len);
  CLEAR(sm->m_bs_wrkbench);
  CLEAR(sm->m_fs_llabels);
  CLEAR(sm->m_forward_trans1);
  CLEAR(sm->m_forward_trans2);
  CLEAR(sm->m_fsid2fs);
  CLEAR(sm->m_backward_trans);
  CLEAR(sm->m_pattern_trans1);
  CLEAR(sm->m_pattern_trans2);
  CLEAR(sm->m_bsid2bs);

  if (sm->m_ring) {
    sm->m_ring->free(sm->m_ring);
    free(sm->m_ring);
    sm->m_ring = NULL;
  }

  RUMAVL_CLEAR(sm->m_forward_states);
  sm->m_num_fs = 0;

  RUMAVL_CLEAR(sm->m_backward_states);
  sm->m_num_bs = 0;

  RUMAVL_CLEAR(sm->m_patterns);
  sm->m_num_patterns = 0;
}

crf1de_semimarkov_t *crf1de_create_semimarkov(void) {
  crf1de_semimarkov_t *sm = calloc(1, sizeof(crf1de_semimarkov_t));
  if (sm == NULL)
    return sm;

  sm->initialize = semimarkov_initialize;
  sm->update = semimarkov_update;
  sm->finalize = semimarkov_finalize;
  sm->clear = semimarkov_clear;

  return sm;
}

/* reset both workbenches */
/* memset(sm->m_fs_wrkbench, -1, sm->m_max_fs_size); */
/* memset(sm->m_bs_wrkbench, -1, sm->m_max_bs_size); */

/* add forward and backward states (a.k.a prefixes and suffixes) */
/* seg_len = 1; */
/* chlink = sm->m_ring->tail->prev; /\* tail points to the one after last element *\/ */
/* last_label = chlink->data; */
/* sm->m_fs_wrkbench[F_START] = sm->m_bs_wrkbench[F_START + 1] = last_label; */


/* /\* append prefixes to forwards states *\/ */
/* for (int i = 1; i < sm->m_ring->num_items; ++i) { */
/*   ++seg_len; */
/*   fs_lbl_idx = F_START + i; */
/*   bs_lbl_idx = fs_lbl_idx + 1; */
/*   chlink = chlink->prev; */
/*   sm->m_fs_wrkbench[F_LEN] = seg_len; */
/*   sm->m_fs_wrkbench[fs_lbl_idx] = sm->m_bs_wrkbench[bs_lbl_idx] = chlink->data; */
/*   if (rumavl_find(sm->m_forward_states, sm->m_fs_wrkbench) == NULL) { */
/*     sm->m_fs_wrkbench[F_ID] = sm->m_num_fs++; */
/*     rumavl_insert(sm->m_forward_states, sm->m_fs_wrkbench); */
/*     /\* increase segment length for backward states *\/ */
/*     sm->m_bs_wrkbench[F_LEN] = seg_len + 1; */
/*     for (j = 0; j < sm->L; ++j) { */
/*       /\* skip equal tags for semi-markov CRFs *\/ */
/*       if (last_label == j && sm->m_seg_len_lim <= 0) */
/* 	continue; */

/*       sm->m_bs_wrkbench[F_ID] = sm->m_num_bs++; */
/*       sm->m_bs_wrkbench[F_START] = j; */
/*       rumavl_insert(sm->m_backward_states, sm->m_bs_wrkbench); */
/*     } */
/*   } */
/*  } */
