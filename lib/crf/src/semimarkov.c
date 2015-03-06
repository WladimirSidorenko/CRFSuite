/*
 * Data structures and auxiliary functions for semi-markov CRFs.
 *
 * Copyright (c) 2014-2015, ANONYMOUS
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

/// function for freeing an item if it is not null
#define CLEAR(a_item)				\
  if ((a_item)) {				\
    free((a_item));				\
    a_item = NULL;				\
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
static void semimarkov_output_state(FILE *a_fstream, const char *a_name, \
				    const crf1de_state_t *a_entry)
{
  int pk_id = a_entry->m_id;
  size_t pk_len = a_entry->m_len;

  if (a_name)
    fprintf(a_fstream, "%s[%d] (%zu) = ", a_name, pk_id, pk_len);

  if (pk_len)
    fprintf(a_fstream, "%d", a_entry->m_seq[0]);

  for (size_t i = 1; i < pk_len; ++i) {
    fprintf(a_fstream, "|%d", a_entry->m_seq[i]);
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
  size_t i = 0;

  for (i = 0; i < sm->L; ++i) {
    fprintf(stderr, "sm->m_max_seg_len[%zu] = %d\n", i, sm->m_max_seg_len[i]);
  }

  /* output patterns */
  for (i = 0; i < sm->m_num_ptrns; ++i) {
    semimarkov_output_state(stderr, "sm->m_ptrns =", &sm->m_ptrns[i]);
  }
  fprintf(stderr, "******************************************************************\n");

  /* output forward states */
  for (i = 0; i < sm->m_num_frw; ++i) {
    semimarkov_output_state(stderr, "sm->m_frw_states", &sm->m_frw_states[i]);
    fprintf(stderr, "sm->m_frw_llabel[%zu] = %d\n", i, sm->m_frw_llabels[i]);
  }
  fprintf(stderr, "******************************************************************\n");

  /* output backward states */
  for (i = 0; i < sm->m_num_bkw; ++i) {
    semimarkov_output_state(stderr, "sm->m_bkw_states", &sm->m_bkw_states[i]);
  }
  fprintf(stderr, "******************************************************************\n");
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
  size_t j = 0;
  const crf1de_state_t *pk_entry = NULL, *pky_entry = NULL;

  /* output forward transition for states */
  for (size_t i = 0; i < sm->m_num_frw; ++i) {
    /* obtain pointer to forward state */
    pk_entry = &sm->m_frw_states[i];
    fprintf(stderr, "forward_transition1[(id = %zu) ", i);
    semimarkov_output_state(stderr, NULL, pk_entry);
    fprintf(stderr, "] =");

    fprintf(stderr, " (%zu):", pk_entry->m_num_affixes);
    for (j = 0; j < pk_entry->m_num_affixes; ++j) {
      fprintf(stderr, " ");
      semimarkov_output_state(stderr, NULL, &sm->m_frw_states[pk_entry->m_trans1[j]]);
      fprintf(stderr, ";");
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "forward_transition2[(id = %zu) ", i);
    semimarkov_output_state(stderr, NULL, pk_entry);
    fprintf(stderr, "] =");

    fprintf(stderr, " (%zu):", pk_entry->m_num_affixes);
    for (j = 0; j < pk_entry->m_num_affixes; ++j) {
      fprintf(stderr, " ");
      semimarkov_output_state(stderr, NULL, &sm->m_bkw_states[pk_entry->m_trans2[j]]);
      fprintf(stderr, ";");
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "******************************************************************\n");

  /* output backward transition for states */
  for (size_t i = 0; i < sm->m_num_bkw; ++i) {
    pky_entry = &sm->m_bkw_states[i];
    for (j = 0; j < sm->L; ++j) {
      fprintf(stderr, "backwardTransition[");
      semimarkov_output_state(stderr, NULL, pky_entry);
      fprintf(stderr, "][%zu] = ", j);

      if (pky_entry->m_bkw_trans[j] >= 0) {
      	semimarkov_output_state(stderr, NULL, &sm->m_bkw_states[pky_entry->m_bkw_trans[j]]);
      }
      fprintf(stderr, "\n");
    }
  }
  fprintf(stderr, "******************************************************************\n");

  /* output suffixes */
  size_t i;
  int ptrn_id = -1, *suffixes = NULL;
  for (size_t pky_id = 0; pky_id < sm->m_num_bkw; ++pky_id) {
    pky_entry = &sm->m_bkw_states[pky_id];
    fprintf(stderr, "suffixes[");
    semimarkov_output_state(stderr, NULL, pky_entry);
    fprintf(stderr, "] =");

    suffixes = &SUFFIXES(sm, pky_id, 0);
    for (i = 0; (ptrn_id = suffixes[i]) >= 0; ++i) {
      fprintf(stderr, " ");
      semimarkov_output_state(stderr, NULL, &sm->m_ptrns[ptrn_id]);
      fprintf(stderr, ";");
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "******************************************************************\n");

  /* output pattern transitions */
  crf1de_state_t *ptrn_entry = NULL;
  for (i = 0; i < sm->m_num_ptrns; ++i) {
    ptrn_entry = &sm->m_ptrns[i];

    /* output pattern transition 1 */
    fprintf(stderr, "patternTransition1[");
    semimarkov_output_state(stderr, NULL, ptrn_entry);
    fprintf(stderr, "] =");
    for (j = 0; j < ptrn_entry->m_num_affixes; ++j) {
      fprintf(stderr, " ");
      semimarkov_output_state(stderr, NULL, &sm->m_frw_states[ptrn_entry->m_trans1[j]]);
      fprintf(stderr, ";");
    }
    fprintf(stderr, "\n");

    /* output pattern transition 2 */
    fprintf(stderr, "patternTransition2[");
    semimarkov_output_state(stderr, NULL, ptrn_entry);
    fprintf(stderr, "] =");
    for (j = 0; j < ptrn_entry->m_num_affixes; ++j) {
      fprintf(stderr, " ");
      semimarkov_output_state(stderr, NULL, &sm->m_bkw_states[ptrn_entry->m_trans2[j]]);
      fprintf(stderr, ";");
    }
    fprintf(stderr, "\n");
  }
}

/**
 * Auxiliary function for deleting label sequences.
 *
 * @param a_tree - pointer to RUMAVL tree
 * @param a_state - state to be inserted
 *
 * @return \c void
 */
static void insert_state(RUMAVL *a_tree, crf1de_state_t *a_state)
{
  int *seq = calloc(a_state->m_len, sizeof(int));
  memcpy((void *) seq, (const void *) a_state->m_seq, a_state->m_len * sizeof(int));

  int *tmp = a_state->m_seq;
  a_state->m_seq = seq;
  rumavl_insert(a_tree, a_state);
  a_state->m_seq = tmp;
}

/**
 * Auxiliary function for copying states.
 *
 * @param a_trg_state - state to which we should copy
 * @param a_src_state - state which should be copied
 *
 * @return \c 0 on success, non-\c 0 otherwise
 */
static int copy_state(crf1de_state_t *a_trg_state, crf1de_state_t *a_src_state)
{
  int *seq = calloc(a_src_state->m_len, sizeof(int));
  if (seq == NULL)
    return -1;

  memcpy((void *) seq, (const void *) a_src_state->m_seq, a_src_state->m_len * sizeof(int));
  memcpy((void *) a_trg_state, (const void *) a_src_state, sizeof(crf1de_state_t));
  a_trg_state->m_seq = seq;
  return 0;
}

/**
 * Auxiliary function for deleting label sequences.
 *
 * @param tree - pointer to RUMAVL tree
 * @param n - pointer to RUMAVL node
 * @param _x - pointer to RUMAVL record
 * @param _y - (not used)
 * @param udata - (not used)
 *
 * @return 0
 */
static int owcb(RUMAVL *tree, RUMAVL_NODE *n, void *_x, const void *_y, void *udata)
{
  crf1de_state_t* x = (crf1de_state_t*)_x;
  free(x->m_seq);
  return 0;
}

/**
 * Auxiliary function for deleting label sequences.
 *
 * @param tree - pointer to RUMAVL tree
 * @param n - pointer to RUMAVL node
 * @param _record - pointer to RUMAVL record
 * @param udata - (not used)
 *
 * @return 0
 */
static int delcb(RUMAVL *tree, RUMAVL_NODE *n, void *_record, void *udata)
{
  crf1de_state_t* record = (crf1de_state_t*)_record;
  free(record->m_seq);
  return 0;
}

/**
 * Function for comparing two label sequences.
 *
 * @param a_lseq1 - first label sequence to compare
 * @param a_lseq2 - second label sequence to compare
 * @param a_size  - record size
 * @param a_udata - unused parameter (needed for compliance)
 *
 * @return > \c 1 if `a_lseq1` is greater than `a_lseq2`, \c 0 if both
 * sequences are the same, and < 0 if `a_lseq1` is smaller than
 * `a_lseq2`
 */
static int sm_cmp_lseq(const void *a_entry1, const void *a_entry2,	\
			   size_t a_size, void *a_udata)
{
  int ret = 0;
  /* pointers to whole workbenches */
  const crf1de_state_t *entry1 = (const crf1de_state_t *) a_entry1;
  const crf1de_state_t *entry2 = (const crf1de_state_t *) a_entry2;

  /* number of labels in both sequences */
  const int n = entry1->m_len;

  if (n > entry2->m_len)
    return RUMAVL_ASC;
  else if (n < entry2->m_len)
    return RUMAVL_DESC;

  /* consecutively check every tag */
  for (int i = 0; i < n; ++i) {
    if ((ret = entry1->m_seq[i] - entry2->m_seq[i]))
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
 * @param a_seq - tag sequence for which we should find the longest suffix
 *
 * @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static crf1de_state_t *semimarkov_find_max_sfx(RUMAVL *a_dic, crf1de_state_t *a_state)
{
  void *ret = NULL;
  int orig_len = a_state->m_len;

  for (int ilen = orig_len; ilen > 0; --ilen) {
    a_state->m_len = ilen;
    if ((ret = rumavl_find(a_dic, a_state)))
      break;
  }
  /* restore original length and return */
  a_state->m_len = orig_len;
  return (crf1de_state_t *) ret;
}

/**
 * Convert set of states to a vector.
 *
 * @param a_set - set of states
 * @param a_vec - vector to be populated
 *
 * @return \c int (\c 0 on SUCCESS and non-\c 0 otherwise)
 */
static int sm_set2vec(RUMAVL *a_set, crf1de_semimarkov_t *a_vec)
{
  int ret = 0;
  RUMAVL_NODE *node = NULL;
  crf1de_state_t *entry = NULL;

  while ((node = rumavl_node_next(a_set, node, 1, (void**) &entry)) != NULL) {
    if ((ret = copy_state(&sm->m_frw_states[entry->m_id], entry)))
      return ret;
  }
  return ret;
}

/**
 * Convert sets of states to vectors.
 *
 * @param sm - pointer to semi-markov model data
 *
 * @return \c int (\c 0 on SUCCESS and non-\c 0 otherwise)
 */
static int semimarkov_generate_state_vecs(crf1de_semimarkov_t *sm)
{
  ret = 0;

  sm->m_frw_states = (crf1de_state_t *) calloc(sm->m_num_frw, sizeof(crf1de_state_t));
  sm->m_bkw_states = (crf1de_state_t *) calloc(sm->m_num_bkw, sizeof(crf1de_state_t));
  sm->m_ptrns = (crf1de_state_t *) calloc(m_num_ptrns, sizeof(crf1de_state_t));

  if (sm->m_frw_states == NULL || sm->m_bkw_states == NULL || sm->m_ptrns == NULL) {
    ret = -1;
    goto final_steps;
  }

  if (sm_set2vec(sm->m_frw_states_set, sm->m_frw_states) || \
      sm_set2vec(sm->m_bkw_states_set, sm->m_bkw_states) || \
      sm_set2vec(sm->m_ptrns_set, sm->m_ptrns)) {
    ret = -2;
    goto final_steps;
  }

  return ret;

 error_exit:
  CLEAR(sm->m_frw_states);
  CLEAR(sm->m_bkw_states);
  CLEAR(sm->m_ptrns);
  return ret;
}

/**
 * Obtain longest suffixes for the given states.
 *
 * @param sm - pointer to semi-markov model
 * @param pky_id2ky_id - pointer to the mapping from suffix id to
 *             state id
 * @param ref_set - reference set for looking up suffixes
 * @param ref_vec - vector of states whose affix counters should be updated
 *
 * @return \c void
 */
void generate_frw_suffixes(crf1de_semimarkov_t *sm, int *pky_id2ky_id, \
			   RUMAVL *ref_set, crf1de_state_t *ref_vec)
{
  int ky_id;
  RUMAVL_NODE *node = NULL;
  crf1de_state_t *pky_entry, *ky_entry;

  for (size_t pky_id = 0; pky_id < sm->m_num_bkw; ++pky_id) {
    pky_entry = &sm->m_bkw_states[pky_id];
    ky_entry = semimarkov_find_max_sfx(ref_set, pky_entry);
    ky_id = ky_entry->m_id;
    pky_id2ky_id[pky_id] = ky_id;
    ++ref_vec[ky_id].m_num_affixes;
  }
}

/**
 * Initialize forward transition tables.
 *
 * @param sm - pointer to semi-markov model data
 *
 * @return \c int (\c 0 on SUCCESS and non-\c 0 otherwise)
 */
static int semimarkov_build_frw_transitions(crf1de_semimarkov_t *sm)
{
  int ret = 0;

  /* allocate necessary memory */
  int *pky_id2ky_id = (int *) calloc(sm->m_num_bkw, sizeof(int));
  sm->m_frw_trans1 = calloc(sm->m_num_bkw, sizeof(int));
  sm->m_frw_trans2 = calloc(sm->m_num_bkw, sizeof(int));

  if (pky_id2ky_id == NULL || sm->m_frw_trans1 == NULL || sm->m_frw_trans2 == NULL) {
    ret = -1;
    goto final_steps;
  }

  /* find longest suffixes for transitions of forwards states */
  generate_frw_suffixes(sm, pky_id2ky_id, sm->m_frw_states_set, sm->m_frw_states);

  /* set forward transitions pointers */
  int *trans1 = sm->m_frw_trans1;
  int *trans2 = sm->m_frw_trans2;
  crf1de_state_t *pk_entry = NULL;
  while (int pk_id = 0; pk_id < sm->m_num_frw; ++pk_id) {
    pk_entry = &sm->m_frw_states[pk_id];

    pk_entry->m_trans1 = trans1;
    pk_entry->m_trans2 = trans2;

    trans1 += pk_entry->m_num_affixes;
    trans2 += pk_entry->m_num_affixes;
    pk_entry->m_num_affixes = 0;
  }

  /* populate forward transitions of the states */
  int ky_id = -1;
  for (size_t pky_id = 0; pky_id < sm->m_num_bkw; ++pky_id) {
    ky_id = pky_id2ky_id[pky_id];
    pk_id = sm->m_bkw_states[pky_id].m_pk_id;

    ky_entry = &sm->m_frw_states[ky_id];
    ky_entry->m_frw_trans1[ky_entry->m_num_affixes] = pk_id;
    ky_entry->m_frw_trans2[ky_entry->m_num_affixes] = pky_id;
    ++ky_entry->m_num_affixes;
  }
  CLEAR(pky_id2ky_id);

  return ret;

 final_steps:
  CLEAR(pky_id2ky_id);
  CLEAR(sm->m_frw_trans2);
  CLEAR(sm->m_frw_trans1);
  return ret;
}

/**
 * Generate suffixes for given state and add them to semi-markov model.
 *
 * @param sm - pointer to semi-markov model
 * @param pky_entry - backward state for which all suffixes should be generated
 *
 * @return \c 0 on success, non-\c 0 otherwise
 */
static void semimarkov_build_suffixes(crf1de_semimarkov_t *sm, crf1de_state_t *pky_entry)
{
  size_t pky_id = pky_entry->m_id;
  size_t pky_len = pky_entry->m_len;
  size_t max_len = pky_len;
  int *sfxp = &SUFFIXES(sm, pky_id, 0);
  crf1de_state_t *ptrnp = NULL;

  for (size_t len = max_len; len > 1; --len) {
    pky_entry->m_len = len;
    if ((ptrnp = rumavl_find(sm->m_ptrns_set, pky_entry))) {
      *sfxp = ptrnp->m_id;
      ++sfxp;			/* increment the suffix pointer for the given `pky_id` */
      ++ptrnp->m_num_affixes;	/* increase the total number of prefixes for given suffix */
    }
  }
  pky_entry->m_len = pky_len;
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
  int ret = 0;

  /* allocate necessary space */
  int ifactor = sm->L - (sm->m_seg_len_lim < 0? 1: 0);
  sm->m_bkw_trans = calloc(sm->m_num_bkw * ifactor, sizeof(int));
  sm->m_ptrnid2bkwid = calloc(sm->m_num_ptrns, sizeof(int)); /* get rif of `m_ptrnid2bkwid`? */
  sm->m_ptrn_llabels = calloc(sm->m_num_ptrns, sizeof(int));

  if (sm->m_bkw_trans == NULL || sm->m_ptrnid2bkwid == NULL || sm->m_ptrn_llabels == NULL) {
    return -1;
    goto error_exit;
  }

  int i;
  size_t pky_len;
  RUMAVL_NODE *node = NULL;
  crf1de_state_t *pky_entry = NULL, *ptrn_entry = NULL;
  int y, pky_id, ptrn_id, last_label, *pky_seq = NULL, *bkw_trans = sm->m_bkw_trans;
  for (int pky_id = 0; pky_id < sm->m_num_bkw; ++pky_id) {
    pky_entry = &sm->m_bkw_states[pky_id];
    last_label = *pky_seq;
    pky_len = pky_entry->m_len;
    pky_seq = pky_entry->m_seq;

    /* obtain pattern correpsonding to the given backward state  */
    if ((ptrn_entry = rumavl_find(sm->m_ptrns_set, pky_entry))) {
      ptrn_id = ptrn_entry->m_id;
      sm->m_ptrn_llabels[ptrn_id] = last_label;
      sm->m_ptrnid2bkwid[ptrn_id] = pky_entry->m_id;
    }

    /* generate backward transitions for all possible labels */
    pky_entry->m_bkw_trans = bkw_trans;
    bkw_trans += ifactor;

    sm->m_wrkbench1.m_len = pky_len + 1;
    memcpy((void *) &sm->m_wrkbench1.m_seq[1], (const void *) pky_seq, pky_len * sizeof(int));

    i = 0;
    for (y = 0; y < sm->L; ++y) {
      if (y == last_label && sm->m_seg_len_lim < 0)
	continue;

      sm->m_wrkbench1.m_seq[0] = y;
      pkyy_entry = semimarkov_find_max_sfx(sm->m_bkw_states_set, &sm->m_wrkbench1);
      pky_entry->m_bkw_trans[i++] = semimarkov_find_max_sfx(sm->m_bkw_states_set, \
							    &sm->m_wrkbench1)->m_id;
    }
  }

  /* allocate space for pattern transitions */
  sm->m_ptrns = calloc(sm->m_num_ptrns, sizeof(crf1de_state_t));
  sm->m_ptrn_trans1 = calloc(sm->m_num_suffixes, sizeof(crf1de_state_t *));
  sm->m_ptrn_trans2 = calloc(sm->m_num_suffixes, sizeof(crf1de_state_t *));
  if (sm->m_ptrns == NULL || sm->m_ptrn_trans1 == NULL || sm->m_ptrn_trans2) {
    ret = -2;
    goto error_exit;
  }

  /* assign addresses of prefix arrays to patterns */
  node = NULL;
  int *ptrn_trans1 = sm->m_ptrn_trans1, *ptrn_trans2 = sm->m_ptrn_trans2;
  while ((node = rumavl_node_next(sm->m_ptrns_set, node, 1, (void**) &ptrn_entry)) != NULL) {
    ptrn_id = ptrn_entry->m_id;
    memcpy((void *) &sm->m_ptrns[ptrn_id], (const void *) ptrn_entry, sizeof(crf1de_state_t));
    ptrn_entry = &sm->m_ptrns[ptrn_id];

    ptrn_entry->m_frw_trans1 = ptrn_trans1;
    ptrn_trans1 += ptrn_entry->m_num_affixes;

    ptrn_entry->m_frw_trans2 = ptrn_trans2;
    ptrn_trans2 += ptrn_entry->m_num_affixes;
  }

  /* populate patterns with prefixes */
  size_t i;
  int *suffixes = NULL;

  for (size_t pky_id = 0; pky_id < sm->m_num_bkw; ++pky_id) {
    pky_entry = &sm->m_bkw_states[pky_id];
    suffixes = &SUFFIXES(sm, pky_id, 0);
    for (i = 0; (ptrn_id = suffixes[i]) >= 0; ++i) {
      ptrn_entry = &sm->m_ptrns[ptrn_id];

      ptrn_entry->m_frw_trans1[ptrn_entry->m_num_affixes] = sm->m_bkwid2frwid[pky_id];
      ptrn_entry->m_frw_trans2[ptrn_entry->m_num_affixes] = pky_id;
      ptrn_entry->m_num_affixes++;
    }
  }
  return ret;

 error_exit:
  CLEAR(sm->m_ptrn_trans1);
  CLEAR(sm->m_ptrns);
  CLEAR(sm->m_ptrn_llabels);
  CLEAR(sm->m_ptrnid2bkwid);
  CLEAR(sm->m_suffixes);
  CLEAR(sm->m_bkw_trans);
  return ret;
}

/**
 * Update frequency counters of detected patterns.
 *
 * @param sm - pointer to semi-markov model
 * @param a_wrkbench - pointer to workbench with pattern
 *
 * @return \c void
 */
static void semimarkov_update_patterns(crf1de_semimarkov_t *sm, crf1de_state_t *a_wrkbench)
{
  crf1de_state_t *ptrn_entry = NULL;
  size_t orig_len = a_wrkbench->m_len;

  while (a_wrkbench->m_len > 0) {
    ptrn_entry = (crf1de_state_t *) rumavl_find(sm->m_ptrns_set, a_wrkbench);
    ptrn_entry->m_freq += 1;
    --a_wrkbench->m_len;
  }

  a_wrkbench->m_len = orig_len;
}

/**
 * Generate all possible patterns and add them to the semimarkov model.
 *
 * @param sm - pointer to semi-markov model
 * @param a_wrkbench - pointer to workbench with pattern
 *
 * @return \c void
 */
static void semimarkov_add_patterns(crf1de_semimarkov_t *sm, crf1de_state_t *a_wrkbench)
{
  size_t orig_len = a_wrkbench->m_len;
  while (a_wrkbench->m_len > 1 && rumavl_find(sm->m_ptrns_set, a_wrkbench) == NULL) {
    a_wrkbench->m_id = sm->m_num_ptrns++;
    insert_state(sm->m_ptrns_set, a_wrkbench);

    /* reduce pattern for the next loop */
    --a_wrkbench->m_len;
  }
  /* increment frequency counters for remaining patterns */
  semimarkov_update_patterns(sm, a_wrkbench);
  /* restore original length */
  a_wrkbench->m_len = orig_len;
}

/**
 * Generate backward states for given prefix.
 *
 * @param sm - pointer to semi-markov model
 * @param a_wrkbench - pointer to workbench with prefix
 *
 * @return \c void
 */
static void semimarkov_add_bkw_states(crf1de_semimarkov_t *sm, crf1de_state_t *a_wrkbench)
{
  int last_label = a_wrkbench->m_seq[1];

  /* append all possible tags to given prefix */
  for (int i = 0; i < sm->L; ++i) {
    if (i == last_label && sm->m_seg_len_lim < 0)
      continue;

    /* we assume that backward state is not known */
    a_wrkbench->m_id = sm->m_num_bkw++;
    a_wrkbench->m_seq[0] = i;
    insert_state(sm->m_bkw_states_set, a_wrkbench);
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
static void semimarkov_add_states(crf1de_semimarkov_t *sm, crf1de_state_t *a_wrkbench)
{
  int wbn_len = a_wrkbench->m_len - 1;
  sm->m_wrkbench2.m_len = wbn_len;
  memcpy((void *) &sm->m_wrkbench2.m_seq, (const void *) &a_wrkbench->m_seq[1], \
	 wbn_len * sizeof(int));

  while (wbn_len-- > 1 && rumavl_find(sm->m_frw_states_set, &sm->m_wrkbench2) == NULL) {
    do {
      a_wrkbench.m_pk_id = sm->m_num_frw;
      sm->m_wrkbench2.m_id = sm->m_num_frw++;
      insert_state(sm->m_frw_states_set, &sm->m_wrkbench2);
      semimarkov_add_bkw_states(sm, a_wrkbench);
      --a_wrkbench->m_len;
    } while (--sm->m_wrkbench2.m_len > 1 && \
	     rumavl_find(sm->m_frw_states_set, &sm->m_wrkbench2) == NULL);

    sm->m_wrkbench2.m_len = wbn_len;
    memmove((void *) &sm->m_wrkbench2.m_seq, (const void *) &sm->m_wrkbench2.m_seq[1], \
	    sm->m_wrkbench2.m_len * sizeof(int));

    a_wrkbench->m_len = wbn_len + 1;
    memmove((void *) &a_wrkbench->m_seq, (const void *) &a_wrkbench->m_seq[1], \
	    a_wrkbench->m_len * sizeof(int));
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
  /* insert zero prefix in the forward state set */
  /* sm->m_wrkbench1.m_len = 0; */
  /* sm->m_wrkbench1.m_id = sm->m_num_frw++; */
  /* rumavl_insert(sm->m_frw_states_set, &sm->m_wrkbench1); */
  sm->m_wrkbench1.m_len = 1;
  sm->m_wrkbench1.m_pk_id = -1;
  sm->m_wrkbench2.m_len = 2;

  for (int i = 0; i < sm->L; ++i) {
    sm->m_wrkbench1.m_seq[0] = i;
    sm->m_wrkbench2.m_seq[1] = i;

    /* add label to the set of prefixes (forward states) */
    sm->m_wrkbench1.m_id = sm->m_num_frw;
    insert_state(sm->m_frw_states_set, &sm->m_wrkbench1);

    /* add label to the set of backward states and patterns */
    sm->m_wrkbench1.m_id = sm->m_num_bkw++;
    insert_state(sm->m_bkw_states_set, &sm->m_wrkbench1);

    sm->m_wrkbench1.m_id = sm->m_num_ptrns++;
    insert_state(sm->m_ptrns_set, &sm->m_wrkbench1);

    /* add two labels tag sequences to the backward state set */
    sm->m_wrkbench2.m_pk_id = sm->m_num_frw++;
    semimarkov_add_bkw_states(sm, &sm->m_wrkbench2);
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
  int ret = 0;

  if (L <= 0)
    return -1;

  sm->L = L;
  sm->m_max_order = ((size_t) a_max_order) + 1;
  sm->m_seg_len_lim = a_seg_len_lim;

  /* allocate memory for storing maximum segment lengths */
  sm->m_max_seg_len = calloc(L, sizeof(int));
  if (sm->m_max_seg_len == NULL)
    return -1;

  /* allocate space for tag sequences */
  sm->m_wrkbench1.m_seq = calloc(sm->m_max_order + 1, sizeof(int));
  sm->m_wrkbench2.m_seq = calloc(sm->m_max_order + 1, sizeof(int));
  if (sm->m_wrkbench1.m_seq == NULL || sm->m_wrkbench2.m_seq == NULL) {
    ret = -2;
    goto error_exit;
  }

  /* allocate sets for patterns, forward and backward states */
  sm->m_num_frw = 0;
  sm->m_frw_states_set = rumavl_new(sizeof(crf1de_state_t), sm_cmp_lseq, NULL, NULL);

  sm->m_num_bkw = 0;
  sm->m_bkw_states_set = rumavl_new(sizeof(crf1de_state_t), sm_cmp_lseq, NULL, NULL);

  sm->m_num_ptrns = 0;
  sm->m_ptrns_set = rumavl_new(sizeof(crf1de_state_t), sm_cmp_lseq, NULL, NULL);

  if (sm->m_frw_states_set == NULL || sm->m_bkw_states_set == NULL || \
      sm->m_ptrns_set == NULL) {
    ret = -1;
    goto error_exit;
  } else {
    *rumavl_delcb(sm->m_frw_states_set) = delcb; *rumavl_owcb(sm->m_frw_states_set) = owcb;
    *rumavl_delcb(sm->m_bkw_states_set) = delcb; *rumavl_owcb(sm->m_bkw_states_set) = owcb;
    *rumavl_delcb(sm->m_ptrns_set) = delcb; *rumavl_owcb(sm->m_ptrns_set) = owcb;
  }

  /* allocate space for auxiliary data structures */
  if (crfsuite_ring_create_instance(&sm->m_ring, sm->m_max_order)) {
    ret = -2;
    goto error_exit;
  }

  /* unconditionally add all labels to the sets of patterns, forward
     and backward states */
  semimarkov_initialize_states(sm);

  return ret;

 error_exit:
  rumavl_destroy(sm->m_frw_states_set);
  rumavl_destroy(sm->m_bkw_states_set);
  rumavl_destroy(sm->m_ptrns_set);
  free(sm->m_wrkbench1.m_seq);
  free(sm->m_wrkbench2.m_seq);
  free(sm->m_max_seg_len);
  return ret;
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

  /* transfer tag sequence from ring to the workbench */
  sm->build_state(&sm->m_wrkbench1, sm->m_ring);
  sm->m_wrkbench1.m_freq = 1;

  /* generate all possible prefixes and multiply them by L */
  if (rumavl_find(sm->m_ptrns_set, &sm->m_wrkbench1) == NULL) {
    semimarkov_add_patterns(sm, &sm->m_wrkbench1);
    semimarkov_add_states(sm, &sm->m_wrkbench1);
  } else {
    /* if pattern is already known, then increment its frequency counters */
    semimarkov_update_patterns(sm, &sm->m_wrkbench1);
  }
}

/**
 * Create all possible transitions (up-to and including max_order)
 *
 * @param sm - pointer to semi-markov model
 * @param a_order - feature order
 * @param a_prev_label - previous sequence label
 * @param a_wrkbench - partially constructed tag sequence
 *
 * @return \c void
 */
static void semimarkov_generate_edges_helper(crf1de_semimarkov_t *sm, size_t a_order, \
					     int a_prev_label, crf1de_state_t *a_wrkbench)
{
  size_t orig_len = a_wrkbench->m_len;
  a_wrkbench->m_len = a_order + 1;
  a_wrkbench->m_freq = 0;

  for (int i = 0; i < sm->L; ++i) {
    if (i == a_prev_label && sm->m_seg_len_lim < 0)
      continue;

    a_wrkbench->m_seq[a_order] = i;

    /* add patterns and states */
    if (! rumavl_find(sm->m_ptrns_set, a_wrkbench)) {
      a_wrkbench->m_id = sm->m_num_ptrns++;
      insert_state(sm->m_ptrns_set, a_wrkbench);
    }

    if (! rumavl_find(sm->m_bkw_states_set, a_wrkbench)) {
      a_wrkbench->m_id = sm->m_num_bkw++;
      insert_state(sm->m_bkw_states_set, a_wrkbench);
    }

    /* recursively invoke function if we did not already reach the maximum order */
    if (a_wrkbench->m_len < sm->m_max_order) {
      if (rumavl_find(sm->m_frw_states_set, a_wrkbench) == NULL) {
	a_wrkbench->m_id = sm->m_num_frw++;
	insert_state(sm->m_frw_states_set, a_wrkbench);
      }
      semimarkov_generate_edges_helper(sm, a_wrkbench->m_len, i, a_wrkbench);
    }
  }
  /* restore original length */
  a_wrkbench->m_len = orig_len;
}

/**
 * Create all possible transitions (up-to and including max_order)
 *
 * @param sm - pointer to semi-markov model
 *
 * @return \c void
 */
static void semimarkov_generate_all_edges(crf1de_semimarkov_t *sm)
{
  semimarkov_generate_edges_helper(sm, 0, -1, &sm->m_wrkbench1);
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
  int ret = 0;

  /* convert sets to vectors */
  if (semimarkov_generate_state_vecs(sm))
    return -1;

  /* generate forward transitions */
  if (semimarkov_build_frw_transitions(sm)) {
    ret = -2;
    goto exit_section;
  }

  /* generate backward transitions */
  if (semimarkov_build_bkw_transitions(sm)) {
    CLEAR(sm->m_frw_states);
    CLEAR(sm->m_bkw_states);
    CLEAR(sm->m_frw_trans1);
    CLEAR(sm->m_frw_trans2);
    CLEAR(sm->m_frw_llabels);
    CLEAR(sm->m_bkwid2frwid);
    ret = -3;
    goto exit_section;
  }

  /* debug states and transitions */
  if (1) {
    semimarkov_debug_states(sm);
    semimarkov_debug_transitions(sm);
    exit(66);
  }

 exit_section:
  return ret;
}

/**
 * Create new state from ring of labels and store it in `a_state`.
 *
 * @param a_state - address of state in which new state should be constructed
 * @param a_ring - ring containing sequence of labels
 *
 * @return \c void
 */
void semimarkov_build_state(crf1de_state_t *a_state, const crfsuite_ring_t *a_ring)
{
  memset(a_state, 0, sizeof(crf1de_state_t));
  a_state->m_len = a_ring->num_items;
  const crfsuite_chain_link_t *chlink = a_ring->tail->prev;

  for (size_t i = 0; i < a_ring->num_items; ++i) {
    a_state->m_seq[i] = chlink->data;
    chlink = chlink->prev;
  }
}

/**
 * Obtain internal `id` of given state.
 *
 * @param a_state - state whose id should be obtained
 * @param a_dic - reference dictionary in which we should search the id
 *
 * @return \c int - id of the state
 */
inline static int semimarkov_get_state_id(crf1de_state_t *a_state,	\
					  RUMAVL *a_dic)
{
  crf1de_state_t *item = (crf1de_state_t *) rumavl_find(a_dic, a_state);

  if (item)
    return item->m_id;

  return -1;
}

/**
 * Deallocate space used by states.
 *
 * @param a_states - pointer to an array of semi-markov states
 * @param a_cnt - number of states to be deallocated
 *
 * @return \c void
 */
static void semimarkov_clear_states(crf1de_state_t **a_states)
{
  if (*a_states == NULL)
    return;

  for (size_t i = 0; i < a_cnt; ++a) {
    free((*a_states)[i].m_seq);
  }

  free(*a_states);
  *a_states = NULL;
}

/**
 * Deallocate space and reset values.
 *
 * @param sm - pointer to semi-markov model data
 *
 * @return \c void
 */
static void semimarkov_clear(crf1de_semimarkov_t *sm)
{
  CLEAR(sm->m_max_seg_len);
  /* clear patterns */
  semimarkov_clear_states(&sm->m_ptrns, sm->m_num_ptrns);
  RUMAVL_CLEAR(sm->m_ptrns_set);
  sm->m_num_ptrns = 0;

  CLEAR(sm->m_ptrn_llabels);
  CLEAR(sm->m_ptrn_trans1);
  CLEAR(sm->m_ptrn_trans2);
  CLEAR(sm->m_ptrnid2bkwid);

  CLEAR(sm->m_suffixes);
  sm->m_num_suffixes = 0;

  /* clear forward states */
  semimarkov_clear_states(&sm->m_frw_states, sm->m_num_frw);
  RUMAVL_CLEAR(sm->m_frw_states_set);
  sm->m_num_frw = 0;

  CLEAR(sm->m_frw_llabels);
  CLEAR(sm->m_frw_trans1);
  CLEAR(sm->m_frw_trans2);

  /* clear backward states */
  semimarkov_clear_states(&sm->m_bkw_states, sm->m_num_bkw);
  RUMAVL_CLEAR(sm->m_bkw_states_set);
  sm->m_num_bkw = 0;

  CLEAR(sm->m_bkw_trans);
  CLEAR(sm->m_bkwid2frwid);

  /* auxiliary data members */
  if (sm->m_ring) {
    sm->m_ring->free(sm->m_ring);
    free(sm->m_ring);
    sm->m_ring = NULL;
  }
}

crf1de_semimarkov_t *crf1de_create_semimarkov(void) {
  crf1de_semimarkov_t *sm = calloc(1, sizeof(crf1de_semimarkov_t));
  if (sm == NULL)
    return sm;
  sm->initialize = semimarkov_initialize;
  sm->update = semimarkov_update;
  sm->generate_all_edges = semimarkov_generate_all_edges;
  sm->finalize = semimarkov_finalize;
  sm->clear = semimarkov_clear;
  sm->output_state = semimarkov_output_state;
  sm->build_state = semimarkov_build_state;
  sm->get_state_id = semimarkov_get_state_id;

  return sm;
}
