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
      semimarkov_output_state(stderr, NULL, &sm->m_frw_states[pk_entry->m_frw_trans1[j]]);
      fprintf(stderr, ";");
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "forward_transition2[(id = %zu) ", i);
    semimarkov_output_state(stderr, NULL, pk_entry);
    fprintf(stderr, "] =");

    fprintf(stderr, " (%zu):", pk_entry->m_num_affixes);
    for (j = 0; j < pk_entry->m_num_affixes; ++j) {
      fprintf(stderr, " ");
      semimarkov_output_state(stderr, NULL, &sm->m_bkw_states[pk_entry->m_frw_trans2[j]]);
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
      semimarkov_output_state(stderr, NULL, &sm->m_frw_states[ptrn_entry->m_frw_trans1[j]]);
      fprintf(stderr, ";");
    }
    fprintf(stderr, "\n");

    /* output pattern transition 2 */
    fprintf(stderr, "patternTransition2[");
    semimarkov_output_state(stderr, NULL, ptrn_entry);
    fprintf(stderr, "] =");
    for (j = 0; j < ptrn_entry->m_num_affixes; ++j) {
      fprintf(stderr, " ");
      semimarkov_output_state(stderr, NULL, &sm->m_bkw_states[ptrn_entry->m_frw_trans2[j]]);
      fprintf(stderr, ";");
    }
    fprintf(stderr, "\n");
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
 * Initialize forward transition tables.
 *
 * @param sm - pointer to semi-markov model data
 *
 * @return \c int (0 on SUCCESS and non-0 otherwise)
 */
static int semimarkov_build_frw_transitions(crf1de_semimarkov_t *sm)
{
  /* allocate memory for storing last labels of the prefixes */
  sm->m_frw_llabels = (int *) calloc(sm->m_num_frw, sizeof(int));
  if (sm->m_frw_llabels == NULL)
    return -1;

  /* allocate memory for mapping of pky states to corresponding pk states */
  sm->m_bkwid2frwid = (int *) calloc(sm->m_num_bkw, sizeof(int));
  if (sm->m_bkwid2frwid == NULL) {
    CLEAR(sm->m_frw_llabels);
    return -2;
  }

  /* create temporary map for storing states to which given affixes
     correspond */
  int *pky_id2ky_id = (int *) malloc(sm->m_num_bkw * sizeof(int));
  if (pky_id2ky_id == NULL) {
    CLEAR(sm->m_bkwid2frwid);
    CLEAR(sm->m_frw_llabels);
    return -3;
  }
  memset((void *) pky_id2ky_id, -1, sm->m_num_bkw * sizeof(int));

  /* in the first run, populate the `m_frw_llabels` array and
     calculate the number of prefixes */
  size_t pk_len;
  int y, pk_id, pky_id, last_label;
  const int *pk_start = NULL;

  RUMAVL_NODE *node = NULL;
  crf1de_state_t *pk_entry = NULL, *pky_entry = NULL, *ky_entry = NULL;

  /* find longest suffixes and states */
  while ((node = rumavl_node_next(sm->m__frw_states_set, node, 1, (void**) &pk_entry)) != NULL) {
    pk_id = pk_entry->m_id;
    pk_len = pk_entry->m_len;
    pk_start = pk_entry->m_seq;

    if (pk_len)
      last_label = *pk_start;
    else
      last_label = -1;

    sm->m_frw_llabels[pk_id] = last_label;
    /* copy label sequence to `sm->m_wrkbench1' and append to it all
       possible tags */
    sm->m_wrkbench1.m_len = pk_len + 1;
    memcpy((void *) &sm->m_wrkbench1.m_seq[1], (const void *) pk_start, pk_len * sizeof(int));
    for (y = 0; y < sm->L; ++y) {
      if (y == last_label && sm->m_seg_len_lim < 0)
	continue;

      /* find which `pky` state corresponds to give `pk` prefix */
      sm->m_wrkbench1.m_seq[0] = y;
      pky_entry = (crf1de_state_t *) rumavl_find(sm->m__bkw_states_set, &sm->m_wrkbench1);
      pky_id = pky_entry->m_id;
      sm->m_bkwid2frwid[pky_id] = pk_id;

      /* increment the number of prefixes corresponding to the given `ky` suffix */
      ky_entry = semimarkov_find_max_sfx(sm->m__frw_states_set, &sm->m_wrkbench1);
      pky_id2ky_id[pky_id] = ky_entry->m_id;
      ++ky_entry->m_num_affixes;
    }
  }

  /* allocate memory for the lists of pk transtitions */
  sm->m_frw_trans1 = calloc(sm->m_num_bkw, sizeof(int));
  if (sm->m_frw_trans1 == NULL) {
    CLEAR(pky_id2ky_id);
    CLEAR(sm->m_bkwid2frwid);
    CLEAR(sm->m_frw_llabels);
    return -4;
  }

  /* allocate memory for the lists of pky states corresponding to pk
     transtitions */
  sm->m_frw_trans2 = calloc(sm->m_num_bkw, sizeof(int));
  if (sm->m_frw_trans2 == NULL) {
    CLEAR(sm->m_frw_trans1);
    CLEAR(pky_id2ky_id);
    CLEAR(sm->m_bkwid2frwid);
    CLEAR(sm->m_frw_llabels);
    return -5;
  }

  /* allocate memory for the vector of forward states */
  sm->m_frw_states = calloc(sm->m_num_frw, sizeof(crf1de_state_t));
  if (sm->m_frw_states == NULL) {
    CLEAR(sm->m_frw_trans2);
    CLEAR(sm->m_frw_trans1);
    CLEAR(pky_id2ky_id);
    CLEAR(sm->m_bkwid2frwid);
    CLEAR(sm->m_frw_llabels);
    return -6;
  }

  /* populate forward transitions on the basis of previously computed
     `pky_id2ky_id' and `pky2pk' */
  node = NULL;
  int *frw_trans1 = sm->m_frw_trans1;
  int *frw_trans2 = sm->m_frw_trans2;
  while ((node = rumavl_node_next(sm->m__frw_states_set, node, 1, (void**) &pk_entry)) != NULL) {
    pk_id = pk_entry->m_id;

    memcpy((void *) &sm->m_frw_states[pk_id], (const void *) pk_entry, sm->m_frw_size);
    pk_entry = &sm->m_frw_states[pk_id];

    pk_entry->m_frw_trans1 = frw_trans1;
    pk_entry->m_frw_trans2 = frw_trans2;

    frw_trans1 += pk_entry->m_num_affixes;
    frw_trans2 += pk_entry->m_num_affixes;
  }
  /* RUMAVL_CLEAR(sm->m__frw_states_set); */

  /* populate forward transitions of the states */
  int ky_id = -1;
  for (size_t pky_id = 0; pky_id < sm->m_num_bkw; ++pky_id) {
    if (pky_id2ky_id[pky_id] < 0)
      continue;

    pk_id = sm->m_bkwid2frwid[pky_id];
    ky_id = pky_id2ky_id[pky_id];

    ky_entry = &sm->m_frw_states[ky_id];
    ky_entry->m_frw_trans1[ky_entry->m__cnt_trans1++] = pk_id;
    ky_entry->m_frw_trans2[ky_entry->m__cnt_trans2++] = pky_id;
  }
  free(pky_id2ky_id);

  return 0;
}

/**
 * Generate suffixes for given state and add them to semi-markov model.
 *
 * @param sm - pointer to semi-markov model
 * @param pky_entry - backward state for which all suffixes should be generated
 *
 * @return \c void
 */
static void semimarkov_build_suffixes(crf1de_semimarkov_t *sm, crf1de_state_t *pky_entry)
{
  size_t pky_id = pky_entry->m_id;
  size_t pky_len = pky_entry->m_len;
  size_t max_len = pky_len;
  int *sfxp = &SUFFIXES(sm, pky_id, 0);
  crf1de_state_t *ptrnp = NULL;
  /* TODO: only suffixes with length greater than 1 should be kept */
  for (size_t len = max_len; len > 0; --len) {
    pky_entry->m_len = len;
    if ((ptrnp = rumavl_find(sm->m__ptrns_set, pky_entry))) {
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
  sm->m_bkw_trans = calloc(sm->m_num_bkw, sm->L * sizeof(int));
  if (sm->m_bkw_trans == NULL)
    return -1;

  /* additional byte in suffix array will serve as a sentinel */
  sm->m_num_suffixes = sm->m_num_bkw * (sm->m_max_order + 1);
  size_t sfx_vec_size = sm->m_num_suffixes * sizeof(int);
  sm->m_suffixes = (int *) malloc(sfx_vec_size);
  if (sm->m_suffixes == NULL) {
    CLEAR(sm->m_bkw_trans);
    return -2;
  }
  memset((void *) sm->m_suffixes, -1, sfx_vec_size);

  sm->m_bkw_states = calloc(sm->m_num_bkw, sizeof(crf1de_state_t));
  if (sm->m_bkw_states == NULL) {
    CLEAR(sm->m_suffixes);
    CLEAR(sm->m_bkw_trans);
    return -3;
  }

  sm->m_ptrnid2bkwid = calloc(sm->m_num_ptrns, sizeof(int));
  if (sm->m_ptrnid2bkwid == NULL) {
    CLEAR(sm->m_bkw_states);
    CLEAR(sm->m_suffixes);
    CLEAR(sm->m_bkw_trans);
    return -4;
  }

  sm->m_ptrn_llabels = calloc(sm->m_num_ptrns, sizeof(int));
  if (sm->m_ptrn_llabels == NULL) {
    CLEAR(sm->m_ptrnid2bkwid);
    CLEAR(sm->m_bkw_states);
    CLEAR(sm->m_suffixes);
    CLEAR(sm->m_bkw_trans);
    return -5;
  }

  size_t pky_len;
  int y, pky_id, ptrn_id, last_label, *pky_seq = NULL, *bkw_trans = sm->m_bkw_trans;

  RUMAVL_NODE *node = NULL;
  crf1de_state_t *pky_entry = NULL, *ptrn_entry = NULL;

  while ((node = rumavl_node_next(sm->m__bkw_states_set, node, 1, (void**) &pky_entry)) != NULL) {
    pky_id = pky_entry->m_id;
    pky_len = pky_entry->m_len;
    pky_seq = pky_entry->m_seq;
    last_label = *pky_seq;

    /* copy backward state to the newly constructed array of backward states */
    memcpy((void *) &sm->m_bkw_states[pky_id], (const void *) pky_entry, sm->m_bkw_size);
    pky_entry = &sm->m_bkw_states[pky_id];

    /* generate suffixes */
    if ((ptrn_entry = rumavl_find(sm->m__ptrns_set, pky_entry))) {
      ptrn_id = ptrn_entry->m_id;
      sm->m_ptrn_llabels[ptrn_id] = last_label;
      sm->m_ptrnid2bkwid[ptrn_id] = pky_entry->m_id;
    }
    semimarkov_build_suffixes(sm, pky_entry);

    /* generate backward transitions for all possible labels */
    pky_entry->m_bkw_trans = bkw_trans;
    bkw_trans += sm->L;

    sm->m_wrkbench1.m_len = pky_len + 1;
    memcpy((void *) &sm->m_wrkbench1.m_seq[1], (const void *) pky_seq, pky_len * sizeof(int));

    for (y = 0; y < sm->L; ++y) {
      if (y == last_label && sm->m_seg_len_lim < 0) {
	pky_entry->m_bkw_trans[y] = -1;
      } else {
	sm->m_wrkbench1.m_seq[0] = y;
	pky_entry->m_bkw_trans[y] = semimarkov_find_max_sfx(sm->m__bkw_states_set, &sm->m_wrkbench1)->m_id;
      }
    }
  }
  /* RUMAVL_CLEAR(sm->m__bkw_states_set); */

  /* allocate space for pattern array */
  sm->m_ptrns = calloc(sm->m_num_ptrns, sizeof(crf1de_state_t));
  if (sm->m_ptrns == NULL) {
    CLEAR(sm->m_ptrn_llabels);
    CLEAR(sm->m_ptrnid2bkwid);
    CLEAR(sm->m_bkw_states);
    CLEAR(sm->m_suffixes);
    CLEAR(sm->m_bkw_trans);
    return -6;
  }

  /* allocate space for pattern transitions */
  sm->m_ptrn_trans1 = calloc(sm->m_num_suffixes, sizeof(crf1de_state_t *));
  if (sm->m_ptrn_trans1 == NULL) {
    CLEAR(sm->m_ptrns);
    CLEAR(sm->m_ptrn_llabels);
    CLEAR(sm->m_ptrnid2bkwid);
    CLEAR(sm->m_suffixes);
    CLEAR(sm->m_bkw_trans);
    return -7;
  }

  sm->m_ptrn_trans2 = calloc(sm->m_num_suffixes, sizeof(crf1de_state_t *));
  if (sm->m_ptrn_trans2 == NULL) {
    CLEAR(sm->m_ptrn_trans1);
    CLEAR(sm->m_ptrns);
    CLEAR(sm->m_ptrn_llabels);
    CLEAR(sm->m_ptrnid2bkwid);
    CLEAR(sm->m_suffixes);
    CLEAR(sm->m_bkw_trans);
    return -8;
  }

  /* assign addresses of prefix arrays to patterns */
  node = NULL;
  int *ptrn_trans1 = sm->m_ptrn_trans1, *ptrn_trans2 = sm->m_ptrn_trans2;
  while ((node = rumavl_node_next(sm->m__ptrns_set, node, 1, (void**) &ptrn_entry)) != NULL) {
    ptrn_id = ptrn_entry->m_id;
    memcpy((void *) &sm->m_ptrns[ptrn_id], (const void *) ptrn_entry, sm->m_bkw_size);
    ptrn_entry = &sm->m_ptrns[ptrn_id];

    ptrn_entry->m_frw_trans1 = ptrn_trans1;
    ptrn_trans1 += ptrn_entry->m_num_affixes;

    ptrn_entry->m_frw_trans2 = ptrn_trans2;
    ptrn_trans2 += ptrn_entry->m_num_affixes;
  }
  /* RUMAVL_CLEAR(sm->m__ptrns_set); */

  /* populate patterns with prefixes */
  size_t i;
  int *suffixes = NULL;

  for (size_t pky_id = 0; pky_id < sm->m_num_bkw; ++pky_id) {
    pky_entry = &sm->m_bkw_states[pky_id];
    suffixes = &SUFFIXES(sm, pky_id, 0);
    for (i = 0; (ptrn_id = suffixes[i]) >= 0; ++i) {
      ptrn_entry = &sm->m_ptrns[ptrn_id];

      ptrn_entry->m_frw_trans1[ptrn_entry->m__cnt_trans1++] = sm->m_bkwid2frwid[pky_id];
      ptrn_entry->m_frw_trans2[ptrn_entry->m__cnt_trans2++] = pky_id;
    }
  }
  return 0;
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
    ptrn_entry = (crf1de_state_t *) rumavl_find(sm->m__ptrns_set, a_wrkbench);
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
  while (a_wrkbench->m_len > 1 && rumavl_find(sm->m__ptrns_set, a_wrkbench) == NULL) {
    a_wrkbench->m_id = sm->m_num_ptrns++;
    rumavl_insert(sm->m__ptrns_set, a_wrkbench);

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
    rumavl_insert(sm->m__bkw_states_set, a_wrkbench);
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

  while (wbn_len-- > 1 && rumavl_find(sm->m__frw_states_set, &sm->m_wrkbench2) == NULL) {
    do {
      sm->m_wrkbench2.m_id = sm->m_num_frw++;
      rumavl_insert(sm->m__frw_states_set, &sm->m_wrkbench2);
      semimarkov_add_bkw_states(sm, a_wrkbench);
      --a_wrkbench->m_len;
    } while (--sm->m_wrkbench2.m_len > 1 && \
	     rumavl_find(sm->m__frw_states_set, &sm->m_wrkbench2) == NULL);

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
  /* insert zero prefix in forward state set */
  sm->m_wrkbench1.m_len = 0;
  sm->m_wrkbench1.m_id = sm->m_num_frw++;
  rumavl_insert(sm->m__frw_states_set, &sm->m_wrkbench1);

  sm->m_wrkbench1.m_len = 1;
  sm->m_wrkbench2.m_len = 2;

  for (int i = 0; i < sm->L; ++i) {
    sm->m_wrkbench2.m_seq[1] = sm->m_wrkbench1.m_seq[0] = i;

    /* add label to the set of prefixes (forward states) */
    sm->m_wrkbench1.m_id = sm->m_num_frw++;
    rumavl_insert(sm->m__frw_states_set, &sm->m_wrkbench1);

    /* add label to the set of backward states and patterns */
    sm->m_wrkbench1.m_id = sm->m_num_ptrns++;
    rumavl_insert(sm->m__ptrns_set, &sm->m_wrkbench1);

    sm->m_wrkbench1.m_id = sm->m_num_bkw++;
    rumavl_insert(sm->m__bkw_states_set, &sm->m_wrkbench1);

    /* add two labels tag sequences to the backward state set */
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
  if (L <= 0)
    return -1;

  sm->L = L;
  sm->m_max_order = ((size_t) a_max_order) + 1;
  if (sm->m_max_order > CRFSUITE_SM_MAX_PTRN_LEN) {
    fprintf(stderr, "Max order (%zu) exceeds the limit (%d). ", \
	    sm->m_max_order, CRFSUITE_SM_MAX_PTRN_LEN);
    fprintf(stderr, "Max order is set to %d. ", CRFSUITE_SM_MAX_PTRN_LEN);
    fprintf(stderr, "To increase this limit, recompile the program with the option "
	    "-DCRFSUITE_SM_MAX_PTRN_LEN=NEW_LIM added to CPPFLAGS.\n");
    sm->m_max_order = CRFSUITE_SM_MAX_PTRN_LEN;
  }
  sm->m_seg_len_lim = a_seg_len_lim;

  /* allocate memory for storing maximum segment lengths */
  sm->m_max_seg_len = calloc(L, sizeof(int));
  if (sm->m_max_seg_len == NULL)
    return -1;

  /* allocate sets for patterns, forward and backward states */
  sm->m_frw_size = sizeof(crf1de_state_t) - sizeof(int) * \
    (CRFSUITE_SM_MAX_PTRN_LEN - sm->m_max_order);
  sm->m_bkw_size = sm->m_frw_size + sizeof(int);

  sm->m_num_frw = 0;
  sm->m__frw_states_set = rumavl_new(sm->m_frw_size, crf1de_cmp_lseq, NULL, NULL);

  sm->m_num_bkw = 0;
  sm->m__bkw_states_set = rumavl_new(sm->m_bkw_size, crf1de_cmp_lseq, NULL, NULL);

  sm->m_num_ptrns = 0;
  sm->m__ptrns_set = rumavl_new(sm->m_bkw_size, crf1de_cmp_lseq, NULL, NULL);

  /* allocate space for auxiliary data structures */
  if (crfsuite_ring_create_instance(&sm->m_ring, sm->m_max_order)) {
    free(sm->m_max_seg_len);
    return -1;
  }

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

  /* transfer tag sequence from ring to the workbench */
  sm->build_state(&sm->m_wrkbench1, sm->m_ring);
  sm->m_wrkbench1.m_freq = 1;

  /* generate all possible prefixes and multiply them by L */
  if (rumavl_find(sm->m__ptrns_set, &sm->m_wrkbench1) == NULL) {
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
    if (! rumavl_find(sm->m__ptrns_set, a_wrkbench)) {
      a_wrkbench->m_id = sm->m_num_ptrns++;
      rumavl_insert(sm->m__ptrns_set, a_wrkbench);
    }

    if (! rumavl_find(sm->m__bkw_states_set, a_wrkbench)) {
      a_wrkbench->m_id = sm->m_num_bkw++;
      rumavl_insert(sm->m__bkw_states_set, a_wrkbench);
    }

    /* recursively invoke function if we did not already reach the maximum order */
    if (a_wrkbench->m_len < sm->m_max_order) {
      if (rumavl_find(sm->m__frw_states_set, a_wrkbench) == NULL) {
	a_wrkbench->m_id = sm->m_num_frw++;
	rumavl_insert(sm->m__frw_states_set, a_wrkbench);
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
  /* generate forward transitions */
  if (semimarkov_build_frw_transitions(sm)) {
    ret = -1;
    goto exit_section;
  }

  /* generate backward transitions */
  if (semimarkov_build_bkw_transitions(sm)) {
    CLEAR(sm->m_frw_states);
    CLEAR(sm->m_frw_trans1);
    CLEAR(sm->m_frw_trans2);
    CLEAR(sm->m_frw_llabels);
    ret = -2;
    goto exit_section;
  }

  /* debug states and transitions */
  if (0) {
    semimarkov_debug_states(sm);
    semimarkov_debug_transitions(sm);
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
 * Deallocate space and reset values.
 *
 * @param sm - pointer to semi-markov model data
 *
 * @return void
 */
static void semimarkov_clear(crf1de_semimarkov_t *sm)
{
  CLEAR(sm->m_max_seg_len);
  /* clear patterns */
  CLEAR(sm->m_ptrns);
  RUMAVL_CLEAR(sm->m__ptrns_set);
  sm->m_num_ptrns = 0;

  CLEAR(sm->m_ptrn_llabels);
  CLEAR(sm->m_ptrn_trans1);
  CLEAR(sm->m_ptrn_trans2);
  CLEAR(sm->m_ptrnid2bkwid);

  CLEAR(sm->m_suffixes);
  sm->m_num_suffixes = 0;

  /* clear forward states */
  CLEAR(sm->m_frw_states);
  RUMAVL_CLEAR(sm->m__frw_states_set);
  sm->m_num_frw = 0;

  CLEAR(sm->m_frw_llabels);
  CLEAR(sm->m_frw_trans1);
  CLEAR(sm->m_frw_trans2);

  /* clear backward states */
  CLEAR(sm->m_bkw_states);
  RUMAVL_CLEAR(sm->m__bkw_states_set);
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
