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

/// function for accessing suffix list
#define SUFFIXES(sm, y, x)				\
  MATRIX(sm->m_suffixes, (sm->m_max_order + 1), x, y)

/// function for freeing an item if it is not null
#define CLEAR(a_item)				\
  if ((a_item)) {				\
    free(a_item);				\
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
  int pk_len = a_entry->m_len;

  if (a_name)
    fprintf(a_fstream, "%s[%d] (%d) = %d", a_name, pk_id, pk_len, a_entry->m_seq[0]);
  else
    fprintf(a_fstream, "%d", a_entry->m_seq[0]);

  for (int i = 1; i < pk_len; ++i) {
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
  RUMAVL_NODE *node = NULL;
  crf1de_state_t *entry = NULL;

  /* output patterns */
  while ((node = rumavl_node_next(sm->m_patterns, node, 1, (void**) &entry)) != NULL) {
    semimarkov_output_state(stderr, "sm->m_patterns", entry);
  }
  fprintf(stderr, "******************************************************************\n");

  /* output forward states */
  node = NULL;
  while ((node = rumavl_node_next(sm->m_frw_states, node, 1, (void**) &entry)) != NULL) {
    semimarkov_output_state(stderr, "sm->m_frw_states", entry);
  }
  fprintf(stderr, "******************************************************************\n");

  /* output backward states */
  node = NULL;
  while ((node = rumavl_node_next(sm->m_bkw_states, node, 1, (void**) &entry)) != NULL) {
    semimarkov_output_state(stderr, "sm->m_bkw_states", entry);
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
    pk_entry = sm->m_frwid2frw[i];
    fprintf(stderr, "forward_transition1[");
    semimarkov_output_state(stderr, NULL, pk_entry);
    fprintf(stderr, "] = ");

    fprintf(stderr, "(%d): ", pk_entry->m_n_prefixes);
    for (j = 0; j < pk_entry->m_n_prefixes; ++j) {
      semimarkov_output_state(stderr, NULL, pk_entry->m_frw_trans1[j]);
      fprintf(stderr, "; ");
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "forward_transition1[(id = %d) ", pk_entry->m_id);
    semimarkov_output_state(stderr, NULL, pk_entry);
    fprintf(stderr, "] = ");

    fprintf(stderr, "(%d): ", pk_entry->m_n_prefixes);
    for (j = 0; j < pk_entry->m_n_prefixes; ++j) {
      semimarkov_output_state(stderr, NULL, pk_entry->m_frw_trans2[j]);
      fprintf(stderr, "; ");
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "******************************************************************\n");

  /* output backward transition for states */
  for (size_t i = 0; i < sm->m_num_bkw; ++i) {
    pky_entry = sm->m_bkwid2bkw[i];
    for (j = 0; j < sm->L; ++j) {
      fprintf(stderr, "backwardTransition[");
      semimarkov_output_state(stderr, NULL, pky_entry);
      fprintf(stderr, "][%d] = ", j);

      if (pky_entry->m_bkw_trans[j])
	semimarkov_output_state(stderr, NULL, pky_entry->m_bkw_trans[j]);

      fprintf(stderr, "\n");
    }
  }
  fprintf(stderr, "******************************************************************\n");

  /* output suffixes */
  size_t i;
  crf1de_state_t *ptrn_entry = NULL, **suffixes = NULL;
  for (size_t pky_id = 0; pky_id < sm->m_num_bkw; ++pky_id) {
    pky_entry = sm->m_bkwid2bkw[pky_id];
    fprintf(stderr, "suffixes[");
    semimarkov_output_state(stderr, NULL, pky_entry);
    fprintf(stderr, "] =");

    suffixes = &SUFFIXES(sm, pky_id, 0);
    for (i = 0; ptrn_entry = suffixes[i]; ++i) {
      fprintf(stderr, " ");
      semimarkov_output_state(stderr, NULL, ptrn_entry);
      fprintf(stderr, ";");
      /* ptrn_entry->m_frw_trans1[ptrn_entry->m__cnt_trans1++] = sm->m_bkwid2frw[pky_entry->m_id]; */
      /* ptrn_entry->m_frw_trans2[ptrn_entry->m__cnt_trans2++] = pky_entry; */
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "******************************************************************\n");
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
    if (ret = entry1->m_seq[i] - entry2->m_seq[i])
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
    if (ret = rumavl_find(a_dic, a_state))
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
  sm->m_frw_llabels = calloc(sm->m_num_frw, sizeof(int));
  if (sm->m_frw_llabels == NULL)
    return -1;

  /* allocate memory for map of prefix id's to prefix objects */
  sm->m_frwid2frw = calloc(sm->m_num_frw, sizeof(crf1de_state_t *));
  if (sm->m_frwid2frw == NULL) {
    CLEAR(sm->m_frw_llabels);
    return -2;
  }

  sm->m_bkwid2bkw  = calloc(sm->m_num_bkw, sizeof(crf1de_state_t *));
  if (sm->m_bkwid2bkw == NULL) {
    CLEAR(sm->m_frwid2frw);
    CLEAR(sm->m_frw_llabels);
    return -3;
  }

  /* allocate memory for the lists of pk transtitions */
  sm->m_frw_trans1 = calloc(sm->m_num_bkw, sizeof(crf1de_state_t *));
  if (sm->m_frw_trans1 == NULL) {
    CLEAR(sm->m_bkwid2bkw);
    CLEAR(sm->m_frwid2frw);
    CLEAR(sm->m_frw_llabels);
    return -4;
  }

  /* allocate memory for the lists of pky states corresponding to pk
     transtitions */
  sm->m_frw_trans2 = calloc(sm->m_num_bkw, sizeof(crf1de_state_t *));
  if (sm->m_frw_trans2 == NULL) {
    CLEAR(sm->m_frw_trans1);
    CLEAR(sm->m_bkwid2bkw);
    CLEAR(sm->m_frwid2frw);
    CLEAR(sm->m_frw_llabels);
    return -5;
  }

  /* create temporary maps for storing states to which given affixes
     correspond */
  crf1de_state_t **pky2ky = calloc(sm->m_num_bkw, sizeof(crf1de_state_t *));
  crf1de_state_t **pky2pk = calloc(sm->m_num_bkw, sizeof(crf1de_state_t *));
  if (pky2ky == NULL || pky2pk == NULL) {
    CLEAR(pky2ky);
    CLEAR(pky2pk);
    CLEAR(sm->m_frw_trans2);
    CLEAR(sm->m_frw_trans1);
    CLEAR(sm->m_bkwid2bkw);
    CLEAR(sm->m_frwid2frw);
    CLEAR(sm->m_frw_llabels);
    return -6;
  }
  /* in the first run, populate the `m_frw_llabels` array and
     calculate the number of prefixes */
  int y;
  const int *pk_start = NULL;
  size_t pk_id, pk_len, sfx_id, last_label;

  RUMAVL_NODE *node = NULL;
  crf1de_state_t *pk_entry = NULL, *pky_entry = NULL, *ky_entry = NULL;

  /* find longest suffixes and states */
  while ((node = rumavl_node_next(sm->m_frw_states, node, 1, (void**) &pk_entry)) != NULL) {
    pk_id = pk_entry->m_id;
    pk_len = pk_entry->m_len;
    pk_start = pk_entry->m_seq;

    sm->m_frw_llabels[pk_id] = last_label = *pk_start;
    sm->m_frwid2frw[pk_id] = pk_entry;

    /* copy label sequence to `sm->m_wrkbench1' and append to it all
       possible tags */
    sm->m_wrkbench1.m_len = pk_len + 1;
    memcpy((void *) &sm->m_wrkbench1.m_seq[1], (const void *) pk_start, pk_len * sizeof(int));
    for (y = 0; y < sm->L; ++y) {
      if (y == last_label && sm->m_seg_len_lim < 0)
	continue;

      sm->m_wrkbench1.m_seq[0] = y;
      pky_entry = (crf1de_state_t *) rumavl_find(sm->m_bkw_states, &sm->m_wrkbench1);
      ky_entry = semimarkov_find_max_sfx(sm->m_frw_states, &sm->m_wrkbench1);

      /* increment the number of prefixes corresponding to the given
	 suffix */
      ++ky_entry->m_n_prefixes;
      /* store which `pk' and `ky' entries correspond to the given
	 `pky' state */
      pky2ky[pky_entry->m_id] = ky_entry;
      pky2pk[pky_entry->m_id] = pk_entry;
    }
  }

  /* populate `m_bkwid2bkw` */
  node = NULL;
  while ((node = rumavl_node_next(sm->m_bkw_states, node, 1, (void**) &pk_entry)) != NULL) {
    sm->m_bkwid2bkw[pk_entry->m_id] = pk_entry;
  }

  /* populate forward transitions on the basis of previously computed
     `pky2ky' and `pky2pk' */
  node = NULL;
  int prev_id = -1;
  crf1de_state_t **frw_trans1 = sm->m_frw_trans1, **frw_trans2 = sm->m_frw_trans2;
  while ((node = rumavl_node_next(sm->m_frw_states, node, 1, (void**) &pk_entry)) != NULL) {
    pk_entry->m_frw_trans1 = frw_trans1;
    pk_entry->m_frw_trans2 = frw_trans2;

    frw_trans1 += pk_entry->m_n_prefixes;
    frw_trans2 += pk_entry->m_n_prefixes;
  }

  /* populate forward transitions of the states */
  for (size_t i = 0; i < sm->m_num_bkw; ++i) {
    if (! pky2ky[i])
      continue;

    ky_entry = pky2ky[i];
    pk_entry = pky2pk[i];
    pky_entry = sm->m_bkwid2bkw[i];

    ky_entry->m_frw_trans1[ky_entry->m__cnt_trans1++] = pk_entry;
    ky_entry->m_frw_trans2[ky_entry->m__cnt_trans2++] = pky_entry;
  }
  free(pky2ky);
  free(pky2pk);

  return 0;
}

/**
 * Generate suffixes for given state and add them to semi-markov model.
 *
 * @param sm - pointer to semi-markov model
 * @param pky_entry - backward state for which all suffixes should be generated
 * @param ptrn_entry - pattern state corresponding to pky state
 *
 * @return \c void
 */
static void semimarkov_build_suffixes(crf1de_semimarkov_t *sm, crf1de_state_t *pky_entry, \
				      crf1de_state_t *ptrn_entry)
{
  size_t pky_id = pky_entry->m_id;
  size_t pky_len = pky_entry->m_len;
  crf1de_state_t **sfxp = &SUFFIXES(sm, pky_id, 0), *ptrnp = NULL;

  for (size_t i = 0; i < pky_len; ++i) {
    pky_entry->m_len = pky_len - i;
    ptrnp = rumavl_find(sm->m_patterns, pky_entry);
    if (ptrnp) {
      *sfxp = ptrnp;
      ++sfxp;
      ++ptrnp->m_n_prefixes;
      ++sm->m_num_ptrn_suffixes;
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
  sm->m_bkw_trans = calloc(sm->m_num_bkw, sm->L * sizeof(crf1de_state_t *));
  if (sm->m_bkw_trans == NULL)
    return -1;

  sm->m_bkwid2frw = calloc(sm->m_num_bkw, sizeof(crf1de_state_t *));
  if (sm->m_bkwid2frw == NULL) {
    CLEAR(sm->m_bkw_trans);
    return -2;
  }

  sm->m_suffixes = calloc(sm->m_num_bkw, (sm->m_max_order + 1) * sizeof(crf1de_state_t *));
  if (sm->m_suffixes == NULL) {
    CLEAR(sm->m_bkwid2frw);
    CLEAR(sm->m_bkw_trans);
    return -3;
  }

  sm->m_ptrnid2bkw = calloc(sm->m_num_ptrns, sizeof(crf1de_state_t *));
  if (sm->m_ptrnid2bkw == NULL) {
    CLEAR(sm->m_suffixes);
    CLEAR(sm->m_bkwid2frw);
    CLEAR(sm->m_bkw_trans);
    return -4;
  }

  sm->m_ptrn_llabels = calloc(sm->m_num_ptrns, sizeof(int));
  if (sm->m_ptrn_llabels == NULL) {
    CLEAR(sm->m_ptrnid2bkw);
    CLEAR(sm->m_suffixes);
    CLEAR(sm->m_bkwid2frw);
    CLEAR(sm->m_bkw_trans);
    return -5;
  }

  size_t pky_len, ptrn_id;
  int y, last_label, *pky_seq = NULL;
  RUMAVL_NODE *node = NULL;
  crf1de_state_t *pky_entry = NULL, *ptrn_entry = NULL, **bkw_trans = sm->m_bkw_trans;

  while ((node = rumavl_node_next(sm->m_bkw_states, node, 1, (void**) &pky_entry)) != NULL) {
    pky_len = pky_entry->m_len;
    pky_seq = pky_entry->m_seq;
    last_label = *pky_seq;

    /* generate suffixes */
    if (ptrn_entry = rumavl_find(sm->m_patterns, pky_entry)) {
      ptrn_id = ptrn_entry->m_id;
      sm->m_ptrnid2bkw[ptrn_id] = pky_entry;
      sm->m_ptrn_llabels[ptrn_id] = last_label;
    }
    semimarkov_build_suffixes(sm, pky_entry, ptrn_entry);

    /* add mapping from backward state to its corresponding forward state */
    sm->m_wrkbench1.m_len = pky_len - 1;
    memcpy((void *) sm->m_wrkbench1.m_seq, (const void *) &pky_seq[1], sm->m_wrkbench1.m_len * sizeof(int));
    sm->m_bkwid2frw[pky_entry->m_id] = rumavl_find(sm->m_patterns, &sm->m_wrkbench1);

    /* generate backward transitions for all possible labels */
    pky_entry->m_bkw_trans = bkw_trans;
    bkw_trans += sm->L;

    sm->m_wrkbench1.m_len = pky_len + 1;
    memcpy((void *) &sm->m_wrkbench1.m_seq[1], (const void *) pky_seq, pky_len * sizeof(int));

    for (y = 0; y < sm->L; ++y) {
      if (y == last_label && sm->m_seg_len_lim < 0) {
	pky_entry->m_bkw_trans[y] = NULL;
      } else {
	sm->m_wrkbench1.m_seq[0] = y;
	pky_entry->m_bkw_trans[y] = semimarkov_find_max_sfx(sm->m_bkw_states, &sm->m_wrkbench1);
      }
    }
  }

  /* allocate space for pattern transitions */
  sm->m_ptrn_trans1 = calloc(sm->m_num_ptrn_suffixes, sizeof(crf1de_state_t *));
  if (sm->m_ptrn_trans1 == NULL) {
    CLEAR(sm->m_ptrn_llabels);
    CLEAR(sm->m_ptrnid2bkw);
    CLEAR(sm->m_suffixes);
    CLEAR(sm->m_bkw_trans);
    return -6;
  }

  sm->m_ptrn_trans2 = calloc(sm->m_num_ptrn_suffixes, sizeof(crf1de_state_t *));
  if (sm->m_ptrn_trans2 == NULL) {
    CLEAR(sm->m_ptrn_trans1);
    CLEAR(sm->m_ptrn_llabels);
    CLEAR(sm->m_ptrnid2bkw);
    CLEAR(sm->m_suffixes);
    CLEAR(sm->m_bkw_trans);
    return -7;
  }

  /* assign addresses of prefix arrays to patterns */
  node = NULL;
  crf1de_state_t **ptrn_trans1 = sm->m_ptrn_trans1, **ptrn_trans2 = sm->m_ptrn_trans2;
  while ((node = rumavl_node_next(sm->m_patterns, node, 1, (void**) &ptrn_entry)) != NULL) {
    ptrn_entry->m_frw_trans1 = ptrn_trans1;
    ptrn_trans1 += ptrn_entry->m_n_prefixes;

    ptrn_entry->m_frw_trans2 = ptrn_trans2;
    ptrn_trans2 += ptrn_entry->m_n_prefixes;
  }

  /* populate patterns with prefixes */
  size_t i;
  node = NULL;
  crf1de_state_t **suffixes;

  for (size_t pky_id = 0; pky_id < sm->m_num_bkw; ++pky_id) {
    pky_entry = sm->m_bkwid2bkw[pky_id];
    suffixes = &SUFFIXES(sm, pky_id, 0);

    for (i = 0; ptrn_entry = suffixes[i]; ++i) {
      ptrn_entry->m_frw_trans1[ptrn_entry->m__cnt_trans1++] = sm->m_bkwid2frw[pky_entry->m_id];
      ptrn_entry->m_frw_trans2[ptrn_entry->m__cnt_trans2++] = pky_entry;
    }
  }
  return 0;
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
  while (a_wrkbench->m_len > 1 && rumavl_find(sm->m_patterns, a_wrkbench) == NULL) {
    a_wrkbench->m_id = sm->m_num_ptrns++;
    rumavl_insert(sm->m_patterns, a_wrkbench);

    /* reduce pattern for the next loop */
    --a_wrkbench->m_len;
  }
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
  semimarkov_output_state(stderr, "semimarkov_add_bkw_states: input state", a_wrkbench);

  int last_label = a_wrkbench->m_seq[1];

  /* append all possible tags to given prefix */
  for (int i = 0; i < sm->L; ++i) {
    if (i == last_label && sm->m_seg_len_lim < 0)
      continue;

    /* we assume that backward state is not known */
    a_wrkbench->m_id = sm->m_num_bkw++;
    a_wrkbench->m_seq[0] = i;
    rumavl_insert(sm->m_bkw_states, a_wrkbench);
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
  sm->m_wrkbench2.m_len = a_wrkbench->m_len - 1;
  memcpy((void *) &sm->m_wrkbench2.m_seq, (const void *) &a_wrkbench->m_seq[1], \
	 sm->m_wrkbench2.m_len * sizeof(int));

  while (sm->m_wrkbench2.m_len > 1 && rumavl_find(sm->m_frw_states, &sm->m_wrkbench2) == NULL) {
    sm->m_wrkbench2.m_id = sm->m_num_frw++;
    rumavl_insert(sm->m_frw_states, &sm->m_wrkbench2);
    semimarkov_add_bkw_states(sm, a_wrkbench);

    --sm->m_wrkbench2.m_len;
    memmove((void *) &sm->m_wrkbench2.m_seq, (const void *) &sm->m_wrkbench2.m_seq[1], \
	    sm->m_wrkbench2.m_len * sizeof(int));

    --a_wrkbench->m_len;
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
  sm->m_wrkbench1.m_len = 1;
  sm->m_wrkbench2.m_len = 2;

  for (int i = 0; i < sm->L; ++i) {
    sm->m_wrkbench2.m_seq[1] = sm->m_wrkbench1.m_seq[0] = i;

    /* add label to the set of prefixes (forward states) */
    sm->m_wrkbench1.m_id = sm->m_num_frw++;
    rumavl_insert(sm->m_frw_states, &sm->m_wrkbench1);

    /* add label to the set of backward states and patterns */
    sm->m_wrkbench1.m_id = sm->m_num_ptrns++;
    rumavl_insert(sm->m_patterns, &sm->m_wrkbench1);

    sm->m_wrkbench1.m_id = sm->m_num_bkw++;
    rumavl_insert(sm->m_bkw_states, &sm->m_wrkbench1);

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
  sm->m_max_order = a_max_order + 1;
  if (sm->m_max_order > CRFSUITE_SM_MAX_PTRN_LEN) {
    fprintf(stderr, "Max order (%z) exceeds limit (%z). ", \
	    sm->m_max_order, CRFSUITE_SM_MAX_PTRN_LEN);
    fprintf(stderr, "Set to %z. ", CRFSUITE_SM_MAX_PTRN_LEN);
    fprintf(stderr, "To increase the limit, recompile the program with the option -DCRFSUITE_SM_MAX_PTRN_LEN=NEW_LIM added to CPPFLAGS.\n");
    sm->m_max_order = CRFSUITE_SM_MAX_PTRN_LEN;
  }
  sm->m_seg_len_lim = a_seg_len_lim;

  /* allocate memory for storing maximum segment lengths */
  sm->m_max_seg_len = calloc(L, sizeof(int));
  if (sm->m_max_seg_len == NULL)
    return -1;

  /* allocate sets for patterns, forward and backward states */
  size_t fs_size = sizeof(crf1de_state_t) - sizeof(int) * \
    (CRFSUITE_SM_MAX_PTRN_LEN - sm->m_max_order);
  size_t bs_size = fs_size + sizeof(int);

  sm->m_num_ptrns = 0;
  sm->m_patterns = rumavl_new(bs_size, crf1de_cmp_lseq, NULL, NULL);

  sm->m_num_frw = 0;
  sm->m_frw_states = rumavl_new(fs_size, crf1de_cmp_lseq, NULL, NULL);

  sm->m_num_bkw = 0;
  sm->m_bkw_states = rumavl_new(bs_size, crf1de_cmp_lseq, NULL, NULL);

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

  /* transfer tag sequence from ring to workbench */
  crfsuite_chain_link_t *chlink = sm->m_ring->tail->prev;
  memset(&sm->m_wrkbench1, -1, sizeof(crf1de_state_t));
  sm->m_wrkbench1.m_len = 0;

  for (size_t i = 0; i < sm->m_ring->num_items; ++i) {
    ++sm->m_wrkbench1.m_len;
    sm->m_wrkbench1.m_seq[i] = chlink->data;
    chlink = chlink->prev;
  }
  /* generate all possible prefixes and multiply them by L */
  if (rumavl_find(sm->m_patterns, &sm->m_wrkbench1) == NULL) {
    semimarkov_add_patterns(sm, &sm->m_wrkbench1);
    semimarkov_add_states(sm, &sm->m_wrkbench1);
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

  /* generate forward transitions */
  if (semimarkov_build_frw_transitions(sm))
    return -1;

  /* generate backward transitions */
  if (semimarkov_build_bkw_transitions(sm)) {
    CLEAR(sm->m_frw_trans1);
    CLEAR(sm->m_frw_trans2);
    CLEAR(sm->m_bkwid2bkw);
    CLEAR(sm->m_frwid2frw);
    CLEAR(sm->m_frw_llabels);
    return -2;
  }

  /* generate pattern transitions */
  semimarkov_debug_transitions(sm);

  /* clear ring */
  if (sm->m_ring) {
    sm->m_ring->free(sm->m_ring);
    free(sm->m_ring);
    sm->m_ring = NULL;
  }
  exit(66);
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
  /* clear patterns */
  RUMAVL_CLEAR(sm->m_patterns);
  sm->m_num_ptrns = 0;

  CLEAR(sm->m_ptrn_llabels);
  CLEAR(sm->m_ptrn_trans1);
  CLEAR(sm->m_ptrn_trans2);
  CLEAR(sm->m_ptrnid2bkw);

  CLEAR(sm->m_suffixes);
  sm->m_num_ptrn_suffixes = 0;

  /* clear forward states */
  RUMAVL_CLEAR(sm->m_frw_states);
  sm->m_num_frw = 0;

  CLEAR(sm->m_frw_llabels);
  CLEAR(sm->m_frw_trans1);
  CLEAR(sm->m_frw_trans2);
  CLEAR(sm->m_frwid2frw);

  /* clear backward states */
  RUMAVL_CLEAR(sm->m_bkw_states);
  sm->m_num_bkw = 0;

  CLEAR(sm->m_bkw_trans);
  CLEAR(sm->m_bkwid2bkw);
  CLEAR(sm->m_bkwid2frw);

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
  sm->finalize = semimarkov_finalize;
  sm->clear = semimarkov_clear;

  return sm;
}
