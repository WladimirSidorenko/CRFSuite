/*
 *      CRF1d context (forward-backward, viterbi, etc).
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

#ifdef    HAVE_CONFIG_H
#include <config.h>
#endif/*HAVE_CONFIG_H*/

#include <assert.h>
#include <os.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <crfsuite.h>

#include "crf1d.h"
#include "vecmath.h"

crf1d_context_t* crf1dc_new(int flag, const int ftype, int L, int T, const crf1de_semimarkov_t *sm)
{
  int ret = 0;
  crf1d_context_t* ctx = NULL;
  int n_src_tags = ftype == FTYPE_SEMIMCRF? sm->m_num_frw: L;

  ctx = (crf1d_context_t*) calloc(1, sizeof(crf1d_context_t));
  if (ctx == NULL)
    return ctx;

  ctx->ftype = ftype;
  ctx->flag = flag;
  ctx->num_labels = L;

  ctx->trans = (floatval_t*) calloc(n_src_tags * L, sizeof(floatval_t));
  if (ctx->trans == NULL) goto error_exit;

  if (ctx->flag & CTXF_MARGINALS) {
    ctx->exp_trans = (floatval_t*)_aligned_malloc((n_src_tags * L + 4) * sizeof(floatval_t), 16);
    if (ctx->exp_trans == NULL) goto error_exit;
    ctx->mexp_trans = (floatval_t*) calloc(n_src_tags * L, sizeof(floatval_t));
    if (ctx->mexp_trans == NULL) goto error_exit;
  }

  if ((ret = crf1dc_set_num_items(ctx, sm, T)))
    goto error_exit;

  /* T gives the 'hint' for maximum length of items. */
  ctx->num_items = 0;

  return ctx;

 error_exit:
  crf1dc_delete(ctx);
  return NULL;
}

int crf1dc_set_num_items(crf1d_context_t* ctx, const crf1de_semimarkov_t *sm, const int T)
{
  const int L = ctx->num_labels;

  ctx->num_items = T;

  if (ctx->cap_items < T) {
    free(ctx->backward_edge);
    free(ctx->mexp_state);
    _aligned_free(ctx->exp_state);
    free(ctx->scale_factor);
    free(ctx->row);
    free(ctx->child_row);
    free(ctx->beta_score);
    free(ctx->alpha_score);
    free(ctx->child_alpha_score);

    /* transition feature vectors will look differently for semimarkov model */
    int n_alpha_states = L, n_beta_states = L;
    if (ctx->ftype == FTYPE_SEMIMCRF) {
      n_alpha_states = sm->m_num_frw;
      n_beta_states = sm->m_num_bkw;
    }

    fprintf(stderr, "crf1dc_set_num_items: allocating T (%d) * n_alpha_states (%d) for alpha\n", \
	    T, n_alpha_states);
    ctx->alpha_score = (floatval_t*)calloc(T * n_alpha_states, sizeof(floatval_t));
    fprintf(stderr, "crf1dc_set_num_items: ctx->alpha_score = %p\n", ctx->alpha_score);
    if (ctx->alpha_score == NULL) return CRFSUITEERR_OUTOFMEMORY;

    ctx->beta_score = (floatval_t*)calloc(T * n_beta_states, sizeof(floatval_t));
    if (ctx->beta_score == NULL) return CRFSUITEERR_OUTOFMEMORY;

    ctx->row = (floatval_t*)calloc(n_beta_states, sizeof(floatval_t));
    if (ctx->row == NULL) return CRFSUITEERR_OUTOFMEMORY;

    if (ctx->ftype == FTYPE_CRF1TREE) {
      ctx->child_alpha_score = (floatval_t*) calloc(T * n_alpha_states, sizeof(floatval_t));
      if (ctx->child_alpha_score == NULL) return CRFSUITEERR_OUTOFMEMORY;

      ctx->child_row = (floatval_t*)calloc(L, sizeof(floatval_t));
      if (ctx->child_row == NULL) return CRFSUITEERR_OUTOFMEMORY;
    }

    if (ctx->flag & CTXF_VITERBI) {
      ctx->backward_edge = (int*)calloc(T * n_beta_states, sizeof(int));
      if (ctx->backward_edge == NULL) return CRFSUITEERR_OUTOFMEMORY;
    }

    ctx->scale_factor = (floatval_t*)calloc(T, sizeof(floatval_t));
    if (ctx->scale_factor == NULL) return CRFSUITEERR_OUTOFMEMORY;

    ctx->state = (floatval_t*)calloc(T * L, sizeof(floatval_t));
    if (ctx->state == NULL) return CRFSUITEERR_OUTOFMEMORY;

    if (ctx->flag & CTXF_MARGINALS) {
      ctx->exp_state = (floatval_t*)_aligned_malloc((T * L + 4) * sizeof(floatval_t), 16);
      if (ctx->exp_state == NULL) return CRFSUITEERR_OUTOFMEMORY;
      ctx->mexp_state = (floatval_t*)calloc(T * L, sizeof(floatval_t));
      if (ctx->mexp_state == NULL) return CRFSUITEERR_OUTOFMEMORY;
    }
    ctx->cap_items = T;
  }
  return 0;
}

void crf1dc_delete(crf1d_context_t* ctx)
{
  if (ctx != NULL) {
    free(ctx->backward_edge);
    free(ctx->mexp_state);
    _aligned_free(ctx->exp_state);
    free(ctx->state);
    free(ctx->scale_factor);
    free(ctx->row);
    free(ctx->child_row);
    free(ctx->beta_score);
    free(ctx->alpha_score);
    free(ctx->child_alpha_score);
    free(ctx->mexp_trans);
    _aligned_free(ctx->exp_trans);
    free(ctx->trans);
  }
  free(ctx);
}

void crf1dc_reset(crf1d_context_t* ctx, int flag, const crf1de_semimarkov_t *sm)
{
  const int ftype = ctx->ftype;
  const int T = ctx->num_items;
  const int L = ctx->num_labels;

  int n_states = ftype == FTYPE_SEMIMCRF? sm->m_num_frw: L;

  if (flag & RF_STATE)
    veczero(ctx->state, n_states * T);

  if (flag & RF_TRANS)
    veczero(ctx->trans, n_states * L);

  if (ctx->flag & CTXF_MARGINALS) {
    veczero(ctx->mexp_state, n_states * T);
    veczero(ctx->mexp_trans, n_states * L);
    ctx->log_norm = 0;
  }
}

void crf1dc_exp_state(crf1d_context_t* ctx)
{
  const int T = ctx->num_items;
  const int L = ctx->num_labels;

  veccopy(ctx->exp_state, ctx->state, L * T);
  vecexp(ctx->exp_state, L * T);
}

void crf1dc_exp_transition(crf1d_context_t* ctx, const crf1de_semimarkov_t *sm)
{
  const int L = sm ? sm->m_num_frw: ctx->num_labels;

  veccopy(ctx->exp_trans, ctx->trans, L * L);
  vecexp(ctx->exp_trans, L * L);
}

void crf1dc_alpha_score(crf1d_context_t* a_ctx,  const void *a_aux)
{
  int i, t;
  floatval_t sum, *cur = NULL, *scale = &a_ctx->scale_factor[0];
  const floatval_t *prev = NULL, *trans = NULL, *state = NULL;
  const int T = a_ctx->num_items;
  const int L = a_ctx->num_labels;
  /* Compute alpha scores on leaves (0, *).
     alpha[0][j] = state[0][j]
  */
  cur = ALPHA_SCORE(a_ctx, 0);
  state = EXP_STATE_SCORE(a_ctx, 0);
  // copy L elements from state to current
  veccopy(cur, state, L);
  for (int i = 0; i < L; ++i) {
    fprintf(stderr, "unscaled alpha[0][%d] = %f\n", i, cur[i]);
  }
  // total sum of L elements in vector
  sum = vecsum(cur, L);
  *scale = (sum != 0.) ? 1. / sum : 1.;
  fprintf(stderr, "scale[0] = %f\n", *scale);
  // multiply L elements in cur by scale factor (i.e. normalize weights)
  vecscale(cur, *scale, L);
  for (int i = 0; i < L; ++i) {
    fprintf(stderr, "scaled alpha[0][%d] = %f\n", i, cur[i]);
  }
  ++scale;

  /* Compute the alpha scores on nodes (t, *).
     alpha[t][j] = state[t][j] * \sum_{i} alpha[t-1][i] * trans[i][j]
  */
  for (t = 1;t < T;++t) {
    prev = ALPHA_SCORE(a_ctx, t-1);
    cur = ALPHA_SCORE(a_ctx, t);
    state = EXP_STATE_SCORE(a_ctx, t);

    veczero(cur, L);
    for (i = 0;i < L;++i) {
      trans = EXP_TRANS_SCORE(a_ctx, i);
      // for each cell j of cur add score prev[i] multiplied by
      // transition score i -> j
      vecaadd(cur, prev[i], trans, L);
    }
    // memberwise multiplication of values in cur and values in state
    vecmul(cur, state, L);
    sum = vecsum(cur, L);
    *scale = (sum != 0.) ? 1. / sum : 1.;
    fprintf(stderr, "scale[%d] = %f\n", t, *scale);
    // normalize weights
    for (int i = 0; i < L; ++i) {
      fprintf(stderr, "unscaled alpha[%d][%d] = %f\n", t, i, cur[i]);
    }
    vecscale(cur, *scale, L);
    for (int i = 0; i < L; ++i) {
      fprintf(stderr, "scaled alpha[%d][%d] = %f\n", t, i, cur[i]);
    }
    ++scale;
  }

  /* Compute the logarithm of the normalization factor here.
     norm = 1. / (C[0] * C[1] ... * C[T-1])
     log(norm) = - \sum_{t = 0}^{T-1} log(C[t]).
  */

  // sum logarithms of all elements in scale factor
  a_ctx->log_norm = -vecsumlog(a_ctx->scale_factor, T);
}

/**
 * Compute alpha score for tree structured CRF.
 *
 * The score will be computed from leaves up to the root.
 *
 * @param a_ctx - gm context for which the score should be computed
 * @param a_aux - auxiliary data-structure (in this case, tree)
 **/
void crf1dc_tree_alpha_score(crf1d_context_t* a_ctx, const void *a_aux)
{
  const crfsuite_node_t *tree = (const crfsuite_node_t *) a_aux;
  assert(tree);

  int i, c, t;
  int item_id, chld_item_id;
  floatval_t sum, *crnt_alpha, *chld_alpha, *chld_alpha_score;
  floatval_t *scale, *row = a_ctx->row;
  const int T = a_ctx->num_items;
  const int L = a_ctx->num_labels;
  const floatval_t *trans, *state;
  const crfsuite_node_t *node, *child;

  /* Since nodes in the tree are ordered topologically, we start at
     the end of the tree where leaf nodes are located.  Once we are
     done with the leaves, we go up and compute the probabilities of
     nodes whose children probabilities have already been computed.

     The probability of a node computed in this way is the exponent of
     its state probablities times the product of the probability
     scores of its children:

     alpha[t][j] = state[t][j] * \sum_{i \in L} trans[i][j] * \sum_{c
     \in Ch_t} alpha[c][i]
  */
  for (t = T - 1; t >= 0; --t) {
    // node corresponding to given item
    node = &tree[t];
    // id of item corresponding to given node
    item_id = node->self_item_id;
    // column for storing transition weights
    crnt_alpha = ALPHA_SCORE(a_ctx, item_id);
    // state weigts for that item
    state = EXP_STATE_SCORE(a_ctx, item_id);
    // slot for scaling factor
    scale = &a_ctx->scale_factor[item_id];

    // if current node is a leaf, set it transition weights to state weights
    if (node->num_children == 0) {
      veccopy(crnt_alpha, state, L);
    } else {
      // set weights for current node to zeros
      vecset(crnt_alpha, 1., L);
      // compute the total sum of all alpha score weights for children
      for (c = 0; c < node->num_children; ++c) {
	// obtain the child item
	child = &tree[node->children[c]];
	chld_item_id = child->self_item_id;
	chld_alpha = ALPHA_SCORE(a_ctx, chld_item_id);
	/* alpha score propagated from this child to parent */
	chld_alpha_score = CHILD_ALPHA_SCORE(a_ctx, chld_item_id);
	// for each label `i', we compute the sum of the weights for i-th labels
	// in all of the children, then we set the weight for each target tag j
	// in the current column to the sum of the i-th weights * the transition
	// weight i --> j
	veczero(row, L);
	for (i = 0; i < L; ++i) {
	  trans = EXP_TRANS_SCORE(a_ctx, i);
	  // for each cell j of cur add score prev[i] multiplied by
	  // transition score i -> j
	  vecaadd(row, chld_alpha[i], trans, L);
	}
	// store computed alpha score in `child_alpha_score`
	veccopy(chld_alpha_score, row, L);
	vecinv(chld_alpha_score, L);
	vecmul(crnt_alpha, row, L);
      }
      // multiply all obtained transition weights with state weights for that
      // node
      vecmul(crnt_alpha, state, L);
    }
    // normalize weights
    sum = vecsum(crnt_alpha, L);
    *scale = (sum != 0.) ? 1. / sum : 1.;
    // multiply L elements in cur by scale factor (i.e. normalize weights)
    vecscale(crnt_alpha, *scale, L);
  }
  // sum logarithms of all elements in scale factor
  a_ctx->log_norm = -vecsumlog(a_ctx->scale_factor, T);
}

/**
 * Compute alpha score for semi-markov CRF (of possibly higher order).
 *
 * The score will be computed for all possible segments.
 *
 * @param a_ctx - gm context for which the score should be computed
 * @param a_aux - auxiliary data structure (in this case semi-markov model)
 **/
void crf1dc_sm_alpha_score(crf1d_context_t* a_ctx, const void *a_aux)
{
  /* Scaling trick will not work here */
  const int T = a_ctx->num_items;
  const crf1de_semimarkov_t *sm = (const crf1de_semimarkov_t *) a_aux;
  /* Compute alpha scores on leaves (0, *).
     alpha[0][j] = state[0][j]
  */
  floatval_t *cur = SM_ALPHA_SCORE(a_ctx, sm, 0);
  veczero(cur, sm->m_num_frw);
  const floatval_t *exp_state_score = EXP_STATE_SCORE(a_ctx, 0);

  int j, y;
  crf1de_state_t *frw_state = NULL;
  for (j = 0; j < sm->m_num_frw; ++j) {
    frw_state = &sm->m_frw_states[j];

    if (frw_state->m_len == 1) {
      y = sm->m_frw_llabels[j];
      cur[j] = exp_state_score[y];
    } else if (frw_state->m_len > 1)
      break;
  }

  /* Compute alpha scores on nodes (t, *).
     alpha[t][j] = \sum_{s = t}^{s - max_seg_len[j]} \prod_{s}^{t} state[s][j] * \
     \sum_{i \in j's prefixes} alpha[s-1][i] * trans[i][j]
  */
  const floatval_t *prev = NULL;
  floatval_t state_score = 0., trans_score = 0.;
  const int *frw_trans1 = NULL, *frw_trans2 = NULL, *suffixes;
  int i, k, seg_start, min_seg_start, prev_seg_end, prev_id1, prev_id2, sfx_id;

  for (int t = 1; t < T;++t) {
    cur = SM_ALPHA_SCORE(a_ctx, sm, t);
    veczero(cur, sm->m_num_frw);

    for (j = 0; j < sm->m_num_frw; ++j) {
      /* obtain semi-markov state, corresponding to #i-th index */
      frw_state = &sm->m_frw_states[j];
      /* obtain possible transitions for that semi-markov state */
      frw_trans1 = frw_state->m_frw_trans1;
      frw_trans2 = frw_state->m_frw_trans2;
      /* get last label and obtain maximum length of a span with that label */
      y = sm->m_frw_llabels[j];
      if (y < 0)
	continue;

      min_seg_start = t - sm->m_max_seg_len[y];
      if (min_seg_start < 0)
	min_seg_start = -1;

      /* iterate over all possible previous states in the range [t -
	 max_seg_len, t) and compute the transition scores */
      state_score = 1.;
      for (seg_start = t; seg_start > min_seg_start; --seg_start) {
	prev_seg_end = seg_start - 1;
	state_score *= (EXP_STATE_SCORE(a_ctx, seg_start))[y];

	if (prev_seg_end < 0) {
	  if (sm->m_frw_states[j].m_len == 1)
	    cur[j] += state_score;
	  /* else if (sm->m_frw_states[j].m_len > 1) */
	  /*   break; */
	} else {
	  prev = SM_ALPHA_SCORE(a_ctx, sm, prev_seg_end);
	  for (i = 0; i < frw_state->m_num_affixes; ++i) {
	    prev_id1 = frw_trans1[i];
	    prev_id2 = frw_trans2[i];
	    suffixes = &SUFFIXES(sm, prev_id2, 0);
	    trans_score = 1.;
	    for (k = 0; (sfx_id = suffixes[k]) >= 0; ++k) {
	      trans_score *= *(EXP_TRANS_SCORE(a_ctx, sm->m_ptrns[sfx_id].m_feat_id));
	    }
	    cur[j] += prev[prev_id1] * trans_score * state_score;
	  }
	}
      }
      fprintf(stderr, "crf1dc_sm_alpha_score: alpha[%d][%d (", t, j);
      sm->output_state(stderr, NULL, &sm->m_frw_states[j]);
      fprintf(stderr, ")] = %f\n", cur[j]);
    }
  }
  // sum logarithms of all elements in scale factor
  /* a_ctx->log_norm = -vecsumlog(a_ctx->scale_factor, T); */
  a_ctx->log_norm = log(vecsum(cur, sm->m_num_frw));
  fprintf(stderr, "crf1dc_sm_alpha_score: a_ctx->log_norm = %f\n", a_ctx->log_norm);
}

void crf1dc_beta_score(crf1d_context_t* a_ctx, const void *a_aux)
{
  int i, t;
  floatval_t *cur = NULL;
  floatval_t *row = a_ctx->row;
  const floatval_t *next = NULL, *state = NULL, *trans = NULL;
  const int T = a_ctx->num_items;
  const int L = a_ctx->num_labels;
  const floatval_t *scale = &a_ctx->scale_factor[T-1];

  /* Compute the beta scores at (T-1, *). */
  cur = BETA_SCORE(a_ctx, T-1);
  // set all elements of cur to *scale
  vecset(cur, *scale, L);
  --scale;

  /* Compute the beta scores at (t, *). */
  for (t = T-2; 0 <= t; --t) {
    cur = BETA_SCORE(a_ctx, t);
    next = BETA_SCORE(a_ctx, t+1);
    state = EXP_STATE_SCORE(a_ctx, t+1);

    veccopy(row, next, L);
    vecmul(row, state, L);

    /* Compute the beta score at (t, i). */
    for (i = 0;i < L;++i) {
      trans = EXP_TRANS_SCORE(a_ctx, i);
      cur[i] = vecdot(trans, row, L);
    }
    vecscale(cur, *scale, L);
    --scale;
  }
}

/**
 * Compute beta score for tree-structured CRF.
 *
 * The score will be computed from root down to the leaves.
 *
 * @param a_ctx - gm context for which the score should be computed
 * @param a_aux - pointer to tree
 **/
void crf1dc_tree_beta_score(crf1d_context_t* a_ctx, const void *a_aux)
{
  const crfsuite_node_t *tree = (const crfsuite_node_t *) a_aux;

  int i, t;
  int item_id, prnt_item_id;
  floatval_t prnt_scale;
  floatval_t *crnt_beta, *row = a_ctx->row;
  const floatval_t *prnt_alpha, *prnt_beta, *prnt_state;
  const floatval_t *chld_alpha_score, *scale, *trans;
  const crfsuite_node_t *node, *prnt_node;
  const int T = a_ctx->num_items;
  const int L = a_ctx->num_labels;

  /* Compute beta score for root noode. */
  node = &tree[0];
  item_id = node->self_item_id;
  crnt_beta = BETA_SCORE(a_ctx, item_id);
  scale = &a_ctx->scale_factor[item_id];
  // set beta score for all root labels to *scale, so that it will be
  // cancelled-out at the end when computing marginals
  vecset(crnt_beta, *scale, L);

  /* Compute beta scores for children. */
  for (t = 1; t < T; ++t) {
    node = &tree[t];
    item_id = node->self_item_id;
    crnt_beta = BETA_SCORE(a_ctx, item_id);
    scale = &a_ctx->scale_factor[item_id];

    assert(node->prnt_node_id >= 0 && node->prnt_node_id < t);
    prnt_node = &tree[node->prnt_node_id];
    prnt_item_id = prnt_node->self_item_id;
    prnt_alpha = ALPHA_SCORE(a_ctx, prnt_item_id);
    prnt_beta = BETA_SCORE(a_ctx, prnt_item_id);
    prnt_state = EXP_STATE_SCORE(a_ctx, prnt_item_id);
    prnt_scale = a_ctx->scale_factor[prnt_item_id];

    if (prnt_node->num_children > 1) {
      veccopy(row, prnt_alpha, L);
      /* cancel-out scale factor */
      if (prnt_scale)
	vecscale(row, 1. / prnt_scale, L);
      /* divide parent's alpha score by the alpha score which came from this
	 child */
      chld_alpha_score = CHILD_ALPHA_SCORE(a_ctx, item_id);
      vecmul(row, chld_alpha_score, L);
    } else {
      veccopy(row, prnt_state, L);
    }
    vecmul(row, prnt_beta, L);
    // sum-out target labels
    for (i = 0; i < L; ++i) {
      trans = EXP_TRANS_SCORE(a_ctx, i);
      crnt_beta[i] = vecdot(trans, row, L);
    }
    vecscale(crnt_beta, *scale, L);
  }
}

/**
 * Compute beta score for semi-markov CRF.
 *
 * The score will be computed for all possible segments.
 *
 * @param a_ctx - gm context for which the score should be computed
 * @param a_aux - auxiliary data structure (semi-markov model in this case)
 **/
void crf1dc_sm_beta_score(crf1d_context_t* a_ctx, const void *a_aux)
{
  const int T = a_ctx->num_items;
  const crf1de_semimarkov_t *sm = (const crf1de_semimarkov_t *) a_aux;

  const floatval_t *nxt = NULL;
  floatval_t state_score = 0., trans_score = 0.;

  /* Compute beta score at nodes (t, *). */
  int j, y, pk_id, pky_id, sfx_id;
  int seg_end, max_seg_end;
  const int *suffixes = NULL, *bkw_trans = NULL;
  floatval_t *cur = SM_BETA_SCORE(a_ctx, sm, T - 1);

  for (int t = T-1; 0 < t; --t) {
    cur = SM_BETA_SCORE(a_ctx, sm, t);

    for (y = 0; y < sm->L; ++y) {
      max_seg_end = t + sm->m_max_seg_len[y];
      if (max_seg_end > T)
	max_seg_end = T;

      state_score = 1.;
      for (seg_end = t; seg_end < max_seg_end; ++seg_end) {
	state_score *= (EXP_STATE_SCORE(a_ctx, seg_end))[y];

	if (seg_end == T - 1) {
	  for (pk_id = 0; pk_id < sm->m_num_bkw; ++pk_id) {
	    bkw_trans = sm->m_bkw_states[pk_id].m_bkw_trans;
	    pky_id = bkw_trans[y];

	    if (pky_id < 0)
	      continue;

	    cur[pk_id] += state_score;
	  }
	} else {
	  nxt = SM_BETA_SCORE(a_ctx, sm, seg_end + 1);

	  for (pk_id = 0; pk_id < sm->m_num_bkw; ++pk_id) {
	    bkw_trans = sm->m_bkw_states[pk_id].m_bkw_trans;
	    pky_id = bkw_trans[y];
	    if (pky_id < 0)
	      continue;

	    trans_score = 1.;
	    suffixes = &SUFFIXES(sm, pky_id, 0);
	    for (j = 0; (sfx_id = suffixes[j]) >= 0; ++j) {
	      trans_score *= *(EXP_TRANS_SCORE(a_ctx, sm->m_ptrns[sfx_id].m_feat_id));
	    }
	    cur[pk_id] += state_score * trans_score * nxt[pky_id];
	  }
	}
      }
    }
  }

  for (int t = T-1; 0 <= t; --t) {
    cur = SM_BETA_SCORE(a_ctx, sm, t);
    for (pk_id = 0; pk_id < sm->m_num_bkw; ++pk_id) {
      fprintf(stderr, "crf1dc_sm_beta_score: beta[%d][", t);
      sm->output_state(stderr, NULL, &sm->m_bkw_states[pk_id]);
      fprintf(stderr, "] = %.4f\n", cur[pk_id]);
    }
  }
}

void crf1dc_marginals(crf1d_context_t* a_ctx, const void *a_aux)
{
  int i, j, t;
  const int T = a_ctx->num_items;
  const int L = a_ctx->num_labels;

  /*
    Compute model expectation of states.
    p(t,i) = fwd[t][i] * bwd[t][i] / norm
    = (1. / C[t]) * fwd'[t][i] * bwd'[t][i][hhh
  */
  for (t = 0;t < T;++t) {
    floatval_t *fwd = ALPHA_SCORE(a_ctx, t);
    floatval_t *bwd = BETA_SCORE(a_ctx, t);
    floatval_t *prob = STATE_MEXP(a_ctx, t);
    veccopy(prob, fwd, L);
    vecmul(prob, bwd, L);
    vecscale(prob, 1. / a_ctx->scale_factor[t], L);
  }

  /*
    Compute the model expectations of transitions.
    p(t,i,t+1,j)
    = fwd[t][i] * edge[i][j] * state[t+1][j] * bwd[t+1][j] / norm
    = (fwd'[t][i] / (C[0] ... C[t])) * edge[i][j] * state[t+1][j] * (bwd'[t+1][j] / (C[t+1] ... C[T-1])) * (C[0] * ... * C[T-1])
    = fwd'[t][i] * edge[i][j] * state[t+1][j] * bwd'[t+1][j]
    The model expectation of a transition (i -> j) is the sum of the marginal
    probabilities p(t,i,t+1,j) over t.
  */
  floatval_t *row = a_ctx->row;
  for (t = 0; t < T-1; ++t) {
    floatval_t *fwd = ALPHA_SCORE(a_ctx, t);
    floatval_t *bwd = BETA_SCORE(a_ctx, t+1);
    floatval_t *state = EXP_STATE_SCORE(a_ctx, t+1);

    /* row[j] = state[t+1][j] * bwd'[t+1][j] */
    veccopy(row, bwd, L);
    vecmul(row, state, L);

    for (i = 0;i < L; ++i) {
      floatval_t *edge = EXP_TRANS_SCORE(a_ctx, i);
      floatval_t *prob = TRANS_MEXP(a_ctx, i);
      for (j = 0; j < L; ++j) {
	prob[j] += fwd[i] * edge[j] * row[j];
      }
    }
  }
}

void crf1dc_tree_marginals(crf1d_context_t* a_ctx, const void *a_aux)
{
  const crfsuite_node_t *tree = (const crfsuite_node_t *) a_aux;

  int i, j, t, c;
  int item_id, prnt_id, chld_item_id;
  floatval_t chck_sum = 0.;
  const int T = a_ctx->num_items;
  const int L = a_ctx->num_labels;
  const crfsuite_node_t *node, *prnt_node, *chld_node;
  const floatval_t *fwd = NULL, *bwd = NULL, *state = NULL, *scale = NULL, *trans = NULL;
  floatval_t *prob = NULL;
  /*
   Compute model expectations of states (this expectation is the same for all
   types of graphical models).

    p(t,i) = fwd[t][i] * bwd[t][i] / norm
    = (1. / C[t]) * fwd'[t][i] * bwd'[t][i]
  */
  for (t = 0; t < T; ++t) {
    fwd = ALPHA_SCORE(a_ctx, t);
    bwd = BETA_SCORE(a_ctx, t);
    prob = STATE_MEXP(a_ctx, t);
    scale = &a_ctx->scale_factor[t];

    veccopy(prob, fwd, L);
    vecmul(prob, bwd, L);
    vecscale(prob, 1. / *scale, L);
    chck_sum = vecsum(prob, L) - 1;
    assert(chck_sum < 0.00001  && chck_sum > -0.00001);
  }

  /*
    Compute model expectations of transitions (these transitions will be
    different for linear and tree structured CRFs).

    p(t,i,t+1,j) = fwd[t][i] * edge[i][j] * state[t+1][j] * bwd[t+1][j] / norm
    = (fwd'[t][i] / (C[0] ... C[t])) * edge[i][j] * state[t+1][j] *
    (bwd'[t+1][j] / (C[t+1] ... C[T-1])) * (C[0] * ... * C[T-1]) = fwd'[t][i]
    * edge[i][j] * state[t+1][j] * bwd'[t+1][j] The model expectation of a
    transition (i -> j) is the sum of the marginal probabilities p(t,i,t+1,j)
    over t.
  */

  floatval_t *row = a_ctx->row, *workbench = a_ctx->child_row;
  const floatval_t *chld_alpha = NULL;

  for (t = T - 1; t > 0; --t) {
    node = &tree[t];
    item_id = node->self_item_id;

    assert(node->prnt_node_id >= 0 && node->prnt_node_id < t);
    prnt_node = &tree[node->prnt_node_id];
    prnt_id = prnt_node->self_item_id;

    fwd = ALPHA_SCORE(a_ctx, item_id);
    bwd = BETA_SCORE(a_ctx, prnt_id);
    state = EXP_STATE_SCORE(a_ctx, prnt_id);

    /* row[j] = state[t+1][j] * bwd'[t+1][j] */

    /* multiply beta score with alpha scores which came to the parent from
       other children */
    if (prnt_node->num_children > 1) {
      vecset(row, 1., L);
      /* compute alpha score which came from other children to this parent */
      for (c = 0; c < prnt_node->num_children; ++c) {
	chld_node = &tree[prnt_node->children[c]];
	chld_item_id = chld_node->self_item_id;
	if (chld_item_id == item_id)
	  continue;
	chld_alpha = ALPHA_SCORE(a_ctx, chld_item_id);
	veczero(workbench, L);
	for (i = 0; i < L; ++i) {
	  trans = EXP_TRANS_SCORE(a_ctx, i);
	  vecaadd(workbench, chld_alpha[i], trans, L);
	}
	vecmul(row, workbench, L);
      }
      vecmul(row, state, L);
    } else {
      veccopy(row, state, L);
    }
    vecmul(row, bwd, L);

    for (i = 0;i < L; ++i) {
      floatval_t *edge = EXP_TRANS_SCORE(a_ctx, i);
      floatval_t *prob = TRANS_MEXP(a_ctx, i);
      for (j = 0; j < L; ++j) {
	prob[j] += fwd[i] * edge[j] * row[j];
      }
    }
  }
}

/**
 * Compute marginals for semi-markov CRF.
 *
 * The score will be computed for all possible segments.
 *
 * @param a_ctx - gm context for which the score should be computed
 * @param a_aux - pointer to semi-markov model
 **/
void crf1dc_sm_marginals(crf1d_context_t* a_ctx, const void *a_aux)
{
  crf1de_semimarkov_t *sm = (crf1de_semimarkov_t *) a_aux;

  const int T = a_ctx->num_items;
  const int L = a_ctx->num_labels;

  /*
   * Compute model expectation of states and transitions.
   */
  const crf1de_state_t *ptrn_entry;
  int i, y, n_affixes, ptrn_id, prfx_id, sfx_id;
  floatval_t *alpha, *beta, *state_mexp;
  floatval_t *trans_mexp;

  for (int t = 0; t < T; ++t) {
    state_mexp = STATE_MEXP(a_ctx, t);
    veczero(state_mexp, L);
    alpha = ALPHA_SCORE(a_ctx, t);

    if (t < (T - 1)) {
      beta = BETA_SCORE(a_ctx, t + 1);
    } else {
      beta = NULL;
    }

    for (ptrn_id = 0; ptrn_id < sm->m_num_ptrns; ++ptrn_id) {
      y = sm->m_ptrn_llabels[ptrn_id];
      ptrn_entry = &sm->m_ptrns[ptrn_id];
      /* iterate over all pattern affixes */
      n_affixes = ptrn_entry->m_num_affixes;
      for (i = 0; i < n_affixes; ++i) {
	prfx_id = ptrn_entry->m_frw_trans1[i];
	sfx_id = ptrn_entry->m_frw_trans2[i];
	if (beta) {
	  state_mexp[y] += alpha[prfx_id] * beta[sfx_id];
	  trans_mexp = TRANS_MEXP(a_ctx, ptrn_id);
	  *trans_mexp += alpha[prfx_id] * beta[sfx_id];
	} else {
	  state_mexp[y] += alpha[prfx_id];
	}
      }
    }
  }
}

floatval_t crf1dc_marginal_point(crf1d_context_t *ctx, int l, int t)
{
  floatval_t *fwd = ALPHA_SCORE(ctx, t);
  floatval_t *bwd = BETA_SCORE(ctx, t);
  return fwd[l] * bwd[l] / ctx->scale_factor[t];
}

floatval_t crf1dc_marginal_path(crf1d_context_t *ctx, const int *path, int begin, int end)
{
  int t;
  /*
    Compute the marginal probability of a (partial) path.
    a = path[begin], b = path[begin+1], ..., y = path[end-2], z = path[end-1]
    fwd[begin][a] = (fwd'[begin][a] / (C[0] ... C[begin])
    bwd[end-1][z] = (bwd'[end-1][z] / (C[end-1] ... C[T-1]))
    norm = 1 / (C[0] * ... * C[T-1])
    p(a, b, ..., z)
    = fwd[begin][a] * edge[a][b] * state[begin+1][b] * ... * edge[y][z] * state[end-1][z] * bwd[end-1][z] / norm
    = fwd'[begin][a] * edge[a][b] * state[begin+1][b] * ... * edge[y][z] * state[end-1][z] * bwd'[end-1][z] * (C[begin+1] * ... * C[end-2])
  */
  floatval_t *fwd = ALPHA_SCORE(ctx, begin);
  floatval_t *bwd = BETA_SCORE(ctx, end-1);
  floatval_t prob = fwd[path[begin]] * bwd[path[end-1]] / ctx->scale_factor[begin];

  for (t = begin;t < end-1;++t) {
    floatval_t *state = EXP_STATE_SCORE(ctx, t+1);
    floatval_t *edge = EXP_TRANS_SCORE(ctx, path[t]);
    prob *= (edge[path[t+1]] * state[path[t+1]] * ctx->scale_factor[t]);
  }

  return prob;
}

floatval_t crf1dc_tree_marginal_path(crf1d_context_t *ctx, const int *path, int begin, int end)
{
  return 0.;
}

floatval_t crf1dc_sm_marginal_path(crf1d_context_t *ctx, const int *path, int begin, int end)
{
  return 0.;
}

#if 0
/* Sigh, this was found to be slower than the forward-backward algorithm. */

#define    ADJACENCY(ctx, i)			\
  (&MATRIX(ctx->adj, ctx->num_labels, 0, i))

void crf1dc_marginal_without_beta(crf1d_context_t* ctx)
{
  int i, j, t;
  floatval_t *prob = NULL;
  floatval_t *row = ctx->row;
  const floatval_t *fwd = NULL;
  const int T = ctx->num_items;
  const int L = ctx->num_labels;

  /*
    Compute marginal probabilities of states at T-1
    p(T-1,j) = fwd'[T-1][j]
  */
  fwd = ALPHA_SCORE(ctx, T-1);
  prob = STATE_MEXP(ctx, T-1);
  veccopy(prob, fwd, L);
  /*
    Repeat the following computation for t = T-1,T-2, ..., 1.
    1) Compute p(t-1,i,t,j) using p(t,j)
    2) Compute p(t,i) using p(t-1,i,t,j)
  */
  for (t = T-1;0 < t;--t) {
    fwd = ALPHA_SCORE(ctx, t-1);
    prob = STATE_MEXP(ctx, t);

    veczero(ctx->adj, L*L);
    veczero(row, L);

    /*
      Compute adj[i][j] and row[j].
      adj[i][j] = fwd'[t-1][i] * edge[i][j]
      row[j] = \sum_{i} adj[i][j]
    */
    for (i = 0;i < L;++i) {
      floatval_t *adj = ADJACENCY(ctx, i);
      floatval_t *edge = EXP_TRANS_SCORE(ctx, i);
      vecaadd(adj, fwd[i], edge, L);
      vecadd(row, adj, L);
    }

    /*
      Find z such that z * \sum_{i] adj[i][j] = p(t,j).
      Thus, z = p(t,j) / row[j]; we overwrite row with z.
    */
    vecinv(row, L);
    vecmul(row, prob, L);

    /*
      Apply the partition factor z (row[j]) to adj[i][j].
    */
    for (i = 0;i < L;++i) {
      floatval_t *adj = ADJACENCY(ctx, i);
      vecmul(adj, row, L);
    }

    /*
      Now that adj[i][j] presents p(t-1,i,t,j),
      accumulate model expectations of transitions.
    */
    for (i = 0;i < L;++i) {
      floatval_t *adj = ADJACENCY(ctx, i);
      floatval_t *prob = TRANS_MEXP(ctx, i);
      vecadd(prob, adj, L);
    }

    /*
      Compute the marginal probability of states at t-1.
      p(t-1,i) = \sum_{j} p(t-1,i,t,j)
    */
    prob = STATE_MEXP(ctx, t-1);
    for (i = 0;i < L;++i) {
      floatval_t *adj = ADJACENCY(ctx, i);
      prob[i] = vecsum(adj, L);
    }
  }
}
#endif

floatval_t crf1dc_score(crf1d_context_t* a_ctx, const int *a_labels, \
			const void *a_aux)
{
  int i, j, t;
  floatval_t ret = 0.;
  const floatval_t *state = NULL, *trans = NULL;
  const int T = a_ctx->num_items;

  /* Stay at (0, labels[0]). */
  i = a_labels[0];
  state = STATE_SCORE(a_ctx, 0);
  ret = state[i];

  /* Loop over the rest of items. */
  for (t = 1; t < T; ++t) {
    j = a_labels[t];
    trans = TRANS_SCORE(a_ctx, i);
    state = STATE_SCORE(a_ctx, t);

    /* Transit from (t-1, i) to (t, j). */
    ret += trans[j];
    ret += state[j];
    i = j;
  }
  return ret;
}

floatval_t crf1dc_tree_score(crf1d_context_t* a_ctx, const int *a_labels, \
			     const void *a_aux)
{
  const crfsuite_node_t *tree = (const crfsuite_node_t *) a_aux;

  int t, c;
  floatval_t score = 0., ret = 0.;
  const floatval_t *state = NULL, *trans = NULL;
  const int T = a_ctx->num_items;

  const crfsuite_node_t *node;
  int item_id, chld_node_id, chld_item_id;
  int label, chld_label;

  /* Loop over items. */
  for (t = 0; t < T; ++t) {
    /* add probability of the state */
    node = &tree[t];
    item_id = node->self_item_id;
    label = a_labels[item_id];
    state = STATE_SCORE(a_ctx, item_id);
    score = state[label];

    /* add transition probabilities from each of the children */
    for (c = 0; c < node->num_children; ++c) {
      chld_node_id = node->children[c];
      chld_item_id = tree[chld_node_id].self_item_id;
      chld_label = a_labels[chld_item_id];
      /* get transition row corresponding to the tag */
      trans = TRANS_SCORE(a_ctx, chld_label);
      /* add transition probability from child tag to parent tag */
      score += trans[label];
    }
    ret += score;
  }
  return ret;
}

floatval_t crf1dc_sm_score(crf1d_context_t* a_ctx, const int *a_labels, \
			     const void *a_aux)
{
  floatval_t ret = 0.;

  fprintf(stderr, "entered crf1dc_sm_score()\n");
  if (a_ctx->num_items == 0)
    return ret;

  crf1de_semimarkov_t *sm = (crf1de_semimarkov_t *) a_aux;
  fprintf(stderr, "obtained sm()\n");
  int semimarkov = sm->m_seg_len_lim < 0;

  /* Obtain label for 0-th element. */
  fprintf(stderr, "obtaining label\n");
  int i_label = a_labels[0];

  /* Add first label to semi-markov ring. */
  fprintf(stderr, "resetting ring\n");
  sm->m_ring->reset(sm->m_ring);
  sm->m_ring->push(sm->m_ring, i_label);

  fprintf(stderr, "ring reset\n");

  /* Obtain state score for 0-th element. */
  const floatval_t *state = STATE_SCORE(a_ctx, 0);
  floatval_t state_score = state[i_label];

  int j_label, ptrn_id, prfx_id, p, n_prefixes;
  int *ptrn_trans = NULL;
  const floatval_t *trans = NULL;
  const int T = a_ctx->num_items;
  fprintf(stderr, "T = %d\n", T);

  /* Loop over the rest of the items. */
  for (int t = 1; t < T; ++t) {
    j_label = a_labels[t];
    state = STATE_SCORE(a_ctx, t);
    fprintf(stderr, "t = %d, j_label = %d, i_label = %d, ret = %.6f\n", t, j_label, i_label, ret);
    fprintf(stderr, "state[%d] = %f\n", j_label, state[j_label]);

    if (semimarkov && j_label == i_label) {
      state_score += state[j_label]; /* TODO: not sure if it shouldn't be a multiplication */
    } else {
      /* add score for segment */
      ret += state_score;
      state_score = state[i_label];

      /* get longest possible transition pattern */
      sm->m_ring->push(sm->m_ring, j_label);
      sm->build_state(&sm->m_wrkbench1, sm->m_ring);
      while ((ptrn_id = sm->get_state_id(&sm->m_wrkbench1, sm->m__ptrns_set)) == -1 && \
	     sm->m_wrkbench1.m_len > 2) {
	/* decrement state length */
	--sm->m_wrkbench1.m_len;
	/* increment sequence pointer */
	memmove((void *) &sm->m_wrkbench1.m_seq[0], (void *) &sm->m_wrkbench1.m_seq[1],
		sm->m_wrkbench1.m_len * sizeof(int));
      }

      /* add scores for all possible transitions */
      if (ptrn_id >= 0) {
	ret += *TRANS_SCORE(a_ctx, sm->m_ptrns[ptrn_id].m_feat_id);

	/* add scores for all possible suffixes */
	n_prefixes = sm->m_ptrns[ptrn_id].m_num_affixes;
	ptrn_trans = sm->m_ptrns[ptrn_id].m_frw_trans1;
	for (p = 0; p < n_prefixes; ++p) {
	  prfx_id = ptrn_trans[p];
	  trans = TRANS_SCORE(a_ctx, prfx_id);
	  ret += trans[j_label];
	}
      }
      i_label = j_label;
    }
  }
  /* add state score of last label */
  ret += state_score;
  return ret;
}

floatval_t crf1dc_lognorm(crf1d_context_t* ctx)
{
  return ctx->log_norm;
}

floatval_t crf1dc_viterbi(crf1d_context_t* ctx, int *labels, const crfsuite_node_t *a_tree)
{
  int i, j, t;
  int *back = NULL;
  floatval_t max_score, score, *cur = NULL;
  const floatval_t *prev = NULL, *state = NULL, *trans = NULL;
  const int T = ctx->num_items;
  const int L = ctx->num_labels;

  /*
    This function assumes state and trans scores to be in the logarithm domain.
  */

  /* Compute scores at (0, *). */
  cur = ALPHA_SCORE(ctx, 0);
  state = STATE_SCORE(ctx, 0);
  for (j = 0;j < L;++j) {
    cur[j] = state[j];
  }

  /* Compute the scores at (t, *). */
  for (t = 1;t < T;++t) {
    prev = ALPHA_SCORE(ctx, t-1);
    cur = ALPHA_SCORE(ctx, t);
    state = STATE_SCORE(ctx, t);
    back = BACKWARD_EDGE_AT(ctx, t);

    /* Compute the score of (t, j). */
    for (j = 0;j < L;++j) {
      max_score = -FLOAT_MAX;

      for (i = 0;i < L;++i) {
	/* Transit from (t-1, i) to (t, j). */
	trans = TRANS_SCORE(ctx, i);
	score = prev[i] + trans[j];

	/* Store this path if it has the maximum score. */
	if (max_score < score) {
	  max_score = score;
	  /* Backward link (#t, #j) -> (#t-1, #i). */
	  back[j] = i;
	}
      }
      /* Add the state score on (t, j). */
      cur[j] = max_score + state[j];
    }
  }

  /* Find the node (#T, #i) that reaches EOS with the maximum score. */
  max_score = -FLOAT_MAX;
  prev = ALPHA_SCORE(ctx, T-1);
  for (i = 0;i < L;++i) {
    if (max_score < prev[i]) {
      max_score = prev[i];
      labels[T-1] = i;        /* Tag the item #T. */
    }
  }

  /* Tag labels by tracing the backward links. */
  for (t = T-2;0 <= t;--t) {
    back = BACKWARD_EDGE_AT(ctx, t+1);
    labels[t] = back[labels[t+1]];
  }

  /* Return the maximum score (without the normalization factor subtracted). */
  return max_score;
}

floatval_t crf1dc_tree_viterbi(crf1d_context_t* ctx, int *labels, const crfsuite_node_t *a_tree)
{
  int i, j, c, t;
  int item_id = -1, chld_item_id, lbl;
  int *back = NULL;
  floatval_t max_score = -FLOAT_MAX, score = -FLOAT_MAX;
  floatval_t *alpha = NULL, *chld_alpha = NULL;
  const floatval_t *state = NULL, *trans = NULL;
  const crfsuite_node_t *node, *child;
  const int T = ctx->num_items;
  const int L = ctx->num_labels;

  /*
    This function assumes state and trans scores to be in the logarithm domain.
  */

  /* Compute scores in bottom-up fashion (i.e. from leaves to the root). */
  for (t = T - 1; t >= 0; --t) {
    node = &a_tree[t];
    item_id = node->self_item_id;
    alpha = ALPHA_SCORE(ctx, item_id);
    state = STATE_SCORE(ctx, item_id);
    /* for leaves, alpha score will be equal to the state score */
    if (node->num_children == 0) {
      veccopy(alpha, state, L);
      /* for nodes other than leaves, we have to sum-in all the
	 ingoing messages to compute the final score for the current
	 node */
    } else {
      veczero(alpha, L);
      /* iterate over all possible target tags */
      for (j = 0; j < L; ++j) {
	/* iterate over all children */
	for (c = 0; c < node->num_children; ++c) {
	  child = &a_tree[node->children[c]];
	  chld_item_id = child->self_item_id;
	  chld_alpha = ALPHA_SCORE(ctx, chld_item_id);
	  back = BACKWARD_EDGE_AT(ctx, chld_item_id);

	  max_score = -FLOAT_MAX;
	  /* iterate over all possible source tags */
	  for (i = 0; i < L; ++i) {
	    trans = TRANS_SCORE(ctx, i);
	    score = chld_alpha[i] + trans[j];
	    if (score > max_score) {
	      max_score = score;
	      back[j] = i;
	    }
	  }
	}
	alpha[j] += max_score;
      }
      vecadd(alpha, state, L);
    }
  }

  if (item_id < 0)
    return 0.;

  /* Find label for root node which has the maximum probability. */
  max_score = -FLOAT_MAX;
  for (i = 0; i < L; ++i) {
    if (max_score < alpha[i]) {
      max_score = alpha[i];
      labels[item_id] = i;        /* Tag the item #T. */
    }
  }

  /* Find remaining labels by tracing-back tag sequence which lead to
     the most probable root tag. */
  for (t = 0; t < T; ++t) {
    node = &a_tree[t];
    item_id = node->self_item_id;
    lbl = labels[item_id];

    /* find most likely labels for children */
    for (c = 0; c < node->num_children; ++c) {
      child = &a_tree[node->children[c]];
      chld_item_id = child->self_item_id;
      back = BACKWARD_EDGE_AT(ctx, chld_item_id);
      labels[chld_item_id] = back[lbl];
    }
  }

  /* Return the maximum score (without the normalization factor subtracted). */
  return max_score;
}

floatval_t crf1dc_sm_viterbi(crf1d_context_t* ctx, int *labels, const crfsuite_node_t *a_tree)
{
  return 0.;
}

static void check_values(FILE *fp, floatval_t cv, floatval_t tv)
{
  if (fabs(cv - tv) < 1e-9) {
    fprintf(fp, "OK (%f)\n", cv);
  } else {
    fprintf(fp, "FAIL: %f (%f)\n", cv, tv);
  }
}

void crf1dc_debug_context(FILE *fp)
{
  int y1, y2, y3;
  floatval_t norm = 0;
  const int L = 3;
  const int T = 3;
  crf1d_context_t *ctx = crf1dc_new(CTXF_MARGINALS, FTYPE_CRF1D, L, T, NULL);
  floatval_t *trans = NULL, *state = NULL;
  floatval_t scores[3][3][3];
  int labels[3];

  /* Initialize the state scores. */
  state = EXP_STATE_SCORE(ctx, 0);
  state[0] = .4;    state[1] = .5;    state[2] = .1;
  state = EXP_STATE_SCORE(ctx, 1);
  state[0] = .4;    state[1] = .1;    state[2] = .5;
  state = EXP_STATE_SCORE(ctx, 2);
  state[0] = .4;    state[1] = .1;    state[2] = .5;

  /* Initialize the transition scores. */
  trans = EXP_TRANS_SCORE(ctx, 0);
  trans[0] = .3;    trans[1] = .1;    trans[2] = .4;
  trans = EXP_TRANS_SCORE(ctx, 1);
  trans[0] = .6;    trans[1] = .2;    trans[2] = .1;
  trans = EXP_TRANS_SCORE(ctx, 2);
  trans[0] = .5;    trans[1] = .2;    trans[2] = .1;

  ctx->num_items = ctx->cap_items;
  crf1dc_alpha_score(ctx, NULL);
  crf1dc_beta_score(ctx, NULL);

  /* Compute the score of every label sequence. */
  for (y1 = 0;y1 < L;++y1) {
    floatval_t s1 = EXP_STATE_SCORE(ctx, 0)[y1];
    for (y2 = 0;y2 < L;++y2) {
      floatval_t s2 = s1;
      s2 *= EXP_TRANS_SCORE(ctx, y1)[y2];
      s2 *= EXP_STATE_SCORE(ctx, 1)[y2];
      for (y3 = 0;y3 < L;++y3) {
	floatval_t s3 = s2;
	s3 *= EXP_TRANS_SCORE(ctx, y2)[y3];
	s3 *= EXP_STATE_SCORE(ctx, 2)[y3];
	scores[y1][y2][y3] = s3;
      }
    }
  }

  /* Compute the partition factor. */
  norm = 0.;
  for (y1 = 0;y1 < L;++y1) {
    for (y2 = 0;y2 < L;++y2) {
      for (y3 = 0;y3 < L;++y3) {
	norm += scores[y1][y2][y3];
      }
    }
  }

  /* Check the partition factor. */
  fprintf(fp, "Check for the partition factor... ");
  check_values(fp, exp(ctx->log_norm), norm);

  /* Compute the sequence probabilities. */
  for (y1 = 0;y1 < L;++y1) {
    for (y2 = 0;y2 < L;++y2) {
      for (y3 = 0;y3 < L;++y3) {
	floatval_t logp;

	labels[0] = y1;
	labels[1] = y2;
	labels[2] = y3;
	/* TODO: provide support for tree structured CRF score */
	logp = crf1dc_score(ctx, labels, NULL) - crf1dc_lognorm(ctx);
	fprintf(fp, "crf1dc_score(ctx, labels, NULL) = %f\n", crf1dc_score(ctx, labels, NULL));

	fprintf(fp, "Check for the sequence %d-%d-%d... ", y1, y2, y3);
	check_values(fp, exp(logp), scores[y1][y2][y3] / norm);
      }
    }
  }

  /* Compute the marginal probability at t=0 */
  for (y1 = 0;y1 < L;++y1) {
    floatval_t a, b, c, s = 0.;
    for (y2 = 0;y2 < L;++y2) {
      for (y3 = 0;y3 < L;++y3) {
	s += scores[y1][y2][y3];
      }
    }

    a = ALPHA_SCORE(ctx, 0)[y1];
    b = BETA_SCORE(ctx, 0)[y1];
    c = 1. / ctx->scale_factor[0];

    fprintf(fp, "Check for the marginal probability (0,%d)... ", y1);
    check_values(fp, a * b * c, s / norm);
  }

  /* Compute the marginal probability at t=1 */
  for (y2 = 0;y2 < L;++y2) {
    floatval_t a, b, c, s = 0.;
    for (y1 = 0;y1 < L;++y1) {
      for (y3 = 0;y3 < L;++y3) {
	s += scores[y1][y2][y3];
      }
    }

    a = ALPHA_SCORE(ctx, 1)[y2];
    b = BETA_SCORE(ctx, 1)[y2];
    c = 1. / ctx->scale_factor[1];

    fprintf(fp, "Check for the marginal probability (1,%d)... ", y2);
    check_values(fp, a * b * c, s / norm);
  }

  /* Compute the marginal probability at t=2 */
  for (y3 = 0;y3 < L;++y3) {
    floatval_t a, b, c, s = 0.;
    for (y1 = 0;y1 < L;++y1) {
      for (y2 = 0;y2 < L;++y2) {
	s += scores[y1][y2][y3];
      }
    }

    a = ALPHA_SCORE(ctx, 2)[y3];
    b = BETA_SCORE(ctx, 2)[y3];
    c = 1. / ctx->scale_factor[2];

    fprintf(fp, "Check for the marginal probability (2,%d)... ", y3);
    check_values(fp, a * b * c, s / norm);
  }

  /* Compute the marginal probabilities of transitions. */
  for (y1 = 0;y1 < L;++y1) {
    for (y2 = 0;y2 < L;++y2) {
      floatval_t a, b, s, t, p = 0.;
      for (y3 = 0;y3 < L;++y3) {
	p += scores[y1][y2][y3];
      }

      a = ALPHA_SCORE(ctx, 0)[y1];
      b = BETA_SCORE(ctx, 1)[y2];
      s = EXP_STATE_SCORE(ctx, 1)[y2];
      t = EXP_TRANS_SCORE(ctx, y1)[y2];

      fprintf(fp, "Check for the marginal probability (0,%d)-(1,%d)... ", y1, y2);
      check_values(fp, a * t * s * b, p / norm);
    }
  }

  for (y2 = 0;y2 < L;++y2) {
    for (y3 = 0;y3 < L;++y3) {
      floatval_t a, b, s, t, p = 0.;
      for (y1 = 0;y1 < L;++y1) {
	p += scores[y1][y2][y3];
      }

      a = ALPHA_SCORE(ctx, 1)[y2];
      b = BETA_SCORE(ctx, 2)[y3];
      s = EXP_STATE_SCORE(ctx, 2)[y3];
      t = EXP_TRANS_SCORE(ctx, y2)[y3];

      fprintf(fp, "Check for the marginal probability (1,%d)-(2,%d)... ", y2, y3);
      check_values(fp, a * t * s * b, p / norm);
    }
  }
}

void crf1dc_debug_tree_context(FILE *fp)
{
  int y1, y2, y3;
  floatval_t norm = 0;
  const int L = 3;
  const int T = 3;
  crf1d_context_t *ctx = crf1dc_new(CTXF_MARGINALS, FTYPE_CRF1TREE, L, T, NULL);
  floatval_t *trans = NULL, *state = NULL;
  floatval_t scores[3][3][3];
  int labels[3];

  /* Initialize tree. */
  crfsuite_node_t tree[3];
  int children[2] = {1, 2};
  crfsuite_node_t *node = &tree[0];
  node->self_item_id = 0; node->prnt_item_id = -1;
  node->prnt_node_id = -1; node->cap_children = node->num_children = 2;
  node->children = children;

  node = &tree[1];
  node->self_item_id = 1; node->prnt_item_id = 0;
  node->prnt_node_id = 0; node->cap_children = node->num_children = 0;
  node->children = NULL;

  node = &tree[2];
  node->self_item_id = 2; node->prnt_item_id = 0;
  node->prnt_node_id = 0; node->cap_children = node->num_children = 0;
  node->children = NULL;

  /* Initialize state scores. */
  state = EXP_STATE_SCORE(ctx, 0);
  state[0] = .5;    state[1] = .4;    state[2] = .1;
  state = EXP_STATE_SCORE(ctx, 1);
  state[0] = .4;    state[1] = .5;    state[2] = .1;
  state = EXP_STATE_SCORE(ctx, 2);
  state[0] = .1;    state[1] = .4;    state[2] = .5;

  /* Initialize transition scores. */
  trans = EXP_TRANS_SCORE(ctx, 0);
  trans[0] = .4;    trans[1] = .1;    trans[2] = .4;
  trans = EXP_TRANS_SCORE(ctx, 1);
  trans[0] = .6;    trans[1] = .3;    trans[2] = .1;
  trans = EXP_TRANS_SCORE(ctx, 2);
  trans[0] = .4;    trans[1] = .4;    trans[2] = .1;

  fprintf(fp, "Computing alpha/beta scores...\n");
  ctx->num_items = ctx->cap_items;
  crf1dc_tree_alpha_score(ctx, tree);
  fprintf(fp, "Alpha score computed...\n");
  crf1dc_tree_beta_score(ctx, tree);
  fprintf(fp, "Beta score computed...\n");

  fprintf(fp, "Computing scores for label sequences...\n");
  /* Compute the score of every label sequence. */
  for (y2 = 0;y2 < L;++y2) {
    floatval_t s2 = EXP_STATE_SCORE(ctx, 1)[y2];

    for (y3 = 0; y3 < L; ++y3) {
      floatval_t s3 = s2 * EXP_STATE_SCORE(ctx, 2)[y3];

      for (y1 = 0; y1 < L; ++y1) {
	floatval_t s1 = s3 * EXP_STATE_SCORE(ctx, 0)[y1];
	s1 *= EXP_TRANS_SCORE(ctx, y2)[y1];
	s1 *= EXP_TRANS_SCORE(ctx, y3)[y1];
	scores[y1][y2][y3] = s1;
      }
    }
  }

  fprintf(fp, "Computing partition factor...\n");
  /* Compute the partition factor. */
  norm = 0.;
  for (y1 = 0;y1 < L;++y1) {
    for (y2 = 0;y2 < L;++y2) {
      for (y3 = 0;y3 < L;++y3) {
	norm += scores[y1][y2][y3];
      }
    }
  }

  /* Check the partition factor. */
  fprintf(fp, "Checking partition factor... ");
  check_values(fp, exp(ctx->log_norm), norm);

  /* Compute the sequence probabilities. */
  for (y1 = 0;y1 < L;++y1) {
    for (y2 = 0;y2 < L;++y2) {
      for (y3 = 0;y3 < L;++y3) {
	floatval_t logp;

	labels[0] = y1;
	labels[1] = y2;
	labels[2] = y3;
	logp = crf1dc_tree_score(ctx, labels, tree) - crf1dc_lognorm(ctx);

	fprintf(fp, "Checking sequence %d-%d-%d... ", y1, y2, y3);
	check_values(fp, exp(logp), scores[y1][y2][y3] / norm);
      }
    }
  }

  /* Check marginal probabilities */
  fprintf(fp, "Starting tree marginals...\n");
  crf1dc_tree_marginals(ctx, tree);
  fprintf(fp, "Tree marginals computed...\n");

  /* Compute marginal probability at t=0 */
  floatval_t true_marginal = 0.;
  for (y1 = 0; y1 < L; ++y1) {
    floatval_t a, b, c, s = 0.;
    for (y2 = 0; y2 < L; ++y2) {
      for (y3 = 0; y3 < L; ++y3) {
	s += scores[y1][y2][y3];
      }
    }

    a = ALPHA_SCORE(ctx, 0)[y1];
    b = BETA_SCORE(ctx, 0)[y1];
    c = 1. / ctx->scale_factor[0];

    true_marginal = s / norm;
    fprintf(fp, "Checking marginal probability (0,%d)...\n", y1);
    check_values(fp, a * b * c, true_marginal);
    fprintf(fp, "Checking model marginal (0,%d)...\n", y1);
  }

  /* Compute the marginal probability at t=1 */
  for (y2 = 0;y2 < L;++y2) {
    floatval_t a, b, c, s = 0.;
    for (y1 = 0;y1 < L;++y1) {
      for (y3 = 0;y3 < L;++y3) {
	s += scores[y1][y2][y3];
      }
    }

    a = ALPHA_SCORE(ctx, 1)[y2];
    b = BETA_SCORE(ctx, 1)[y2];
    c = 1. / ctx->scale_factor[1];

    fprintf(fp, "Check for the marginal probability (1,%d)... ", y2);
    check_values(fp, a * b * c, s / norm);
  }

  /* Compute the marginal probability at t=2 */
  for (y3 = 0;y3 < L;++y3) {
    floatval_t a, b, c, s = 0.;
    for (y1 = 0;y1 < L;++y1) {
      for (y2 = 0;y2 < L;++y2) {
	s += scores[y1][y2][y3];
      }
    }

    a = ALPHA_SCORE(ctx, 2)[y3];
    b = BETA_SCORE(ctx, 2)[y3];
    c = 1. / ctx->scale_factor[2];

    fprintf(fp, "Check for the marginal probability (2,%d)... ", y3);
    check_values(fp, a * b * c, s / norm);
  }

  /* Compute marginal probabilities of transitions y2 -- y1. */
  for (y2 = 0;y2 < L;++y2) {
    for (y1 = 0;y1 < L;++y1) {
      floatval_t a, *a_c, b, s, t, t_c, p = 0.;
      for (y3 = 0;y3 < L;++y3) {
  	p += scores[y1][y2][y3];
      }

      a = ALPHA_SCORE(ctx, 1)[y2];
      s = EXP_STATE_SCORE(ctx, 0)[y1];
      t = EXP_TRANS_SCORE(ctx, y2)[y1];

      b = 0;
      a_c = ALPHA_SCORE(ctx, 2);
      for (int i = 0; i < L; ++i) {
	t_c = EXP_TRANS_SCORE(ctx, i)[y1];
	b += a_c[i] * t_c;
      }
      b *= BETA_SCORE(ctx, 0)[y1];

      fprintf(fp, "Check for marginal probability (1,%d)-(0,%d)... ", y2, y1);
      check_values(fp, a * t * s * b, p / norm);
    }
  }

  for (y3 = 0; y3 < L; ++y3) {
    for (y1 = 0; y1 < L; ++y1) {
      floatval_t a, *a_c, b, s, t, t_c, p = 0.;
      for (y2 = 0; y2 < L; ++y2) {
  	p += scores[y1][y2][y3];
      }

      a = ALPHA_SCORE(ctx, 2)[y3];
      b = BETA_SCORE(ctx, 0)[y1];
      s = EXP_STATE_SCORE(ctx, 0)[y1];
      t = EXP_TRANS_SCORE(ctx, y3)[y1];

      b = 0;
      a_c = ALPHA_SCORE(ctx, 1);
      for (int i = 0; i < L; ++i) {
	t_c = EXP_TRANS_SCORE(ctx, i)[y1];
	b += a_c[i] * t_c;
      }
      b *= BETA_SCORE(ctx, 0)[y1];

      fprintf(fp, "Check for marginal probability (2,%d)-(0,%d)... ", y2, y3);
      check_values(fp, a * t * s * b, p / norm);
    }
  }
}

void crf1dc_debug_sm_context(FILE *fp)
{}
