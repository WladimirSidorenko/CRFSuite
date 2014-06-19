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



crf1d_context_t* crf1dc_new(int flag, int L, int T)
{
    int ret = 0;
    crf1d_context_t* ctx = NULL;

    ctx = (crf1d_context_t*)calloc(1, sizeof(crf1d_context_t));
    if (ctx != NULL) {
        ctx->flag = flag;
        ctx->num_labels = L;

        ctx->trans = (floatval_t*)calloc(L * L, sizeof(floatval_t));
        if (ctx->trans == NULL) goto error_exit;

        if (ctx->flag & CTXF_MARGINALS) {
            ctx->exp_trans = (floatval_t*)_aligned_malloc((L * L + 4) * sizeof(floatval_t), 16);
            if (ctx->exp_trans == NULL) goto error_exit;
            ctx->mexp_trans = (floatval_t*)calloc(L * L, sizeof(floatval_t));
            if (ctx->mexp_trans == NULL) goto error_exit;
        }

        if (ret = crf1dc_set_num_items(ctx, T)) {
            goto error_exit;
        }

        /* T gives the 'hint' for maximum length of items. */
        ctx->num_items = 0;
    }

    return ctx;

error_exit:
    crf1dc_delete(ctx);
    return NULL;
}

int crf1dc_set_num_items(crf1d_context_t* ctx, int T)
{
    const int L = ctx->num_labels;

    ctx->num_items = T;

    if (ctx->cap_items < T) {
        free(ctx->backward_edge);
        free(ctx->mexp_state);
        _aligned_free(ctx->exp_state);
        free(ctx->scale_factor);
        free(ctx->row);
        free(ctx->beta_score);
        free(ctx->alpha_score);

        ctx->alpha_score = (floatval_t*)calloc(T * L, sizeof(floatval_t));
        if (ctx->alpha_score == NULL) return CRFSUITEERR_OUTOFMEMORY;
        ctx->beta_score = (floatval_t*)calloc(T * L, sizeof(floatval_t));
        if (ctx->beta_score == NULL) return CRFSUITEERR_OUTOFMEMORY;
        ctx->scale_factor = (floatval_t*)calloc(T, sizeof(floatval_t));
        if (ctx->scale_factor == NULL) return CRFSUITEERR_OUTOFMEMORY;
        ctx->row = (floatval_t*)calloc(L, sizeof(floatval_t));
        if (ctx->row == NULL) return CRFSUITEERR_OUTOFMEMORY;

        if (ctx->flag & CTXF_VITERBI) {
            ctx->backward_edge = (int*)calloc(T * L, sizeof(int));
            if (ctx->backward_edge == NULL) return CRFSUITEERR_OUTOFMEMORY;
        }

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
        free(ctx->beta_score);
        free(ctx->alpha_score);
        free(ctx->mexp_trans);
        _aligned_free(ctx->exp_trans);
        free(ctx->trans);
    }
    free(ctx);
}

void crf1dc_reset(crf1d_context_t* ctx, int flag)
{
    const int T = ctx->num_items;
    const int L = ctx->num_labels;

    if (flag & RF_STATE)
        veczero(ctx->state, T*L);

    if (flag & RF_TRANS)
        veczero(ctx->trans, L*L);

    if (ctx->flag & CTXF_MARGINALS) {
        veczero(ctx->mexp_state, T*L);
        veczero(ctx->mexp_trans, L*L);
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

void crf1dc_exp_transition(crf1d_context_t* ctx)
{
    const int L = ctx->num_labels;

    veccopy(ctx->exp_trans, ctx->trans, L * L);
    vecexp(ctx->exp_trans, L * L);
}

void crf1dc_alpha_score(crf1d_context_t* a_ctx,  const crfsuite_node_t *a_tree)
{
    int i, t;
    floatval_t sum, *cur = NULL, *scale = &a_ctx->scale_factor[0];
    const floatval_t *prev = NULL, *trans = NULL, *state = NULL;
    const int T = a_ctx->num_items;
    const int L = a_ctx->num_labels;
    /* Compute the alpha scores on leaves (0, *).
        alpha[0][j] = state[0][j]
     */
    cur = ALPHA_SCORE(a_ctx, 0);
    state = EXP_STATE_SCORE(a_ctx, 0);
    // copy L elements from state to current
    veccopy(cur, state, L);
    for (i = 0; i < L; ++i) {
      fprintf(stderr, "alpha[%d][%d] = %f (unscaled)\n", 0, i, cur[i]);
    }
    // total sum of L elements in vector
    sum = vecsum(cur, L);
    *scale = (sum != 0.) ? 1. / sum : 1.;
    // multiply L elements in cur by scale factor (i.e. normalize weights)
    vecscale(cur, *scale, L);
    for (i = 0; i < L; ++i) {
      fprintf(stderr, "alpha[%d][%d] = %f (scaled)\n", 0, i, cur[i]);
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
    	// normalize the weights
	for (i = 0; i < L; ++i) {
	  fprintf(stderr, "alpha[%d][%d] = %f (unscaled)\n", t, i, cur[i]);
	}
        vecscale(cur, *scale, L);
	for (i = 0; i < L; ++i) {
	  fprintf(stderr, "alpha[%d][%d] = %f (scaled)\n", t, i, cur[i]);
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
 * @param a_tree - tree corresponding to given instance
 **/
void crf1dc_tree_alpha_score(crf1d_context_t* a_ctx, const crfsuite_node_t *a_tree)
{
  assert(a_tree);

  int i, t, c;
  floatval_t sum, *cur_alpha, *chld_alpha;
  floatval_t *scale, *row = a_ctx->row;
  const floatval_t *trans, *state;
  const int T = a_ctx->num_items;
  const int L = a_ctx->num_labels;
  const crfsuite_node_t *node, *child; // current and child nodes
  int item_id, ch_item_id;	       // indices of current and child items

  /* Since nodes in the tree are ordered topologically, we start at the end
     of the tree array since that is where the leaf nodes are located.  Once
     we are done with the leaves, we go up the tree and compute the
     probabilities of nodes whose child probabilities have already been
     computed.

     We define the probabilities of the leaves as the exponent of state
     probablities, and probabilities of parent nodes as:

     alpha[t][j] = state[t][j] * \sum_{i \in L} trans[i][j] * \sum_{c \in
     Ch_j} alpha[c][i]
  */
  for (t = T - 1; t >= 0; --t) {
    // node corresponding to given item
    node = &a_tree[t];
    // id of item corresponding to given node
    item_id = node->self_item_id;
    // column for storing transition weights
    cur_alpha = ALPHA_SCORE(a_ctx, item_id);
    // state weigts for that item
    state = EXP_STATE_SCORE(a_ctx, item_id);
    // slot for scaling factor
    scale = &a_ctx->scale_factor[item_id];

    // if current node is a leaf, set it transition weights to state weights
    if (node->num_children == 0) {
      veccopy(cur_alpha, state, L);
    } else {
      // set weights for cirrent node to zeros
      veczero(cur_alpha, L);
      // clear the workspace
      vecset(row, 1., L);
      // compute the total sum of all alpha score weights for children
      for (c = 0; c < node->num_children; ++c) {
	// obtain the child item
	child = &a_tree[node->children[c]];
	ch_item_id = child->self_item_id;
	chld_alpha = ALPHA_SCORE(a_ctx, ch_item_id);
	vecmul(row, chld_alpha, L); /* row will serve as an auxiliary
				       container for total sum of child
				       scores */
      }
      // for each label `i', we compute the sum of the weights for i-th labels
      // in all of the children, then we set the weight for each target tag j
      // in the current column to the sum of the i-th weights * the transition
      // weight i --> j
      for (i = 0; i < L; ++i) {
	trans = EXP_TRANS_SCORE(a_ctx, i);
	// for each cell j of cur add score prev[i] multiplied by
	// transition score i -> j
	vecaadd(cur_alpha, row[i], trans, L);
      }
      // multiply all obtained transition weights by state weights for that
      // labels
      vecmul(cur_alpha, state, L);
    }
    // normalize weights
    sum = vecsum(cur_alpha, L);
    *scale = (sum != 0.) ? 1. / sum : 1.;
    // multiply L elements in cur by scale factor (i.e. normalize weights)
    vecscale(cur_alpha, *scale, L);
    for (i = 0; i < L; ++i) {
      fprintf(stderr, "alpha[%d][%d] = %f (scaled)\n", t, i, cur_alpha[i]);
    }
  }
  // sum logarithms of all elements in scale factor
  /* a_ctx->log_norm = -vecsumlog(a_ctx->scale_factor, T); */
}

void crf1dc_beta_score(crf1d_context_t* a_ctx, const crfsuite_node_t *a_tree)
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
    for (i = 0; i < L; ++i) {
      fprintf(stderr, "beta[%d][%d] = %f (unscaled)\n", T - 1, i, cur[i]);
    }
    for (i = 0; i < L; ++i) {
      fprintf(stderr, "beta[%d][%d] = %f (scaled)\n", T - 1, i, cur[i]);
    }

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
	for (i = 0; i < L; ++i) {
	  fprintf(stderr, "beta[%d][%d] = %f (unscaled)\n", t, i, cur[i]);
	}
        vecscale(cur, *scale, L);
	for (i = 0; i < L; ++i) {
	  fprintf(stderr, "beta[%d][%d] = %f (scaled)\n", t, i, cur[i]);
	}
        --scale;
    }
}

void crf1dc_tree_beta_score(crf1d_context_t* a_ctx, const crfsuite_node_t *a_tree)
{
  /* In tree-structured CRF, beta score will be computed from parents to children */
  int i, t;
  floatval_t sum;
  floatval_t *scale = NULL;
  floatval_t *crnt_beta = NULL;
  floatval_t *row = a_ctx->row;
  const floatval_t *prnt_beta = NULL, *crnt_state = NULL, *prnt_state = NULL, *trans = NULL;
  const int T = a_ctx->num_items;
  const int L = a_ctx->num_labels;
  const crfsuite_node_t *node;	// current node
  int item_id, prnt_item_id;	// indices of current and child items

  /* Compute beta score for root noode. */
  node = &a_tree[0];
  item_id = node->self_item_id;
  crnt_beta = BETA_SCORE(a_ctx, item_id);
  // set beta score for all root labels to *scale
  scale = &a_ctx->scale_factor[item_id]; // ? a_ctx->scale_factor[item_id]: 1.;
  if (*scale != 0)
    vecset(crnt_beta, 1. / *scale, L);
  else
    vecset(crnt_beta, 1., L);

  /* sum = vecsum(crnt_beta, L); */
  /* *scale = (sum != 0.) ? 1. / sum : 1.; */
  /* vecscale(crnt_beta, *scale, L); */
  for (i = 0; i < L; ++i) {
    fprintf(stderr, "beta[0][%d] = %f (scaled)\n", i, crnt_beta[i]);
  }

  /* Compute beta scores for children. */
  for (t = 1; t < T; ++t) {
    node = &a_tree[t];
    item_id = node->self_item_id;
    crnt_beta = BETA_SCORE(a_ctx, item_id);
    crnt_state = EXP_STATE_SCORE(a_ctx, item_id);
    scale = &a_ctx->scale_factor[item_id]; // ? a_ctx->scale_factor[item_id]: 1.;

    prnt_item_id = node->prnt_item_id;
    prnt_beta = BETA_SCORE(a_ctx, prnt_item_id);
    prnt_state = EXP_STATE_SCORE(a_ctx, prnt_item_id);

    veccopy(row, prnt_beta, L);
    vecmul(row, prnt_state, L);

    /* Compute the beta score at (crnt, i). */
    for (i = 0; i < L; ++i) {
      trans = EXP_TRANS_SCORE(a_ctx, i);
      crnt_beta[i] = vecdot(trans, row, L);
    }
    /* scale the vector */
    sum = vecsum(crnt_beta, L);
    *scale = (sum != 0.) ? 1. / sum : 1.;
    vecscale(crnt_beta, *scale, L);
    for (i = 0; i < L; ++i) {
      fprintf(stderr, "beta[%d][%d] = %f (scaled)\n", t, i, crnt_beta[i]);
    }
  }
  // sum logarithms of all elements in scale factor
  a_ctx->log_norm = -vecsumlog(a_ctx->scale_factor, T);
}

void crf1dc_marginals(crf1d_context_t* a_ctx, const crfsuite_node_t *a_tree)
{
  int i, j, t;
  const int T = a_ctx->num_items;
  const int L = a_ctx->num_labels;

  /*
    Compute the model expectations of states.
    p(t,i) = fwd[t][i] * bwd[t][i] / norm
    = (1. / C[t]) * fwd'[t][i] * bwd'[t][i]
  */
  for (t = 0;t < T;++t) {
    floatval_t *fwd = ALPHA_SCORE(a_ctx, t);
    floatval_t *bwd = BETA_SCORE(a_ctx, t);
    floatval_t *prob = STATE_MEXP(a_ctx, t);
    veccopy(prob, fwd, L);
    vecmul(prob, bwd, L);
    vecscale(prob, 1. / a_ctx->scale_factor[t], L);
    fprintf(stderr, "vecsum(prob, L) = %f\n", vecsum(prob, L));
    assert(vecsum(prob, L) - 1. < 0.00001);
    for (i = 0; i < L; ++i) {
      fprintf(stderr, "state_expectation[%d][%d] = %f\n", t, i, prob[i]);
    }
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
  for (t = 0;t < T-1;++t) {
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
	/* fprintf(stderr, "prob[%d][%d] = %f\n", i, j, prob[j]); */
      }
    }
  }
}

void crf1dc_tree_marginals(crf1d_context_t* a_ctx, const crfsuite_node_t *a_tree)
{
  assert(a_tree);

  int i, j, t, item_id, prnt_id;
  const int T = a_ctx->num_items;
  const int L = a_ctx->num_labels;
  const crfsuite_node_t *node, *prnt;

  /*
    Compute the model expectations of states (this expectation is the same for
    all types of graphical models).

    p(t,i) = fwd[t][i] * bwd[t][i] / norm
    = (1. / C[t]) * fwd'[t][i] * bwd'[t][i]
  */
  for (t = 0; t < T; ++t) {
    floatval_t *fwd = ALPHA_SCORE(a_ctx, t);
    floatval_t *bwd = BETA_SCORE(a_ctx, t);
    floatval_t *prob = STATE_MEXP(a_ctx, t);
    veccopy(prob, fwd, L);
    vecmul(prob, bwd, L);
    vecscale(prob, 1. / a_ctx->scale_factor[t], L);
    fprintf(stderr, "vecsum(prob, L) = %f\n", vecsum(prob, L));
    assert(vecsum(prob, L) - 1. < 0.00001);
    for (i = 0; i < L; ++i) {
      fprintf(stderr, "state_expectation[%d][%d] = %f\n", t, i, prob[i]);
    }
  }

  /*
    Compute the model expectations of transitions (these transitions will be
    different for linear and tree structured CRFs).

    p(t,i,t+1,j)
    = fwd[t][i] * edge[i][j] * state[t+1][j] * bwd[t+1][j] / norm
    = (fwd'[t][i] / (C[0] ... C[t])) * edge[i][j] * state[t+1][j] * (bwd'[t+1][j] / (C[t+1] ... C[T-1])) * (C[0] * ... * C[T-1])
    = fwd'[t][i] * edge[i][j] * state[t+1][j] * bwd'[t+1][j]
    The model expectation of a transition (i -> j) is the sum of the marginal
    probabilities p(t,i,t+1,j) over t.
  */
  floatval_t *row = a_ctx->row;
  for (t = T - 1; t > 0; --t) {
    node = &a_tree[t];
    item_id = node->self_item_id;
    prnt_id = node->prnt_item_id;

    floatval_t *fwd = ALPHA_SCORE(a_ctx, item_id);
    floatval_t *bwd = BETA_SCORE(a_ctx, prnt_id);
    floatval_t *state = EXP_STATE_SCORE(a_ctx, prnt_id);

    /* row[j] = state[t+1][j] * bwd'[t+1][j] */
    veccopy(row, bwd, L);
    vecmul(row, state, L);

    for (i = 0;i < L; ++i) {
      floatval_t *edge = EXP_TRANS_SCORE(a_ctx, i);
      floatval_t *prob = TRANS_MEXP(a_ctx, i);
      for (j = 0; j < L; ++j) {
	prob[j] += fwd[i] * edge[j] * row[j];
	fprintf(stderr, "prob[%d][%d] = %f\n", i, j, prob[j]);
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

#if 0
/* Sigh, this was found to be slower than the forward-backward algorithm. */

#define    ADJACENCY(ctx, i) \
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
			const crfsuite_node_t *a_tree)
{
  int i, j, k, t;
    floatval_t ret = 0;
    const floatval_t *state = NULL, *cur = NULL, *trans = NULL;
    const int T = a_ctx->num_items;
    const int L = a_ctx->num_labels;

    /* Stay at (0, labels[0]). */
    i = a_labels[0];
    state = STATE_SCORE(a_ctx, 0);
    ret = state[i];
    /* for (k = 0; k < L; ++k) { */
    /*   fprintf(stderr, "STATE_SCORE[%d][%d] = %f\n", 0, k, state[k]); */
    /* } */

    /* Loop over the rest of items. */
    for (t = 1; t < T; ++t) {
        j = a_labels[t];
	fprintf(stderr, "j = %d\n", j);
        trans = TRANS_SCORE(a_ctx, i);
	/* fprintf(stderr, "TRANS_SCORE[%d][%d][%d] = %f\n", t, i, j, trans[j]); */
        state = STATE_SCORE(a_ctx, t);
	/* for (k = 0; k < L; ++k) { */
	/*   fprintf(stderr, "STATE_SCORE[%d][%d] = %f\n", t, k, state[k]); */
	/* } */

        /* Transit from (t-1, i) to (t, j). */
        ret += trans[j];
        ret += state[j];
	/* fprintf(stderr, "ret = %f\n", ret); */
	i = j;
    }
    /* fprintf(stderr, "ret = %f\n", ret); */
    return ret;
}

floatval_t crf1dc_tree_score(crf1d_context_t* a_ctx, const int *a_labels, \
			     const crfsuite_node_t *a_tree)
{
  assert(a_tree);

  int i, j, t, c;
  floatval_t score = 0., ret = 0.;
  const floatval_t *state = NULL, *cur = NULL, *trans = NULL;
  const int T = a_ctx->num_items;
  const int L = a_ctx->num_labels;

  const crfsuite_node_t *node;
  int item_id, chld_node_id, chld_item_id;
  int label, chld_label;

  /* Loop over items. */
  for (t = 0; t < T; ++t) {
    fprintf(stderr, "t = %d\n", t);
    /* add probability of the state */
    node = &a_tree[t];
    item_id = node->self_item_id;
    fprintf(stderr, "t = %d; item_id = %d\n", t, item_id);
    label = a_labels[item_id];
    state = STATE_SCORE(a_ctx, item_id);
    score = state[label];
    fprintf(stderr, "label = %d\n", label);
    for (i = 0; i < L; ++i) {
      fprintf(stderr, "STATE_SCORE[%d][%d] = %f\n", t, i, state[i]);
    }

    /* add transition probabilities from each of the children */
    for (c = 0; c < node->num_children; ++c) {
      chld_node_id = node->children[c];
      chld_item_id = a_tree[chld_node_id].self_item_id;
      chld_label = a_labels[chld_item_id];
      /* get transition row corresponding to the tag */
      trans = TRANS_SCORE(a_ctx, chld_label);
      fprintf(stderr, "TRANS_SCORE[%d][%d][%d] = %f\n", t, chld_label, label, \
	      trans[label]);
      /* add transition probability from child tag to parent tag */
      score += trans[label];
    }
    ret += score;
    fprintf(stderr, "score = %f; ret = %f\n", score, ret);
  }
  fprintf(stderr, "Finished: score = %f; ret = %f\n", score, ret);
  return ret;
}

floatval_t crf1dc_lognorm(crf1d_context_t* ctx)
{
    return ctx->log_norm;
}

floatval_t crf1dc_viterbi(crf1d_context_t* ctx, int *labels)
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

    /* Compute the scores at (0, *). */
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
    crf1d_context_t *ctx = crf1dc_new(CTXF_MARGINALS, L, T);
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
    crf1dc_alpha_score(ctx, NULL); /* TODO: re-write to support tree_alpha_score */
    crf1dc_beta_score(ctx, NULL); /* TODO: re-write to support tree_beta_score */

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
