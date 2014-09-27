/*
 *      CRF1d feature generator (dyad features).
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <crfsuite.h>

#include "logging.h"
#include "crf1d.h"
#include "rumavl.h"    /* AVL tree library necessary for feature generation. */

/**
 * Feature set.
 */
typedef struct {
  RUMAVL* avl;    /**< Root node of the AVL tree. */
  int num;        /**< Number of features in the AVL tree. */
} featureset_t;


#define    COMP(a, b)    ((a)>(b))-((a)<(b))

static int featureset_comp(const void *x, const void *y, size_t n, void *udata)
{
  int ret = 0;
  const crf1df_feature_t* f1 = (const crf1df_feature_t*)x;
  const crf1df_feature_t* f2 = (const crf1df_feature_t*)y;

  ret = COMP(f1->type, f2->type);
  if (ret == 0) {
    ret = COMP(f1->src, f2->src);
    if (ret == 0) {
      ret = COMP(f1->dst, f2->dst);
    }
  }
  return ret;
}

static featureset_t* featureset_new()
{
  featureset_t* set = NULL;
  set = (featureset_t*)calloc(1, sizeof(featureset_t));
  if (set != NULL) {
    set->num = 0;
    set->avl = rumavl_new(
			  sizeof(crf1df_feature_t), featureset_comp, NULL, NULL);
    if (set->avl == NULL) {
      free(set);
      set = NULL;
    }
  }
  return set;
}

static void featureset_delete(featureset_t* set)
{
  if (set != NULL) {
    rumavl_destroy(set->avl);
    free(set);
  }
}

static int featureset_add(featureset_t* set, const crf1df_feature_t* f)
{
  /* Check whether if the feature already exists. */
  crf1df_feature_t *p = (crf1df_feature_t*)rumavl_find(set->avl, f);
  if (p == NULL) {
    /* Insert the feature to the feature set. */
    rumavl_insert(set->avl, f);
    ++set->num;
  } else {
    /* An existing feature: add the observation expectation. */
    p->freq += f->freq;
  }
  return 0;
}

static crf1df_feature_t* featureset_generate(int *ptr_num_features,
					     featureset_t* set,
					     floatval_t minfreq)
{
  int n = 0, k = 0;
  RUMAVL_NODE *node = NULL;
  crf1df_feature_t *f = NULL;
  crf1df_feature_t *features = NULL;

  /* The first pass: count the number of valid features. */
  while ((node = rumavl_node_next(set->avl, node, 1, (void**)&f)) != NULL) {
    if (minfreq <= f->freq) {
      ++n;
    }
  }

  /* The second path: copy the valid features to the feature array. */
  features = (crf1df_feature_t*)calloc(n, sizeof(crf1df_feature_t));
  if (features != NULL) {
    node = NULL;
    while ((node = rumavl_node_next(set->avl, node, 1, (void**)&f)) != NULL) {
      if (minfreq <= f->freq) {
	memcpy(&features[k], f, sizeof(crf1df_feature_t));
	++k;
      }
    }
    *ptr_num_features = n;
    return features;
  } else {
    *ptr_num_features = 0;
    return NULL;
  }
}

/* Custom function for comparing label sequences.

   @param a_lseq1 - first label sequence to compare
   @param a_lseq2 - second label sequence to compare
   @param a_size - maximum sequence size
   @param a_udata - unused parameter (needed for compliance)

   @return \c void
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

crf1df_feature_t* crf1df_generate(int *ptr_num_features,		\
				  crf1de_semimarkov_t *sm,		\
				  int *max_items,			\
				  const crf1de_option_t *opt,		\
				  const dataset_t *ds,			\
				  const int ftype,			\
				  const int num_labels,			\
				  const crfsuite_logging_callback func,	\
				  void *instance)
{
  int c, i, j, s, t;
  int prev, cur, seg_len;

  crf1df_feature_t f;
  crf1df_feature_t *features = NULL;
  featureset_t* set = NULL;

  const int N = ds->num_instances;
  const int L = num_labels;
  const int connect_all_attrs = opt->feature_possible_states ? 1 : 0;
  const int connect_all_edges = opt->feature_possible_transitions ? 1 : 0;
  const int minfreq = opt->feature_minfreq;
  const int max_order = opt->feature_max_order;
  const int restr_seg_len = opt->feature_max_seg_len > 0;

  logging_t lg;
  lg.func = func;
  lg.instance = instance;
  lg.percent = 0;

  /* Create an instance of feature set. */
  set = featureset_new();

  /* Initialize workbench and ring for storing previous labels for
     higher-order features */
  int *wb = NULL;
  crfsuite_ring_t *labelseq = NULL;
  if (ftype == FTYPE_SEMIMCRF) {
    if (crfsuite_ring_create_instance(&labelseq, max_order)) {
      featureset_delete(set);
      return NULL;
    }
    /* Initialize workbench for constructing affixes */
    wb = (int *) malloc((max_order + 1) * sizeof(int));
    if (wb == NULL) {
      featureset_delete(set);
      return NULL;
    }
  }

  /* Loop over the sequences in the training data. */
  logging_progress_start(&lg);

  // auxiliary variables for iteration
  const crfsuite_item_t* item = NULL;
  const crfsuite_node_t *node_p = NULL, *chld_node_p = NULL;
  // iterate over instances
  for (s = 0; s < N; ++s) {
    const crfsuite_instance_t* seq = dataset_get(ds, s);
    const int T = seq->num_items;
    if (T > *max_items)
      *max_items = T;

    /* reset counters */
    prev = L; cur = 0; seg_len = 1;
    if (ftype == FTYPE_SEMIMCRF)
      labelseq->reset(labelseq);

    /* Loop over items in the sequence. */
    for (t = 0; t < T; ++t) {
      item = &seq->items[t];
      cur = seq->labels[t];

      // Transition feature: label #prev -> label #(item->yid)

      /* If feature type is tree, add all edges from children to the
	 current node as features (nothing is added for leaves). */
      if (ftype == FTYPE_CRF1TREE) {
	// obtain pointer to node which corresponds to this item
	assert(item->id >= 0 && item->id < seq->num_items);
	node_p = &seq->tree[item->id];
	// iterate over all children of that node
	for (int i = 0; i < node_p->num_children; ++i) {
	  f.type = FT_TRANS;
	  chld_node_p = &seq->tree[node_p->children[i]];
	  f.src = seq->labels[chld_node_p->self_item_id];
	  f.dst = cur;
	  f.freq = 1;
	  featureset_add(set, &f);
	}
	/* In semi-markov model, we generate all possible prefixes and affixes
	   up to and including max order. */
      } else if (prev != L) {
	if (ftype == FTYPE_SEMIMCRF) {
	  if (prev != cur) {
	    labelseq->push(labelseq, cur);
	    sm->update(sm, prev, seg_len, labelseq);
	    seg_len = 1;
	  } else {
	    ++seg_len;
	  }
	}
	f.type = FT_TRANS;
	f.src = prev;
	f.dst = cur;
	f.freq = 1;
	featureset_add(set, &f);
      }

      /* Iterate over state features. */
      for (c = 0; c < item->num_contents; ++c) {
	/* State feature: attribute #a -> state #(item->yid). */
	f.type = FT_STATE;
	f.src = item->contents[c].aid;
	f.dst = cur;
	f.freq = item->contents[c].value;
	featureset_add(set, &f);

	/* Generate state features connecting attributes with all
	   output labels. These features are not unobserved in the
	   training data (zero expexctations). */
	if (connect_all_attrs) {
	  for (i = 0;i < L;++i) {
	    f.type = FT_STATE;
	    f.src = item->contents[c].aid;
	    f.dst = i;
	    f.freq = 0;
	    featureset_add(set, &f);
	  }
	}
      }
      prev = cur;
    }
    if (ftype == FTYPE_SEMIMCRF)
      sm->update(sm, prev, seg_len, labelseq);

    logging_progress(&lg, s * 100 / N);
  }
  if (ftype == FTYPE_SEMIMCRF) {
    sm->finalize(sm);
    labelseq->free(labelseq);
  }
  logging_progress_end(&lg);

  /* Generate edge features representing all pairs of labels.
     These features are not unobserved in the training data
     (zero expexcations). */
  if (connect_all_edges) {
    for (i = 0; i < L; ++i) {
      for (j = 0; j < L; ++j) {
	f.type = FT_TRANS;
	f.src = i;
	f.dst = j;
	f.freq = 0;
	featureset_add(set, &f);
      }
    }
  }

  /* Convert the feature set to an feature array. */
  features = featureset_generate(ptr_num_features, set, minfreq);

  /* Delete the feature set. */
  featureset_delete(set);

  return features;
}

int crf1df_init_references(feature_refs_t **ptr_attributes,
			   feature_refs_t **ptr_trans,
			   const crf1df_feature_t *features,
			   const int K,
			   const int A,
			   const int L)
{
  int i, k;
  feature_refs_t *fl = NULL;
  feature_refs_t *attributes = NULL;
  feature_refs_t *trans = NULL;

  /*
    The purpose of this routine is to collect references (indices) of:
    - state features fired by each attribute (attributes)
    - transition features pointing from each label (trans)
  */

  /* Allocate arrays for feature references. */
  attributes = (feature_refs_t*)calloc(A, sizeof(feature_refs_t));
  if (attributes == NULL) goto error_exit;
  trans = (feature_refs_t*) calloc(L, sizeof(feature_refs_t));
  if (trans == NULL) goto error_exit;

  /*
    First, loop over the features to count the number of references.  We
    don't use realloc() to avoid memory fragmentation.
  */
  for (k = 0; k < K; ++k) {
    const crf1df_feature_t *f = &features[k];
    switch (f->type) {
    case FT_STATE:
      ++attributes[f->src].num_features;
      break;
    case FT_TRANS:
      ++trans[f->src].num_features;
      break;
    }
  }
  /*
    Second, allocate memory blocks to store feature references.  We also
    clear fl->num_features fields, which will be used as indices in the
    next phase.
  */
  for (i = 0;i < A;++i) {
    fl = &attributes[i];
    fl->fids = (int*)calloc(fl->num_features, sizeof(int));
    if (fl->fids == NULL) goto error_exit;
    fl->num_features = 0;
  }

  for (i = 0;i < L;++i) {
    fl = &trans[i];
    fl->fids = (int*)calloc(fl->num_features, sizeof(int));
    if (fl->fids == NULL) goto error_exit;
    fl->num_features = 0;
  }

  /*
    Finally, store feature indices.
  */
  for (k = 0; k < K;++k) {
    const crf1df_feature_t *f = &features[k];
    switch (f->type) {
    case FT_STATE:
      fl = &attributes[f->src];
      fl->fids[fl->num_features++] = k;
      break;
    case FT_TRANS:
      fl = &trans[f->src];
      fl->fids[fl->num_features++] = k;
      break;
    }
  }

  *ptr_attributes = attributes;
  *ptr_trans = trans;
  return 0;

 error_exit:
  if (attributes != NULL) {
    for (i = 0;i < A;++i) free(attributes[i].fids);
    free(attributes);
  }
  if (trans != NULL) {
    for (i = 0;i < L;++i) free(trans[i].fids);
    free(trans);
  }
  *ptr_attributes = NULL;
  *ptr_trans = NULL;
  return -1;
}
