/*
 * Data structures and auxiliary functions for semi-markov CRF.
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
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* $Id$ */

#ifndef    CRFSUITE_SEMIMARKOV_H_
# define   CRFSUITE_SEMIMARKOV_H_

/* Libraries */
# include "ring.h"
# include "rumavl.h"

/* Macros */

/// macro for accessing suffix list
#define SUFFIXES(sm, y, x)				\
  MATRIX(sm->m_suffixes, (sm->m_max_order + 1), x, y)

/**
 * \addtogroup crfsuite_object Object interfaces and utilities.
 * @{
 */

/**
 * Synonym for semi-markov state struct.
 */
typedef struct crf1de_state crf1de_state_t;

/**
 * Auxiliary structure with information about single forward or backward
 * state.
 */
struct crf1de_state {
  int m_id;		/**< id of the label sequence */
  int m_feat_id;	/**< id of the label corresponding feature */
  size_t m_len;		/**< length of the label sequence */
  size_t m_num_affixes;	/**< number of prefixes for given label sequence */

  int *m_frw_trans1; /**< array of prefix indices of forward transitions (pk states) */
  int *m_frw_trans2; /**< array of indices of forward transitions (pky states) */
  int *m_bkw_trans;  /**< array of backward transition indices */

  floatval_t m_freq;		/**< frequency of label pattern */
  int *m_seq;			/**< label sequence */
};

/**
 * Synonym for semi-markov struct.
 */
typedef struct crf1de_semimarkov crf1de_semimarkov_t;

/**
 * Auxiliary structure for holding data specific to semi-markov model.
 */
/* Interface */
struct crf1de_semimarkov {
  /* General data */
  int L;			/**< Number of distinct labels.  */
  size_t m_max_order;		/**< Maximum order of label sequences, this
				   value is one more then the order specified
				   by `-p feature.max_order` parameter). */
  int m_seg_len_lim; /**< Limit of the maximum segment length (value < 0 means
			unconstrained (semi-markov), value >= 0 implies
			standard CRF. */
  int *m_max_seg_len;  /**< Array holding maximum observed segment lengths for
			  spans with given labels. */

  /* Label patterns */
  size_t m_num_ptrns;	    /**< Number of possible tag patterns. */
  crf1de_state_t *m_ptrns;  /**< Array of possible tag sequences. */
  RUMAVL *m_ptrns_set;	    /**< Auxiliary set of possible tag sequences (used
			       during construction). */

  int *m_ptrn_llabels;		/**< Array of last labels of tag patterns. */
  int *m_ptrn_trans1;		/**< Array holding frw state id's of possible
				   pattern transitions. */
  int *m_ptrn_trans2;	 /**< Array holding bkw state id's of possible pattern
			       transitions. */
  int *m_ptrnid2bkwid;	    /**< Array representing mapping from pattern id to
			       backward state id. */

  /* Pattern suffixes */
  int *m_suffixes;		/**< Array of pattern suffixes. */
  size_t m_num_suffixes;	/**< Number of possible pattern suffixes. */
  size_t m_cap_suffixes;	/**< Capacity for storing suffixes. */

  /* Forward states */
  size_t m_num_frw;		/**< Number of forward states. */
  crf1de_state_t *m_frw_states;	/**< Array of forward states (`pk` states). */
  RUMAVL *m_frw_states_set; /**< Auxiliary set of possible forward
			       states (used during construction). */

  int *m_frw_llabels;	       /**< Array of last labels of forward states. */
  int *m_frw_trans1; /**< Array holding possible prefixes for given states. */
  int *m_frw_trans2; /**< Array holding full form of the former prefixes. */

  /* Backward states */
  size_t m_num_bkw;	  /**< Number of backward states. */
  crf1de_state_t *m_bkw_states;	/**< Array of backward states (`pky` states). */
  RUMAVL *m_bkw_states_set;	  /**< Set of backward states. */

  int *m_bkw_trans;   /**< Array holding possible backward transitions. */
  int *m_bkwid2frwid; /**< Mapping from backward state id to forward state id */

  /* Auxiliary data members */
  crf1de_state_t m_wrkbench1;  /**< Auxiliary array for constructing
				  states and transitions. */
  crf1de_state_t m_wrkbench2;  /**< Auxiliary array for constructing
				  states and transitions. */
  crfsuite_ring_t *m_ring; /**< Circular buffer for storing tag sequences. */

  /* Methods */
  /** Allocate memory for necessary data. */
  int (*initialize)(crf1de_semimarkov_t *sm, const int a_max_order, \
		    const int a_seg_len_lim, const int L);
  /** Update relevent information for semi-markov model. */
  void (*update)(crf1de_semimarkov_t *sm, int a_lbl, int a_seg_len);
  /** Generate all possible label sequences for the given order. */
  void (*generate_all_edges)(crf1de_semimarkov_t *sm);
  /** Update relevent information for semi-markov model. */
  int (*finalize)(crf1de_semimarkov_t *sm);
  /** Clear data stored in semi-markov model. */
  void (*clear)(crf1de_semimarkov_t *sm);
  /** Create state from circular buffer of labels. */
  void (*build_state)(crf1de_state_t *a_state, const crfsuite_ring_t *a_ring);
  /** Obtain id of state. */
  int (*get_state_id)(crf1de_state_t *a_state, RUMAVL *a_dic);
  /** Output state. */
  void (*output_state)(FILE *a_fstream, const char *a_name, const crf1de_state_t *a_entry);
};

/**
 * Allocate and initialize semi-markov data.
 *
 * @return pointer to initialized semi-markov data.
 */
crf1de_semimarkov_t *crf1de_create_semimarkov(void);

#endif	/* CRFSUITE_SEMIMARKOV_H_ */
