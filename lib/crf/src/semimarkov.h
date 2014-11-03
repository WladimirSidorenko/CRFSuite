/*
 *      Data structures and auxiliary functions for semi-markov CRF.
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

#ifndef    CRFSUITE_SEMIMARKOV_H_
# define   CRFSUITE_SEMIMARKOV_H_

/* Libraries */
# include "ring.h"
# include "rumavl.h"

/* Macros */
# ifndef CRFSUITE_SM_MAX_PTRN_LEN
#  define CRFSUITE_SM_MAX_PTRN_LEN 64
# endif

/**
 * \addtogroup crfsuite_object Object interfaces and utilities.
 * @{
 */

/**
 * Synonym for semi-markov state struct.
 */
typedef struct crf1de_state crf1de_state_t;

/**
 * Auxiliary structure for holding information about single forward or
 * backward state.
 */
struct crf1de_state {
  size_t m_id;	       /**< id of the label sequence */
  size_t m_len;	       /**< length of the label sequence */
  size_t m_n_prefixes; /**< number of prefixes for given label sequence */
  size_t m__cnt_trans1; /**< internal counter of transitions */
  size_t m__cnt_trans2; /**< internal counter of transitions */

  crf1de_state_t **m_frw_trans1; /**< array of prefixes (pk states) */
  crf1de_state_t **m_frw_trans2; /**< array of prefixes (pky states) */
  crf1de_state_t **m_bkw_trans;	 /**< array of backward states */

  floatval_t m_freq;		       /**< frequency of label pattern */
  int m_seq[CRFSUITE_SM_MAX_PTRN_LEN]; /**< label sequence */
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
  int L;	     /**< Number of distinct labels.  */
  int m_seg_len_lim; /**< Limit on the maximum segment length (value < 0 means
			unconstrained (semi-markov), value >= 0 implies standard CRF. */

  size_t m_max_order;   /**< Maximum order of the label sequence. */
  int *m_max_seg_len;  /**< Array holding maximum observed segment lengths for
			 spans with given labels. */

  /* Label patterns */
  size_t m_num_ptrns;		/**< Number of possible tag patterns. */
  RUMAVL *m_patterns;		/**< Set of possible tag sequences. */
  int *m_ptrn_llabels;		/**< Array of patterns' last labels. */
  crf1de_state_t **m_ptrn_trans1; /**< Array holding possible pattern
				   transitions. */
  crf1de_state_t **m_ptrn_trans2; /**< Array holding possible pattern
				     transitions. */
  crf1de_state_t **m_ptrnid2ptrn; /**< Mapping from pattern id to pattern state. */
  crf1de_state_t **m_ptrnid2bkw; /**< Mapping from pattern id to backward state. */

  /* Pattern suffixes */
  size_t m_num_ptrn_suffixes;	/**< Number of possible pattern suffixes. */
  crf1de_state_t **m_suffixes;	/**< Array of pattern suffixes. */

  /* Forward states */
  size_t m_num_frw;	/**< Number of forward states. */
  RUMAVL *m_frw_states;	/**< Set of possible forward states. */
  int *m_frw_llabels;	/**< Array of last labels of forward states. */
  crf1de_state_t **m_frw_trans1; /**< Array holding possible prefixes for given states. */
  crf1de_state_t **m_frw_trans2; /**< Array holding full form of the former prefixes. */
  crf1de_state_t **m_frwid2frw;	/**< Mapping from forward state id to forward state */

  /* Backward states */
  size_t m_num_bkw;	  /**< Number of backward states. */
  RUMAVL *m_bkw_states;	  /**< Set of backward states. */
  crf1de_state_t **m_bkw_trans;	  /**< Array holding possible backward transitions. */
  crf1de_state_t **m_bkwid2bkw;   /**< Mapping from backward state id to backward state */
  crf1de_state_t **m_bkwid2frw;   /**< Mapping from backward state id to forward state */

  /* Auxiliary data members */
  crf1de_state_t m_wrkbench1;  /**< Auxiliary array for constructing
				  states and transitions. */
  crf1de_state_t m_wrkbench2;  /**< Auxiliary array for constructing
				  states and transitions. */
  crfsuite_ring_t *m_ring; /**< Circular buffer for storing tag sequences. */

  /* Functions */
  /** Allocate memory for necessary data. */
  int (*initialize)(crf1de_semimarkov_t *sm, const int a_max_order, \
		    const int a_seg_len_lim, const int L);
  /** Update relevent information for semi-markov model. */
  void (*update)(crf1de_semimarkov_t *sm, int a_lbl, int a_seg_len);
  /** Update relevent information for semi-markov model. */
  int (*finalize)(crf1de_semimarkov_t *sm);
  /** Clear data stored in semi-markov model. */
  void (*clear)(crf1de_semimarkov_t *sm);
};

/**
 * Allocate and initialize semi-markov data.
 *
 * @return pointer to initialized semi-markov data.
 */
crf1de_semimarkov_t *crf1de_create_semimarkov(void);
#endif	/* CRFSUITE_SEMIMARKOV_H_ */
