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

#ifndef    __SEMIMARKOV_H__
# define   __SEMIMARKOV_H__

/* Libraries */
# include "ring.h"
# include "rumavl.h"

/**
 * \addtogroup crfsuite_object Object interfaces and utilities.
 * @{
 */

struct crf1de_semimarkov;
/**
 * Synonym for semi-markov struct.
 */
typedef struct crf1de_semimarkov crf1de_semimarkov_t;

/**
 * Structure for holding data specific to semi-markov model.
 */
/* Interface */
struct crf1de_semimarkov {
  int L;		  /**< Number of distinct labels.  */
  int num_fs;		     /**< Number of forward state prefixes. */
  RUMAVL *forward_states;    /**< Dictionary of possible forward state prefixes. */
  int **forward_trans1; /**< Array holding possible forward transitions. */
  int **forward_trans2; /**< Array holding possible forward transitions. */
  int *fs_llabels;  /**< Array of last labels of forward states. */


  int num_bs;		/**< Number of backward states (prefixes * labels). */
  RUMAVL *backward_states;	/**< Array of possible backward states. */
  int *backward_trans; /**< Array holding possible backward transitions. */

  int *pattern_trans1; /**< Array holding possible patterns. */
  int *pattern_trans2;  /**< Array holding possible patterns. */

  int *max_seg_len; /**< Array holding maximum lengths of spans with same label. */
  int *wrkbench;    /**< Auxiliary array for constructing prefixes. */

  /** Allocate memory for necessary data. */
  int (*initialize)(crf1de_semimarkov_t *sm, int max_order, int L);
  /** Update relevent information for semi-markov model. */
  void (*update)(crf1de_semimarkov_t *sm, int a_prev_lbl, int a_seg_len, \
		 crfsuite_ring_t *a_labelseq);
  /** Update relevent information for semi-markov model. */
  void (*finalize)(crf1de_semimarkov_t *sm);
  /** Clear data stored in semi-markov model. */
  void (*clear)(crf1de_semimarkov_t *sm);
};

/**
 * Allocate and initialize semi-markov data.
 *
 * @return pointer to initialized semi-markov data.
 */
crf1de_semimarkov_t *crf1de_create_semimarkov(void);
#endif	/* __SEMIMARKOV_H__ */
