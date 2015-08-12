/*
 * Declaration of circular buffer interface.
 *
 * Copyright (c) 2015, Uladzimir Sidarenka
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

#ifndef    CRFSUITE_RING_H_
# define   CRFSUITE_RING_H_

# include <stdlib.h>		/* size_t */

/**
 * \addtogroup crfsuite_object Object interfaces and utilities.
 * @{
 */

/**
 * CRFSuite single chain link of a ring.
 */
struct tag_crfsuite_chain_link;
typedef struct tag_crfsuite_chain_link crfsuite_chain_link_t;

struct tag_crfsuite_chain_link {
  /** Stored element.*/
  int data;
  /** Link to the previous element.*/
  crfsuite_chain_link_t *prev;
  /** Link to the next element.*/
  crfsuite_chain_link_t *next;
};

/**
 * CRFSuite ring interface.
 */
struct tag_crfsuite_ring;
typedef struct tag_crfsuite_ring crfsuite_ring_t;

/**
 * Create an instance of \ref crfsuite_ring_t.
 *
 * @param a_ring - address of pointer which should point to the newly created ring
 * @param a_size - maximum number of elements to store in the ring
 *
 * @return \c 0 on success and non-\c 0 otherwise
 */
int crfsuite_ring_create_instance(crfsuite_ring_t **a_ring, const size_t a_size);

struct tag_crfsuite_ring {
  /**
   * Pointer to instance data (internal use only).
   */
  void *internal;

  /**
   * Pointer to the initial element of the ring.
   */
  crfsuite_chain_link_t *head;

  /**
   * Pointer to the last element of the ring.
   */
  crfsuite_chain_link_t *tail;

  /**
   * Clear elements in the ring.
   *
   * @param a_ring - pointer to queue instance to be freed
   */
  void (*free)(crfsuite_ring_t *a_ring);

  /**
   * Add element to the ring.
   *
   * @param a_ring - queue instance to which element should be added
   * @param a_el - element to push in the ring
   */
  void (*push)(crfsuite_ring_t *a_ring, int a_el);

  /**
   * Remove last element from the ring.
   *
   * @param a_ring - queue instance from which element should be removed
   */
  void (*pop)(crfsuite_ring_t *a_ring);

  /**
   * Reset counters of elements in the ring (without deallocation).
   *
   * @param a_ring - pointer to ring instance to be reset
   */
  void (*reset)(crfsuite_ring_t *a_ring);

  /**
   * Maximum number of items to store in the ring.
   */
  size_t max_items;

  /**
   * Actual number of stored items.
   */
  size_t num_items;
};
/**@}*/
#endif	/* CRFSUITE_RING_H_ */
