/*
 *      Implementation of dictionary.
 *
 * Copyright (c) 2007-2010, Uladzimir Sidarenka
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

///////////////
// Libraries //
///////////////

#include "ring.h"

#include <stdlib.h>

/////////////
// Methods //
/////////////

// Push an element in the ring.
static int crfsuite_ring_push(crfsuite_ring_t *a_ring, int a_el)
{
  if (a_ring->max_items == 0)
    return -1;

  /* if we overflow the buffer, move the head of ringed queue */
  if (a_ring->n_items >= a_ring->max_items)
    ++a_ring->head;
  else
    ++a_ring->n_items;

  /* if tail or head go beyond array boundaries, reset them to the beginning
     of the array */
  if (++a_ring->tail >= a_ring->end)
    a_ring->tail = a_ring->internal;

  if (++a_ring->head >= a_ring->end)
    a_ring->head = a_ring->internal;

  *a_ring->tail = a_el;
  return 0;
}

// Clear ring and reset pointers.
static void crfsuite_ring_free(crfsuite_ring_t *a_ring)
{
  free(a_ring->internal);
  a_ring->internal = a_ring->end = NULL;
  a_ring->head = a_ring->tail = NULL;
  a_ring->max_items = a_ring->n_items = 0;
}

int crfsuite_ring_create_instance(crfsuite_ring_t **a_ring, int a_size)
{
  crfsuite_ring_t* iring = (crfsuite_ring_t*) calloc(1, sizeof(crfsuite_ring_t));

  if (iring != NULL) {
    iring->n_items = 0;
    iring->internal = (int *) calloc(a_size, sizeof(int));
    if (iring->internal == NULL)
      return 1;

    iring->max_items = a_size;
    iring->end = iring->internal + iring->max_items + 1;
    iring->head = iring->tail = iring->internal;
    iring->push = crfsuite_ring_push;
    iring->free = crfsuite_ring_free;

    *a_ring = iring;
    return 0;
  } else {
    return -1;
  }
}
