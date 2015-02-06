/*
 *      CRFsuite library.
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

#include <os.h>

#include <assert.h>
#include <stdarg.h>
#include <stddef.h>		/* ptrdiff_t, size_t */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <crfsuite.h>
#include "logging.h"

int crf1de_create_instance(const char *iid, void **ptr);
int crfsuite_dictionary_create_instance(const char *interface, void **ptr);
int crf1m_create_instance_from_file(const char *filename, void **ptr, const int ftype);

static void swap_vars(int *a, int *b, int *tmp)
{
  *tmp = *a;
  *a = *b;
  *b = *tmp;
}


int crfsuite_create_instance(const char *iid, void **ptr)
{
  int ret = \
    crf1de_create_instance(iid, ptr) == 0 ||
    crfsuite_dictionary_create_instance(iid, ptr) == 0;
  return ret;
}

int crfsuite_create_instance_from_file(const char *filename, void **ptr, const int ftype)
{
  int ret = crf1m_create_instance_from_file(filename, ptr, ftype);
  return ret;
}

void crfsuite_attribute_init(crfsuite_attribute_t* cont)
{
  memset(cont, 0, sizeof(*cont));
  cont->value = 1;
}

void crfsuite_attribute_set(crfsuite_attribute_t* cont, int aid, floatval_t value)
{
  crfsuite_attribute_init(cont);
  cont->aid = aid;
  cont->value = value;
}

void crfsuite_attribute_copy(crfsuite_attribute_t* dst, const crfsuite_attribute_t* src)
{
  dst->aid = src->aid;
  dst->value = src->value;
}

void crfsuite_attribute_swap(crfsuite_attribute_t* x, crfsuite_attribute_t* y)
{
  crfsuite_attribute_t tmp = *x;
  x->aid = y->aid;
  x->value = y->value;
  y->aid = tmp.aid;
  y->value = tmp.value;
}


void crfsuite_item_init(crfsuite_item_t* item)
{
  memset(item, 0, sizeof(*item));
}

void crfsuite_item_init_n(crfsuite_item_t* item, int num_contents)
{
  crfsuite_item_init(item);
  item->num_contents = num_contents;
  item->cap_contents = num_contents;
  item->contents = (crfsuite_attribute_t*)calloc(num_contents, sizeof(crfsuite_attribute_t));
}

void crfsuite_item_finish(crfsuite_item_t* item)
{
  free(item->contents);
  free(item->node_label);
  crfsuite_item_init(item);
}

void crfsuite_item_copy(crfsuite_item_t* dst, const crfsuite_item_t* src)
{
  int i;

  dst->id = src->id;
  dst->prnt = src->prnt;
  dst->num_contents = src->num_contents;
  dst->cap_contents = src->cap_contents;
  dst->contents = (crfsuite_attribute_t*)calloc(dst->num_contents, sizeof(crfsuite_attribute_t));

  for (i = 0;i < dst->num_contents;++i)
    crfsuite_attribute_copy(&dst->contents[i], &src->contents[i]);

  free(dst->node_label);
  if (src->node_label) {
    dst->node_label = (char *) malloc(sizeof(char) * (strlen(src->node_label) + 1));

    if (dst->node_label == NULL)
      fprintf(stderr, "ERROR: Could not allocate memory for node label copy.\n");
    else
      strcpy(dst->node_label, src->node_label);

  } else {
    dst->node_label = NULL;
  }
}

void crfsuite_item_swap(crfsuite_item_t* x, crfsuite_item_t* y)
{
  crfsuite_item_t tmp = *x;
  /* memberwise copy */
  x->num_contents = y->num_contents;
  x->cap_contents = y->cap_contents;
  x->contents = y->contents;
  /* memberwise copy */
  y->num_contents = tmp.num_contents;
  y->cap_contents = tmp.cap_contents;
  y->contents = tmp.contents;
}

int crfsuite_item_append_attribute(crfsuite_item_t* item, const crfsuite_attribute_t* cont)
{
  if (item->cap_contents <= item->num_contents) {
    item->cap_contents = (item->cap_contents + 1) * 2;
    item->contents = (crfsuite_attribute_t*) realloc(
						     item->contents, sizeof(crfsuite_attribute_t) * item->cap_contents);
  }
  crfsuite_attribute_copy(&item->contents[item->num_contents++], cont);
  return 0;
}

int  crfsuite_item_empty(crfsuite_item_t* item)
{
  return (item->num_contents == 0);
}

static int crfsuite_node_add_child(crfsuite_node_t *a_node, int a_chld)
{
  if (a_node->num_children >= a_node->cap_children) {
    a_node->cap_children = (a_node->cap_children + 1) * 2;
    a_node->children = (int *) realloc(a_node->children, \
				       a_node->cap_children * sizeof(int));
    if (a_node->children == NULL) {
      fprintf(stderr, "ERROR: Could not allocate space for children.\n");
      return -1;
    }
  }
  a_node->children[a_node->num_children++] = a_chld;
  return 0;
}

static int crfsuite_node_copy(crfsuite_node_t* a_dst, const crfsuite_node_t* const a_src)
{
  // perform memberwise copy
  a_dst->self_item_id = a_src->self_item_id;
  a_dst->prnt_item_id = a_src->prnt_item_id;
  a_dst->prnt_node_id = a_src->prnt_node_id;

  a_dst->num_children = a_src->num_children;
  a_dst->cap_children = a_src->cap_children;
  // delete destination children if necessary
  if (a_dst->children)
    free(a_dst->children);

  // populate destination's children with children from source
  a_dst->children = (int *) calloc(a_src->num_children, sizeof(int));
  if (a_dst->children == NULL) {
    fprintf(stderr, "ERROR: Could not allocate space for children copy.\n");
    return -2;
  }
  memcpy(a_dst->children, a_src->children, a_src->num_children * sizeof(int));
  return 0;
}

static void crfsuite_node_swap(crfsuite_node_t *a_trg, crfsuite_node_t *a_src)
{
  if (a_trg == a_src)
    return;

  int aux, *auxp;
  swap_vars(&a_trg->self_item_id, &a_src->self_item_id, &aux);
  swap_vars(&a_trg->prnt_item_id, &a_src->prnt_item_id, &aux);
  swap_vars(&a_trg->prnt_node_id, &a_src->prnt_node_id, &aux);
  swap_vars(&a_trg->cap_children, &a_src->cap_children, &aux);
  swap_vars(&a_trg->num_children, &a_src->num_children, &aux);

  auxp = a_trg->children;
  a_trg->children = a_src->children;
  a_src->children = auxp;
}

/**
 * Delete tree with all its nodes.
 *  @param  a_tree      Tree's address.
 *  @param  a_n_nodes   Number of nodes in tree.
 */
static void crfsuite_tree_finish(crfsuite_node_t **a_tree, const int n_nodes)
{
  if (*a_tree == NULL)
    return;

  for (int i = 0; i < n_nodes; ++i)
    free((*a_tree)[i].children);

  free(*a_tree);
  *a_tree = NULL;
}

static int crfsuite_tree_copy(crfsuite_instance_t* dst, const crfsuite_instance_t* const src)
{
  int ret, item_cnt = 0;
  if (src->tree == NULL) {
    dst->tree = NULL;
    return 0;
  }

  if (dst->tree == NULL) {
    dst->tree = (crfsuite_node_t *) calloc(src->num_items, sizeof(crfsuite_node_t));
    if (dst->tree == NULL) {
      fprintf(stderr, "ERROR: Could not allocate space for tree copy.\n");
      goto error_exit;
    }
    memset(dst->tree, 0, src->num_items * sizeof(crfsuite_node_t *));
  }

  for (int i = 0; i < src->num_items; ++i) {
    if ((ret = crfsuite_node_copy(&dst->tree[i], &src->tree[i]) != 0)) {
      item_cnt = i;
      goto error_exit;
    }
  }
  return 0;

 error_exit:
  fprintf(stderr, "ERROR: could not copy tree\n");
  crfsuite_tree_finish(&dst->tree, item_cnt);
  return -1;
}

static int crfsuite_tree_populate_children(crfsuite_node_t *a_tree,	\
					   const int a_root_id,		\
					   const int a_n_items)
{
  assert(a_tree);

  int ret = 0;
  crfsuite_node_t *node_p = NULL, *prnt_node_p = NULL;

  // iterate over all nodes and populate their parents and children
  for (int i = 0; i < a_n_items; ++i) {
    node_p = &a_tree[i];

    if (i == a_root_id) {
      node_p->prnt_item_id = -1;
    } else {
      assert(node_p->prnt_node_id >= 0);
      prnt_node_p = &a_tree[node_p->prnt_node_id];

      node_p->prnt_item_id = prnt_node_p->self_item_id;
      // add this node as child to parent
      if ((ret = crfsuite_node_add_child(prnt_node_p, i)))
	return ret;
    }
  }
  return ret;
}

static void crfsuite_tree_get_order(const crfsuite_node_t *a_tree,	\
				    int a_root_id,			\
				    const int a_n_items,		\
				    int **old2new,			\
				    int **new2old)
{
  assert(a_root_id >= 0 && a_root_id < a_n_items);

  int old_id, n2o_id = 0;
  int n_children;
  const crfsuite_node_t *crnt_node_p;

  *old2new = (int *) calloc(a_n_items, sizeof(int));
  *new2old = (int *) calloc(a_n_items, sizeof(int));
  if (*old2new == NULL || *old2new == NULL)
    goto error_exit;

  // populate base case
  (*new2old)[n2o_id++] = a_root_id;
  // traverse tree in breadth first search manner and sequentially assign next
  // available indices to node's children
  for (int i = 0; i < a_n_items; ++i) {
    old_id = (*new2old)[i];
    (*old2new)[old_id] = i;

    crnt_node_p = &a_tree[old_id];
    n_children = crnt_node_p->num_children;

    // iterate over children and place them at the next available indices
    for (int j = 0; j < n_children; ++j) {
      if (n2o_id >= a_n_items) {
	fprintf(stderr, "ERROR: Tree has a loop.\n");
	goto error_exit;
      }
      (*new2old)[n2o_id++] = crnt_node_p->children[j];
    }
  }
  return;

 error_exit:
  free(*old2new); *old2new = NULL;
  free(*new2old); *new2old = NULL;
  return;
}

static int crfsuite_tree_reorder(crfsuite_node_t *a_tree,		\
				 const int a_root_id,			\
				 const int a_n_items)
{
  assert(a_root_id >= 0 && a_root_id < a_n_items);
  int ret = 0, old_i, n_children;
  int *old2new = NULL, *new2old = NULL;
  crfsuite_node_t *crnt_node = NULL;
  // get mappings from old node indices to new node indices and vice versa
  crfsuite_tree_get_order(a_tree, a_root_id, a_n_items, &old2new, &new2old);
  if (old2new == NULL || new2old == NULL) {
    ret = -1;
    goto final_steps;
  }

  // swap nodes according to their topological order and update their children
  int aux_i;
  for (int i = 0; i < a_n_items; ++i) {
    // get index of the node wich should be placed at i-th position in
    // topological order
    old_i = new2old[i];
    // swap i-th node with the node at old_i position
    crnt_node = &a_tree[i];
    crfsuite_node_swap(crnt_node, &a_tree[old_i]);
    crnt_node->prnt_node_id = old2new[crnt_node->prnt_node_id];
    // update child indices
    n_children = crnt_node->num_children;
    for (int j = 0; j < n_children; ++j) {
      crnt_node->children[j] = old2new[crnt_node->children[j]];
    }
    // since we have changed the node at the `i'-th and `old_i'-th positions,
    // we need to update `new2old' mapping
    aux_i = i;
    do {
      aux_i = old2new[aux_i];
    } while (aux_i < i);
    new2old[aux_i] = old_i;
  }

 final_steps:
  free(old2new);
  free(new2old);
  return ret;

}

int crfsuite_tree_init(crfsuite_instance_t* const a_inst)
{
  int i = 0, root_id = -1, node_id, prnt_id;
  int n_items = a_inst->num_items;
  crfsuite_node_t *node_p = NULL;
  crfsuite_item_t *item_p = NULL;

  crfsuite_node_t *tmp_tree = (crfsuite_node_t *) calloc(n_items, sizeof(crfsuite_node_t));
  char *active_nodes = (char *) calloc(n_items, sizeof(char));
  if (tmp_tree == NULL || active_nodes == NULL) {
      fprintf(stderr, "ERROR: Could not allocate memory for tree.\n");
      goto error_exit;
  }
  // iterate over instance items and populate their corresponding nodes
  for (i = 0; i < n_items; ++i) {
    item_p = &a_inst->items[i];
    node_id = item_p->id;
    prnt_id = item_p->prnt;

    if (node_id < 0 || node_id >= n_items) {
      fprintf(stderr, "ERROR: Node id '%d' for label '%s' is out of range "
"(perhaps more labels than tree nodes are present).\n", node_id, item_p->node_label);
      goto error_exit;
    }
    node_p = &tmp_tree[node_id];
    if (active_nodes[node_id]) {
      fprintf(stderr, "ERROR: Duplicate node with label %s\n", item_p->node_label);
      goto error_exit;
    }
    active_nodes[node_id] = 1;
    node_p->self_item_id = i;
    node_p->prnt_node_id = prnt_id;
    if (prnt_id < 0) {
      if (root_id >= 0)
	goto error_exit;
      else
	root_id = node_id;
    }
  }
  /* do nothing for empty instances */
  if (n_items == 0) {
    free(active_nodes);
    return 0;
  }
  /* check tree and reorder its nodes in topological order */
  if (root_id < 0) {
    fprintf(stderr, "ERROR: No root found in tree.  Root node should have parent specified as '_'.\n");
    goto error_exit;
  }
  // once all tree nodes have been populated, we can start populating their
  // parent and children fields
  if (crfsuite_tree_populate_children(tmp_tree, root_id, n_items) != 0)
    goto error_exit;

  // after that, we reorder nodes in topological order
  if (crfsuite_tree_reorder(tmp_tree, root_id, n_items) != 0)
    goto error_exit;

  // remember new id's of tree nodes in items
  for (i = 0; i < n_items; ++i) {
    item_p = &a_inst->items[tmp_tree[i].self_item_id];
    item_p->id = i;
    item_p->prnt = tmp_tree[i].prnt_node_id;
  }
  // assign address of temporary tree to item and return
  a_inst->tree = tmp_tree;
  free(active_nodes);
  return 0;

  // take clean-up actions and return
 error_exit:
  crfsuite_tree_finish(&tmp_tree, i);
  crfsuite_instance_finish(a_inst);
  free(active_nodes);
  return -2;
}


void crfsuite_instance_init(crfsuite_instance_t* inst)
{
  memset(inst, 0, sizeof(*inst));
}

void crfsuite_instance_init_n(crfsuite_instance_t* inst, int num_items)
{
  crfsuite_instance_init(inst);
  inst->num_items = num_items;
  inst->cap_items = num_items;
  inst->items = (crfsuite_item_t*) calloc(num_items, sizeof(crfsuite_item_t));
  inst->labels = (int*) calloc(num_items, sizeof(int));
  /* TODO: What to do with trees? */
}

void crfsuite_instance_finish(crfsuite_instance_t* inst)
{
  int i;

  for (i = 0;i < inst->num_items; ++i)
    crfsuite_item_finish(&inst->items[i]);

  crfsuite_tree_finish(&inst->tree, inst->num_items);
  free(inst->labels);
  free(inst->items);
  crfsuite_instance_init(inst);
}

void crfsuite_instance_copy(crfsuite_instance_t* dst, const crfsuite_instance_t* const src)
{
  int i;
  dst->group = src->group;
  dst->num_items = src->num_items;
  dst->cap_items = src->cap_items;
  dst->items = (crfsuite_item_t*) calloc(dst->num_items, sizeof(crfsuite_item_t));
  if (dst->items == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for items copy.\n");
    goto error_exit;
  }
  dst->labels = (int*) calloc(dst->num_items, sizeof(int));
  if (dst->items == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for labels copy.\n");
    goto error_exit;
  }

  // copy items and labels
  for (i = 0; i < src->num_items; ++i) {
    dst->labels[i] = src->labels[i];
    crfsuite_item_copy(&dst->items[i], &src->items[i]);
  }

  // copy tree if necessary
  if (src->tree) {
    if (crfsuite_tree_copy(dst, src) != 0) {
	fprintf(stderr, "ERROR: Failed to copy the tree.\n");
	goto error_exit;
      }
  } else {
    dst->tree = NULL;
  }
  return;

 error_exit:
  if (dst->items) {free(dst->items); dst->items = NULL;}
  if (dst->labels) {free(dst->labels); dst->labels = NULL;}
  if (dst->tree) crfsuite_tree_finish(&dst->tree, dst->num_items);
}

void crfsuite_instance_swap(crfsuite_instance_t* x, crfsuite_instance_t* y)
{
  crfsuite_instance_t tmp = *x;
  x->num_items = y->num_items;
  x->cap_items = y->cap_items;
  x->items = y->items;
  x->labels = y->labels;
  x->group = y->group;
  x->tree = y->tree;

  y->num_items = tmp.num_items;
  y->cap_items = tmp.cap_items;
  y->items = tmp.items;
  y->labels = tmp.labels;
  y->group = tmp.group;
  y->tree = tmp.tree;
}

int crfsuite_instance_append(crfsuite_instance_t* inst, const crfsuite_item_t* item, int label)
{
  if (inst->cap_items <= inst->num_items) {
    inst->cap_items = (inst->cap_items + 1) * 2;
    inst->labels = (int*) realloc(inst->labels, sizeof(int) * inst->cap_items);
    inst->items = (crfsuite_item_t*) realloc(inst->items, sizeof(crfsuite_item_t) * inst->cap_items);
    /* zero-out new memory */
    memset(inst->items + inst->num_items, 0, sizeof(crfsuite_item_t) * \
	   (inst->cap_items - inst->num_items));
  }
  crfsuite_item_copy(&inst->items[inst->num_items], item);
  inst->labels[inst->num_items] = label;
  ++inst->num_items;
  return 0;
}

int crfsuite_instance_empty(crfsuite_instance_t* inst)
{
  return (inst->num_items == 0);
}


void crfsuite_data_init(crfsuite_data_t* data)
{
  memset(data, 0, sizeof(*data));
}

void crfsuite_data_init_n(crfsuite_data_t* data, int n)
{
  crfsuite_data_init(data);
  data->num_instances = n;
  data->cap_instances = n;
  data->instances = (crfsuite_instance_t*)calloc(n, sizeof(crfsuite_instance_t));
}

void crfsuite_data_finish(crfsuite_data_t* data)
{
  int i;

  for (i = 0;i < data->num_instances;++i) {
    crfsuite_instance_finish(&data->instances[i]);
  }
  free(data->instances);
  crfsuite_data_init(data);
}

void crfsuite_data_copy(crfsuite_data_t* dst, const crfsuite_data_t* src)
{
  int i;

  dst->num_instances = src->num_instances;
  dst->cap_instances = src->cap_instances;
  dst->instances = (crfsuite_instance_t*) calloc(dst->num_instances, sizeof(crfsuite_instance_t));
  for (i = 0;i < dst->num_instances;++i)
    crfsuite_instance_copy(&dst->instances[i], &src->instances[i]);
}

void crfsuite_data_swap(crfsuite_data_t* x, crfsuite_data_t* y)
{
  crfsuite_data_t tmp = *x;
  x->num_instances = y->num_instances;
  x->cap_instances = y->cap_instances;
  x->instances = y->instances;
  y->num_instances = tmp.num_instances;
  y->cap_instances = tmp.cap_instances;
  y->instances = tmp.instances;
}

int crfsuite_data_append(crfsuite_data_t* data, const crfsuite_instance_t* inst)
{
  if (0 < inst->num_items) {
    if (data->cap_instances <= data->num_instances) {
      data->cap_instances = (data->cap_instances + 1) * 2;
      data->instances = (crfsuite_instance_t*)realloc(
						      data->instances, sizeof(crfsuite_instance_t) * data->cap_instances);
      memset(&data->instances[data->num_instances], 0, sizeof(crfsuite_instance_t) * \
	     (data->cap_instances - data->num_instances));
    }
    crfsuite_instance_copy(&data->instances[data->num_instances++], inst);
  }
  return 0;
}

int crfsuite_data_maxlength(crfsuite_data_t* data)
{
  int i, T = 0;
  for (i = 0;i < data->num_instances;++i) {
    if (T < data->instances[i].num_items) {
      T = data->instances[i].num_items;
    }
  }
  return T;
}

int  crfsuite_data_totalitems(crfsuite_data_t* data)
{
  int i, n = 0;
  for (i = 0;i < data->num_instances;++i) {
    n += data->instances[i].num_items;
  }
  return n;
}

void crfsuite_evaluation_init(crfsuite_evaluation_t* eval, int n)
{
  memset(eval, 0, sizeof(*eval));
  eval->tbl = (crfsuite_label_evaluation_t*)calloc(n+1, sizeof(crfsuite_label_evaluation_t));
  if (eval->tbl != NULL) {
    eval->num_labels = n;
  }
}

void crfsuite_evaluation_clear(crfsuite_evaluation_t* eval)
{
  int i;
  for (i = 0;i <= eval->num_labels;++i) {
    memset(&eval->tbl[i], 0, sizeof(eval->tbl[i]));
  }

  eval->item_total_correct = 0;
  eval->item_total_num = 0;
  eval->item_total_model = 0;
  eval->item_total_observation = 0;
  eval->item_accuracy = 0;

  eval->inst_total_correct = 0;
  eval->inst_total_num = 0;
  eval->inst_accuracy = 0;

  eval->macro_precision = 0;
  eval->macro_recall = 0;
  eval->macro_fmeasure = 0;
}

void crfsuite_evaluation_finish(crfsuite_evaluation_t* eval)
{
  free(eval->tbl);
  memset(eval, 0, sizeof(*eval));
}

int crfsuite_evaluation_accmulate(crfsuite_evaluation_t* eval, const int* reference, const int* prediction, int T)
{
  int t, nc = 0;

  for (t = 0;t < T;++t) {
    int lr = reference[t];
    int lt = prediction[t];

    if (eval->num_labels <= lr || eval->num_labels <= lt) {
      return 1;
    }

    ++eval->tbl[lr].num_observation;
    ++eval->tbl[lt].num_model;
    if (lr == lt) {
      ++eval->tbl[lr].num_correct;
      ++nc;
    }
    ++eval->item_total_num;
  }

  if (nc == T) {
    ++eval->inst_total_correct;
  }
  ++eval->inst_total_num;

  return 0;
}

void crfsuite_evaluation_finalize(crfsuite_evaluation_t* eval)
{
  int i;

  for (i = 0;i <= eval->num_labels;++i) {
    crfsuite_label_evaluation_t* lev = &eval->tbl[i];

    /* Do not evaluate labels that does not in the test data. */
    if (lev->num_observation == 0) {
      continue;
    }

    /* Sum the number of correct labels for accuracy calculation. */
    eval->item_total_correct += lev->num_correct;
    eval->item_total_model += lev->num_model;
    eval->item_total_observation += lev->num_observation;

    /* Initialize the precision, recall, and f1-measure values. */
    lev->precision = 0;
    lev->recall = 0;
    lev->fmeasure = 0;

    /* Compute the precision, recall, and f1-measure values. */
    if (lev->num_model > 0) {
      lev->precision = lev->num_correct / (double)lev->num_model;
    }
    if (lev->num_observation > 0) {
      lev->recall = lev->num_correct / (double)lev->num_observation;
    }
    if (lev->precision + lev->recall > 0) {
      lev->fmeasure = lev->precision * lev->recall * 2 / (lev->precision + lev->recall);
    }

    /* Exclude unknown labels from calculation of macro-average values. */
    if (i != eval->num_labels) {
      eval->macro_precision += lev->precision;
      eval->macro_recall += lev->recall;
      eval->macro_fmeasure += lev->fmeasure;
    }
  }

  /* Copute the macro precision, recall, and f1-measure values. */
  eval->macro_precision /= eval->num_labels;
  eval->macro_recall /= eval->num_labels;
  eval->macro_fmeasure /= eval->num_labels;

  /* Compute the item accuracy. */
  eval->item_accuracy = 0;
  if (0 < eval->item_total_num) {
    eval->item_accuracy = eval->item_total_correct / (double)eval->item_total_num;
  }

  /* Compute the instance accuracy. */
  eval->inst_accuracy = 0;
  if (0 < eval->inst_total_num) {
    eval->inst_accuracy = eval->inst_total_correct / (double)eval->inst_total_num;
  }
}

void crfsuite_evaluation_output(crfsuite_evaluation_t* eval, crfsuite_dictionary_t* labels, crfsuite_logging_callback cbm, void *instance)
{
  int i;
  const char *lstr = NULL;
  logging_t lg;

  lg.func = cbm;
  lg.instance = instance;

  logging(&lg, "Performance by label (#match, #model, #ref) (precision, recall, F1):\n");

  for (i = 0;i < eval->num_labels;++i) {
    const crfsuite_label_evaluation_t* lev = &eval->tbl[i];

    labels->to_string(labels, i, &lstr);
    if (lstr == NULL) lstr = "[UNKNOWN]";

    if (lev->num_observation == 0) {
      logging(&lg, "    %s: (%d, %d, %d) (******, ******, ******)\n",
	      lstr, lev->num_correct, lev->num_model, lev->num_observation
	      );
    } else {
      logging(&lg, "    %s: (%d, %d, %d) (%1.4f, %1.4f, %1.4f)\n",
	      lstr, lev->num_correct, lev->num_model, lev->num_observation,
	      lev->precision, lev->recall, lev->fmeasure
	      );
    }
    labels->free(labels, lstr);
  }
  logging(&lg, "Macro-average precision, recall, F1: (%f, %f, %f)\n",
	  eval->macro_precision, eval->macro_recall, eval->macro_fmeasure
	  );
  logging(&lg, "Item accuracy: %d / %d (%1.4f)\n",
	  eval->item_total_correct, eval->item_total_num, eval->item_accuracy
	  );
  logging(&lg, "Instance accuracy: %d / %d (%1.4f)\n",
	  eval->inst_total_correct, eval->inst_total_num, eval->inst_accuracy
	  );
}

int crfsuite_interlocked_increment(int *count)
{
  return ++(*count);
}

int crfsuite_interlocked_decrement(int *count)
{
  return --(*count);
}
