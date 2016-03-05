/*
 *      CRFsuite C++/SWIG API wrapper.
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

#ifndef __CRFSUITE_HPP__
#define __CRFSUITE_HPP__

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>

#include "crfsuite_api.hpp"

namespace CRFSuite
{

  Trainer::Trainer()
  {
    data = new crfsuite_data_t;
    if (data != NULL) {
      crfsuite_data_init(data);
    }
    tr = NULL;
  }

  Trainer::~Trainer()
  {
    if (data) {
      clear();
      delete data;
      data = NULL;
    }
    if (tr) {
      tr->release(tr);
      tr = NULL;
    }
  }

  void Trainer::init()
  {
    // Create an instance of attribute dictionary.
    if (data->attrs == NULL) {
      int ret = crfsuite_create_instance("dictionary", (void**)&data->attrs);
      if (!ret) {
	throw std::runtime_error("Failed to create a dictionary instance for attributes.");
      }
    }

    // Create an instance of label dictionary.
    if (data->labels == NULL) {
      int ret = crfsuite_create_instance("dictionary", (void**)&data->labels);
      if (!ret) {
	throw std::runtime_error("Failed to create a dictionary instance for labels.");
      }
    }

    // Create an instance of node label dictionary.
    if (tr && tr->ftype == FTYPE_CRF1TREE && data->node_labels == NULL) {
      int ret = crfsuite_create_instance("dictionary", (void**)&data->node_labels);
      if (!ret) {
	throw std::runtime_error("Failed to create a dictionary instance for node labels.");
      }
    }
  }

  void Trainer::clear()
  {
    if (data != NULL) {
      if (data->labels != NULL) {
	data->labels->release(data->labels);
	data->labels = NULL;
      }

      if (data->attrs != NULL) {
	data->attrs->release(data->attrs);
	data->attrs = NULL;
      }

      if (data->node_labels != NULL) {
	data->node_labels->release(data->node_labels);
	data->node_labels = NULL;
      }
      // crfsuite_data_init() is automatically called from `finish()`
      crfsuite_data_finish(data);
    }
  }

  void Trainer::append(const ItemSequence& xseq, const StringList& yseq, int group)
  {
    // Create dictionary objects if necessary.
    if (data->attrs == NULL || data->labels == NULL || \
	(tr && tr->ftype == FTYPE_CRF1TREE && data->node_labels == NULL))
      init();

    // Make sure |y| == |x|.
    if (xseq.size() != yseq.size()) {
      std::stringstream ss;
      ss << "The numbers of items and labels differ: |x| = " << \
	xseq.size() << ", |y| = " << yseq.size();
      throw std::invalid_argument(ss.str());
    }

    // Convert instance_type to crfsuite_instance_t.
    int i, n_items;
    crfsuite_instance_t _inst;
    crfsuite_instance_init_n(&_inst, xseq.size());
    for (size_t t = 0; t < xseq.size(); ++t) {
      const Item& item = xseq[t];
      crfsuite_item_t* _item = &_inst.items[t];

      // Set the attributes in the item.
      i = 0;
      n_items = (int) item.size();

      if (tr && tr->ftype == FTYPE_CRF1TREE) {
      	if (n_items < 1)
      	  throw std::runtime_error("Invalid tree format: node label should be given as attribute.");
      	else if (n_items < 2)
      	  throw std::runtime_error("Invalid tree format: parent label should be given as attribute.");

      	// Allocate memory for attributes
      	i = 2;
      	crfsuite_item_init_n(_item, n_items - i);

      	// remember string label of this node
      	_item->id = data->node_labels->get(data->node_labels, item[0].attr.c_str());
      	_item->node_label = (char *) calloc((item[0].attr.length() + 1), sizeof(char));

      	if (_item->node_label)
      	  strcpy(_item->node_label, item[0].attr.c_str());
      	else
      	  throw std::runtime_error("ERROR: Could not allocate memory for storing node label.\n");

      	// remember parent of this node
      	if (item[1].attr == "_")
      	  _item->prnt = -1;
      	else
      	  _item->prnt = data->node_labels->get(data->node_labels, item[1].attr.c_str());
      } else {
      	crfsuite_item_init_n(_item, n_items);
      }

      // add attributes
      for (; i < n_items; ++i) {
	if (item[i].attr.empty())
	  continue;

	_item->contents[i].aid = data->attrs->get(data->attrs, item[i].attr.c_str());
	_item->contents[i].value = (floatval_t)item[i].value;
      }

      // Set the label of the item.
      _inst.labels[t] = data->labels->get(data->labels, yseq[t].c_str());
    }

    // initialize instance tree
    if (tr && tr->ftype == FTYPE_CRF1TREE && crfsuite_tree_init(&_inst) != 0)
      throw std::runtime_error("ERROR: Could not create tree for training instance.\n");

    // assign instance to a group
    _inst.group = group;

    // Append the instance to the training set.
    crfsuite_data_append(data, &_inst);

    // Finish the instance.
    crfsuite_instance_finish(&_inst);

    /* clear dictionary of node labels so that new instances will
       have dense representation of node ids again */
    if (data->node_labels)
      data->node_labels->reset(data->node_labels);
  }

  bool Trainer::select(const std::string& algorithm, const std::string& type)
  {
    int ret = 0;

    // Release the trainer if it is already initialized.
    if (tr != NULL) {
      tr->release(tr);
      tr = NULL;
    }

    if (algorithm != "lbfgs" && type == "semim") {
      std::stringstream ss;
      ss << "ERROR: Training algorithm '" << algorithm << \
	"' is not supported for this type of graphical model.  Try `lbfgs' instead";
      throw std::invalid_argument(ss.str());
    }

    // Build the trainer string ID.
    std::string tid = "train/";
    tid += type;
    tid += '/';
    tid += algorithm;

    // Create an instance of a trainer.
    ret = crfsuite_create_instance(tid.c_str(), (void**)&tr);
    if (!ret)
      return false;

    // Set the callback function for receiving messages.
    tr->set_message_callback(tr, this, __logging_callback);
    return true;
  }

  int Trainer::train(const std::string& model, int holdout)
  {
    // Run the training algorithm.
    return tr->train(tr, data, model.empty()? NULL: model.c_str(), holdout);
  }

  StringList Trainer::params()
  {
    StringList pars;
    if (!tr)
      return pars;

    crfsuite_params_t* params = tr->params(tr);
    int n = params->num(params);
    for (int i = 0; i < n; ++i) {
      char *name = NULL;
      params->name(params, i, &name);
      pars.push_back(name);
      params->free(params, name);
    }
    return pars;
  }

  void Trainer::set(const std::string& name, const std::string& value)
  {
    crfsuite_params_t* params = tr->params(tr);
    if (params->set(params, name.c_str(), value.c_str()) != 0) {
      std::stringstream ss;
      ss << "Parameter not found: " << name << " = " << value;
      params->release(params);
      throw std::invalid_argument(ss.str());
    }
    params->release(params);
  }

  std::string Trainer::get(const std::string& name)
  {
    std::string value;
    char *_value = NULL;
    crfsuite_params_t* params = tr->params(tr);
    if (params->get(params, name.c_str(), &_value) != 0) {
      std::stringstream ss;
      ss << "Parameter not found: " << name << " = " << value;
      params->release(params);
      throw std::invalid_argument(ss.str());
    }
    value = _value;
    params->free(params, _value);
    params->release(params);
    return value;
  }

  std::string Trainer::help(const std::string& name)
  {
    std::string str;
    crfsuite_params_t* params = tr->params(tr);
    char *_str = NULL;
    params->help(params, name.c_str(), NULL, &_str);
    str = _str;
    params->free(params, _str);
    params->release(params);
    return str;
  }

  void Trainer::message(const std::string& msg)
  {
  }

  int Trainer::__logging_callback(void *instance, const char *format, va_list args)
  {
    char buffer[65536];
    vsnprintf(buffer, sizeof(buffer)-1, format, args);
    reinterpret_cast<Trainer*>(instance)->message(buffer);
    return 0;
  }


  // all members are deault initialized in class
  Tagger::Tagger():
    m_ftype(FTYPE_NONE)
  {}

  Tagger::~Tagger()
  {
    this->close();
  }

  bool Tagger::open(const std::string& name, const int ftype)
  {
    m_ftype = ftype;
    int ret;

    // Close the model if it is already opened.
    this->close();

    // Open the model file.
    if ((ret = crfsuite_create_instance_from_file(name.c_str(), (void**)&model, m_ftype))) {
      return false;
    }

    // Obtain the tagger interface.
    if ((ret = model->get_tagger(model, &tagger))) {
      throw std::runtime_error("Failed to obtain the tagger interface");
    }

    // Obtain the dictionary interface representing the attributes in the model.
    if ((ret = model->get_attrs(model, &m_attrs))) {
      throw std::runtime_error("Failed to obtain the dictionary interface for attributes");
    }

    // Obtain the dictionary interface representing the labels in the model.
    if ((ret = model->get_labels(model, &m_labels))) {
      throw std::runtime_error("Failed to obtain the dictionary interface for labels");
    }

    // initialize auxiliary member
    if (m_ftype == FTYPE_CRF1TREE) {
      if (!(ret = crfsuite_create_instance("dictionary", (void**) &m_node_labels)))
	throw std::runtime_error("Failed to create a dictionary instance for tree labels.");
    } else if (m_ftype == FTYPE_SEMIMCRF) {
      model->get_sm(model, &m_aux);
    }

    return true;
  }

  void Tagger::close()
  {
    if (tagger != NULL) {
      tagger->release(tagger);
      tagger = NULL;
    }
    if (model != NULL) {
      model->release(model);
      model = NULL;
    }
    if (m_labels != NULL) {
      m_labels->release(m_labels);
      m_labels = NULL;
    }
    if (m_node_labels != NULL) {
      m_node_labels->release(m_node_labels);
      m_node_labels = NULL;
    }
    if (m_attrs != NULL) {
      m_attrs->release(m_attrs);
      m_attrs = NULL;
    }
  }

  StringList Tagger::tag(ItemSequence& xseq)
  {
    set(xseq);
    return viterbi();
  }

  void Tagger::set(ItemSequence& xseq)
  {
    int ret;
    StringList yseq;
    crfsuite_instance_t _inst;
    const bool itree = (m_ftype == FTYPE_CRF1TREE);

    if (model == NULL || tagger == NULL || m_attrs == NULL) {
      throw std::invalid_argument("Tagger is not opened.");
    }


    // Build an instance.
    crfsuite_instance_init_n(&_inst, xseq.size());
    for (size_t t = 0; t < xseq.size(); ++t) {
      const Item& item = xseq[t];
      crfsuite_item_t* _item = &_inst.items[t];

      // Set the attributes in the item.
      crfsuite_item_init(_item);
      size_t i = 0;

      if (itree) {
	if (item.size() < 3)
	  throw std::runtime_error(\
				   "Invalid input format (2-nd and 3-rd" \
				   " fields should denote current and parent labels).");
	// set node id of the current item
	_item->id = m_node_labels->get(m_node_labels, item[1].attr.c_str());
	// set parent id of the current item
	if (item[2].attr.compare("_") == 0)
	  _item->prnt = -1;
	else
	  _item->prnt = m_node_labels->get(m_node_labels, item[2].attr.c_str());

	i = 3;
      }

      for (; i < item.size(); ++i) {
	int aid = m_attrs->to_id(m_attrs, item[i].attr.c_str());
	if (0 <= aid) {
	  crfsuite_attribute_t cont;
	  crfsuite_attribute_set(&cont, aid, item[i].value);
	  crfsuite_item_append_attribute(_item, &cont);
	}
      }
    }
    // create tree instance for tree-structired model
    if (itree) {
      if ((ret = crfsuite_tree_init(&_inst)) != 0)
        throw std::runtime_error("Could not create tree instance for sequence.");

      m_aux = (const void *) _inst.tree;
    }

    // Set the instance to the tagger.
    if ((ret = tagger->set(tagger, &_inst))) {
      crfsuite_instance_finish(&_inst);
      throw std::runtime_error("Failed to set the instance to the tagger.");
    }
    crfsuite_instance_finish(&_inst);

    // reset labels of tree node
    if (m_node_labels)
      m_node_labels->reset(m_node_labels);
  }

  StringList Tagger::viterbi()
  {
    if (model == NULL || tagger == NULL)
      throw std::invalid_argument("The tagger is not opened");

    int ret;
    StringList yseq;
    // Make sure that the current instance is not empty.
    const size_t T = (size_t)tagger->length(tagger);
    if (T <= 0)
      return yseq;

    // Run the Viterbi algorithm.
    floatval_t score;
    int *path = new int[T];
    if ((ret = tagger->viterbi(tagger, path, &score, m_aux))) {
      delete[] path;
      throw std::runtime_error("Failed to find the Viterbi path.");
    }

    // Convert the Viterbi path to a label sequence.
    yseq.resize(T);
    for (size_t t = 0; t < T; ++t) {
      const char *label = NULL;
      if (m_labels->to_string(m_labels, path[t], &label) != 0) {
	delete[] path;
	throw std::runtime_error("Failed to convert a label identifier to string.");
      }
      yseq[t] = label;
      m_labels->free(m_labels, label);
    }
    delete[] path;
    return yseq;
  }

  double Tagger::probability(StringList& yseq)
  {
    int ret;
    size_t T;
    int *path = NULL;
    std::stringstream msg;
    floatval_t score, lognorm;
    crfsuite_dictionary_t *labels = NULL;

    if (model == NULL || tagger == NULL) {
      msg << "The tagger is not opened";
      throw std::invalid_argument(msg.str());
    }

    // Make sure that the current instance is not empty.
    T = (size_t)tagger->length(tagger);
    if (T <= 0)
      return 0.;

    // Make sure that |y| == |x|.
    if (yseq.size() != T) {
      msg << "The numbers of items and labels differ: |x| = " << T << ", |y| = " << yseq.size();
      throw std::invalid_argument(msg.str());
    }

    // Obtain the dictionary interface representing the labels in the model.
    if ((ret = model->get_labels(model, &labels))) {
      msg << "Failed to obtain the dictionary interface for labels";
      goto error_exit;
    }

    // Convert string labels into label IDs.
    path = new int[T];
    for (size_t t = 0;t < T;++t) {
      int l = labels->to_id(labels, yseq[t].c_str());
      if (l < 0) {
	msg << "Failed to convert into label identifier: " << yseq[t];
	goto error_exit;
      }
      path[t] = l;
    }

    // Compute the score of the path.
    if ((ret = tagger->score(tagger, path, &score))) {
      msg << "Failed to score the label sequence";
      goto error_exit;
    }

    // Compute the partition factor.
    if ((ret = tagger->lognorm(tagger, &lognorm, m_aux))) {
      msg << "Failed to compute the partition factor";
      goto error_exit;
    }

    labels->release(labels);
    delete[] path;
    return std::exp((double)(score - lognorm));

  error_exit:
    if (labels != NULL) {
      labels->release(labels);
      labels = NULL;
    }
    delete[] path;
    throw std::runtime_error(msg.str());
  }

  double Tagger::marginal(const std::string& y, const int t)
  {
    int l, ret, T;
    floatval_t prob;
    std::stringstream msg;
    crfsuite_dictionary_t *labels = NULL;

    if (model == NULL || tagger == NULL) {
      msg << "The tagger is not opened";
      throw std::invalid_argument(msg.str());
    }

    // Make sure that the current instance is not empty.
    T = tagger->length(tagger);
    if (T <= 0) {
      return 0.;
    }

    // Make sure that 0 <= t < |x|.
    if (t < 0 || T <= t) {
      msg << "The position, " << t << "is out of range of " << T;
      throw std::invalid_argument(msg.str());
    }

    // Obtain the dictionary interface representing the labels in the model.
    if ((ret = model->get_labels(model, &labels))) {
      msg << "Failed to obtain the dictionary interface for labels";
      goto error_exit;
    }

    // Convert string labels into label IDs.
    l = labels->to_id(labels, y.c_str());
    if (l < 0) {
      msg << "Failed to convert into label identifier: " << y;
      goto error_exit;
    }

    // Compute the score of the path.
    if ((ret = tagger->marginal_point(tagger, l, t, &prob, m_aux))) {
      msg << "Failed to compute the marginal probability of '" << y << "' at " << t;
      goto error_exit;
    }

    labels->release(labels);
    return prob;

  error_exit:
    if (labels != NULL) {
      labels->release(labels);
      labels = NULL;
    }
    throw std::runtime_error(msg.str());
  }


  std::string version()
  {
    return CRFSUITE_VERSION;
  }

};

#endif/*__CRFSUITE_HPP__*/

