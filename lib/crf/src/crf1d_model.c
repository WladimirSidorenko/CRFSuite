/*
 *      CRF1d model.
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

#include <string.h>
#include <crfsuite.h>

#include "crf1d.h"
#include "os.h"

#define FILEMAGIC       "lCRF"
#define MODELTYPE_TREE  "TREE"
#define MODELTYPE_CRF1D "FOMC"
#define MODELTYPE_SEMIM "HOSM"
#define VERSION_NUMBER  (101)
#define CHUNK_LABELREF  "LFRF"
#define CHUNK_ATTRREF   "AFRF"
#define CHUNK_FEATURE   "FEAT"
#define CHUNK_HSMM      "HSMM"
#define CHUNK_HSMM_SIZE 4
#define HEADER_SIZE     52
#define CHUNK_SIZE      12
#define FEATURE_SIZE    20
#define HSM_CHUNK_SIZE  32

enum {
  WSTATE_NONE,
  WSTATE_LABELS,
  WSTATE_ATTRS,
  WSTATE_LABELREFS,
  WSTATE_ATTRREFS,
  WSTATE_FEATURES,
  WSTATE_SM,
};

enum {
  KT_GLOBAL = 'A',
  KT_NUMATTRS,
  KT_NUMLABELS,
  KT_STR2LID,
  KT_LID2STR,
  KT_STR2AID,
  KT_FEATURE,
};

static int write_uint8(FILE *fp, uint8_t value)
{
  return fwrite(&value, sizeof(value), 1, fp) == 1 ? 0 : 1;
}

static int read_uint8(uint8_t* buffer, uint8_t* value)
{
  *value = *buffer;
  return sizeof(*value);
}

static int write_uint32(FILE *fp, uint32_t value)
{
  uint8_t buffer[4];
  buffer[0] = (uint8_t)(value & 0xFF);
  buffer[1] = (uint8_t)(value >> 8);
  buffer[2] = (uint8_t)(value >> 16);
  buffer[3] = (uint8_t)(value >> 24);
  return fwrite(buffer, sizeof(uint8_t), 4, fp) == 4 ? 0 : 1;
}

static int read_uint32(uint8_t* buffer, uint32_t* value)
{
  *value  = ((uint32_t)buffer[0]);
  *value |= ((uint32_t)buffer[1] << 8);
  *value |= ((uint32_t)buffer[2] << 16);
  *value |= ((uint32_t)buffer[3] << 24);
  return sizeof(*value);
}

static int write_uint8_array(FILE *fp, uint8_t *array, size_t n)
{
  size_t i;
  int ret = 0;
  for (i = 0;i < n;++i) {
    ret |= write_uint8(fp, array[i]);
  }
  return ret;
}

static int read_uint8_array(uint8_t* buffer, uint8_t *array, size_t n)
{
  size_t i;
  int ret = 0;
  for (i = 0;i < n;++i) {
    int size = read_uint8(buffer, &array[i]);
    buffer += size;
    ret += size;
  }
  return ret;
}

static void write_float(FILE *fp, floatval_t value)
{
  /*
    We assume:
    - sizeof(floatval_t) = sizeof(double) = sizeof(uint64_t)
    - the byte order of floatval_t and uint64_t is the same
    - ARM's mixed-endian is not supported
  */
  uint64_t iv;
  uint8_t buffer[8];

  /* Copy the memory image of floatval_t value to uint64_t. */
  memcpy(&iv, &value, sizeof(iv));

  buffer[0] = (uint8_t)(iv & 0xFF);
  buffer[1] = (uint8_t)(iv >> 8);
  buffer[2] = (uint8_t)(iv >> 16);
  buffer[3] = (uint8_t)(iv >> 24);
  buffer[4] = (uint8_t)(iv >> 32);
  buffer[5] = (uint8_t)(iv >> 40);
  buffer[6] = (uint8_t)(iv >> 48);
  buffer[7] = (uint8_t)(iv >> 56);
  fwrite(buffer, sizeof(uint8_t), 8, fp);
}

static int read_float(uint8_t* buffer, floatval_t* value)
{
  uint64_t iv;
  iv  = ((uint64_t)buffer[0]);
  iv |= ((uint64_t)buffer[1] << 8);
  iv |= ((uint64_t)buffer[2] << 16);
  iv |= ((uint64_t)buffer[3] << 24);
  iv |= ((uint64_t)buffer[4] << 32);
  iv |= ((uint64_t)buffer[5] << 40);
  iv |= ((uint64_t)buffer[6] << 48);
  iv |= ((uint64_t)buffer[7] << 56);
  memcpy(value, &iv, sizeof(*value));
  return sizeof(*value);
}

crf1dmw_t* crf1mmw(const char *filename, const int ftype)
{
  header_t *header = NULL;
  crf1dmw_t *writer = NULL;

  /* Create a writer instance. */
  writer = (crf1dmw_t*) calloc(1, sizeof(crf1dmw_t));
  if (writer == NULL) {
    goto error_exit;
  }

  /* Open the file for writing. */
  writer->fp = fopen(filename, "wb");
  if (writer->fp == NULL) {
    goto error_exit;
  }

  /* Fill the members in the header. */
  header = &writer->header;
  strncpy(header->magic, FILEMAGIC, 4);
  if (ftype == FTYPE_CRF1TREE)
    strncpy(header->type, MODELTYPE_TREE, 4);
  else if (ftype == FTYPE_SEMIMCRF)
    strncpy(header->type, MODELTYPE_SEMIM, 4);
  else
    strncpy(header->type, MODELTYPE_CRF1D, 4);
  header->version = VERSION_NUMBER;

  /* Advance the file position to skip the file header. */
  if (fseek(writer->fp, HEADER_SIZE, SEEK_CUR) != 0) {
    goto error_exit;
  }

  return writer;

 error_exit:
  if (writer != NULL) {
    if (writer->fp != NULL) {
      fclose(writer->fp);
    }
    free(writer);
  }
  return NULL;
}

int crf1dmw_close(crf1dmw_t* writer)
{
  FILE *fp = writer->fp;
  header_t *header = &writer->header;

  /* Store the file size. */
  header->size = (uint32_t)ftell(fp);

  /* Move the file position to the head. */
  if (fseek(fp, 0, SEEK_SET) != 0) {
    goto error_exit;
  }

  /* Write the file header. */
  write_uint8_array(fp, header->magic, sizeof(header->magic));
  write_uint32(fp, header->size);
  write_uint8_array(fp, header->type, sizeof(header->type));
  write_uint32(fp, header->version);
  write_uint32(fp, header->num_features);
  write_uint32(fp, header->num_labels);
  write_uint32(fp, header->num_attrs);
  write_uint32(fp, header->off_features);
  write_uint32(fp, header->off_labels);
  write_uint32(fp, header->off_attrs);
  write_uint32(fp, header->off_labelrefs);
  write_uint32(fp, header->off_attrrefs);
  write_uint32(fp, header->off_sm);

  /* Check for any error occurrence. */
  if (ferror(fp)) {
    goto error_exit;
  }

  /* Close the writer. */
  fclose(fp);
  free(writer);
  return 0;

 error_exit:
  if (writer != NULL) {
    if (writer->fp != NULL) {
      fclose(writer->fp);
    }
    free(writer);
  }
  return 1;
}

int crf1dmw_open_labels(crf1dmw_t* writer, int num_labels)
{
  /* Check if we aren't writing anything at this moment. */
  if (writer->state != WSTATE_NONE) {
    return 1;
  }

  /* Store the current offset. */
  writer->header.off_labels = (uint32_t) ftell(writer->fp);

  /* Open a CQDB chunk for writing. */
  writer->dbw = cqdb_writer(writer->fp, 0);
  if (writer->dbw == NULL) {
    writer->header.off_labels = 0;
    return 1;
  }

  writer->state = WSTATE_LABELS;
  writer->header.num_labels = num_labels;
  return 0;
}

int crf1dmw_close_labels(crf1dmw_t* writer)
{
  /* Make sure that we are writing labels. */
  if (writer->state != WSTATE_LABELS) {
    return 1;
  }

  /* Close the CQDB chunk. */
  if (cqdb_writer_close(writer->dbw)) {
    return 1;
  }

  writer->dbw = NULL;
  writer->state = WSTATE_NONE;
  return 0;
}

int crf1dmw_put_label(crf1dmw_t* writer, int lid, const char *value)
{
  /* Make sure that we are writing labels. */
  if (writer->state != WSTATE_LABELS) {
    return 1;
  }

  /* Put the label. */
  if (cqdb_writer_put(writer->dbw, value, lid)) {
    return 1;
  }

  return 0;
}

int crf1dmw_open_attrs(crf1dmw_t* writer, int num_attrs)
{
  /* Check if we aren't writing anything at this moment. */
  if (writer->state != WSTATE_NONE) {
    return 1;
  }

  /* Store the current offset. */
  writer->header.off_attrs = (uint32_t)ftell(writer->fp);

  /* Open a CQDB chunk for writing. */
  writer->dbw = cqdb_writer(writer->fp, 0);
  if (writer->dbw == NULL) {
    writer->header.off_attrs = 0;
    return 1;
  }

  writer->state = WSTATE_ATTRS;
  writer->header.num_attrs = num_attrs;
  return 0;
}

int crf1dmw_close_attrs(crf1dmw_t* writer)
{
  /* Make sure that we are writing attributes. */
  if (writer->state != WSTATE_ATTRS) {
    return 1;
  }

  /* Close the CQDB chunk. */
  if (cqdb_writer_close(writer->dbw)) {
    return 1;
  }

  writer->dbw = NULL;
  writer->state = WSTATE_NONE;
  return 0;
}

int crf1dmw_put_attr(crf1dmw_t* writer, int aid, const char *value)
{
  /* Make sure that we are writing labels. */
  if (writer->state != WSTATE_ATTRS) {
    return 1;
  }

  /* Put the attribute. */
  if (cqdb_writer_put(writer->dbw, value, aid)) {
    return 1;
  }

  return 0;
}

int crf1dmw_open_labelrefs(crf1dmw_t* writer, int num_labels)
{
  uint32_t offset;
  FILE *fp = writer->fp;
  featureref_header_t* href = NULL;
  size_t size = CHUNK_SIZE + sizeof(uint32_t) * num_labels;

  /* Check if we aren't writing anything at this moment. */
  if (writer->state != WSTATE_NONE) {
    return CRFSUITEERR_INTERNAL_LOGIC;
  }

  /* Allocate a feature reference array. */
  href = (featureref_header_t*) calloc(size, 1);
  if (href == NULL) {
    return CRFSUITEERR_OUTOFMEMORY;
  }

  /* Align the offset to a DWORD boundary. */
  offset = (uint32_t)ftell(fp);
  while (offset % 4 != 0) {
    uint8_t c = 0;
    fwrite(&c, sizeof(uint8_t), 1, fp);
    ++offset;
  }

  /* Store the current offset position to the file header. */
  writer->header.off_labelrefs = offset;
  fseek(fp, size, SEEK_CUR);

  /* Fill members in the feature reference header. */
  strncpy(href->chunk, CHUNK_LABELREF, 4);
  href->size = 0;
  href->num = num_labels;

  writer->href = href;
  writer->state = WSTATE_LABELREFS;
  return 0;
}

int crf1dmw_close_labelrefs(crf1dmw_t* writer)
{
  uint32_t i;
  FILE *fp = writer->fp;
  featureref_header_t* href = writer->href;
  uint32_t begin = writer->header.off_labelrefs, end = 0;

  /* Make sure that we are writing label feature references. */
  if (writer->state != WSTATE_LABELREFS) {
    return CRFSUITEERR_INTERNAL_LOGIC;
  }

  /* Store the current offset position. */
  end = (uint32_t)ftell(fp);

  /* Compute the size of this chunk. */
  href->size = (end - begin);

  /* Write the chunk header and offset array. */
  fseek(fp, begin, SEEK_SET);
  write_uint8_array(fp, href->chunk, 4);
  write_uint32(fp, href->size);
  write_uint32(fp, href->num);
  for (i = 0;i < href->num;++i) {
    write_uint32(fp, href->offsets[i]);
  }

  /* Move the file pointer to the tail. */
  fseek(fp, end, SEEK_SET);

  /* Uninitialize. */
  free(href);
  writer->href = NULL;
  writer->state = WSTATE_NONE;
  return 0;
}

int crf1dmw_put_labelref(crf1dmw_t* writer, int lid, const feature_refs_t* ref, int *map)
{
  int i, fid;
  uint32_t n = 0;
  FILE *fp = writer->fp;
  featureref_header_t* href = writer->href;

  /* Make sure that we are writing label feature references. */
  if (writer->state != WSTATE_LABELREFS) {
    return CRFSUITEERR_INTERNAL_LOGIC;
  }

  /* Store the current offset to the offset array. */
  href->offsets[lid] = ftell(fp);

  /* Count the number of references to active features. */
  for (i = 0;i < ref->num_features;++i) {
    if (0 <= map[ref->fids[i]]) ++n;
  }

  /* Write the feature reference. */
  write_uint32(fp, (uint32_t)n);
  for (i = 0;i < ref->num_features;++i) {
    fid = map[ref->fids[i]];
    if (0 <= fid) write_uint32(fp, (uint32_t)fid);
  }

  return 0;
}

int crf1dmw_open_attrrefs(crf1dmw_t* writer, int num_attrs)
{
  uint32_t offset;
  FILE *fp = writer->fp;
  featureref_header_t* href = NULL;
  size_t size = CHUNK_SIZE + sizeof(uint32_t) * num_attrs;

  /* Check if we aren't writing anything at this moment. */
  if (writer->state != WSTATE_NONE) {
    return CRFSUITEERR_INTERNAL_LOGIC;
  }

  /* Allocate a feature reference array. */
  href = (featureref_header_t*)calloc(size, 1);
  if (href == NULL) {
    return CRFSUITEERR_OUTOFMEMORY;
  }

  /* Align the offset to a DWORD boundary. */
  offset = (uint32_t)ftell(fp);
  while (offset % 4 != 0) {
    uint8_t c = 0;
    fwrite(&c, sizeof(uint8_t), 1, fp);
    ++offset;
  }

  /* Store the current offset position to the file header. */
  writer->header.off_attrrefs = offset;
  fseek(fp, size, SEEK_CUR);

  /* Fill members in the feature reference header. */
  strncpy(href->chunk, CHUNK_ATTRREF, 4);
  href->size = 0;
  href->num = num_attrs;

  writer->href = href;
  writer->state = WSTATE_ATTRREFS;
  return 0;
}

int crf1dmw_close_attrrefs(crf1dmw_t* writer)
{
  uint32_t i;
  FILE *fp = writer->fp;
  featureref_header_t* href = writer->href;
  uint32_t begin = writer->header.off_attrrefs, end = 0;

  /* Make sure that we are writing attribute feature references. */
  if (writer->state != WSTATE_ATTRREFS) {
    return CRFSUITEERR_INTERNAL_LOGIC;
  }

  /* Store the current offset position. */
  end = (uint32_t)ftell(fp);

  /* Compute the size of this chunk. */
  href->size = (end - begin);

  /* Write the chunk header and offset array. */
  fseek(fp, begin, SEEK_SET);
  write_uint8_array(fp, href->chunk, 4);
  write_uint32(fp, href->size);
  write_uint32(fp, href->num);
  for (i = 0;i < href->num;++i) {
    write_uint32(fp, href->offsets[i]);
  }

  /* Move the file pointer to the tail. */
  fseek(fp, end, SEEK_SET);

  /* Uninitialize. */
  free(href);
  writer->href = NULL;
  writer->state = WSTATE_NONE;
  return 0;
}

int crf1dmw_put_attrref(crf1dmw_t* writer, int aid, const feature_refs_t* ref, int *map)
{
  int i, fid;
  uint32_t n = 0;
  FILE *fp = writer->fp;
  featureref_header_t* href = writer->href;

  /* Make sure that we are writing attribute feature references. */
  if (writer->state != WSTATE_ATTRREFS) {
    return CRFSUITEERR_INTERNAL_LOGIC;
  }

  /* Store the current offset to the offset array. */
  href->offsets[aid] = ftell(fp);

  /* Count the number of references to active features. */
  for (i = 0;i < ref->num_features;++i) {
    if (0 <= map[ref->fids[i]]) ++n;
  }

  /* Write the feature reference. */
  write_uint32(fp, (uint32_t)n);
  for (i = 0;i < ref->num_features;++i) {
    fid = map[ref->fids[i]];
    if (0 <= fid) write_uint32(fp, (uint32_t)fid);
  }

  return 0;
}

int crf1dmw_open_features(crf1dmw_t* writer)
{
  FILE *fp = writer->fp;
  feature_header_t* hfeat = NULL;

  /* Check if we aren't writing anything at this moment. */
  if (writer->state != WSTATE_NONE) {
    return CRFSUITEERR_INTERNAL_LOGIC;
  }

  /* Allocate a feature chunk header. */
  hfeat = (feature_header_t*) calloc(sizeof(feature_header_t), 1);
  if (hfeat == NULL) {
    return CRFSUITEERR_OUTOFMEMORY;
  }

  writer->header.off_features = (uint32_t)ftell(fp);
  fseek(fp, CHUNK_SIZE, SEEK_CUR);

  strncpy(hfeat->chunk, CHUNK_FEATURE, 4);
  writer->hfeat = hfeat;

  writer->state = WSTATE_FEATURES;
  return 0;
}

int crf1dmw_close_features(crf1dmw_t* writer)
{
  FILE *fp = writer->fp;
  feature_header_t* hfeat = writer->hfeat;
  uint32_t begin = writer->header.off_features, end = 0;

  /* Make sure that we are writing attribute feature references. */
  if (writer->state != WSTATE_FEATURES) {
    return CRFSUITEERR_INTERNAL_LOGIC;
  }

  /* Store the current offset position. */
  end = (uint32_t)ftell(fp);

  /* Compute the size of this chunk. */
  hfeat->size = (end - begin);

  /* Write the chunk header and offset array. */
  fseek(fp, begin, SEEK_SET);
  write_uint8_array(fp, hfeat->chunk, 4);
  write_uint32(fp, hfeat->size);
  write_uint32(fp, hfeat->num);

  /* Move the file pointer to the tail. */
  fseek(fp, end, SEEK_SET);

  /* Uninitialize. */
  free(hfeat);
  writer->hfeat = NULL;
  writer->state = WSTATE_NONE;
  return 0;
}

int crf1dmw_put_feature(crf1dmw_t* writer, int fid, const crf1dm_feature_t* f)
{
  FILE *fp = writer->fp;
  feature_header_t* hfeat = writer->hfeat;

  /* Make sure that we are writing attribute feature references. */
  if (writer->state != WSTATE_FEATURES) {
    return CRFSUITEERR_INTERNAL_LOGIC;
  }

  /* We must put features #0, #1, ..., #(K-1) in this order. */
  if (fid != hfeat->num) {
    return CRFSUITEERR_INTERNAL_LOGIC;
  }

  write_uint32(fp, f->type);
  write_uint32(fp, f->src);
  write_uint32(fp, f->dst);
  write_float(fp, f->weight);
  ++hfeat->num;
  return 0;
}

int crf1dmw_open_sm(crf1dmw_t *writer, const crf1de_semimarkov_t *a_sm)
{
  /* Check that we are not writing anything at this moment. */
  if (writer->state != WSTATE_NONE) {
    return CRFSUITEERR_INTERNAL_LOGIC;
  }
  size_t size = HSM_CHUNK_SIZE + sizeof(uint32_t) * a_sm->m_num_frw;

  sm_header_t* hsm = (sm_header_t *) calloc(1, size);
  if (hsm == NULL) {
    return CRFSUITEERR_OUTOFMEMORY;
  }

  /* Align the offset to a DWORD boundary. */
  FILE *fp = writer->fp;
  uint8_t c = 0;
  uint32_t offset = (uint32_t) ftell(fp);
  while (offset % 4 != 0) {
    fwrite(&c, sizeof(uint8_t), 1, fp);
    ++offset;
  }

  /* Store current offset position in file header. */
  writer->header.off_sm = offset;
  /* Skip bytes that are needed for sm header. */
  if (fseek(fp, size, SEEK_CUR) != 0) {
    goto error_exit;
  }

  /* Fill members in the header that we already know. */
  strncpy(hsm->chunk, CHUNK_HSMM, CHUNK_HSMM_SIZE);
  hsm->num_labels = 0;
  hsm->num_suffixes = 0;
  hsm->num_states = 0;
  hsm->num_bkw_states = (uint32_t) a_sm->m_num_bkw;
  hsm->max_order = (uint32_t) a_sm->m_max_order;

  /* Store reference to the header. */
  writer->hsm = hsm;
  writer->state = WSTATE_SM;

  /* Write maximum segment length for each label. */
  hsm->off_max_seg_len = (uint32_t) ftell(fp);
  uint32_t max_seg_len = 0;
  for (size_t i = 0; i < a_sm->L; ++i) {
    max_seg_len = (uint32_t) a_sm->m_max_seg_len[i];
    fwrite(&max_seg_len, sizeof(uint32_t), 1, fp);
    ++hsm->num_labels;
  }

  /* Write array of sufixes. */
  int ptrn_id = 0;
  uint32_t feat_id = 0;
  hsm->off_suffixes = (uint32_t) ftell(fp);
  for (size_t i = 0; i < a_sm->m_num_suffixes; ++i) {
    ptrn_id = a_sm->m_suffixes[i];
    if (ptrn_id < 0 || a_sm->m_ptrns[ptrn_id].m_len < 2)
      feat_id = (uint32_t) -1;
    else
      feat_id = (uint32_t) a_sm->m_ptrns[ptrn_id].m_feat_id;
    /* skip suffixes with length lesser than two */
    fwrite(&feat_id, sizeof(uint32_t), 1, fp);
    ++hsm->num_suffixes;
  }

  return 0;

 error_exit:
  free(hsm);
  return 1;
}

int crf1dmw_close_sm(crf1dmw_t* writer)
{
  int ret = 0;
  FILE *fp = writer->fp;
  sm_header_t *hsm = writer->hsm;
  uint32_t begin = writer->header.off_sm;
  uint32_t end = (uint32_t)ftell(fp);

  /* Make sure that we are writing the semi-markov model. */
  if (writer->state != WSTATE_SM) {
    ret = CRFSUITEERR_INTERNAL_LOGIC;
    goto error_exit;
  }

  /* Write the chunk header and offset array. */
  if (fseek(fp, begin, SEEK_SET) == -1) {
    ret = 1;
    goto error_exit;
  }
  write_uint8_array(fp, hsm->chunk, 4);
  write_uint32(fp, hsm->max_order);
  write_uint32(fp, hsm->num_labels);
  write_uint32(fp, hsm->num_states);
  write_uint32(fp, hsm->num_bkw_states);
  write_uint32(fp, hsm->num_suffixes);
  write_uint32(fp, hsm->off_max_seg_len);
  write_uint32(fp, hsm->off_suffixes);
  for (uint32_t i = 0; i < hsm->num_states; ++i) {
    fprintf(stderr, "crf1dmw_close_sm: hsm->off_states[%u] = %u\n", i, hsm->off_states[i]);
    write_uint32(fp, hsm->off_states[i]);
  }

  if (ferror(fp)) {
    ret = 1;
    goto error_exit;
  }
 /* Move the file pointer to the end of file. */
  fseek(fp, end, SEEK_SET);

  /* Uninitialize. */
 error_exit:
  free(hsm);
  writer->hsm = NULL;
  writer->state = WSTATE_NONE;
  return ret;
}

int crf1dmw_put_sm_state(crf1dmw_t* writer, int sid, \
			 const crf1de_state_t *state, \
			 const crf1de_semimarkov_t* sm)
{
  FILE *fp = writer->fp;
  sm_header_t *hsm = writer->hsm;

  /* Make sure that we are writing semi-markov states. */
  if (writer->state != WSTATE_SM) {
    return CRFSUITEERR_INTERNAL_LOGIC;
  }

  /* We must put states #0, #1, ..., #(K-1) in order, rebuke if this is not
     the case. */
  if ((uint32_t) sid != hsm->num_states) {
    return CRFSUITEERR_INTERNAL_LOGIC;
  }

  /* Store the current offset to the offset array. */
  hsm->off_states[sid] = (uint32_t) ftell(fp);
  fprintf(stderr, "crf1dmw_put_sm_state: hsm->off_states[%d] = %u\n", sid, hsm->off_states[sid]);

  /* Write information about state. */
  /* id of corresponding feature */
  fprintf(stderr, "crf1dmw_put_sm_state: state->m_feat_id = %d\n", state->m_feat_id);
  /* write_uint32(fp, (uint32_t) state->m_feat_id); */
  /* length of underlying label sequence */
  fprintf(stderr, "crf1dmw_put_sm_state: state->m_len = %d\n", state->m_len);
  write_uint32(fp, (uint32_t) state->m_len);
  /* write the label sequence */
  for (size_t i = 0; i < state->m_len; ++i) {
    write_uint32(fp, (uint32_t) state->m_seq[i]);
  }
  /* write the number of affixes */
  write_uint32(fp, (uint32_t) state->m_num_affixes);
  /* write prefixes (as id's of corresponding states) */
  uint32_t affx_id = 0;
  for (size_t i = 0; i < state->m_num_affixes; ++i) {
    affx_id = (uint32_t) state->m_frw_trans1[i];
    write_uint32(fp, affx_id);
  }

  /* write the offsets in the corresponding array of suffixes */
  for (size_t i = 0; i < state->m_num_affixes; ++i) {
    affx_id = (uint32_t) state->m_frw_trans2[i];
    write_uint32(fp, affx_id);
  }
  ++hsm->num_states;
  return 0;
}

static crf1de_semimarkov_t *crf1dm_get_sm(void *buffer, void *sm_buffer, size_t size)
{
  /* The minimum size of a valid semi-markov model is the size of its
     header */
  if (size < sizeof(sm_header_t)) {
    return NULL;
  }
  /* Check chunk id */
  uint8_t *model_buffer = (uint8_t *) buffer;
  uint8_t *p = (uint8_t *) sm_buffer;
  if (memcmp(CHUNK_HSMM, p, CHUNK_HSMM_SIZE) != 0)
    return NULL;

  p += CHUNK_HSMM_SIZE;

  fprintf(stderr, "crf1dm_get_sm: CHUNK_HSMM checked\n");
  /* Populate semi-markov header */
  sm_header_t hsm;
  p += read_uint32(p, &hsm.max_order);
  p += read_uint32(p, &hsm.num_labels);
  p += read_uint32(p, &hsm.num_states);
  p += read_uint32(p, &hsm.num_bkw_states);
  p += read_uint32(p, &hsm.num_suffixes);
  p += read_uint32(p, &hsm.off_max_seg_len);
  p += read_uint32(p, &hsm.off_suffixes);
  fprintf(stderr, "crf1dm_get_sm: hsm read\n");

  if (hsm.max_order >= CRFSUITE_SM_MAX_PTRN_LEN) {
    fprintf(stderr, "Maximum order of stored states (%zu) exceeds maximum allowed length.\n\
 To increase this limit, consider recompiling the program with the option\n\
-DCRFSUITE_SM_MAX_PTRN_LEN=NEW_LIM added to CPPFLAGS.\n", hsm.max_order);
    return NULL;
  }

  /* Create and initially populate an instance of semi-markov model */
  crf1de_semimarkov_t *sm = crf1de_create_semimarkov();
  if (!sm)
    return sm;

  fprintf(stderr, "crf1dm_get_sm: sm created\n");
  /* allocate memory for members of semi-markov model */
  sm->m_max_seg_len = (int *) calloc(hsm.num_labels, sizeof(int));
  sm->m_suffixes = (int *) calloc(hsm.num_suffixes, sizeof(int));
  sm->m_frw_states = (crf1de_state_t *) calloc(hsm.num_states, sizeof(crf1de_state_t));
  sm->m_frw_llabels = (int *) calloc(hsm.num_states, sizeof(int));
  sm->m_frw_trans1 = (int *) calloc(hsm.num_bkw_states, sizeof(int));
  sm->m_frw_trans2 = (int *) calloc(hsm.num_bkw_states, sizeof(int));
  fprintf(stderr, "crf1dm_get_sm: sm members allocated\n");

  if (!sm->m_max_seg_len || !sm->m_suffixes || !sm->m_frw_states || \
      !sm->m_frw_trans1 || !sm->m_frw_trans2 || !sm->m_frw_llabels)
    goto error_exit;

  /* populate forward states of semi-markov model */
  uint32_t val;
  uint8_t *saved_state = NULL;
  size_t i = 0, trans1_c = 0, trans2_c = 0;
  crf1de_state_t *sm_state = NULL;
  fprintf(stderr, "crf1dm_get_sm: populating sm states\n");
  for (sm->m_num_frw = 0; sm->m_num_frw < hsm.num_states; ++sm->m_num_frw) {
    /* obtain addresses of semi-markov and saved state */
    sm_state = &sm->m_frw_states[sm->m_num_frw];
    fprintf(stderr, "crf1dm_get_sm: sm_state = %p\n", sm_state);
    p += read_uint32(p, &val);
    fprintf(stderr, "crf1dm_get_sm: val = %u\n", val);
    saved_state = model_buffer + val;
    fprintf(stderr, "crf1dm_get_sm: saved_state = %p\n", saved_state);
    /* obtain feature id and value */
    /* saved_state += read_uint32(saved_state, &val); */
    /* sm_state->m_feat_id = (int) val; */
    /* fprintf(stderr, "crf1dm_get_sm: sm_state->m_feat_id = %d\n", sm_state->m_feat_id); */
    saved_state += read_uint32(saved_state, &val);
    sm_state->m_len = (int) val;
    fprintf(stderr, "crf1dm_get_sm: sm_state->m_len = %d\n", val);
    /* populate label sequence */
    for (i = 0; i < sm_state->m_len; ++i) {
      saved_state += read_uint32(saved_state, &val);
      sm_state->m_seq[i] = (int) val;
      fprintf(stderr, "crf1dm_get_sm: sm_state->m_seq[%d] = %d\n", i, sm_state->m_seq[i]);
    }
    /* populate prefixes */
    saved_state += read_uint32(saved_state, &val);
    sm_state->m_num_affixes = (int) val;
    fprintf(stderr, "crf1dm_get_sm: sm_state->m_num_affixes = %d\n", sm_state->m_num_affixes);

    sm_state->m_frw_trans1 = &sm->m_frw_trans1[trans1_c];
    fprintf(stderr, "crf1dm_get_sm: sm_state->m_frw_trans1 = %p\n", sm_state->m_frw_trans1);
    for (i = 0; i < sm_state->m_num_affixes; ++i) {
      saved_state += read_uint32(saved_state, &val);
      sm->m_frw_trans1[trans1_c++] = (int) val;
      fprintf(stderr, "crf1dm_get_sm: sm->m_frw_trans1[%d] = %d\n", trans1_c - 1, (int) val);
    }
    /* populate indices of suffixes */
    sm_state->m_frw_trans2 = &sm->m_frw_trans2[trans2_c];
    fprintf(stderr, "crf1dm_get_sm: sm_state->m_frw_trans2 = %p\n", sm_state->m_frw_trans2);
    for (i = 0; i < sm_state->m_num_affixes; ++i) {
      saved_state += read_uint32(saved_state, &val);
      sm->m_frw_trans2[trans2_c++] = (int) val;
      fprintf(stderr, "crf1dm_get_sm: sm->m_frw_trans2[%d] = %d\n", trans2_c - 1, (int) val);
    }
  }
  fprintf(stderr, "crf1dm_get_sm: sm states populated\n");

  if (trans1_c != hsm.num_bkw_states || trans2_c != hsm.num_bkw_states) {
    fprintf(stderr, "Unmatching number of read prefixes and suffixes: \
num_bkw = %zu, prefixes = %zu, suffixes = %zu.\n", hsm.num_bkw_states, trans1_c, \
	    trans2_c);
    goto error_exit;
  }

  /* read and populate maximum segment lengths */
  p = model_buffer + hsm.off_max_seg_len;
  for (sm->L = 0; sm->L < hsm.num_labels; ++sm->L) {
    p += read_uint32(p, &val);
    sm->m_max_seg_len[sm->L] = (int) val;
  }

  /* read and populate suffixes */
  p = model_buffer + hsm.off_suffixes;
  for (sm->m_num_suffixes = 0; sm->m_num_suffixes < hsm.num_suffixes; \
       ++sm->m_num_suffixes) {
    p += read_uint32(p, &val);
    sm->m_suffixes[sm->m_num_suffixes] = (int) val;
  }

  /* return constructed model */
  return sm;

 error_exit:
  sm->clear(sm);
  free(sm);
  return NULL;
}

crf1dm_t* crf1dm_new(const char *filename, const int ftype)
{
  FILE *fp = NULL;
  uint8_t* p = NULL;
  crf1dm_t *model = NULL;
  header_t *header = NULL;

  model = (crf1dm_t*)calloc(1, sizeof(crf1dm_t));
  if (model == NULL) {
    goto error_exit;
  }

  fp = fopen(filename, "rb");
  if (fp == NULL) {
    goto error_exit;
  }

  fseek(fp, 0, SEEK_END);
  model->size = (uint32_t)ftell(fp);
  fseek(fp, 0, SEEK_SET);

  model->buffer = (uint8_t*)malloc(model->size + 16);
  model->buffer_orig = model->buffer;
  while ((uintptr_t)model->buffer % 16 != 0) {
    ++model->buffer;
  }

  if (fread(model->buffer, 1, model->size, fp) != model->size) {
    free(model->buffer_orig);
    goto error_exit;
  }
  fclose(fp);
  fp = NULL;

  /* Write the file header. */
  header = (header_t*)calloc(1, sizeof(header_t));

  p = model->buffer;
  p += read_uint8_array(p, header->magic, sizeof(header->magic));
  p += read_uint32(p, &header->size);
  p += read_uint8_array(p, header->type, sizeof(header->type));
  if ((ftype == FTYPE_CRF1TREE &&					\
       strncmp(header->type, MODELTYPE_TREE, sizeof(header->type)) != 0) || \
      (ftype == FTYPE_SEMIMCRF &&					\
       strncmp(header->type, MODELTYPE_SEMIM, sizeof(header->type)) != 0) || \
      (ftype == FTYPE_CRF1D &&					\
       strncmp(header->type, MODELTYPE_CRF1D, sizeof(header->type)) != 0)) {
    fprintf(stderr, "ERROR: Stored model has another type than the model you are loading.\n");
    free(model->buffer_orig);
    goto error_exit;
  }
  p += read_uint32(p, &header->version);
  p += read_uint32(p, &header->num_features);
  p += read_uint32(p, &header->num_labels);
  p += read_uint32(p, &header->num_attrs);
  p += read_uint32(p, &header->off_features);
  p += read_uint32(p, &header->off_labels);
  p += read_uint32(p, &header->off_attrs);
  p += read_uint32(p, &header->off_labelrefs);
  p += read_uint32(p, &header->off_attrrefs);
  p += read_uint32(p, &header->off_sm);
  model->header = header;

  model->labels = cqdb_reader(
			      model->buffer + header->off_labels,
			      model->size - header->off_labels
			      );

  model->attrs = cqdb_reader(
			     model->buffer + header->off_attrs,
			     model->size - header->off_attrs
			     );
  if (header->off_sm) {
    fprintf(stderr, "crf1dm_new: calling crf1dm_get_sm()\n");
    model->sm = crf1dm_get_sm(model->buffer, model->buffer + header->off_sm, \
			      model->size - header->off_sm);
    fprintf(stderr, "crf1dm_new: crf1dm_get_sm() finished\n");
  } else {
    model->sm = NULL;
  }
  return model;

 error_exit:
  if (model != NULL) {
    free(model);
  }
  if (fp != NULL) {
    fclose(fp);
  }
  return NULL;
}

void crf1dm_close(crf1dm_t* model)
{
  if (model->labels != NULL) {
    cqdb_delete(model->labels);
  }
  if (model->attrs != NULL) {
    cqdb_delete(model->attrs);
  }
  if (model->header != NULL) {
    free(model->header);
    model->header = NULL;
  }
  if (model->buffer_orig != NULL) {
    free(model->buffer_orig);
    model->buffer_orig = model->buffer = NULL;
  }
  free(model);
}

int crf1dm_get_num_attrs(crf1dm_t* model)
{
  return model->header->num_attrs;
}

int crf1dm_get_num_labels(crf1dm_t* model)
{
  return model->header->num_labels;
}

const char *crf1dm_to_label(crf1dm_t* model, int lid)
{
  if (model->labels != NULL) {
    return cqdb_to_string(model->labels, lid);
  } else {
    return NULL;
  }
}

int crf1dm_to_lid(crf1dm_t* model, const char *value)
{
  if (model->labels != NULL) {
    return cqdb_to_id(model->labels, value);
  } else {
    return -1;
  }
}

int crf1dm_to_aid(crf1dm_t* model, const char *value)
{
  if (model->attrs != NULL) {
    return cqdb_to_id(model->attrs, value);
  } else {
    return -1;
  }
}

const char *crf1dm_to_attr(crf1dm_t* model, int aid)
{
  if (model->attrs != NULL) {
    return cqdb_to_string(model->attrs, aid);
  } else {
    return NULL;
  }
}

int crf1dm_get_labelref(crf1dm_t* model, int lid, feature_refs_t* ref)
{
  uint8_t *p = model->buffer;
  uint32_t offset;

  p += model->header->off_labelrefs;
  p += CHUNK_SIZE;
  p += sizeof(uint32_t) * lid;
  read_uint32(p, &offset);

  p = model->buffer + offset;
  p += read_uint32(p, &ref->num_features);
  ref->fids = (int*)p;
  return 0;
}

int crf1dm_get_attrref(crf1dm_t* model, int aid, feature_refs_t* ref)
{
  uint8_t *p = model->buffer;
  uint32_t offset;

  p += model->header->off_attrrefs;
  p += CHUNK_SIZE;
  p += sizeof(uint32_t) * aid;
  read_uint32(p, &offset);

  p = model->buffer + offset;
  p += read_uint32(p, &ref->num_features);
  ref->fids = (int*)p;
  return 0;
}

int crf1dm_get_featureid(feature_refs_t* ref, int i)
{
  uint32_t fid;
  uint8_t* p = (uint8_t*)ref->fids;
  p += sizeof(uint32_t) * i;
  read_uint32(p, &fid);
  return (int)fid;
}

int crf1dm_get_feature(crf1dm_t* model, int fid, crf1dm_feature_t* f)
{
  uint8_t *p = NULL;
  uint32_t val = 0;
  uint32_t offset = model->header->off_features + CHUNK_SIZE;
  offset += FEATURE_SIZE * fid;
  p = model->buffer + offset;
  p += read_uint32(p, &val);
  f->type = val;
  p += read_uint32(p, &val);
  f->src = val;
  p += read_uint32(p, &val);
  f->dst = val;
  p += read_float(p, &f->weight);
  return 0;
}

void crf1dm_dump(crf1dm_t* crf1dm, FILE *fp)
{
  int j;
  uint32_t i;
  feature_refs_t refs;
  const header_t* hfile = crf1dm->header;

  /* Dump the file header. */
  fprintf(fp, "FILEHEADER = {\n");
  fprintf(fp, "  magic: %c%c%c%c\n",
	  hfile->magic[0], hfile->magic[1], hfile->magic[2], hfile->magic[3]);
  fprintf(fp, "  size: %d\n", hfile->size);
  fprintf(fp, "  type: %c%c%c%c\n",
	  hfile->type[0], hfile->type[1], hfile->type[2], hfile->type[3]);
  fprintf(fp, "  version: %d\n", hfile->version);
  fprintf(fp, "  num_features: %d\n", hfile->num_features);
  fprintf(fp, "  num_labels: %d\n", hfile->num_labels);
  fprintf(fp, "  num_attrs: %d\n", hfile->num_attrs);
  fprintf(fp, "  off_features: 0x%X\n", hfile->off_features);
  fprintf(fp, "  off_labels: 0x%X\n", hfile->off_labels);
  fprintf(fp, "  off_attrs: 0x%X\n", hfile->off_attrs);
  fprintf(fp, "  off_labelrefs: 0x%X\n", hfile->off_labelrefs);
  fprintf(fp, "  off_attrrefs: 0x%X\n", hfile->off_attrrefs);
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  /* Dump the labels. */
  fprintf(fp, "LABELS = {\n");
  for (i = 0;i < hfile->num_labels;++i) {
    const char *str = crf1dm_to_label(crf1dm, i);
#if 0
    int check = crf1dm_to_lid(crf1dm, str);
    if (i != check) {
      fprintf(fp, "WARNING: inconsistent label CQDB\n");
    }
#endif
    fprintf(fp, "  %5d: %s\n", i, str);
  }
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  /* Dump the attributes. */
  fprintf(fp, "ATTRIBUTES = {\n");
  for (i = 0;i < hfile->num_attrs;++i) {
    const char *str = crf1dm_to_attr(crf1dm, i);
#if 0
    int check = crf1dm_to_aid(crf1dm, str);
    if (i != check) {
      fprintf(fp, "WARNING: inconsistent attribute CQDB\n");
    }
#endif
    fprintf(fp, "  %5d: %s\n", i, str);
  }
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  /* Dump the transition features. */
  fprintf(fp, "TRANSITIONS = {\n");
  for (i = 0;i < hfile->num_labels;++i) {
    crf1dm_get_labelref(crf1dm, i, &refs);
    for (j = 0;j < refs.num_features;++j) {
      crf1dm_feature_t f;
      int fid = crf1dm_get_featureid(&refs, j);
      const char *from = NULL, *to = NULL;

      crf1dm_get_feature(crf1dm, fid, &f);
      from = crf1dm_to_label(crf1dm, f.src);
      to = crf1dm_to_label(crf1dm, f.dst);
      fprintf(fp, "  (%d) %s --> %s: %f\n", f.type, from, to, f.weight);
    }
  }
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  /* Dump the transition features. */
  fprintf(fp, "STATE_FEATURES = {\n");
  for (i = 0;i < hfile->num_attrs;++i) {
    crf1dm_get_attrref(crf1dm, i, &refs);
    for (j = 0;j < refs.num_features;++j) {
      crf1dm_feature_t f;
      int fid = crf1dm_get_featureid(&refs, j);
      const char *attr = NULL, *to = NULL;

      crf1dm_get_feature(crf1dm, fid, &f);
#if 0
      if (f.src != i) {
	fprintf(fp, "WARNING: an inconsistent attribute reference.\n");
      }
#endif
      attr = crf1dm_to_attr(crf1dm, f.src);
      to = crf1dm_to_label(crf1dm, f.dst);
      fprintf(fp, "  (%d) %s --> %s: %f\n", f.type, attr, to, f.weight);
    }
  }
  fprintf(fp, "}\n");
  fprintf(fp, "\n");
}
