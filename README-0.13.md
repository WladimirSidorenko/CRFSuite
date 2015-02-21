CRFSuite 0.13
=============

This version of CRFSuite has been extended with the following variants of CRFs:

* tree-structured CRFs;
* semi-Markov CRFs of arbitrary orders;
* linear-chain CRFs  of arbitrary orders.

To invoke tree-structured CRFs, you should provide the option `--type=tree` when
running `crfsuite learn` and also specify this option when you later envoke
`crfsuite tag` with the trained model.

To use higher-order linear-chain and semi-markov CRFs, you should specify the
option `--type=semim` both during the training and during the tagging, e.g.:

`crfsuite learn --type=semim -p feature.max_seg_len=-1 -m semim.model tests/test_sm_1.input`

`crfsuite tag --type=semim -m semim.model tests/test_sm_1.input`

Setting the option `-p feature.max_seg_len=` to a negative value will envoke the
semi-Markov variant of CRF, providing a non-negative integer will activate the
linear-chain model.  To specify the maximum order of transition features, use
the option `-p feature.max_order=`, e.g.:

`crfsuite learn --type=semim -p feature.max_seg_len=1 -p feature.max_order=4 \
-m semim.model tests/test_sm_1.input`

will train a 4-th order linear-chain model.

FORMAT
======
The format for the first and higher order linear-chain and semi-Markov
looks as follows:
`label1 \t feat_name1:value1 \t feat_name2:value2 \t feat_name3:value3`
`label2 \t feat_name4:value4 \t feat_name5:value5 \t feat_name6:value6`
`label3 \t feat_name7:value7 \t feat_name8:value8`
`label4 \t feat_name9:value9 \t feat_name10:value10 \t feat_name11:value11`
``
For testing, you can either specify a valid label as the first field or
put any value (e.g `_' underscore ) which does not coincide with any known
tagset label.  This first field is skipped during the testing.  Empty lines
delimit the sequences.

For the tree-structured CRFs, you should specify the id of the node as the
second field and the id of the parent node as the third field, e.g:
`label1 \t node_id1 \t node_id2 \t feat_name1:value1`
`label2 \t node_id2 \t _ \t feat_name2:value2`
`label3 \t node_id3 \t node_id4`
`label4 \t node_id4 \t node_id2 \t feat_name3:value3 \t feat_name4:value4`
``
The parent of the root node should be specified as `_' (underscore) (see file
`tests/test_tree_2.input` for an example).

WARNINGS
========
1) Training algorithms other than `l-BFGS` seem to work and converge but have
   not been thoroughly tested yet;
2) Semi-Markov and higher-order linear-chain models do not support the options
`-i` and `-p' for tagging yet;
3) No speed optimization was done for the higher-order semi-Markov and linear-
chain models.
