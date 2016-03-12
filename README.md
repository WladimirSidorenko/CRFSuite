CRFSuite (0.13)
==============

[![Build Status](https://travis-ci.org/WladimirSidorenko/CRFSuite.svg?branch=master)](https://travis-ci.org/WladimirSidorenko/CRFSuite)

Table of Contents
-----------------

  * [CRFSuite (0.13)](#crfsuite-013)
	* [Introduction](#introduction)
	* [Version 0.13](#version-013)
	* [Format](#format)
	* [Warnings](#warnings)
	* [Copyright and Licensing](#copyright-and-licensing)
	* [Acknowledgment](#acknowledgment)

Introduction
------------

CRFSuite 0.13 is a fork of the
[Naoaki Okazaki's](http://www.chokkan.org/) implementation of the
Conditional Random Fields (CRFs) for labeling sequential data.  Please
refer to the [web site](http://www.chokkan.org/software/crfsuite/) for
more information about the original software.

Version 0.13
------------

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

`crfsuite learn --type=semim -p feature.max_seg_len=1 -p feature.max_order=4
-m semim.model tests/test_sm_1.input`

will train a 4-th order linear-chain model.

Format
------
The format for the first and higher order linear-chain and semi-Markov
looks as follows:

`label1 \t feat_name1:value1 \t feat_name2:value2 \t feat_name3:value3`

`label2 \t feat_name4:value4 \t feat_name5:value5 \t feat_name6:value6`

`label3 \t feat_name7:value7 \t feat_name8:value8`

`label4 \t feat_name9:value9 \t feat_name10:value10 \t feat_name11:value11`

` `

For testing, you can either specify a valid label as the first field or
put any value (e.g `_` underscore ) which does not coincide with any known
tagset label.  This first field is skipped during the testing.  Empty lines
delimit the sequences.

For the tree-structured CRFs, you should specify the id of the node as the
second field and the id of the parent node as the third field, e.g:

`label1 \t node_id1 \t node_id2 \t feat_name1:value1`

`label2 \t node_id2 \t _ \t feat_name2:value2`

`label3 \t node_id3 \t node_id4`

`label4 \t node_id4 \t node_id2 \t feat_name3:value3 \t feat_name4:value4`

` `

The parent of the root node should be specified as `_` (underscore) (see file
`tests/test_tree_2.input` for an example).

Warnings
--------
1. Only `l-BFGS` is supported so far for the higher-order and
   semi-markov models;
2. Semi-Markov and higher-order linear-chain models do not support the options
`-i` and `-p` for tagging yet;
3. No speed optimization was done for the higher-order semi-Markov and linear-
chain models;
4. C++ interface has not been updated to support the new types.

Copyright and Licensing
-----------------------

This program is distributed under the modified BSD license. Refer to
COPYING file for the precise description of the license.


Portions of this software are based on libLBFGS.

The MIT License

Copyright (c) 1990 Jorge Nocedal
Copyright (c) 2007 Naoaki Okazaki

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


Portions of this software are based on Constant Quark Database (CQDB).

The BSD license.

Copyright (c) 2007, Naoaki Okazaki
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Northwestern University, University of Tokyo,
      nor the names of its contributors may be used to endorse or promote
      products derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Portions of this software are based on RumAVL.

MIT/X Consortium License.

Copyright (c) 2005-2007 Jesse Long <jpl@unknown.za.net>
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

   1. The above copyright notice and this permission notice shall be
      included in all copies or substantial portions of the Software.
   2. The origin of the Software must not be misrepresented; you must not
      claim that you wrote the original Software.
   3. Altered source versions of the Software must be plainly marked as
      such, and must not be misrepresented as being the original Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.


Portions of this software are based on a portable stdint.h (for MSVC).

Copyright (c) 2005-2007 Paul Hsieh

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

    Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

    Redistributions in binary form must not misrepresent the orignal
    source in the documentation and/or other materials provided
    with the distribution.

    The names of the authors nor its contributors may be used to
    endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.


Portions of this software are based on Mersenne Twister.

Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  1. Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.

  3. The names of its contributors may not be used to endorse or
     promote products derived from this software without specific
     prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Acknowledgment
--------------

Special thanks goes to:

* Olivier Grisel
* Andreas Holzbach
* Baoli Li
* Yoshimasa Tsuruoka
* Hiroshi Manabe
* Riza Theresa B. Batista-Navarro
