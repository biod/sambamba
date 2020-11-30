module bio.std.hts.sam.utils.recordparser;

#line 1 "sam_alignment.rl"
/*
    This file is part of BioD.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

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
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

*/

#line 28 "sam_alignment.d"
static byte[] _sam_alignment_actions = [
	0, 1, 0, 1, 2, 1, 4, 1, 
	6, 1, 7, 1, 8, 1, 9, 1, 
	10, 1, 11, 1, 12, 1, 13, 1, 
	14, 1, 15, 1, 16, 1, 17, 1, 
	18, 1, 19, 1, 21, 1, 22, 1, 
	23, 1, 27, 1, 28, 1, 29, 1, 
	30, 1, 31, 1, 32, 1, 33, 1, 
	34, 1, 35, 1, 36, 1, 37, 1, 
	39, 1, 40, 1, 41, 1, 42, 1, 
	43, 1, 44, 1, 45, 1, 46, 1, 
	48, 1, 49, 1, 51, 1, 53, 1, 
	57, 1, 60, 1, 61, 1, 62, 1, 
	63, 1, 64, 2, 1, 2, 2, 3, 
	38, 2, 3, 58, 2, 5, 59, 2, 
	6, 7, 2, 20, 23, 2, 24, 25, 
	2, 26, 29, 2, 47, 50, 2, 55, 
	64, 2, 56, 64, 3, 3, 52, 64, 
	3, 3, 58, 64, 3, 5, 54, 64, 
	3, 5, 59, 64, 3, 26, 1, 2
	
];

static short[] _sam_alignment_key_offsets = [
	0, 0, 5, 7, 10, 15, 18, 20, 
	23, 25, 28, 31, 32, 36, 39, 41, 
	44, 48, 50, 53, 60, 61, 63, 67, 
	73, 74, 80, 81, 83, 84, 91, 92, 
	96, 98, 99, 106, 110, 112, 116, 118, 
	119, 120, 121, 122, 123, 129, 130, 132, 
	133, 140, 144, 146, 150, 152, 153, 154, 
	155, 156, 157, 161, 163, 170, 173, 176, 
	179, 182, 185, 188, 191, 194, 197, 200, 
	203, 206, 209, 212, 215, 218, 219, 222, 
	225, 228, 231, 234, 237, 240, 243, 246, 
	249, 252, 255, 258, 261, 264, 267, 268, 
	269, 270, 281, 292, 303, 314, 325, 336, 
	347, 358, 369, 380, 391, 402, 413, 424, 
	435, 446, 457, 466, 469, 472, 475, 478, 
	481, 484, 487, 490, 493, 496, 499, 502, 
	505, 508, 511, 514, 517, 518, 521, 524, 
	527, 530, 533, 536, 539, 542, 545, 548, 
	551, 554, 557, 560, 563, 566, 567, 568, 
	571, 574, 577, 580, 583, 586, 589, 592, 
	595, 598, 601, 604, 607, 610, 613, 616, 
	617, 622, 623, 624, 625, 626, 627, 628, 
	629, 630, 631, 632, 633, 634, 635, 636, 
	637, 638, 639, 640, 641, 642, 643, 644, 
	647, 648, 652, 656, 660, 664, 668, 672, 
	676, 680, 684, 688, 692, 696, 700, 704, 
	708, 712, 716, 718, 724, 728, 735, 737, 
	744, 747, 752, 755, 761, 762, 765, 768, 
	771, 774, 777, 780, 783, 786, 789, 792, 
	795, 798, 801, 804, 807, 810, 813, 814, 
	814, 814, 814, 814, 814, 814, 814, 814, 
	814, 814, 814, 814
];

static char[] _sam_alignment_trans_keys = [
	9u, 33u, 63u, 65u, 126u, 48u, 57u, 9u, 
	48u, 57u, 42u, 33u, 60u, 62u, 126u, 9u, 
	33u, 126u, 48u, 57u, 9u, 48u, 57u, 48u, 
	57u, 9u, 48u, 57u, 42u, 48u, 57u, 9u, 
	42u, 61u, 33u, 126u, 9u, 33u, 126u, 48u, 
	57u, 9u, 48u, 57u, 43u, 45u, 48u, 57u, 
	48u, 57u, 9u, 48u, 57u, 42u, 46u, 61u, 
	65u, 90u, 97u, 122u, 9u, 33u, 126u, 65u, 
	90u, 97u, 122u, 48u, 57u, 65u, 90u, 97u, 
	122u, 58u, 65u, 66u, 72u, 90u, 102u, 105u, 
	58u, 33u, 126u, 58u, 67u, 73u, 83u, 99u, 
	102u, 105u, 115u, 44u, 43u, 45u, 48u, 57u, 
	48u, 57u, 44u, 43u, 45u, 46u, 105u, 110u, 
	48u, 57u, 46u, 105u, 48u, 57u, 48u, 57u, 
	43u, 45u, 48u, 57u, 48u, 57u, 110u, 102u, 
	97u, 110u, 58u, 48u, 57u, 65u, 70u, 97u, 
	102u, 58u, 32u, 126u, 58u, 43u, 45u, 46u, 
	105u, 110u, 48u, 57u, 46u, 105u, 48u, 57u, 
	48u, 57u, 43u, 45u, 48u, 57u, 48u, 57u, 
	110u, 102u, 97u, 110u, 58u, 43u, 45u, 48u, 
	57u, 48u, 57u, 9u, 46u, 61u, 65u, 90u, 
	97u, 122u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 9u, 9u, 61u, 68u, 
	80u, 83u, 88u, 48u, 57u, 72u, 73u, 77u, 
	78u, 61u, 68u, 80u, 83u, 88u, 48u, 57u, 
	72u, 73u, 77u, 78u, 61u, 68u, 80u, 83u, 
	88u, 48u, 57u, 72u, 73u, 77u, 78u, 61u, 
	68u, 80u, 83u, 88u, 48u, 57u, 72u, 73u, 
	77u, 78u, 61u, 68u, 80u, 83u, 88u, 48u, 
	57u, 72u, 73u, 77u, 78u, 61u, 68u, 80u, 
	83u, 88u, 48u, 57u, 72u, 73u, 77u, 78u, 
	61u, 68u, 80u, 83u, 88u, 48u, 57u, 72u, 
	73u, 77u, 78u, 61u, 68u, 80u, 83u, 88u, 
	48u, 57u, 72u, 73u, 77u, 78u, 61u, 68u, 
	80u, 83u, 88u, 48u, 57u, 72u, 73u, 77u, 
	78u, 61u, 68u, 80u, 83u, 88u, 48u, 57u, 
	72u, 73u, 77u, 78u, 61u, 68u, 80u, 83u, 
	88u, 48u, 57u, 72u, 73u, 77u, 78u, 61u, 
	68u, 80u, 83u, 88u, 48u, 57u, 72u, 73u, 
	77u, 78u, 61u, 68u, 80u, 83u, 88u, 48u, 
	57u, 72u, 73u, 77u, 78u, 61u, 68u, 80u, 
	83u, 88u, 48u, 57u, 72u, 73u, 77u, 78u, 
	61u, 68u, 80u, 83u, 88u, 48u, 57u, 72u, 
	73u, 77u, 78u, 61u, 68u, 80u, 83u, 88u, 
	48u, 57u, 72u, 73u, 77u, 78u, 61u, 68u, 
	80u, 83u, 88u, 48u, 57u, 72u, 73u, 77u, 
	78u, 61u, 68u, 80u, 83u, 88u, 72u, 73u, 
	77u, 78u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 9u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 9u, 33u, 63u, 65u, 126u, 9u, 9u, 
	9u, 9u, 9u, 9u, 9u, 9u, 9u, 9u, 
	9u, 9u, 9u, 9u, 9u, 9u, 9u, 9u, 
	9u, 9u, 9u, 9u, 9u, 33u, 126u, 9u, 
	9u, 44u, 48u, 57u, 9u, 44u, 48u, 57u, 
	9u, 44u, 48u, 57u, 9u, 44u, 48u, 57u, 
	9u, 44u, 48u, 57u, 9u, 44u, 48u, 57u, 
	9u, 44u, 48u, 57u, 9u, 44u, 48u, 57u, 
	9u, 44u, 48u, 57u, 9u, 44u, 48u, 57u, 
	9u, 44u, 48u, 57u, 9u, 44u, 48u, 57u, 
	9u, 44u, 48u, 57u, 9u, 44u, 48u, 57u, 
	9u, 44u, 48u, 57u, 9u, 44u, 48u, 57u, 
	9u, 44u, 48u, 57u, 9u, 44u, 9u, 44u, 
	69u, 101u, 48u, 57u, 9u, 44u, 48u, 57u, 
	9u, 44u, 46u, 69u, 101u, 48u, 57u, 9u, 
	44u, 9u, 48u, 57u, 65u, 70u, 97u, 102u, 
	9u, 32u, 126u, 9u, 69u, 101u, 48u, 57u, 
	9u, 48u, 57u, 9u, 46u, 69u, 101u, 48u, 
	57u, 9u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 9u, 0
];

static byte[] _sam_alignment_single_lengths = [
	0, 1, 0, 1, 1, 1, 0, 1, 
	0, 1, 1, 1, 2, 1, 0, 1, 
	2, 0, 1, 3, 1, 0, 0, 0, 
	1, 6, 1, 0, 1, 7, 1, 2, 
	0, 1, 5, 2, 0, 2, 0, 1, 
	1, 1, 1, 1, 0, 1, 0, 1, 
	5, 2, 0, 2, 0, 1, 1, 1, 
	1, 1, 2, 0, 3, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 5, 5, 5, 5, 5, 5, 5, 
	5, 5, 5, 5, 5, 5, 5, 5, 
	5, 5, 5, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 4, 2, 5, 2, 1, 
	1, 3, 1, 4, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 1
];

static byte[] _sam_alignment_range_lengths = [
	0, 2, 1, 1, 2, 1, 1, 1, 
	1, 1, 1, 0, 1, 1, 1, 1, 
	1, 1, 1, 2, 0, 1, 2, 3, 
	0, 0, 0, 1, 0, 0, 0, 1, 
	1, 0, 1, 1, 1, 1, 1, 0, 
	0, 0, 0, 0, 3, 0, 1, 0, 
	1, 1, 1, 1, 1, 0, 0, 0, 
	0, 0, 1, 1, 2, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 0, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 0, 0, 
	0, 3, 3, 3, 3, 3, 3, 3, 
	3, 3, 3, 3, 3, 3, 3, 3, 
	3, 3, 2, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 0, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 0, 0, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 0, 
	2, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 1, 
	0, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 0, 1, 1, 1, 0, 3, 
	1, 1, 1, 1, 0, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0
];

static short[] _sam_alignment_index_offsets = [
	0, 0, 4, 6, 9, 13, 16, 18, 
	21, 23, 26, 29, 31, 35, 38, 40, 
	43, 47, 49, 52, 58, 60, 62, 65, 
	69, 71, 78, 80, 82, 84, 92, 94, 
	98, 100, 102, 109, 113, 115, 119, 121, 
	123, 125, 127, 129, 131, 135, 137, 139, 
	141, 148, 152, 154, 158, 160, 162, 164, 
	166, 168, 170, 174, 176, 182, 185, 188, 
	191, 194, 197, 200, 203, 206, 209, 212, 
	215, 218, 221, 224, 227, 230, 232, 235, 
	238, 241, 244, 247, 250, 253, 256, 259, 
	262, 265, 268, 271, 274, 277, 280, 282, 
	284, 286, 295, 304, 313, 322, 331, 340, 
	349, 358, 367, 376, 385, 394, 403, 412, 
	421, 430, 439, 447, 450, 453, 456, 459, 
	462, 465, 468, 471, 474, 477, 480, 483, 
	486, 489, 492, 495, 498, 500, 503, 506, 
	509, 512, 515, 518, 521, 524, 527, 530, 
	533, 536, 539, 542, 545, 548, 550, 552, 
	555, 558, 561, 564, 567, 570, 573, 576, 
	579, 582, 585, 588, 591, 594, 597, 600, 
	602, 606, 608, 610, 612, 614, 616, 618, 
	620, 622, 624, 626, 628, 630, 632, 634, 
	636, 638, 640, 642, 644, 646, 648, 650, 
	653, 655, 659, 663, 667, 671, 675, 679, 
	683, 687, 691, 695, 699, 703, 707, 711, 
	715, 719, 723, 726, 732, 736, 743, 746, 
	751, 754, 759, 762, 768, 770, 773, 776, 
	779, 782, 785, 788, 791, 794, 797, 800, 
	803, 806, 809, 812, 815, 818, 821, 823, 
	824, 825, 826, 827, 828, 829, 830, 831, 
	832, 833, 834, 835
];

static ubyte[] _sam_alignment_trans_targs = [
	2, 168, 168, 0, 3, 0, 4, 151, 
	0, 150, 5, 5, 0, 6, 5, 0, 
	7, 0, 8, 133, 0, 9, 0, 10, 
	116, 0, 11, 97, 0, 12, 0, 95, 
	96, 13, 0, 14, 13, 0, 15, 0, 
	16, 78, 0, 17, 17, 18, 0, 18, 
	0, 19, 61, 0, 20, 60, 60, 60, 
	60, 0, 21, 0, 191, 0, 23, 23, 
	0, 24, 24, 24, 0, 25, 0, 26, 
	28, 43, 45, 47, 57, 0, 27, 0, 
	192, 0, 29, 0, 30, 30, 30, 30, 
	33, 30, 30, 0, 31, 0, 32, 32, 
	193, 0, 193, 0, 34, 0, 35, 35, 
	36, 39, 41, 213, 0, 36, 39, 213, 
	0, 211, 0, 38, 38, 212, 0, 212, 
	0, 40, 0, 214, 0, 42, 0, 214, 
	0, 44, 0, 215, 215, 215, 0, 46, 
	0, 216, 0, 48, 0, 49, 49, 50, 
	53, 55, 219, 0, 50, 53, 219, 0, 
	217, 0, 52, 52, 218, 0, 218, 0, 
	54, 0, 220, 0, 56, 0, 220, 0, 
	58, 0, 59, 59, 221, 0, 221, 0, 
	21, 60, 60, 60, 60, 0, 19, 62, 
	0, 19, 63, 0, 19, 64, 0, 19, 
	65, 0, 19, 66, 0, 19, 67, 0, 
	19, 68, 0, 19, 69, 0, 19, 70, 
	0, 19, 71, 0, 19, 72, 0, 19, 
	73, 0, 19, 74, 0, 19, 75, 0, 
	19, 76, 0, 19, 77, 0, 19, 0, 
	16, 79, 0, 16, 80, 0, 16, 81, 
	0, 16, 82, 0, 16, 83, 0, 16, 
	84, 0, 16, 85, 0, 16, 86, 0, 
	16, 87, 0, 16, 88, 0, 16, 89, 
	0, 16, 90, 0, 16, 91, 0, 16, 
	92, 0, 16, 93, 0, 16, 94, 0, 
	16, 0, 14, 0, 14, 0, 115, 115, 
	115, 115, 115, 98, 115, 115, 0, 115, 
	115, 115, 115, 115, 99, 115, 115, 0, 
	115, 115, 115, 115, 115, 100, 115, 115, 
	0, 115, 115, 115, 115, 115, 101, 115, 
	115, 0, 115, 115, 115, 115, 115, 102, 
	115, 115, 0, 115, 115, 115, 115, 115, 
	103, 115, 115, 0, 115, 115, 115, 115, 
	115, 104, 115, 115, 0, 115, 115, 115, 
	115, 115, 105, 115, 115, 0, 115, 115, 
	115, 115, 115, 106, 115, 115, 0, 115, 
	115, 115, 115, 115, 107, 115, 115, 0, 
	115, 115, 115, 115, 115, 108, 115, 115, 
	0, 115, 115, 115, 115, 115, 109, 115, 
	115, 0, 115, 115, 115, 115, 115, 110, 
	115, 115, 0, 115, 115, 115, 115, 115, 
	111, 115, 115, 0, 115, 115, 115, 115, 
	115, 112, 115, 115, 0, 115, 115, 115, 
	115, 115, 113, 115, 115, 0, 115, 115, 
	115, 115, 115, 114, 115, 115, 0, 115, 
	115, 115, 115, 115, 115, 115, 0, 12, 
	97, 0, 10, 117, 0, 10, 118, 0, 
	10, 119, 0, 10, 120, 0, 10, 121, 
	0, 10, 122, 0, 10, 123, 0, 10, 
	124, 0, 10, 125, 0, 10, 126, 0, 
	10, 127, 0, 10, 128, 0, 10, 129, 
	0, 10, 130, 0, 10, 131, 0, 10, 
	132, 0, 10, 0, 8, 134, 0, 8, 
	135, 0, 8, 136, 0, 8, 137, 0, 
	8, 138, 0, 8, 139, 0, 8, 140, 
	0, 8, 141, 0, 8, 142, 0, 8, 
	143, 0, 8, 144, 0, 8, 145, 0, 
	8, 146, 0, 8, 147, 0, 8, 148, 
	0, 8, 149, 0, 8, 0, 6, 0, 
	4, 152, 0, 4, 153, 0, 4, 154, 
	0, 4, 155, 0, 4, 156, 0, 4, 
	157, 0, 4, 158, 0, 4, 159, 0, 
	4, 160, 0, 4, 161, 0, 4, 162, 
	0, 4, 163, 0, 4, 164, 0, 4, 
	165, 0, 4, 166, 0, 4, 167, 0, 
	4, 0, 2, 168, 168, 0, 239, 169, 
	240, 170, 241, 171, 242, 172, 243, 173, 
	244, 174, 245, 175, 246, 176, 247, 177, 
	248, 178, 249, 179, 250, 180, 2, 0, 
	4, 0, 6, 0, 8, 0, 10, 0, 
	12, 0, 14, 0, 16, 0, 19, 0, 
	21, 0, 22, 191, 0, 22, 0, 22, 
	31, 194, 0, 22, 31, 195, 0, 22, 
	31, 196, 0, 22, 31, 197, 0, 22, 
	31, 198, 0, 22, 31, 199, 0, 22, 
	31, 200, 0, 22, 31, 201, 0, 22, 
	31, 202, 0, 22, 31, 203, 0, 22, 
	31, 204, 0, 22, 31, 205, 0, 22, 
	31, 206, 0, 22, 31, 207, 0, 22, 
	31, 208, 0, 22, 31, 209, 0, 22, 
	31, 210, 0, 22, 31, 0, 22, 34, 
	37, 37, 211, 0, 22, 34, 212, 0, 
	22, 34, 36, 37, 37, 213, 0, 22, 
	34, 0, 22, 215, 215, 215, 0, 22, 
	216, 0, 22, 51, 51, 217, 0, 22, 
	218, 0, 22, 50, 51, 51, 219, 0, 
	22, 0, 22, 222, 0, 22, 223, 0, 
	22, 224, 0, 22, 225, 0, 22, 226, 
	0, 22, 227, 0, 22, 228, 0, 22, 
	229, 0, 22, 230, 0, 22, 231, 0, 
	22, 232, 0, 22, 233, 0, 22, 234, 
	0, 22, 235, 0, 22, 236, 0, 22, 
	237, 0, 22, 238, 0, 22, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 22, 0, 0
];

static ubyte[] _sam_alignment_trans_actions = [
	111, 7, 7, 11, 99, 17, 15, 3, 
	17, 0, 21, 21, 25, 23, 0, 25, 
	99, 31, 29, 3, 31, 99, 35, 114, 
	3, 35, 0, 99, 41, 45, 41, 0, 
	0, 49, 53, 51, 0, 53, 99, 59, 
	57, 3, 59, 1, 1, 99, 63, 99, 
	63, 102, 3, 63, 0, 67, 67, 67, 
	67, 71, 75, 71, 77, 79, 89, 89, 
	0, 0, 0, 0, 0, 91, 0, 0, 
	0, 0, 0, 0, 0, 93, 0, 93, 
	83, 93, 0, 93, 87, 87, 87, 87, 
	87, 87, 87, 93, 0, 93, 1, 1, 
	99, 93, 99, 93, 0, 93, 5, 5, 
	5, 5, 5, 5, 93, 0, 0, 0, 
	93, 0, 93, 0, 0, 0, 93, 0, 
	93, 0, 93, 0, 93, 0, 93, 0, 
	93, 0, 93, 85, 85, 85, 93, 0, 
	93, 85, 93, 0, 93, 5, 5, 5, 
	5, 5, 5, 93, 0, 0, 0, 93, 
	0, 93, 0, 0, 0, 93, 0, 93, 
	0, 93, 0, 93, 0, 93, 0, 93, 
	0, 93, 1, 1, 99, 93, 99, 93, 
	69, 0, 0, 0, 0, 71, 102, 3, 
	63, 102, 3, 63, 102, 3, 63, 102, 
	3, 63, 102, 3, 63, 102, 3, 63, 
	102, 3, 63, 102, 3, 63, 102, 3, 
	63, 102, 3, 63, 102, 3, 63, 102, 
	3, 63, 102, 3, 63, 102, 3, 63, 
	102, 3, 63, 102, 3, 63, 102, 63, 
	57, 3, 59, 57, 3, 59, 57, 3, 
	59, 57, 3, 59, 57, 3, 59, 57, 
	3, 59, 57, 3, 59, 57, 3, 59, 
	57, 3, 59, 57, 3, 59, 57, 3, 
	59, 57, 3, 59, 57, 3, 59, 57, 
	3, 59, 57, 3, 59, 57, 3, 59, 
	57, 59, 0, 53, 47, 53, 117, 117, 
	117, 117, 117, 3, 117, 117, 41, 117, 
	117, 117, 117, 117, 3, 117, 117, 41, 
	117, 117, 117, 117, 117, 3, 117, 117, 
	41, 117, 117, 117, 117, 117, 3, 117, 
	117, 41, 117, 117, 117, 117, 117, 3, 
	117, 117, 41, 117, 117, 117, 117, 117, 
	3, 117, 117, 41, 117, 117, 117, 117, 
	117, 3, 117, 117, 41, 117, 117, 117, 
	117, 117, 3, 117, 117, 41, 117, 117, 
	117, 117, 117, 3, 117, 117, 41, 117, 
	117, 117, 117, 117, 3, 117, 117, 41, 
	117, 117, 117, 117, 117, 3, 117, 117, 
	41, 117, 117, 117, 117, 117, 3, 117, 
	117, 41, 117, 117, 117, 117, 117, 3, 
	117, 117, 41, 117, 117, 117, 117, 117, 
	3, 117, 117, 41, 117, 117, 117, 117, 
	117, 3, 117, 117, 41, 117, 117, 117, 
	117, 117, 3, 117, 117, 41, 117, 117, 
	117, 117, 117, 3, 117, 117, 41, 117, 
	117, 117, 117, 117, 117, 117, 41, 120, 
	148, 41, 114, 3, 35, 114, 3, 35, 
	114, 3, 35, 114, 3, 35, 114, 3, 
	35, 114, 3, 35, 114, 3, 35, 114, 
	3, 35, 114, 3, 35, 114, 3, 35, 
	114, 3, 35, 114, 3, 35, 114, 3, 
	35, 114, 3, 35, 114, 3, 35, 114, 
	3, 35, 114, 35, 29, 3, 31, 29, 
	3, 31, 29, 3, 31, 29, 3, 31, 
	29, 3, 31, 29, 3, 31, 29, 3, 
	31, 29, 3, 31, 29, 3, 31, 29, 
	3, 31, 29, 3, 31, 29, 3, 31, 
	29, 3, 31, 29, 3, 31, 29, 3, 
	31, 29, 3, 31, 29, 31, 0, 25, 
	15, 3, 17, 15, 3, 17, 15, 3, 
	17, 15, 3, 17, 15, 3, 17, 15, 
	3, 17, 15, 3, 17, 15, 3, 17, 
	15, 3, 17, 15, 3, 17, 15, 3, 
	17, 15, 3, 17, 15, 3, 17, 15, 
	3, 17, 15, 3, 17, 15, 3, 17, 
	15, 17, 9, 0, 0, 11, 13, 0, 
	19, 0, 27, 0, 33, 0, 37, 0, 
	43, 0, 55, 0, 61, 0, 65, 0, 
	73, 0, 81, 0, 95, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 39, 0, 
	45, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 123, 77, 79, 97, 93, 136, 
	105, 3, 93, 136, 105, 3, 93, 136, 
	105, 3, 93, 136, 105, 3, 93, 136, 
	105, 3, 93, 136, 105, 3, 93, 136, 
	105, 3, 93, 136, 105, 3, 93, 136, 
	105, 3, 93, 136, 105, 3, 93, 136, 
	105, 3, 93, 136, 105, 3, 93, 136, 
	105, 3, 93, 136, 105, 3, 93, 136, 
	105, 3, 93, 136, 105, 3, 93, 136, 
	105, 3, 93, 136, 105, 93, 144, 108, 
	0, 0, 0, 93, 144, 108, 0, 93, 
	144, 108, 0, 0, 0, 0, 93, 144, 
	108, 93, 129, 0, 0, 0, 93, 126, 
	0, 93, 140, 0, 0, 0, 93, 140, 
	0, 93, 140, 0, 0, 0, 0, 93, 
	140, 93, 132, 3, 93, 132, 3, 93, 
	132, 3, 93, 132, 3, 93, 132, 3, 
	93, 132, 3, 93, 132, 3, 93, 132, 
	3, 93, 132, 3, 93, 132, 3, 93, 
	132, 3, 93, 132, 3, 93, 132, 3, 
	93, 132, 3, 93, 132, 3, 93, 132, 
	3, 93, 132, 3, 93, 132, 93, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0
];

static ubyte[] _sam_alignment_eof_actions = [
	0, 11, 17, 17, 25, 25, 31, 31, 
	35, 35, 41, 41, 53, 53, 59, 59, 
	63, 63, 63, 71, 71, 79, 0, 0, 
	0, 93, 93, 93, 93, 93, 93, 93, 
	93, 93, 93, 93, 93, 93, 93, 93, 
	93, 93, 93, 93, 93, 93, 93, 93, 
	93, 93, 93, 93, 93, 93, 93, 93, 
	93, 93, 93, 93, 71, 63, 63, 63, 
	63, 63, 63, 63, 63, 63, 63, 63, 
	63, 63, 63, 63, 63, 63, 59, 59, 
	59, 59, 59, 59, 59, 59, 59, 59, 
	59, 59, 59, 59, 59, 59, 59, 53, 
	53, 41, 41, 41, 41, 41, 41, 41, 
	41, 41, 41, 41, 41, 41, 41, 41, 
	41, 41, 41, 41, 35, 35, 35, 35, 
	35, 35, 35, 35, 35, 35, 35, 35, 
	35, 35, 35, 35, 35, 31, 31, 31, 
	31, 31, 31, 31, 31, 31, 31, 31, 
	31, 31, 31, 31, 31, 31, 25, 17, 
	17, 17, 17, 17, 17, 17, 17, 17, 
	17, 17, 17, 17, 17, 17, 17, 17, 
	11, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 123, 
	97, 136, 136, 136, 136, 136, 136, 136, 
	136, 136, 136, 136, 136, 136, 136, 136, 
	136, 136, 136, 144, 144, 144, 144, 129, 
	126, 140, 140, 140, 140, 132, 132, 132, 
	132, 132, 132, 132, 132, 132, 132, 132, 
	132, 132, 132, 132, 132, 132, 132, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0
];

static int sam_alignment_start = 1;
static int sam_alignment_first_final = 191;
static int sam_alignment_error = 0;

static int sam_alignment_en_recover_from_invalid_qname = 169;
static int sam_alignment_en_recover_from_invalid_flag = 170;
static int sam_alignment_en_recover_from_invalid_rname = 171;
static int sam_alignment_en_recover_from_invalid_pos = 172;
static int sam_alignment_en_recover_from_invalid_mapq = 173;
static int sam_alignment_en_recover_from_invalid_cigar = 174;
static int sam_alignment_en_recover_from_invalid_rnext = 175;
static int sam_alignment_en_recover_from_invalid_pnext = 176;
static int sam_alignment_en_recover_from_invalid_tlen = 177;
static int sam_alignment_en_recover_from_invalid_seq = 178;
static int sam_alignment_en_recover_from_invalid_qual = 179;
static int sam_alignment_en_recover_from_invalid_tag = 180;
static int sam_alignment_en_alignment = 1;
static int sam_alignment_en_alignment_field_parsing_mandatoryfields_flag_parsing = 181;
static int sam_alignment_en_alignment_field_parsing_mandatoryfields_rname_parsing = 182;
static int sam_alignment_en_alignment_field_parsing_mandatoryfields_pos_parsing = 183;
static int sam_alignment_en_alignment_field_parsing_mandatoryfields_mapq_parsing = 184;
static int sam_alignment_en_alignment_field_parsing_mandatoryfields_cigar_parsing = 185;
static int sam_alignment_en_alignment_field_parsing_mandatoryfields_rnext_parsing = 186;
static int sam_alignment_en_alignment_field_parsing_mandatoryfields_pnext_parsing = 187;
static int sam_alignment_en_alignment_field_parsing_mandatoryfields_tlen_parsing = 188;
static int sam_alignment_en_alignment_field_parsing_mandatoryfields_seq_parsing = 189;
static int sam_alignment_en_alignment_field_parsing_mandatoryfields_qual_parsing = 190;
static int sam_alignment_en_alignment_tag_parsing = 251;


#line 419 "sam_alignment.rl"


import bio.std.hts.sam.header;
import bio.std.hts.bam.cigar;
import bio.std.hts.bam.read;
import bio.std.hts.bam.bai.bin;
import bio.core.utils.outbuffer;
import bio.core.base;
import std.conv;
import std.array;
import std.exception;

BamRead parseAlignmentLine(string line, SamHeader header, OutBuffer buffer=null) {
    char* p = cast(char*)line.ptr;
    char* pe = p + line.length;
    char* eof = pe;
    int cs;

    if (buffer is null)
        buffer = new OutBuffer(8192);
    else
        buffer.clear();

    size_t rollback_size; // needed in case of invalid data

    byte current_sign = 1;

    size_t read_name_beg; // position of beginning of QNAME

    size_t sequence_beg; // position of SEQ start
    int l_seq;           // sequence length

    uint cigar_op_len;   // length of CIGAR operation
    char cigar_op_chr;   // CIGAR operation

    size_t quals_length;  // number of QUAL characters
    char quals_last_char; // needed in order to handle '*' correctly

    size_t cigar_op_len_start; // position of start of CIGAR operation

    long int_value;                      // for storing temporary integers
    float float_value;                   // for storing temporary floats
    size_t float_beg;                    // position of start of current float
    char arraytype;                      // type of last array tag value
    size_t tag_array_length_offset;      // where the length is stored in the buffer

    string read_name;
    ushort flag;
    int pos = -1;
    int end_pos; // for bin calculation
    int mate_pos = -1;
    ubyte mapping_quality = 255;
    int template_length = 0;

    size_t tag_key_beg, tagvalue_beg;
    ubyte[] tag_key;
    size_t rname_beg, rnext_beg;

    int ref_id = -1;

    
#line 640 "sam_alignment.d"
	{
	cs = sam_alignment_start;
	}

#line 480 "sam_alignment.rl"
    
#line 647 "sam_alignment.d"
	{
	int _klen;
	uint _trans;
	byte* _acts;
	uint _nacts;
	char* _keys;

	if ( p == pe )
		goto _test_eof;
	if ( cs == 0 )
		goto _out;
_resume:
	_keys = &_sam_alignment_trans_keys[_sam_alignment_key_offsets[cs]];
	_trans = _sam_alignment_index_offsets[cs];

	_klen = _sam_alignment_single_lengths[cs];
	if ( _klen > 0 ) {
		char* _lower = _keys;
		char* _mid;
		char* _upper = _keys + _klen - 1;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + ((_upper-_lower) >> 1);
			if ( (*p) < *_mid )
				_upper = _mid - 1;
			else if ( (*p) > *_mid )
				_lower = _mid + 1;
			else {
				_trans += cast(uint)(_mid - _keys);
				goto _match;
			}
		}
		_keys += _klen;
		_trans += _klen;
	}

	_klen = _sam_alignment_range_lengths[cs];
	if ( _klen > 0 ) {
		char* _lower = _keys;
		char* _mid;
		char* _upper = _keys + (_klen<<1) - 2;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + (((_upper-_lower) >> 1) & ~1);
			if ( (*p) < _mid[0] )
				_upper = _mid - 2;
			else if ( (*p) > _mid[1] )
				_lower = _mid + 2;
			else {
				_trans += cast(uint)((_mid - _keys)>>1);
				goto _match;
			}
		}
		_trans += _klen;
	}

_match:
	cs = _sam_alignment_trans_targs[_trans];

	if ( _sam_alignment_trans_actions[_trans] == 0 )
		goto _again;

	_acts = &_sam_alignment_actions[_sam_alignment_trans_actions[_trans]];
	_nacts = cast(uint) *_acts++;
	while ( _nacts-- > 0 )
	{
		switch ( *_acts++ )
		{
	case 0:
#line 27 "sam_alignment.rl"
	{ current_sign = (*p) == '-' ? -1 : 1; }
	break;
	case 1:
#line 28 "sam_alignment.rl"
	{ int_value = 0; }
	break;
	case 2:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	break;
	case 3:
#line 30 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
	break;
	case 4:
#line 37 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	break;
	case 5:
#line 38 "sam_alignment.rl"
	{
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }
	break;
	case 6:
#line 48 "sam_alignment.rl"
	{ read_name_beg = p - line.ptr; }
	break;
	case 7:
#line 49 "sam_alignment.rl"
	{ read_name = line[read_name_beg .. p - line.ptr]; }
	break;
	case 8:
#line 50 "sam_alignment.rl"
	{ p--; {cs = 169; if (true) goto _again;} }
	break;
	case 9:
#line 51 "sam_alignment.rl"
	{ p--; {cs = 181; if (true) goto _again;} }
	break;
	case 10:
#line 56 "sam_alignment.rl"
	{ flag = to!ushort(int_value); }
	break;
	case 11:
#line 58 "sam_alignment.rl"
	{ p--; {cs = 170; if (true) goto _again;} }
	break;
	case 12:
#line 59 "sam_alignment.rl"
	{ p--; {cs = 182; if (true) goto _again;} }
	break;
	case 13:
#line 62 "sam_alignment.rl"
	{ rname_beg = p - line.ptr; }
	break;
	case 14:
#line 63 "sam_alignment.rl"
	{
        ref_id = header.getSequenceIndex(line[rname_beg .. p - line.ptr]);
    }
	break;
	case 15:
#line 67 "sam_alignment.rl"
	{ p--; {cs = 171; if (true) goto _again;} }
	break;
	case 16:
#line 68 "sam_alignment.rl"
	{ p--; {cs = 183; if (true) goto _again;} }
	break;
	case 17:
#line 73 "sam_alignment.rl"
	{ end_pos = pos = to!uint(int_value); }
	break;
	case 18:
#line 75 "sam_alignment.rl"
	{ p--; {cs = 172; if (true) goto _again;} }
	break;
	case 19:
#line 76 "sam_alignment.rl"
	{ p--; {cs = 184; if (true) goto _again;} }
	break;
	case 20:
#line 79 "sam_alignment.rl"
	{ mapping_quality = to!ubyte(int_value); }
	break;
	case 21:
#line 81 "sam_alignment.rl"
	{ p--; {cs = 173; if (true) goto _again;} }
	break;
	case 22:
#line 82 "sam_alignment.rl"
	{ p--; {cs = 185; if (true) goto _again;} }
	break;
	case 23:
#line 85 "sam_alignment.rl"
	{
        buffer.capacity = 32 + read_name.length + 1;
        buffer.putUnsafe!int(ref_id);
        buffer.putUnsafe!int(pos - 1);

        enforce(read_name.length + 1 <= 255, "Read name " ~ read_name ~ " is too long!");

        // bin will be set later
        auto bin_mq_nl = ((cast(uint)mapping_quality) << 8) | (read_name.length + 1);
        buffer.putUnsafe(cast(uint)bin_mq_nl);

        // number of CIGAR operations will be set later
        buffer.putUnsafe!uint(flag << 16);

        buffer.putUnsafe!int(0);
        buffer.putUnsafe!int(-1); // mate ref. id
        buffer.putUnsafe!int(-1); // mate pos
        buffer.putUnsafe!int(0);  // tlen

        buffer.putUnsafe(cast(ubyte[])read_name);
        buffer.putUnsafe!ubyte(0);

        rollback_size = buffer.length;
    }
	break;
	case 24:
#line 111 "sam_alignment.rl"
	{ cigar_op_len = to!uint(int_value); }
	break;
	case 25:
#line 112 "sam_alignment.rl"
	{ cigar_op_chr = (*p); }
	break;
	case 26:
#line 113 "sam_alignment.rl"
	{
        auto op = CigarOperation(cigar_op_len, cigar_op_chr);
        if (op.is_reference_consuming)
            end_pos += op.length;
        buffer.put!CigarOperation(op);
        {
        auto ptr = cast(uint*)(buffer.data.ptr + 3 * uint.sizeof);
        *ptr = (*ptr) + 1;
        }
    }
	break;
	case 27:
#line 124 "sam_alignment.rl"
	{
        auto ptr = cast(uint*)(buffer.data.ptr + 3 * uint.sizeof);
        *ptr = (*ptr) & 0xFFFF0000;
        buffer.shrink(rollback_size);
        end_pos = pos + 1;
        p--; {cs = 174; if (true) goto _again;}
    }
	break;
	case 28:
#line 131 "sam_alignment.rl"
	{ p--; {cs = 186; if (true) goto _again;} }
	break;
	case 29:
#line 137 "sam_alignment.rl"
	{
        if (end_pos == pos)
            ++end_pos;
        {
        auto bin = reg2bin(pos - 1, end_pos - 1); // 0-based [) interval
        auto ptr = cast(uint*)(buffer.data.ptr + 2 * uint.sizeof);
        *ptr = (*ptr) | ((cast(uint)bin) << 16);
        }
    }
	break;
	case 30:
#line 148 "sam_alignment.rl"
	{
        {
        auto ptr = cast(int*)(buffer.data.ptr + 5 * int.sizeof);
        *ptr = ref_id;
        }
    }
	break;
	case 31:
#line 155 "sam_alignment.rl"
	{ rnext_beg = p - line.ptr; }
	break;
	case 32:
#line 156 "sam_alignment.rl"
	{
        {
        auto ptr = cast(int*)(buffer.data.ptr + 5 * int.sizeof);
        *ptr = header.getSequenceIndex(line[rnext_beg .. p - line.ptr]);
        }
    }
	break;
	case 33:
#line 162 "sam_alignment.rl"
	{ p--; {cs = 175; if (true) goto _again;} }
	break;
	case 34:
#line 163 "sam_alignment.rl"
	{ p--; {cs = 187; if (true) goto _again;} }
	break;
	case 35:
#line 169 "sam_alignment.rl"
	{
        {
        auto ptr = cast(int*)(buffer.data.ptr + 6 * int.sizeof);
        *ptr = to!int(int_value) - 1;
        }
    }
	break;
	case 36:
#line 175 "sam_alignment.rl"
	{ p--; {cs = 176; if (true) goto _again;} }
	break;
	case 37:
#line 176 "sam_alignment.rl"
	{ p--; {cs = 188; if (true) goto _again;} }
	break;
	case 38:
#line 181 "sam_alignment.rl"
	{
        {
        auto ptr = cast(int*)(buffer.data.ptr + 7 * int.sizeof);
        *ptr = to!int(int_value);
        }
    }
	break;
	case 39:
#line 187 "sam_alignment.rl"
	{ p--; {cs = 177; if (true) goto _again;} }
	break;
	case 40:
#line 188 "sam_alignment.rl"
	{ p--; {cs = 189; if (true) goto _again;} }
	break;
	case 41:
#line 193 "sam_alignment.rl"
	{ sequence_beg = p - line.ptr; }
	break;
	case 42:
#line 194 "sam_alignment.rl"
	{
        auto data = cast(ubyte[])line[sequence_beg .. p - line.ptr];
        l_seq = cast(int)data.length;
        auto raw_len = (l_seq + 1) / 2;

        // reserve space for base qualities, too
        buffer.capacity = buffer.length + raw_len + l_seq;

        for (size_t i = 0; i < raw_len; ++i) {
            auto b = cast(ubyte)(Base(data[2 * i]).internal_code << 4);
            if (2 * i + 1 < l_seq)
                b |= cast(ubyte)(Base(data[2 * i + 1]).internal_code);
            buffer.putUnsafe!ubyte(b);
        }

        // set l_seq
        {
        auto ptr = cast(int*)(buffer.data.ptr + 4 * int.sizeof);
        *ptr = l_seq;
        }

        rollback_size = buffer.length;
    }
	break;
	case 43:
#line 217 "sam_alignment.rl"
	{
        rollback_size = buffer.length;
        p--; {cs = 178; if (true) goto _again;}
    }
	break;
	case 44:
#line 221 "sam_alignment.rl"
	{ p--; {cs = 190; if (true) goto _again;} }
	break;
	case 45:
#line 223 "sam_alignment.rl"
	{
        rollback_size = buffer.length;
    }
	break;
	case 46:
#line 230 "sam_alignment.rl"
	{
        ++quals_length;
        quals_last_char = (*p);
        buffer.putUnsafe!ubyte(cast(ubyte)((*p) - 33));
    }
	break;
	case 47:
#line 236 "sam_alignment.rl"
	{
        // '*' may correspond either to a one-base long sequence
        // or to absence of information
        if (quals_length == 1 && quals_last_char == '*' && l_seq == 0)
            buffer.shrink(rollback_size);
    }
	break;
	case 48:
#line 243 "sam_alignment.rl"
	{
        buffer.shrink(rollback_size);
        for (size_t i = 0; i < l_seq; ++i)
            buffer.putUnsafe!ubyte(0xFF);
        rollback_size = buffer.length;
        p--; {cs = 179; if (true) goto _again;}
    }
	break;
	case 49:
#line 251 "sam_alignment.rl"
	{ p--; {cs = 251; if (true) goto _again;} }
	break;
	case 50:
#line 253 "sam_alignment.rl"
	{
        if (buffer.length - rollback_size != l_seq) {
            buffer.shrink(rollback_size);
            for (size_t i = 0; i < l_seq; ++i)
                buffer.putUnsafe!ubyte(0xFF);
        }
        rollback_size = buffer.length;
    }
	break;
	case 51:
#line 278 "sam_alignment.rl"
	{
        buffer.capacity = buffer.length + 4;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('A');
        buffer.putUnsafe!char((*p));
    }
	break;
	case 52:
#line 285 "sam_alignment.rl"
	{
        buffer.capacity = buffer.length + 7;
        buffer.putUnsafe(tag_key);
        if (int_value < 0) {
            if (int_value >= byte.min) {
                buffer.putUnsafe!char('c');
                buffer.putUnsafe(cast(byte)int_value);
            } else if (int_value >= short.min) {
                buffer.putUnsafe!char('s');
                buffer.putUnsafe(cast(short)int_value);
            } else if (int_value >= int.min) {
                buffer.putUnsafe!char('i');
                buffer.putUnsafe(cast(int)int_value);
            } else {
                throw new Exception("integer out of range");
            }
        } else {
            if (int_value <= ubyte.max) {
                buffer.putUnsafe!char('C');
                buffer.putUnsafe(cast(ubyte)int_value);
            } else if (int_value <= ushort.max) {
                buffer.putUnsafe!char('S');
                buffer.putUnsafe(cast(ushort)int_value);
            } else if (int_value <= uint.max) {
                buffer.putUnsafe!char('I');
                buffer.putUnsafe(cast(uint)int_value);
            } else {
                throw new Exception("integer out of range");
            }
        }
    }
	break;
	case 53:
#line 317 "sam_alignment.rl"
	{ tagvalue_beg = p - line.ptr; }
	break;
	case 54:
#line 319 "sam_alignment.rl"
	{
        buffer.capacity = buffer.length + 7;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('f');
        buffer.putUnsafe!float(float_value);
    }
	break;
	case 55:
#line 326 "sam_alignment.rl"
	{
        {
        auto data = cast(ubyte[])(line[tagvalue_beg .. p - line.ptr]);
        buffer.capacity = buffer.length + 4 + data.length;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('Z');
        buffer.putUnsafe(data);
        buffer.putUnsafe!ubyte(0);
        }
    }
	break;
	case 56:
#line 337 "sam_alignment.rl"
	{
        {
        auto data = cast(ubyte[])(line[tagvalue_beg .. p - line.ptr]);
        buffer.capacity = buffer.length + 4 + data.length;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('H');
        buffer.putUnsafe(data);
        buffer.putUnsafe!ubyte(0);
        }
    }
	break;
	case 57:
#line 352 "sam_alignment.rl"
	{
        arraytype = (*p);
        buffer.capacity = buffer.length + 8;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('B');
        buffer.putUnsafe!char(arraytype);
        buffer.putUnsafe!uint(0);
        tag_array_length_offset = buffer.length - uint.sizeof;
    }
	break;
	case 58:
#line 362 "sam_alignment.rl"
	{
        // here, we assume that compiler is smart enough to move switch out of loop.
        switch (arraytype) {
            case 'c': buffer.put(to!byte(int_value)); break;
            case 'C': buffer.put(to!ubyte(int_value)); break;
            case 's': buffer.put(to!short(int_value)); break;
            case 'S': buffer.put(to!ushort(int_value)); break;
            case 'i': buffer.put(to!int(int_value)); break;
            case 'I': buffer.put(to!uint(int_value)); break;
            default: assert(0);
        }
        {
            auto ptr = cast(uint*)(buffer.data.ptr + tag_array_length_offset);
            ++*ptr;
        }
    }
	break;
	case 59:
#line 379 "sam_alignment.rl"
	{
        buffer.put!float(float_value);
        {
            auto ptr = cast(uint*)(buffer.data.ptr + tag_array_length_offset);
            ++*ptr;
        }
    }
	break;
	case 60:
#line 400 "sam_alignment.rl"
	{ tag_key_beg = p - line.ptr; }
	break;
	case 61:
#line 401 "sam_alignment.rl"
	{ tag_key = cast(ubyte[])(line[tag_key_beg .. p - line.ptr]); }
	break;
	case 62:
#line 403 "sam_alignment.rl"
	{
        buffer.shrink(rollback_size);
        p--; {cs = 180; if (true) goto _again;}
    }
	break;
	case 63:
#line 408 "sam_alignment.rl"
	{ p--; {cs = 251; if (true) goto _again;} }
	break;
	case 64:
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	break;
#line 1189 "sam_alignment.d"
		default: break;
		}
	}

_again:
	if ( cs == 0 )
		goto _out;
	if ( ++p != pe )
		goto _resume;
	_test_eof: {}
	if ( p == eof )
	{
	byte* __acts = &_sam_alignment_actions[_sam_alignment_eof_actions[cs]];
	uint __nacts = cast(uint) *__acts++;
	while ( __nacts-- > 0 ) {
		switch ( *__acts++ ) {
	case 3:
#line 30 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
	break;
	case 5:
#line 38 "sam_alignment.rl"
	{
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }
	break;
	case 8:
#line 50 "sam_alignment.rl"
	{ p--; {cs = 169; if (true) goto _again;} }
	break;
	case 11:
#line 58 "sam_alignment.rl"
	{ p--; {cs = 170; if (true) goto _again;} }
	break;
	case 15:
#line 67 "sam_alignment.rl"
	{ p--; {cs = 171; if (true) goto _again;} }
	break;
	case 18:
#line 75 "sam_alignment.rl"
	{ p--; {cs = 172; if (true) goto _again;} }
	break;
	case 21:
#line 81 "sam_alignment.rl"
	{ p--; {cs = 173; if (true) goto _again;} }
	break;
	case 27:
#line 124 "sam_alignment.rl"
	{
        auto ptr = cast(uint*)(buffer.data.ptr + 3 * uint.sizeof);
        *ptr = (*ptr) & 0xFFFF0000;
        buffer.shrink(rollback_size);
        end_pos = pos + 1;
        p--; {cs = 174; if (true) goto _again;}
    }
	break;
	case 33:
#line 162 "sam_alignment.rl"
	{ p--; {cs = 175; if (true) goto _again;} }
	break;
	case 36:
#line 175 "sam_alignment.rl"
	{ p--; {cs = 176; if (true) goto _again;} }
	break;
	case 39:
#line 187 "sam_alignment.rl"
	{ p--; {cs = 177; if (true) goto _again;} }
	break;
	case 43:
#line 217 "sam_alignment.rl"
	{
        rollback_size = buffer.length;
        p--; {cs = 178; if (true) goto _again;}
    }
	break;
	case 47:
#line 236 "sam_alignment.rl"
	{
        // '*' may correspond either to a one-base long sequence
        // or to absence of information
        if (quals_length == 1 && quals_last_char == '*' && l_seq == 0)
            buffer.shrink(rollback_size);
    }
	break;
	case 48:
#line 243 "sam_alignment.rl"
	{
        buffer.shrink(rollback_size);
        for (size_t i = 0; i < l_seq; ++i)
            buffer.putUnsafe!ubyte(0xFF);
        rollback_size = buffer.length;
        p--; {cs = 179; if (true) goto _again;}
    }
	break;
	case 50:
#line 253 "sam_alignment.rl"
	{
        if (buffer.length - rollback_size != l_seq) {
            buffer.shrink(rollback_size);
            for (size_t i = 0; i < l_seq; ++i)
                buffer.putUnsafe!ubyte(0xFF);
        }
        rollback_size = buffer.length;
    }
	break;
	case 52:
#line 285 "sam_alignment.rl"
	{
        buffer.capacity = buffer.length + 7;
        buffer.putUnsafe(tag_key);
        if (int_value < 0) {
            if (int_value >= byte.min) {
                buffer.putUnsafe!char('c');
                buffer.putUnsafe(cast(byte)int_value);
            } else if (int_value >= short.min) {
                buffer.putUnsafe!char('s');
                buffer.putUnsafe(cast(short)int_value);
            } else if (int_value >= int.min) {
                buffer.putUnsafe!char('i');
                buffer.putUnsafe(cast(int)int_value);
            } else {
                throw new Exception("integer out of range");
            }
        } else {
            if (int_value <= ubyte.max) {
                buffer.putUnsafe!char('C');
                buffer.putUnsafe(cast(ubyte)int_value);
            } else if (int_value <= ushort.max) {
                buffer.putUnsafe!char('S');
                buffer.putUnsafe(cast(ushort)int_value);
            } else if (int_value <= uint.max) {
                buffer.putUnsafe!char('I');
                buffer.putUnsafe(cast(uint)int_value);
            } else {
                throw new Exception("integer out of range");
            }
        }
    }
	break;
	case 54:
#line 319 "sam_alignment.rl"
	{
        buffer.capacity = buffer.length + 7;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('f');
        buffer.putUnsafe!float(float_value);
    }
	break;
	case 55:
#line 326 "sam_alignment.rl"
	{
        {
        auto data = cast(ubyte[])(line[tagvalue_beg .. p - line.ptr]);
        buffer.capacity = buffer.length + 4 + data.length;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('Z');
        buffer.putUnsafe(data);
        buffer.putUnsafe!ubyte(0);
        }
    }
	break;
	case 56:
#line 337 "sam_alignment.rl"
	{
        {
        auto data = cast(ubyte[])(line[tagvalue_beg .. p - line.ptr]);
        buffer.capacity = buffer.length + 4 + data.length;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('H');
        buffer.putUnsafe(data);
        buffer.putUnsafe!ubyte(0);
        }
    }
	break;
	case 58:
#line 362 "sam_alignment.rl"
	{
        // here, we assume that compiler is smart enough to move switch out of loop.
        switch (arraytype) {
            case 'c': buffer.put(to!byte(int_value)); break;
            case 'C': buffer.put(to!ubyte(int_value)); break;
            case 's': buffer.put(to!short(int_value)); break;
            case 'S': buffer.put(to!ushort(int_value)); break;
            case 'i': buffer.put(to!int(int_value)); break;
            case 'I': buffer.put(to!uint(int_value)); break;
            default: assert(0);
        }
        {
            auto ptr = cast(uint*)(buffer.data.ptr + tag_array_length_offset);
            ++*ptr;
        }
    }
	break;
	case 59:
#line 379 "sam_alignment.rl"
	{
        buffer.put!float(float_value);
        {
            auto ptr = cast(uint*)(buffer.data.ptr + tag_array_length_offset);
            ++*ptr;
        }
    }
	break;
	case 62:
#line 403 "sam_alignment.rl"
	{
        buffer.shrink(rollback_size);
        p--; {cs = 180; if (true) goto _again;}
    }
	break;
	case 64:
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	break;
#line 1404 "sam_alignment.d"
		default: break;
		}
	}
	}

	_out: {}
	}

#line 481 "sam_alignment.rl"

    BamRead read;
    read.raw_data = buffer.data[];
    return read;
}

unittest {
    import std.algorithm;
    import std.math;

    auto line = "ERR016155.15021091\t185\t20\t60033\t25\t66S35M\t=\t60033\t0\tAGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC\t#####################################################################################################\tX0:i:1\tX1:i:0\tXC:i:35\tMD:Z:17A8A8\tRG:Z:ERR016155\tAM:i:0\tNM:i:2\tSM:i:25\tXT:A:U\tBQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\tY0:B:c,1,2,3\tY1:B:f,13.263,-3.1415,52.63461";

    auto header = new SamHeader("@SQ\tSN:20\tLN:1234567");
    auto alignment = parseAlignmentLine(line, header);
    assert(alignment.name == "ERR016155.15021091");
    assert(equal(alignment.sequence(), "AGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC"));
    assert(alignment.cigarString() == "66S35M");
    assert(alignment.flag == 185);
    assert(alignment.position == 60032);
    assert(alignment.mapping_quality == 25);
    assert(alignment.mate_position == 60032);
    assert(alignment.ref_id == 0);
    assert(alignment.mate_ref_id == 0);
    assert(to!ubyte(alignment["AM"]) == 0);
    assert(to!ubyte(alignment["SM"]) == 25);
    assert(to!string(alignment["MD"]) == "17A8A8");
    assert(equal(to!(byte[])(alignment["Y0"]), [1, 2, 3]));
    assert(equal!approxEqual(to!(float[])(alignment["Y1"]), [13.263, -3.1415, 52.63461]));
    assert(to!char(alignment["XT"]) == 'U');

    import bio.std.hts.bam.reference;

    auto info = ReferenceSequenceInfo("20", 1234567);

    auto invalid_cigar_string = "1\t100\t20\t50000\t30\tMZABC\t=\t50000\t0\tACGT\t####";
    alignment = parseAlignmentLine(invalid_cigar_string, header);
    assert(equal(alignment.sequence(), "ACGT"));

    auto invalid_tag_and_qual = "2\t100\t20\t5\t40\t27M30X5D\t=\t3\t10\tACT\t !\n\tX1:i:7\tX3:i:zzz\tX4:i:5";
    alignment = parseAlignmentLine(invalid_tag_and_qual, header);
    assert(alignment.base_qualities == [255, 255, 255]); // i.e. invalid
    assert(to!ubyte(alignment["X1"]) == 7);
    assert(alignment["X3"].is_nothing);
    assert(to!ubyte(alignment["X4"]) == 5);
}
