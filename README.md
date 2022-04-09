## JPEG Decoder


WIP toy JPEG decoder.

#### Progress so far
```
00002 :START OF IMAGE

00004 :QUANTIZATION TABLES 67
00006 :QUANTIZATION TABLE 0 0
[  4,   3,   3,   4,   6,  10,  13,  16]
[  3,   3,   4,   5,   7,  15,  16,  14]
[  4,   3,   4,   6,  10,  15,  18,  15]
[  4,   4,   6,   8,  13,  23,  21,  16]
[  5,   6,  10,  15,  18,  28,  27,  20]
[  6,   9,  14,  17,  21,  27,  29,  24]
[ 13,  17,  20,  23,  27,  31,  31,  26]
[ 19,  24,  25,  25,  29,  26,  27,  26]

00049 :QUANTIZATION TABLES 67
0004b :QUANTIZATION TABLE 0 1
[  4,   5,   6,  12,  26,  26,  26,  26]
[  5,   5,   7,  17,  26,  26,  26,  26]
[  6,   7,  15,  26,  26,  26,  26,  26]
[ 12,  17,  26,  26,  26,  26,  26,  26]
[ 26,  26,  26,  26,  26,  26,  26,  26]
[ 26,  26,  26,  26,  26,  26,  26,  26]
[ 26,  26,  26,  26,  26,  26,  26,  26]
[ 26,  26,  26,  26,  26,  26,  26,  26]

0008e :START OF FRAME 0
00090 :PRECISION 8
00091 :NUMBER OF LINES 58
00093 :NUMBER OF SAMPLES PER LINE 128
00095 :NUMBER OF COMPONENTS 3
00096 :PARAMETERS 2 2 0
00099 :PARAMETERS 1 1 1
0009c :PARAMETERS 1 1 1

000a1 :HUFFMAN TABLES 28
000a3 :HUFFMAN TABLE DC 0
[0, 1, 4, 5, 6, 14, 30, 62, 126]

000bf :HUFFMAN TABLES 53
000c1 :HUFFMAN TABLE AC 0
[0, 1, 4, 10, 11, 12, 26, 27, 28, 29, 120, 121, 122, 123, 248, 249, 250, 251, 252, 506, 507, 508, 1018, 1019, 1020, 2042, 2043, 4088, 4089, 4090, 4091, 4092, 4093, 4094]

000f6 :HUFFMAN TABLES 24
000f8 :HUFFMAN TABLE DC 1
[0, 1, 2, 6, 14]

00110 :HUFFMAN TABLES 44
00112 :HUFFMAN TABLE AC 1
[0, 1, 4, 10, 11, 24, 25, 26, 27, 56, 57, 58, 118, 119, 120, 121, 122, 123, 124, 250, 251, 252, 253, 254, 510]

0013e :START OF SCAN 12
00140 :NUMBER OF COMPONENTS 3
00141 :SCAN COMPONENT 1
00142 :DC - AC ENTROPY 0 0
00143 :SCAN COMPONENT 2
00144 :DC - AC ENTROPY 1 1
00145 :SCAN COMPONENT 3
00146 :DC - AC ENTROPY 1 1
00147 :START OF SELECTION 0
00148 :END OF SELECTION 63
00149 :HIGH - LOW BIT TRANSFORM 0 0

01060 :END OF IMAGE
```