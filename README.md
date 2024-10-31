## Overview
To begin the height, width and quality are pushed to the stream for the decompressor. Then each frame is proceeded by a single bit flag, 1 if the frame exists and 0 if there is no more frames to read.

Each frame is read as a YUVframe420 object already having the 4:2:0 subsampling applied. Then the information for each channel is separated and stored in a 2D vector object. The Y_matrix hold a Y_value for each pixel in the frame and the Cb_matrix and Cr_matrix hold the Cb and Cr values respectively for the scaled frame.
Next the matrices are "partitioned" into smaller blocks. The Y_matrix is first partitioned into 16x16 blocks; then each 16x16 block is partitioned into 4 8x8 sub-blocks. The blocks are then pushed to the Y_blocks vector in the order seen in the figure below.

The Cb_matrix and Cr_matrix are each partitioned into 8x8 blocks in row-major-order, as can be seen in the figure below.

With the blocks organized as such each macro-block(16x16) of the frame can be obtained by looking at 4 blocks of the Y_blocks vector, and first block in the Cb_blocks and Cr_blocks in order.

For each macro-block the program first looks for a "good motion vector" using the AAD(v) calculation for all blocks in the vicinity. If a good enough motion vector is found, it is used and the macro-block is encoded as a P-block. Otherwise the macro-block will be encoded as an I-block. In each case a 1-bit flag is pushed to the flags list to represent the decision (0 for I-block and 1 for P-block).

The helper functions take the 6 blocks pertaining to the macro-block (4 Y, 1 Cb and 1 Cr) and encode them accordingly. The previous frame and the vector is used for P-blocks, to calculate delta values. Then both frames undergo DCT and quantization steps. Then the encoded (compressed) blocks are stored in the compressed_blocks list in order (4 Y, 1 Cb and 1 Cr) and the decompressed versions are stored in the uncompressed_blocks list in the same order.

Once an entire frame has been processed as such, the information is pushed to the stream as the program prepares for the next frame. First the set of used motion vectors are send over with delta encoding and unary. The each macro-block is pushed: first a 1-bit flag (I-block or P-block), then 4 Y blocks, 1 Cb block and 1 Cr block. Lastly, the uncompressed_blocks list is used to reproduce the block as seen by the decompressor for the next frame.

Note, the quantized blocks are pushed to the stream using delta encoding and static Huffman codes.

On the decompressor side, the motion vectors are saved into a queue-like data structure. If a macro-block is proceeded with a 1-bit flag for P-block the first vector is popped off the queue and used to decode the block.

## Video Examples
The "videos" directory contains example files showing how the video quality degrades after undergoing compression/decompression.
Both files where compressed with the "low" quality mode.

## Features Implemented
- Generates temporally compressed frames (P-frames)
- Implements motion compensation (Generates different motion vectors for different parts of the image)
	- A vector is assigned to a 16x16 macro-block region of the frame
	- Not every macro-block in a given frame is a P-block
- Achieves a compression ratio of at least 12 on medium quality compared to the YCbCr input (on most inputs)
	- news_cif.y4m        compression ratio: 13.2844
	- ice_cif.y4m         compression ratio: 12.6426
	- pamphlet_cif.y4m    compression ratio: 13.3937
	- flower_cif.y4m      compression ratio: 10.7215
- Achieves real-time decompression on video samples with resolution at leas 640x480

## Architecture
The logic of the compressor and decompressor is largely within the uvid_compress.cpp and uvid_decompress.cpp files respectively.

The functions used to compute the DCT on a block by block basis or carry out basic block operations (add, multiply, transpose) are in the "dct" namespace. They are declared in discrete_cosine_transform.hpp and defined in discrete_cosine_transform.cpp. 
Any hard-coded matrices used to compute the DCT (the C matrix, the quantization matrices, etc.) are also defined in the discrete_cosine_transform.hpp file.

The functions in charge of pushing content to the input or output streams are declared in stream.hpp and defined in the stream.cpp file. To improve readability, these functions are within the namespace "stream".

Lastly, the functions that are specific to the video compression logic (compressing P-blocks or handling motion vectors) are declared in helper.hpp and defined in the helper.cpp file. These functions are within the "helper" namespace.

Please note that the majority of the code in the "dct" and "stream" namespaces have been carried over from Assignment 3 with few modifications. Also, each file has been organized to group functions used by the compressor together, followed by the functions used by the decompressor.

### Data Structures
The program largely relies on 8x8 blocks of Y, Cb or Cr values.
This data structure is defined with a using directive in discrete_cosine_transform.hpp
```
using Block8x8 = std::array<std::array<double, 8>, 8>;
```
When pushing the quantized DCT values to the stream the 8x8 block is converted into an array in the optimal ordering for RLE to be applied. This array is also defined in discrete_cosine_transform.hpp
```
using Array64 = std::array<int, 64>;
```
The conversion between Block8x8 and Array64 is delegated to the functions below within discrete_cosine_transform.hpp
```
Array64 block_to_array(const Block8x8& block);
Block8x8 array_to_block(const Array64& array);
```

When searching for motion vectors the program constructs a 16x16 block (often called a macroblock) from 4 8x8 blocks. 
```
Block16x16 create_macroblock(const Block8x8& b1, const Block8x8& b2, const Block8x8& b3, const Block8x8& b4);
```
The 16x16 block object is also defined in the discrete_cosine_transform.hpp file
```
using Block16x16 = std::array<std::array<double, 16>, 16>;
```

The program read the frame and partitions the frame into 8x8blocks of Y, Cb and Cr values stored with the Y_blocks, Cb_blocks and Cr_blocks vector objects.

The compressor and decompressor each consider one macroblock of the frame as a unit consisting of 4 8x8 Y blocks, 1 8x8 Cb block and 1 8x8 Cr block.
Motion vectors are found based on the Y block values and extended to the Cb and Cr blocks.
Then each of the 6 8x8 blocks is compressed and added to the compressed list in the order of:
$$Y_{0,0}\ Y_{0,1}\ Y_{1,0}\ Y_{1,1}\ Cb\ Cr$$


## Bitstream
File header:
- 2-bit quality flag (0=low 1=medium and 2=high)
- 16-bit height
- 16-bit width

For each frame:
- 1-bit flag (0=no frame and 1=frame coming)
- motion vectors $v = (v_x,v_y)$
	- 16-bits number of motion vectors used in the frame
	- 2 x 4 bits the x and u components of the first motion vector
	- Every other motion vector is sent as a delta value in unary where $$ \delta_x^i  = v_x^i - v_x^{i-1}\ and\ \delta_y^i  = v_y^i - v_y^{i-1} $$
- the encoded blocks
	- 1-bit flag (0=I-block and 1=P-block)
	- 4 Y blocks (8x8), 1 Cb block and 1 Cr block

For each 8x8 block, the block is first converted into an array of size 64 in sig-zag order which is ideal for delta compression. The first 2 values (DC and AC) are pushed to the stream with 1-bit flag (1=negative and 0=positive), followed by a 16-bit representation of the absolute value of DC/AC. Every other value is pushed as a delta value using a set of static Huffman codes. They Huffman symbols used, their lengths and encodings can be found in the stream.hpp file. 
### Huffman Codes
- The symbol -100 is a negative-escape symbol for delta values less than -5
	(ie. push symbol -100 then push the absolute value of the delta value in unary)
- The symbol 100 is a positive-escape symbol for delta values greater than 5
	(ie. push symbol 100 then push the delta value in unary)
- The symbols -5 to 5 correspond to a delta value of -5 to 5 respectively
- The symbol 120 corresponds to a run of 8 zeros in the delta values
- The symbol 150 corresponds to an "End of Block"
	(ie. the rest of the delta values for the block are all zero)

## Bibliography
only the lecture slides were used


## Useful commands
Step 1: Create the raw file (y4m -> raw)
```
ffmpeg -i original.y4m -f rawvideo -pixel_format yuv420p - > input.raw
```

Step 2: Compress the raw file (raw -> uvi)
```
./uvid_compress 720 480 <low/medium/high> < input.raw > compresssed.uvi
```

Step 3: Decompress the (uvi -> raw)
```
./uvid_decompress < compressed.uvi > decompressed.raw
```

Step 4: Convert to playable format (raw -> y4m)
```
ffmpeg -f rawvideo -pixel_format yuv420 -framerate 30 -videosize 720x480 -i - -f yuv4mpegpipe playable.y4m < decompressed.raw
```

Step 5: Play y4m file
```
ffplay playable.y4m
```

Note: due to changes to ffmpeg Steps 1 and 4 are no longer valid