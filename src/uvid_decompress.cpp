/* uvid_decompress.cpp
   CSC 485B/578B - Data Compression - Summer 2023

   Starter code for Assignment 4
   
   This placeholder code reads the (basically uncompressed) data produced by
   the uvid_compress starter code and outputs it in the uncompressed 
   YCbCr (YUV) format used for the sample video input files. To play the 
   the decompressed data stream directly, you can pipe the output of this
   program to the ffplay program, with a command like 

     ffplay -f rawvideo -pixel_format yuv420p -framerate 30 -video_size 352x288 - 2>/dev/null
   (where the resolution is explicitly given as an argument to ffplay).

   B. Bird - 2023-07-08
*/

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <queue>
#include <list>
#include <cstdint>
#include <tuple>
#include "input_stream.hpp"
#include "yuv_stream.hpp"
#include "discrete_cosine_transform.hpp"
#include "stream.hpp"
#include "helper.hpp"


int main(int argc, char** argv){

    //Note: This program must not take any command line arguments. (Anything
    //      it needs to know about the data must be encoded into the bitstream)
    
    InputBitStream input_stream {std::cin};

    dct::Quality quality;
    u16 height, width;
    stream::read_header(input_stream, quality, height, width);

    // calculate number of macro blocks expected
    u16 scaled_height = height/2;
    u16 scaled_width = width/2; 
    u16 C_blocks_wide = (scaled_width%8 == 0) ? scaled_width/8 : (scaled_width/8)+1;
    u16 C_blocks_high = (scaled_height%8 == 0) ? scaled_height/8 : (scaled_height/8)+1;
    u16 num_macro_blocks = C_blocks_wide * C_blocks_high;

    YUVStreamWriter writer {std::cout, width, height};

    // To store uncompressed blocks 
    YUVFrame420 previous_frame {width, height};

    while (input_stream.read_bit()){

        // Read the motion vectors
        std::list<std::pair<int, int>> motion_vectors;
        helper::read_motion_vectors(motion_vectors, input_stream);
        
        // read blocks for each color channel in row major order
        std::vector<Block8x8> Y_blocks, Cb_blocks, Cr_blocks;

        for(u32 macro_idx = 0; macro_idx < num_macro_blocks; macro_idx++){
            bool block_type = input_stream.read_bit();
            if(block_type == 0){
                // I-block
                helper::decompress_I_block(Y_blocks, Cb_blocks, Cr_blocks, quality, input_stream);
            }else{
                //P-block
                helper::decompress_P_block(Y_blocks, Cb_blocks, Cr_blocks, quality, input_stream, macro_idx, motion_vectors.front(), previous_frame);
                motion_vectors.pop_front();
            }
        }

        auto Y_matrix = helper::create_2d_vector<unsigned char>(height,width);
        auto Cb_matrix = helper::create_2d_vector<unsigned char>(height/2 ,width/2);
        auto Cr_matrix = helper::create_2d_vector<unsigned char>(height/2 ,width/2);
        dct::undo_partition_Y_channel(Y_blocks, height, width, Y_matrix);
        dct::undo_partition_C_channel(Cb_blocks, height/2 ,width/2, Cb_matrix);
        dct::undo_partition_C_channel(Cr_blocks, height/2 ,width/2, Cr_matrix);

        // Create and write into frame
        YUVFrame420& active_frame = writer.frame();
        writer.write_frame();
        for (u32 y = 0; y < height; y++)
            for (u32 x = 0; x < width; x++)
                active_frame.Y(x,y) = Y_matrix.at(y).at(x);
        for (u32 y = 0; y < height/2; y++)
            for (u32 x = 0; x < width/2; x++){
                active_frame.Cb(x,y) = Cb_matrix.at(y).at(x);
                active_frame.Cr(x,y) = Cr_matrix.at(y).at(x);
            }

        previous_frame = active_frame;
    }

    return 0;
}