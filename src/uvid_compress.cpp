/* uvid_compress.cpp
   CSC 485B/578B - Data Compression - Summer 2023

   Starter code for Assignment 4

   Reads video data from stdin in uncompresed YCbCr (YUV) format 
   (With 4:2:0 subsampling). To produce this format from 
   arbitrary video data in a popular format, use the ffmpeg
   tool and a command like 

     ffmpeg -i videofile.mp4 -f rawvideo -pixel_format yuv420p - 2>/dev/null | ./this_program <width> height>

   Note that since the width/height of each frame is not encoded into the raw
   video stream, those values must be provided to the program as arguments.

   B. Bird - 2023-07-08
*/

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <cstdint>
#include <tuple>
#include <queue>
#include <map>
#include "output_stream.hpp"
#include "stream.hpp"
#include "yuv_stream.hpp"
#include "discrete_cosine_transform.hpp"
#include "helper.hpp"


int main(int argc, char** argv){

    if (argc < 4){
        std::cerr << "Usage: " << argv[0] << " <width> <height> <low/medium/high>" << std::endl;
        return 1;
    }

    // Parse command line arguments
    u16 width = std::stoi(argv[1]);
    u16 height = std::stoi(argv[2]);
    dct::Quality quality = helper::get_quality(argv[3]);
    if(quality == dct::Quality::ERROR){
        std::cerr << "Usage: " << argv[0] << " <width> <height> <low/medium/high>" << std::endl;
        return 1;
    }

    // calculate number of macro blocks expected
    u16 scaled_height = height/2;
    u16 scaled_width = width/2; 
    u16 C_blocks_wide = (scaled_width%8 == 0) ? scaled_width/8 : (scaled_width/8)+1;
    u16 C_blocks_high = (scaled_height%8 == 0) ? scaled_height/8 : (scaled_height/8)+1;
    u16 num_macro_blocks = C_blocks_wide * C_blocks_high;

    YUVStreamReader reader {std::cin, width, height};
    OutputBitStream output_stream {std::cout};

    stream::push_header(output_stream, quality, height, width);

    // To manage previous frame 
    YUVFrame420 previous_frame {width, height};
    u32 frame_number {0};

    while (reader.read_next_frame()){
        // Get the active frame
        YUVFrame420& active_frame = reader.frame();
        output_stream.push_bit(1);

        // Separate Y Cb and Cr channels
        auto Y_matrix = helper::create_2d_vector<unsigned char>(height, width);
        auto Cb_matrix = helper::create_2d_vector<unsigned char>(height/2, width/2);
        auto Cr_matrix = helper::create_2d_vector<unsigned char>(height/2, width/2);
        for (u32 y = 0; y < height; y++)
            for (u32 x = 0; x < width; x++)
                Y_matrix.at(y).at(x) = active_frame.Y(x,y);
        for (u32 y = 0; y < height/2; y++)
            for (u32 x = 0; x < width/2; x++){
                Cb_matrix.at(y).at(x) = active_frame.Cb(x,y);
                Cr_matrix.at(y).at(x) = active_frame.Cr(x,y);
            }
        
        // Partition color channels into 8x8 blocks
        std::vector<Block8x8> Y_blocks, Cb_blocks, Cr_blocks;
        dct::partition_Y_channel(Y_blocks, height, width, Y_matrix);
        dct::partition_C_channel(Cb_blocks, height/2, width/2, Cb_matrix);
        dct::partition_C_channel(Cr_blocks, height/2, width/2, Cr_matrix);

        // To manage previous frame
        std::list<Block8x8> uncompressed_blocks;

        // To manage active frame
        std::list<bool> flags;
        std::list<Block8x8> compressed_blocks;
        std::list<std::pair<int, int>> motion_vectors;
        double num_bad_motion_vectors {0};
        for(u32 macro_idx = 0; macro_idx < num_macro_blocks; macro_idx++){
            // create 16x16 Y-block
            u32 Y_idx = 4 * macro_idx;
            Block16x16 macroblock = dct::create_macroblock(Y_blocks.at(Y_idx), Y_blocks.at(Y_idx+1), Y_blocks.at(Y_idx+2), Y_blocks.at(Y_idx+3));

            // Look for motion vector (assume non found)
            std::pair<int, int> vector {0, 0};
            bool good_motion_vector = helper::find_motion_vector(macroblock, previous_frame, macro_idx, vector);
            if(!good_motion_vector)
                num_bad_motion_vectors++;
            if (frame_number && good_motion_vector){
                flags.push_back(1);
                motion_vectors.push_back(vector);
                helper::compress_P_block(compressed_blocks, uncompressed_blocks, macro_idx, Y_blocks, Cb_blocks, Cr_blocks, quality, previous_frame, vector);

            }else{
                flags.push_back(0);
                helper::compress_I_block(compressed_blocks, uncompressed_blocks, macro_idx, Y_blocks, Cb_blocks, Cr_blocks, quality);
            }
        }
        // Begin to push the frame
        helper::push_motion_vectors(motion_vectors, output_stream);
        // send compressed blocks
        helper::push_compressed_blocks(flags, compressed_blocks, output_stream);
        // reconstruct prev frame
        previous_frame = helper::reconstruct_prev_frame(uncompressed_blocks, num_macro_blocks, height, width);

        // Send an I-frame every 120 frames or if too many bad motion vectors
        if(frame_number > 175 && (num_bad_motion_vectors/num_macro_blocks) >= 0.35)
            frame_number = 0;
        else    
            frame_number++;
    }

    output_stream.push_bit(0); //Flag to indicate end of data
    output_stream.flush_to_byte();
    return 0;
}