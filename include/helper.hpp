#include <queue>
#include <vector>
#include <string>
#include <list>
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include "discrete_cosine_transform.hpp"
#include "yuv_stream.hpp"
#include "output_stream.hpp"
#include "stream.hpp"

namespace helper{
    
    /* ----- Helper Functions - written by Bill ----- */
    //Convenience function to wrap around the nasty notation for 2d vectors
    //Convenience function to wrap around the nasty notation for 2d vectors
    template<typename T>
    std::vector<std::vector<T> > create_2d_vector(unsigned int outer, unsigned int inner){
        std::vector<std::vector<T> > V {outer, std::vector<T>(inner,T() )};
        return V;
    }

    /* ----- My Helpers ----- */

    // Returns the quality enum corresponding to the input string
    // If invalid returns Quality::ERROR
    dct::Quality get_quality(std::string input_quality){
        if(input_quality == "low")
            return dct::Quality::low;
        else if(input_quality == "medium")
            return dct::Quality::medium;
        else if(input_quality == "high")
            return dct::Quality::high;
        else
            return dct::Quality::ERROR;
    }

    /* ----- Compressor Code ----- */

    bool find_motion_vector(const Block16x16& block, YUVFrame420& prev_frame, u32 macro_idx, std::pair<int,int>& vector, int radius = 8){
        u32 macroblocks_wide = prev_frame.get_Width() / 16;

        // (0,0) coordinate of active block in the frame
        int B_x = (macro_idx % macroblocks_wide) * 16;
        int B_y = (macro_idx / macroblocks_wide) * 16;

        // Search region boundaries (radius of 8)
        int v_x_min = (B_x-radius < 0)? 0 : B_x-radius;
        int v_x_max = (B_x+radius < prev_frame.get_Width()) ? B_x+radius : prev_frame.get_Width();
        int v_y_min = (B_y-radius < 0)? 0 : B_y-radius;
        int v_y_max = (B_y+radius < prev_frame.get_Height()) ? B_y+radius : prev_frame.get_Height();

        double min_avg_difference {INT32_MAX};
        // Look for motion vectors
        for(int v_x = v_x_min; v_x < v_x_max; v_x++){
            for(int v_y = v_y_min; v_y < v_y_max; v_y++){
                // calculate the AAD(v)
                double avg_difference {};

                for(u32 r = 0; r < 16 && r+v_y < prev_frame.get_Height(); r++){
                    for(u32 c = 0; c < 16 && c+v_x < prev_frame.get_Width(); c++){
                        int value = int(block.at(r).at(c)) - int(prev_frame.Y(c+v_x, r+v_y));
                        avg_difference += std::abs(value);
                    }
                }
                avg_difference/=256;
                // update the minimum value
                if(avg_difference < min_avg_difference){
                    min_avg_difference = avg_difference;
                    vector.first = v_x - B_x;
                    vector.second= v_y - B_y;
                }
            }
        }
        // if a good enough vector is found return
        if(min_avg_difference <= 50){
            return true;
        }
        return false;
    }

    void compress_I_block(std::list<Block8x8>& compressed_blocks, std::list<Block8x8>& uncompressed_blocks, u32 C_idx, 
    const std::vector<Block8x8>& Y_blocks, const std::vector<Block8x8>& Cb_blocks, const std::vector<Block8x8>& Cr_blocks, dct::Quality quality){
        u32 Y_idx = 4 * C_idx;
        for(u32 count = 0; count < 4; count++){
            // Take the DCT and quantize
            Block8x8 quantized_block = dct::quantize_block(dct::get_dct(Y_blocks.at(Y_idx+count)), quality, true, false);
            // Push in array format
            compressed_blocks.push_back(quantized_block);
            // Unquantize and take the inverse DCT
            uncompressed_blocks.push_back(dct::get_inverse_dct(dct::unquantize_block(quantized_block, quality, true, false)));
        }

        Block8x8 quantized_Cb_block = dct::quantize_block(dct::get_dct(Cb_blocks.at(C_idx)), quality, false, false);
        compressed_blocks.push_back(quantized_Cb_block);
        uncompressed_blocks.push_back(dct::get_inverse_dct(dct::unquantize_block(quantized_Cb_block, quality, false, false)));

        Block8x8 quantized_Cr_block = dct::quantize_block(dct::get_dct(Cr_blocks.at(C_idx)), quality, false, false);
        compressed_blocks.push_back(quantized_Cr_block);
        uncompressed_blocks.push_back(dct::get_inverse_dct(dct::unquantize_block(quantized_Cr_block, quality, false, false)));
    }

    void compress_P_block(std::list<Block8x8>& compressed_blocks, std::list<Block8x8>& uncompressed_blocks, u32 macro_idx, 
    const std::vector<Block8x8>& Y_blocks, const std::vector<Block8x8>& Cb_blocks, const std::vector<Block8x8>& Cr_blocks, dct::Quality quality,
    YUVFrame420& prev_frame, const std::pair<int, int>& vector){

        std::vector<Block8x8> prev_blocks;
        dct::get_prev_blocks(macro_idx, prev_frame, vector, prev_blocks);
        // std::cerr<< "prev_blocks " << prev_blocks.size() << std::endl;
        u32 Y_idx = 4 * macro_idx;
        for(u32 count = 0; count < 4; count++){
            //Get the delta values 
            Block8x8 delta_block = dct::get_delta_block(Y_blocks.at(Y_idx+count), prev_blocks.at(count));
            // Take the DCT and quantize the delta values
            Block8x8 quantized_block = dct::quantize_block(dct::get_dct(delta_block), quality, true, true);
            // Push in array format
            compressed_blocks.push_back(quantized_block);
            // Unquantize and take the inverse DCT of the delta values 
            Block8x8 uncompressed_delta = dct::get_inverse_dct(dct::unquantize_block(quantized_block, quality, true, true));
            // Unquantize and take the inverse DCT
            uncompressed_blocks.push_back(dct::add_delta_block(prev_blocks.at(count), uncompressed_delta));
        }

        Block8x8 delta_block = dct::get_delta_block(Cb_blocks.at(macro_idx), prev_blocks.at(4));
        Block8x8 quantized_block = dct::quantize_block(dct::get_dct(delta_block), quality, false, true);
        compressed_blocks.push_back(quantized_block);
        Block8x8 uncompressed_delta = dct::get_inverse_dct(dct::unquantize_block(quantized_block, quality, false, true));
        uncompressed_blocks.push_back(dct::add_delta_block(prev_blocks.at(4), uncompressed_delta));

        delta_block = dct::get_delta_block(Cr_blocks.at(macro_idx), prev_blocks.at(5));
        quantized_block = dct::quantize_block(dct::get_dct(delta_block), quality, false, true);
        compressed_blocks.push_back(quantized_block);
        uncompressed_delta = dct::get_inverse_dct(dct::unquantize_block(quantized_block, quality, false, true));
        uncompressed_blocks.push_back(dct::add_delta_block(prev_blocks.at(5), uncompressed_delta));
    }

    void push_motion_vectors(std::list<std::pair<int, int>>& motion_vectors, OutputBitStream& output_stream){
        // Push number of motion vectors
        output_stream.push_u16(motion_vectors.size());

        // No motion vectors to push so return 
        if(motion_vectors.size() == 0)
            return;

        // Push the first motion vector as u16
        std::pair<int, int>& first_vector = motion_vectors.front();
        stream::push_value_n(output_stream, first_vector.first, 4);
        stream::push_value_n(output_stream, first_vector.second, 4);
        std::pair<int, int> prev_vector = first_vector;
        motion_vectors.pop_front();

        // Push the rest of the motion vectors as delta values 
        while(!motion_vectors.empty()){
            std::pair<int, int>& curr_vector = motion_vectors.front();
            stream::push_delta_value(output_stream, curr_vector.first-prev_vector.first);
            stream::push_delta_value(output_stream, curr_vector.second-prev_vector.second);
            prev_vector = curr_vector;
            motion_vectors.pop_front();
        }
    }

    void push_compressed_blocks(const std::list<bool>& flags, std::list<Block8x8>& compressed_blocks, OutputBitStream& output_stream){
        for(bool block_type : flags){
            // Push block-type bit (0=I-block and 1=P-block)
            output_stream.push_bit(block_type);
            // Push the macro block (in Y Cb Cr order)
            for(u32 count = 0; count < 6; count++){
                stream::push_quantized_array_delta(output_stream, dct::block_to_array(compressed_blocks.front()));
                compressed_blocks.pop_front();
            }
        }
    }

    YUVFrame420 reconstruct_prev_frame(std::list<Block8x8>& compressed_blocks, u32 num_macro_blocks, u32 height, u32 width){  
        std::vector<Block8x8> Y_blocks, Cb_blocks, Cr_blocks;
        for(u32 macro_count = 0; macro_count < num_macro_blocks; macro_count++){
            for(u32 y_count = 0; y_count < 4; y_count++){
                Y_blocks.push_back(compressed_blocks.front());
                compressed_blocks.pop_front();
            }
            Cb_blocks.push_back(compressed_blocks.front());
            compressed_blocks.pop_front();

            Cr_blocks.push_back(compressed_blocks.front());
            compressed_blocks.pop_front();
        }

        auto Y_matrix = helper::create_2d_vector<unsigned char>(height,width);
        auto Cb_matrix = helper::create_2d_vector<unsigned char>(height/2 ,width/2);
        auto Cr_matrix = helper::create_2d_vector<unsigned char>(height/2 ,width/2);
        dct::undo_partition_Y_channel(Y_blocks, height, width, Y_matrix);
        dct::undo_partition_C_channel(Cb_blocks, height/2 ,width/2, Cb_matrix);
        dct::undo_partition_C_channel(Cr_blocks, height/2 ,width/2, Cr_matrix);

        YUVFrame420 previous_frame {width, height};
        for (u32 y = 0; y < height; y++)
            for (u32 x = 0; x < width; x++)
                previous_frame.Y(x,y) = Y_matrix.at(y).at(x);
        for (u32 y = 0; y < height/2; y++)
            for (u32 x = 0; x < width/2; x++){
                previous_frame.Cb(x,y) = Cb_matrix.at(y).at(x);
                previous_frame.Cr(x,y) = Cr_matrix.at(y).at(x);
            }
        return previous_frame;

    }

    /*----- Decompressor Code -----*/
    void decompress_I_block(std::vector<Block8x8>& Y_blocks, std::vector<Block8x8>& Cb_blocks, std::vector<Block8x8>& Cr_blocks, dct::Quality quality, InputBitStream& input_stream){
        for(u32 count = 0; count < 4; count++){
            Block8x8 Y_block = dct::array_to_block(stream::read_quantized_array_delta(input_stream));
            // Unquantize and take the inverse dct
            Y_blocks.push_back(dct::get_inverse_dct(dct::unquantize_block(Y_block, quality, true, false)));
        }
        Block8x8 Cb_block = dct::array_to_block(stream::read_quantized_array_delta(input_stream));
        Cb_blocks.push_back(dct::get_inverse_dct(dct::unquantize_block(Cb_block, quality, false, false)));

        Block8x8 Cr_block = dct::array_to_block(stream::read_quantized_array_delta(input_stream));
        Cr_blocks.push_back(dct::get_inverse_dct(dct::unquantize_block(Cr_block, quality, false, false)));
    }

    void decompress_P_block(std::vector<Block8x8>& Y_blocks, std::vector<Block8x8>& Cb_blocks, std::vector<Block8x8>& Cr_blocks, dct::Quality quality, InputBitStream& input_stream, 
    u32 macro_idx, std::pair<int, int>& motion_vector, YUVFrame420& prev_frame){

        std::vector<Block8x8> prev_blocks;
        dct::get_prev_blocks(macro_idx, prev_frame, motion_vector, prev_blocks);

        for(u32 count = 0; count < 4; count++){
            // Create a block of delta values
            Block8x8 delta_block = dct::array_to_block(stream::read_quantized_array_delta(input_stream));
            // Unquantize and take the inverse dct
            delta_block = dct::get_inverse_dct(dct::unquantize_block(delta_block, quality, true, true));
            // Add delta_values to previous block
            Y_blocks.push_back(dct::add_delta_block(prev_blocks.at(count), delta_block));
        }

        Block8x8 delta_block = dct::array_to_block(stream::read_quantized_array_delta(input_stream));
        delta_block = dct::get_inverse_dct(dct::unquantize_block(delta_block, quality, false, true));
        Cb_blocks.push_back(dct::add_delta_block(prev_blocks.at(4), delta_block));

        delta_block = dct::array_to_block(stream::read_quantized_array_delta(input_stream));
        delta_block = dct::get_inverse_dct(dct::unquantize_block(delta_block, quality, false, true));
        Cr_blocks.push_back(dct::add_delta_block(prev_blocks.at(5), delta_block));
    }

    void read_motion_vectors(std::list<std::pair<int, int>>& motion_vectors, InputBitStream& input_stream){
        // push number of motiocln vectors
        int num_vectors = input_stream.read_u16();

        std::pair<int, int> first_vector;
        if (num_vectors > 0){
            first_vector.first = stream::read_value_n(input_stream, 4);
            first_vector.second = stream::read_value_n(input_stream, 4);
            motion_vectors.push_back(first_vector);
            num_vectors--;
        }

        std::pair<int, int> prev_vector = first_vector;
        while(num_vectors > 0){
            std::pair<int, int> curr_vector;
            curr_vector.first = stream::read_delta_value(input_stream) + prev_vector.first;
            curr_vector.second = stream::read_delta_value(input_stream) + prev_vector.second;
            motion_vectors.push_back(curr_vector);
            prev_vector = curr_vector;
            num_vectors--;
        }
    }
}
