#include "discrete_cosine_transform.hpp"

#include <iostream>

namespace dct{

    /* ----- Block Operations ----- */
    // returns the c_matrix for n = 8
    Block8x8 create_c_matrix(){
        Block8x8 c_matrix {};
        double n = 8;
        double root_1_over_n = std::sqrt(1/n);
        double root_2_over_n = std::sqrt(2/n);

        for(u32 r = 0; r < n; r++){
            for(u32 c = 0; c < n; c++){
                if(r == 0){
                    c_matrix.at(r).at(c) = root_1_over_n;
                }else{
                    double angle = ((2*c+1) * r * M_PI)/(2*n);
                    c_matrix.at(r).at(c) = root_2_over_n * std::cos(angle);
                }
            }
        }
        return c_matrix;
    }
    
    // prints array of 64 to standard out
    void print_array(const Array64& array){
        std::cout << "-------" << std::endl;
        for(u32 idx = 0; idx < 64; idx++){
            if(idx > 0 && idx%8 == 0)
                std::cout << std::endl;
            std::cout << array.at(idx) << " ";
        }
        std::cout << std::endl;
    }
    // prints 8x8 block to standard out
    void print_block(const Block8x8& block){
        for(u32 r = 0; r < 8; r++){
            std::cout << "{ ";
            for(u32 c = 0; c < 8; c++){
                std::cout << block.at(r).at(c) << " ";
            }
            std::cout << "}" << std::endl;
        }
    }

    // prints all blocks in a block vector
    void print_blocks(const std::vector<Block8x8>& blocks){
        for(const Block8x8& b : blocks){
            std::cout << "-------" << std::endl;
            print_block(b);
        }
    }

    // returns the 8x8 block result of multiplying blockA by blockB
    Block8x8 multiply_block(const Block8x8& blockA, const Block8x8& blockB){
        Block8x8 result;
        for(u32 r = 0; r < 8; r++){
            for(u32 c = 0; c < 8; c++){
                double sum = 0;
                for(u32 idx = 0; idx < 8; idx++){
                    sum += (blockA.at(r).at(idx) * blockB.at(idx).at(c));
                }
                // std::cout << " = " << sum << std::endl;
                result.at(r).at(c) = sum;
            }
        }
        return result;
    }

    // returns the 8x8 block result of transposing the block
    Block8x8 transpose_block(const Block8x8& block){
        Block8x8 transpose;
        for (u32 r = 0; r < 8; r++)
            for(u32 c = 0; c < 8; c++)
                transpose.at(c).at(r) = block.at(r).at(c);
        return transpose;
    }

    // Returns the 8x8 block of delta values of block1 - block2
    Block8x8 get_delta_block(const Block8x8& block1, const Block8x8& block2){
        Block8x8 delta;
        for (u32 r = 0; r < 8; r++)
            for(u32 c = 0; c < 8; c++)
                delta.at(r).at(c) = block1.at(r).at(c) - block2.at(r).at(c);
        return delta;
    }

    // Returns the 8x8 block of adding the delta values to the block
    Block8x8 add_delta_block(const Block8x8& block, const Block8x8& delta){
        Block8x8 result;
        for (u32 r = 0; r < 8; r++)
            for(u32 c = 0; c < 8; c++)
                result.at(r).at(c) = block.at(r).at(c) + delta.at(r).at(c);
        return result;
    }

    Block16x16 create_macroblock(const Block8x8& b1, const Block8x8& b2, const Block8x8& b3, const Block8x8& b4){
        Block16x16 macroblock;
        for(u32 b1_r = 0; b1_r < 8; b1_r++){
            for(u32 b1_c = 0; b1_c < 8; b1_c++){
                macroblock.at(b1_r).at(b1_c) = b1.at(b1_r).at(b1_c);
            }
        }
        for(u32 b2_r = 0; b2_r < 8; b2_r++){
            for(u32 b2_c = 0; b2_c < 8; b2_c++){
                macroblock.at(b2_r).at(b2_c+8) = b2.at(b2_r).at(b2_c);
            }
        }
        for(u32 b3_r = 0; b3_r < 8; b3_r++){
            for(u32 b3_c = 0; b3_c < 8; b3_c++){
                macroblock.at(b3_r+8).at(b3_c) = b3.at(b3_r).at(b3_c);
            }
        }
        for(u32 b4_r = 0; b4_r < 8; b4_r++){
            for(u32 b4_c = 0; b4_c < 8; b4_c++){
                macroblock.at(b4_r+8).at(b4_c+8) = b4.at(b4_r).at(b4_c);
            }
        }
        return macroblock;
    }

    void get_prev_blocks(u32 macro_idx, YUVFrame420& prev_frame, const std::pair<u32, u32>& vector, std::vector<Block8x8>& prev_blocks){
        
        u32 macroblocks_wide = prev_frame.get_Width() / 16;

        // (0,0) coordinate of active block in the frame
        u32 B_x = (macro_idx % macroblocks_wide) * 16;
        u32 B_y = (macro_idx / macroblocks_wide) * 16;

        // (0,0) coordinate of compare block
        u32 P_x = B_x + vector.first;
        u32 P_y = B_y + vector.second;
       
        Block8x8 block;
        // Push back Y blocks
        // Top-left Y 8x8 sub-block
        for(u32 r = 0; r < 8; r++)
            for(u32 c = 0; c < 8; c++)
                block.at(r).at(c) = prev_frame.Y(P_x+c, P_y+r);
        prev_blocks.push_back(block);
        
        // Top-right 8x8 sub-block
        for(u32 r = 0; r < 8; r++)
            for(u32 c = 0; c < 8; c++)
                block.at(r).at(c) = prev_frame.Y(P_x+8+c, P_y+r);
        prev_blocks.push_back(block);

        // Bottom-left 8x8 sub-block
        for(u32 r = 0; r < 8; r++)
            for(u32 c = 0; c < 8; c++)
                block.at(r).at(c) = prev_frame.Y(P_x+c, P_y+8+r);
        prev_blocks.push_back(block);

        // Bottom-right 8x8 sub-block
        for(u32 r = 0; r < 8; r++)
            for(u32 c = 0; c < 8; c++)
                block.at(r).at(c) = prev_frame.Y(P_x+8+c, P_y+8+r);
        prev_blocks.push_back(block);

        //Push back Cb block
        for(u32 r = 0; r < 8; r++)
            for(u32 c = 0; c < 8; c++)
                block.at(r).at(c) = prev_frame.Cb((P_x+c)/2,(P_y+r)/2);
        prev_blocks.push_back(block);

        // Push back Cr block
        for(u32 r = 0; r < 8; r++)
            for(u32 c = 0; c < 8; c++)
                block.at(r).at(c) = prev_frame.Cr((P_x+c)/2,(P_y+r)/2);
        prev_blocks.push_back(block);
    }

    /* ----- Compressor Functions ----- */

    // given a color channel partitions into 8x8 blocks and adds blocks to vector in row major order
    void partition_C_channel(std::vector<Block8x8>& blocks, u32 height, u32 width, const std::vector<std::vector<unsigned char>>& channel){
        for(u32 r = 0; r < height; r+=8){
            for(u32 c = 0; c < width; c+=8){
                Block8x8 current_block;
                // index into 8x8 sub block
                for(u32 sub_r = 0; sub_r < 8; sub_r++){
                    for(u32 sub_c = 0; sub_c < 8; sub_c++){
                        // copy element 
                        if((r+sub_r) >= height && (c+sub_c) < width){
                            current_block.at(sub_r).at(sub_c) = double(channel.at(height-1).at(c+sub_c));
                        }else if((c+sub_c) >= width && (r+sub_r) < height){
                            current_block.at(sub_r).at(sub_c) = double(channel.at(r+sub_r).at(width-1));
                        }else if((r+sub_r) >= height && (c+sub_c) >= width){
                            current_block.at(sub_r).at(sub_c) = double(channel.at(height-1).at(width-1));
                        }else{
                            current_block.at(sub_r).at(sub_c) = double(channel.at(r+sub_r).at(c+sub_c));
                        }
                    }
                }
                blocks.push_back(current_block);
            }
        }
    }

    // given a Y channel partitions into 8x8 blocks and adds blocks to vector in macroblock row major order
    void partition_Y_channel(std::vector<Block8x8>& blocks, u32 height, u32 width, const std::vector<std::vector<unsigned char>>& channel){
        // break up into 16x 16 blocks
        for(u32 r = 0; r < height; r+=16){
            for(u32 c = 0; c < width; c+=16){
                // break up each 16x16 into 8x8
                for(u32 sub_r = 0; sub_r < 16; sub_r+=8){
                    for(u32 sub_c = 0; sub_c < 16; sub_c+=8){
                        // create the block
                        Block8x8 current_block;
                        // index into 8x8 sub block
                        for(u32 block_r = 0; block_r < 8; block_r++){
                            for(u32 block_c = 0; block_c < 8; block_c++){
                                u32 r_offset = r+sub_r;
                                u32 c_offset = c+sub_c;
                                // copy element 
                                if((r_offset+block_r) >= height && (c_offset+block_c) < width)
                                    current_block.at(block_r).at(block_c) = double(channel.at(height-1).at(c_offset+block_c));
                                else if((c_offset+block_c) >= width && (r_offset+block_r) < height)
                                    current_block.at(block_r).at(block_c) = double(channel.at(r_offset+block_r).at(width-1));
                                else if((r_offset+block_r) >= height && (c_offset+block_c) >= width)
                                    current_block.at(block_r).at(block_c) = double(channel.at(height-1).at(width-1));
                                else
                                    current_block.at(block_r).at(block_c) = double(channel.at(r_offset+block_r).at(c_offset+block_c));
                            }
                        } 
                        blocks.push_back(current_block);
                    }
                }
                
            }
        }
    }
    
    // returns the dct of block A by computing [C][A][C]_transpose
    Block8x8 get_dct(const Block8x8& block){
        Block8x8 result = multiply_block(c_matrix, block);
        return multiply_block(result, c_matrix_transpose);
    }

    double get_multiplier(Quality quality, bool is_luminance, bool is_P_block){
        if(is_luminance && is_P_block){
            if(quality == low)
                return 6;
            else if(quality == medium)
                return 5;
            else
                return 2;
        }else if(is_luminance && !is_P_block){
            if(quality == low)
                return 4;
            else if(quality == medium)
                return 3;
            else
                return 1;
        }else if(!is_luminance && is_P_block){
            if(quality == low)
                return 10;
            else if(quality == medium)
                return 8;
            else
                return 3;
        }else{
            if(quality == low)
                return 6;
            else if(quality == medium)
                return 5;
            else
                return 2;
        }
    }

    // returns the quantized block calculated using the provided quantization matrix at the provided quality 
    Block8x8 quantize_block(const Block8x8& block, Quality quality, bool is_luminance, bool is_P_block){
        double multiplier = get_multiplier(quality, is_luminance, is_P_block);

        Block8x8 result;
        if(is_luminance){
            for(u32 r = 0; r < 8; r++)
                for(u32 c = 0; c < 8; c++)
                    result.at(r).at(c) = std::round(block.at(r).at(c) / (multiplier * luminance.at(r).at(c)) );
        }else{
            for(u32 r = 0; r < 8; r++)
                for(u32 c = 0; c < 8; c++)
                    result.at(r).at(c) = std::round(block.at(r).at(c) / (multiplier * chrominance.at(r).at(c)) );
        }
        return result;
    }

    Direction get_direction(u32 r, u32 c, Direction curr){
        u32 first = 0;
        u32 last = 7;

        if((r == first || r == last) && c%2 == 0){
            return right;
        }else if((c == first || c == last) && r%2 == 1){
            return down;
        }else if((r == first && c%2 == 1) || (c == last && r%2 == 0)){
            return down_left;
        }else if((c == first && r%2 == 0) || (r == last && c%2 == 1)){
            return up_right;
        }
        return curr;
    }

    // converts an 8x8 block to an array of 64 elements in "ideal" order
    Array64 block_to_array(const Block8x8& block){
        Direction dir = right;
        Array64 result;   

        u32 r = 0, c = 0, count = 0;
        while(r < 8 && c < 8){
            result.at(count++) = block.at(r).at(c);
            dir = get_direction(r, c, dir);

            if(dir == right){
                c++;
            }else if(dir == down){
                r++;
            }else if(dir == down_left){
                r++; c--;
            }else{
                r--; c++;
            }
        }
        return result;
    }

    /* ----- Decompressor Functions ----- */

    // converts an array of 64 elements in "ideal" order to an 8x8 block
    Block8x8 array_to_block(const Array64& array){
        Direction dir = right;
        Block8x8 result;   

        u32 r = 0, c = 0, count = 0;
        while(r < 8 && c < 8){
            result.at(r).at(c) = array.at(count++);
            dir = get_direction(r, c, dir);

            if(dir == right){
                c++;
            }else if(dir == down){
                r++;
            }else if(dir == down_left){
                r++; c--;
            }else{
                r--; c++;
            }
        }
        return result;
    }

    // returns the unquantized block calculated using the provided quantization matrix at the provided quality 
    Block8x8 unquantize_block(const Block8x8& block, Quality quality, bool is_luminance, bool is_P_block){
        double multiplier = get_multiplier(quality, is_luminance, is_P_block);

        Block8x8 result;
        if (is_luminance){
            for(u32 r = 0; r < 8; r++)
                for(u32 c = 0; c < 8; c++)
                    result.at(r).at(c) = block.at(r).at(c) * (multiplier * luminance.at(r).at(c));
        }else{
            for(u32 r = 0; r < 8; r++)
                for(u32 c = 0; c < 8; c++)
                    result.at(r).at(c) = block.at(r).at(c) * (multiplier * chrominance.at(r).at(c));
        }
        return result;
    }

    // returns the dct of block A by computing [C][A][C]_transpose
    Block8x8 get_inverse_dct(const Block8x8& block){
        Block8x8 result = multiply_block(c_matrix_transpose, block);
        return multiply_block(result, c_matrix);
    }

    // given a vector of blocks in row major order color reconstructs the channel matrix
    void undo_partition_C_channel(const std::vector<Block8x8>& blocks, u32 height, u32 width, std::vector<std::vector<unsigned char>>& channel){
        u32 idx = 0;
        for(u32 r = 0; r < height; r+=8){
            for(u32 c = 0; c < width; c+=8){
                const Block8x8& current_block = blocks.at(idx++);
                // index into 8x8 sub block
                for(u32 sub_r = 0; sub_r < 8; sub_r++){
                    for(u32 sub_c = 0; sub_c < 8; sub_c++){
                        // copy element 
                        if( (r+sub_r) < height && (c+sub_c) < width ){
                            channel.at(r+sub_r).at(c+sub_c) = round_and_clamp_to_char(current_block.at(sub_r).at(sub_c));
                        }
                    }
                }
            }
        }
    }

    // given a vector of blocks in row major order color reconstructs the channel matrix
    void undo_partition_Y_channel(const std::vector<Block8x8>& blocks, u32 height, u32 width, std::vector<std::vector<unsigned char>>& channel){
        u32 idx = 0;
        for(u32 r = 0; r < height; r+=16){
            for(u32 c = 0; c < width; c+=16){
                for(u32 sub_r = 0; sub_r < 16; sub_r+=8){
                    for(u32 sub_c = 0; sub_c < 16; sub_c+=8){
                        const Block8x8& current_block = blocks.at(idx++);
                        u32 r_offset = r + sub_r;
                        u32 c_offset = c + sub_c;
                        // index into 8x8 sub block
                        for(u32 sub_r = 0; sub_r < 8; sub_r++){
                            for(u32 sub_c = 0; sub_c < 8; sub_c++){
                                // copy element 
                                if( (r_offset+sub_r) < height && (c_offset+sub_c) < width ){
                                    channel.at(r_offset+sub_r).at(c_offset+sub_c) = round_and_clamp_to_char(current_block.at(sub_r).at(sub_c));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}