#include <vector>
#include <array>
#include "stream.hpp"

namespace stream{

    std::map<int,int> delta_frequency {};
    std::map<int,int> RLE_frequency {}; 

    std::map<int, u32> symbol_length {
        {-100, 9},  // negative escape symbol
        {-5, 9},
        {-4, 8},
        {-3, 7},
        {-2, 5},
        {-1, 2},
        {0, 2},
        {1, 3},
        {2, 6},
        {3, 7},
        {4, 8},
        {5, 9},
        {100, 9},   // positive escape symbol
        {120, 5},   // 8 zeros
        {150, 2}    // EOB - the rest of the block is zeros
    };

    std::map<int, u32> symbol_encoding {
        {-100, 500},  // negative escape symbol
        {-5, 501},
        {-4, 248},
        {-3, 122},
        {-2, 28},
        {-1, 1},
        {0, 0},
        {1, 6},
        {2, 60},
        {3, 123},
        {4, 249},
        {5, 502},
        {100, 503},   // positive escape symbol
        {120, 29},   // 8 zeros
        {150, 2}    // EOB - the rest of the block is zeros
    };

    std::map<u32, int> encoding_symbol {
        {500, -100},  // negative escape symbol
        {501, -5},
        {248, -4},
        {122, -3},
        {28, -2},
        {1, -1},
        {0, 0},
        {6, 1},
        {60, 2},
        {123, 3},
        {249, 4},
        {502, 5},
        {503, 100},   // positive escape symbol
        {29, 120},   // 8 zeros
        {2, 150}    // EOB - the rest of the block is zeros
    };

    void print_histograms(){
        std::cerr << "delta histogram" << std::endl;
        int sum_delta {0};
        int neg_x {};
        int pos_x {};
        for (const auto& [value, frequency] : delta_frequency){
            if(value < -5)
                neg_x += frequency;
            else if(value > 5)
                pos_x += frequency;
            sum_delta+=frequency;
        }
        for(int value = -5; value <= 5; value++){
            std::cerr << value << " " << delta_frequency[value] << std::endl;
        }
        std::cerr << "neg_x" << neg_x << std::endl;
        std::cerr << "neg_y" << pos_x << std::endl;
        std::cerr << "sum delta " << sum_delta << std::endl;
        std::cerr << "------------------------------------------------" << std::endl;
        int sum_RLE {0};
        std::cerr << "RLE histogram" << std::endl;
        for (const auto& [value, frequency] : RLE_frequency){
            std::cerr << "RLE: " << value+1 << " with frequency " << frequency << std::endl;
            sum_RLE+=frequency;
        }
        std::cerr << "sum RLE " << sum_RLE << std::endl;
    }

    void push_header(OutputBitStream& stream, dct::Quality q, u16 height, u16 width){
        stream.push_bits(q, 2);
        stream.push_u16(height);
        stream.push_u16(width);
    }

    void push_value(OutputBitStream& stream, int num){
        if(num < 0){
            // negative value push (1) value
            stream.push_bit(1);
            stream.push_u16(-1*num);
        }else{
            // positive or zero value push (0) value
            stream.push_bit(0);
            stream.push_u16(num);
        }
    }

    void push_value_n(OutputBitStream& stream, int value, u16 num_bits){
        if(value < 0){
            // negative value push (1) value
            stream.push_bit(1);
            stream.push_bits((-1*value), num_bits);
        }else{
            // positive or zero value push (0) value
            stream.push_bit(0);
            stream.push_bits((value), num_bits);
        }
    }

    void push_delta_value(OutputBitStream& stream, int num){
        if(num > 0){
            // positive start with 10
            stream.push_bit(1);
            stream.push_bit(0);
        }else if(num < 0){
            //negative start with 11
            stream.push_bit(1);
            stream.push_bit(1);
            num *= -1;
        }
        for(u16 i = 0; i < num-1; i++){
            stream.push_bit(1);
        }
        stream.push_bit(0);
    }

    void push_quantized_array(OutputBitStream& stream, const Array64& array){
        for(u32 idx = 0; idx < 64; idx++){
            push_value(stream, array.at(idx));
        }
    }

    Array64 quantized_to_delta(const Array64& quantized){
        Array64 delta_values;
        delta_values.at(0) = quantized.at(0);
        delta_values.at(1) = quantized.at(1);

        for(u32 idx = 2; idx < 64; idx++){
            delta_values.at(idx) = quantized.at(idx) - quantized.at(idx-1);
        }
        return delta_values;
    }

    void push_symbol_huffman(OutputBitStream& stream, int symbol){
        u32 num_bits = symbol_length[symbol];
        u32 code_bits = symbol_encoding[symbol];

        for(u32 idx = 0; idx < num_bits; idx++){
            stream.push_bit(code_bits >> (num_bits - idx - 1) & 1);
        }
    }

    void push_unary(OutputBitStream& stream, u32 value){
        for(u32 idx = 0; idx < value; idx++)
            stream.push_bit(1);
        stream.push_bit(0);
    }

    u32 count_RLE_zeros(const Array64& array, u32 start){
        u32 idx = start; 
        u32 count = 0;
        while(idx < 64 && array.at(idx) == 0){
            count++;
            idx++;
        }
        return count;
    }

    void push_quantized_array_delta(OutputBitStream& stream, const Array64& array){

        Array64 delta_values = quantized_to_delta(array);
        for(double delta : delta_values)
            delta_frequency[delta]++;

        // Send first 2 values as normal
        push_value(stream, delta_values.at(0));
        push_value(stream, delta_values.at(1));

        // Use huffman codes to send over delta values
        u32 idx = 2;
        while(idx < 64){
            if(delta_values.at(idx) < -5){
                push_symbol_huffman(stream, -100);
                push_unary(stream, -1 * delta_values.at(idx++));
            }else if(delta_values.at(idx) > 5){
                push_symbol_huffman(stream, 100);
                push_unary(stream, delta_values.at(idx++));
            }else if(delta_values.at(idx) != 0){
                push_symbol_huffman(stream, delta_values.at(idx++));
            }else{
                // count the run length of zero
                u32 num_zeros = count_RLE_zeros(delta_values, idx);
                idx += num_zeros;

                if(idx == 64){  
                    push_symbol_huffman(stream, 150);
                    return;
                }
                while(num_zeros >= 8){
                    push_symbol_huffman(stream, 120);   // 8 zeros
                    num_zeros -= 8;
                }
                while(num_zeros > 0){
                    push_symbol_huffman(stream, 0);     // single zero
                    num_zeros --;
                }
            }
        }
    }

    /* ----- Decompressor code -----*/

    void read_header(InputBitStream& stream, dct::Quality& quality, u16& height, u16& width){
        u32 q = stream.read_bits(2);
        if(q == 0){
            quality = dct::Quality::low;
        }else if (q == 1){
            quality = dct::Quality::medium;
        }else{
            quality = dct::Quality::high;
        }

       height = stream.read_u16();
       width = stream.read_u16();
    }

    int read_value(InputBitStream& stream){
        // (+) sign = 0    (-) sign = 1
        bool sign = stream.read_bit();
        int num  = stream.read_u16();
        num = (sign == 1) ? (-1*num) : num ;
        return num;
    }

    int read_value_n(InputBitStream& stream, u16 num_bits){
        // (+) sign = 0    (-) sign = 1
        bool sign = stream.read_bit();
        int num  = stream.read_bits(num_bits);
        num = (sign == 1) ? (-1*num) : num ;
        return num;
    }

    int read_delta_value(InputBitStream& stream){
        if(stream.read_bit() == 0){
            return 0;
        }
        // 0 --> (+)   and 1 --> (-)
        bool sign = stream.read_bit();
        int num = 1;
        while(stream.read_bit()){
            num++;
        }
        num = (sign == 1) ? -1*(num) : num;
        return num;
    }

    Array64 read_quantized_array(InputBitStream& stream){
        Array64 block;
        for(u32 idx = 0; idx < 64; idx++){
            block.at(idx) = read_value(stream);
        }
        return block;
    }

    Array64 delta_to_quantized(const Array64& delta){
        Array64 quantized;
        quantized.at(0) = delta.at(0);
        quantized.at(1) = delta.at(1);

        for(u32 idx = 2; idx < 64; idx++){
            quantized.at(idx) = quantized.at(idx-1) + delta.at(idx);
        }
        return quantized;
    }

    int read_symbol_huffman(InputBitStream& stream){
        int value = 0;
        u32 iter = 0;
        while(iter < 10){
            bool bit_read = stream.read_bit();
            value = (value << 1) | bit_read;
            auto symbol_ref = encoding_symbol.find(value);
            if (symbol_ref != encoding_symbol.end() && iter >= 1)
                return symbol_ref->second;      
            iter++;
        }
    }

    int read_unary(InputBitStream& stream){
        int value = 0;
        while(stream.read_bit())
            value++;
        return value;
    }

    Array64 read_quantized_array_delta(InputBitStream& stream){
        Array64 delta_values;

        // Read first 2 as normal
        delta_values.at(0) = read_value(stream);
        delta_values.at(1) = read_value(stream);

        // Read the rest with huffman codes 
        u32 idx = 2;
        while(idx < 64){
            int curr_symbol = read_symbol_huffman(stream);
            if(curr_symbol == -100){
                delta_values.at(idx++) = -1 * read_unary(stream);
            }else if(curr_symbol == 100){
                delta_values.at(idx++) = 1 * read_unary(stream);
            }else if(curr_symbol == 120){
                for(u32 i = 0; i < 8; i++)
                    delta_values.at(idx++) = 0;
            }else if(curr_symbol == 150){
                while(idx < 64)
                    delta_values.at(idx++) = 0;
            }else{
                delta_values.at(idx++) = curr_symbol;
            }
        }

        Array64 quantized = delta_to_quantized(delta_values);
        return quantized;
    }

}