#ifndef STREAM
#define STREAM

#include <vector>
#include <array>
#include <map>
#include "output_stream.hpp"
#include "input_stream.hpp"
#include "discrete_cosine_transform.hpp"

// extern std::map<int,int> delta_frequency;
// extern std::map<int,int> RLE_frequency;


namespace stream{

    void print_histograms();
    void huffman_print();

    /* ----- Compressor code -----*/

    void push_header(OutputBitStream& stream, dct::Quality q, u16 height, u16 width);
    void push_value(OutputBitStream& stream, int num);
    void push_value_n(OutputBitStream& stream, int value, u16 num_bits);
    void push_delta_value(OutputBitStream& stream, int num);
    void push_quantized_array(OutputBitStream& stream, const Array64& array);
    u32 push_RLE_zeros(OutputBitStream& stream, const Array64& array, u32 start);
    void push_motion_vector_RLE(OutputBitStream& stream, const std::vector<int>& mv);
    Array64 quantized_to_delta(const Array64& quantized);
    void push_quantized_array_delta(OutputBitStream& stream, const Array64& array);

    /* ----- Decompressor code -----*/

    void read_header(InputBitStream& stream, dct::Quality& quality, u16& height, u16& width);
    int read_value(InputBitStream& stream);
    int read_value_n(InputBitStream& stream, u16 num_bits);
    int read_delta_value(InputBitStream& stream);
    Array64 read_quantized_array(InputBitStream& stream);
    Array64 delta_to_quantized(const Array64& delta);
    void add_RLE_zeros(Array64& delta_values, u32 start, u32 count);
    std::vector<int> read_motion_vector_RLE(InputBitStream& stream, int num_vectors);
    Array64 read_quantized_array_delta(InputBitStream& stream);
  
}

#endif