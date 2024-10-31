#ifndef DISCRETE_COSINE_TRANSFORM
#define DISCRETE_COSINE_TRANSFORM

#include <vector>
#include <cmath>
#include <array>
#include <cassert>
#include "yuv_stream.hpp"

using Block8x8 = std::array<std::array<double, 8>, 8>;
using Block16x16 = std::array<std::array<double, 16>, 16>;
using Array64 = std::array<int, 64>;
using u32 = std::uint32_t;
using u16 = std::uint16_t;
using u8 = std::uint8_t;

namespace dct {

    enum Quality {
        low = 0,
        medium,
        high,
        ERROR
    };

    enum Direction {
        right = 0,
        down,
        down_left,
        up_right
    };

    // the result of running create_c_matrix()  
    const Block8x8 c_matrix {{
        {0.353553,  0.353553,   0.353553,   0.353553,   0.353553,   0.353553,   0.353553,   0.353553    },
        {0.490393,  0.415735,   0.277785,   0.0975452,  -0.0975452, -0.277785,  -0.415735,  -0.490393   },
        {0.46194,   0.191342,   -0.191342,  -0.46194,   -0.46194,   -0.191342,  0.191342,   0.46194     },
        {0.415735,  -0.0975452, -0.490393,  -0.277785,  0.277785,   0.490393,   0.0975452,  -0.415735   },
        {0.353553,  -0.353553,  -0.353553,  0.353553,   0.353553,   -0.353553,  -0.353553,  0.353553    },
        {0.277785,  -0.490393,  0.0975452,  0.415735,   -0.415735,  -0.0975452, 0.490393,   -0.277785   },
        {0.191342,  -0.46194,   0.46194,    -0.191342,  -0.191342,  0.46194,    -0.46194,   0.191342    },
        {0.0975452, -0.277785,  0.415735,   -0.490393,  0.490393,   -0.415735,  0.277785,   -0.0975452  }
    }};

    // the result of running transpose_block() from discrete_cosine_transfrom.cpp 
    const Block8x8 c_matrix_transpose {{
        {0.353553, 0.490393,   0.46194,   0.415735,   0.353553,  0.277785,   0.191342,  0.0975452   },
        {0.353553, 0.415735,   0.191342,  -0.0975452, -0.353553, -0.490393,  -0.46194,  -0.277785   },
        {0.353553, 0.277785,   -0.191342, -0.490393,  -0.353553, 0.0975452,  0.46194,   0.415735    },
        {0.353553, 0.0975452,  -0.46194,  -0.277785,  0.353553,  0.415735,   -0.191342, -0.490393   },
        {0.353553, -0.0975452, -0.46194,  0.277785,   0.353553,  -0.415735,  -0.191342, 0.490393    },
        {0.353553, -0.277785,  -0.191342, 0.490393,   -0.353553, -0.0975452, 0.46194,   -0.415735   },
        {0.353553, -0.415735,  0.191342,  0.0975452,  -0.353553, 0.490393,   -0.46194,  0.277785    },
        {0.353553, -0.490393,  0.46194,   -0.415735,  0.353553,  -0.277785,  0.191342,  -0.0975452  }
    }};

    // quantization matrix used by JPEG - from lecture slides
    const Block8x8 luminance {{
        {16, 11, 10, 16, 24,  40,  51,  61},
        {12, 12, 14, 19, 26,  58,  60,  55},
        {14, 13, 16, 24, 40,  57,  69,  56},
        {14, 17, 22, 29, 51,  87,  80,  62},
        {18, 22, 37, 56, 68,  109, 103, 77},
        {24, 35, 55, 64, 81,  104, 113, 92},
        {49, 64, 78, 87, 103, 121, 120, 101},
        {72, 92, 95, 98, 112, 100, 103, 99}
    }};

    // quantization matrix used by JPEG - from lecture slides
    const Block8x8 chrominance {{
        {17, 18, 24, 47, 99, 99, 99, 99},
        {18, 21, 26, 66, 99, 99, 99, 99},
        {24, 26, 56, 99, 99, 99, 99, 99},
        {47, 99, 99, 99, 99, 99, 99, 99},
        {99, 99, 99, 99, 99, 99, 99, 99},
        {99, 99, 99, 99, 99, 99, 99, 99},
        {99, 99, 99, 99, 99, 99, 99, 99},
        {99, 99, 99, 99, 99, 99, 99, 99},
    }};

    const Block8x8 quantization_order {{
        {0,  1,  5,  6,  14, 15, 27, 28},
        {2,  4,  7,  13, 16, 26, 29, 42},
        {3,  8,  12, 17, 25, 30, 41, 43},
        {9,  11, 18, 24, 31, 40, 44, 53},
        {10, 19, 23, 32, 39, 45, 52, 54},
        {20, 22, 33, 38, 46, 51, 55, 60},
        {21, 34, 37, 47, 50, 56, 59, 61},
        {35, 36, 48, 49, 57, 58, 62, 63}
    }};

    /* ----- Written by Bill ------ */
    inline unsigned char round_and_clamp_to_char(double v){
        //Round to int 
        int i = (int)(v+0.5);
        //Clamp to the range [0,255]
        if (i < 0)
            return 0;
        else if (i > 255)
            return 255;
        return i;
    }

    /* ----- Block Operations ----- */
    Block8x8 create_c_matrix();
    void print_array(const Array64& array);
    void print_block(const Block8x8& block);
    void print_blocks(const std::vector<Block8x8>& blocks);
    Block8x8 multiply_block(const Block8x8& blockA, const Block8x8& blockB);
    Block8x8 transpose_block(const Block8x8& block);
    Block8x8 get_delta_block(const Block8x8& block1, const Block8x8& block2);
    Block8x8 add_delta_block(const Block8x8& block, const Block8x8& delta);
    Block16x16 create_macroblock(const Block8x8& b1, const Block8x8& b2, const Block8x8& b3, const Block8x8& b4);
    void get_prev_blocks(u32 macro_idx, YUVFrame420& prev_frame, const std::pair<u32, u32>& vector, std::vector<Block8x8>& prev_blocks);

    /* ----- Compressor Functions ----- */
    void partition_Y_channel(std::vector<Block8x8>& blocks, u32 height, u32 width, const std::vector<std::vector<unsigned char>>& channel);
    void partition_C_channel(std::vector<Block8x8>& blocks, u32 height, u32 width, const std::vector<std::vector<unsigned char>>& channel);
    Block8x8 get_dct(const Block8x8 &block);
    Block8x8 quantize_block(const Block8x8& block, Quality quality, bool is_luminance, bool is_P_block);
    Direction get_direction(u32 r, u32 c, Direction curr);
    Array64 block_to_array(const Block8x8& block);

    /* ----- Decompressor Functions ----- */
    Block8x8 array_to_block(const Array64& array);
    Block8x8 unquantize_block(const Block8x8& block, Quality quality, bool is_luminance, bool is_P_block);
    Block8x8 get_inverse_dct(const Block8x8& block);
    void undo_partition_C_channel(const std::vector<Block8x8>& blocks, u32 height, u32 width, std::vector<std::vector<unsigned char>>& channel);
    void undo_partition_Y_channel(const std::vector<Block8x8>& blocks, u32 height, u32 width, std::vector<std::vector<unsigned char>>& channel);

} // namespace dct

#endif
