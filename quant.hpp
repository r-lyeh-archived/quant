// Quantization suite: float <-> 16-bit half, s/unorm <-> 8/10-bit short, quaternion <-> 32-bit int, position <-> 32/48-bit int, scale <-> 32/16-bit int
// r-lyeh, public domain

#pragma once
#include <math.h>
#include <stdint.h>

#define QUANT_VERSION "1.0.2" /* (2015/06/05) expanded api: 48-bit and 16-bit vector support; improved numerical stability
#define QUANT_VERSION "1.0.1" /* (2015/06/02) improved quat precision; fixed rotation bug
#define QUANT_VERSION "1.0.0" // (2015/05/29) initial revision */

// [usage]
// - Quantize animations (from standard 64-bytes (matrix4x4f) to 12-bytes (uint32_t pos,rot,sca) per bone).
// - Quantize colors.
// - Quantize sounds.
// - Quantize user input.
// - Quantize network packets.
// - ...
// 
// [todo]
// - Research: Integrate (and finish) curve simplification and hermite splines that I have somewhere around.
// - Research: Lossless/lossy animation format proposal:
//   - demultiplex vertex streams to a single giant (mono) stream
//   - quantize whole stream `(iter - min) / ( max - min ) -> [0..1]`
//   - apply lossless audio (FLAC) or lossy audio codec (mp3/ogg)
//   - decode, upscale stream `iter * ( max - min ) + min`, and multiplex vertex
//   - profit (?)
//
// [refs]
// - https://gist.github.com/rygorous/2156668
// - http://zeuxcg.org/2010/12/14/quantizing-floats/
// - http://en.wikipedia.org/wiki/Fast_inverse_square_root
// - http://bitsquid.blogspot.com.es/2009/11/bitsquid-low-level-animation-system.html 

// [api: generics, per component]
namespace quant {
/*/ float to 16-bit half
/*/ static uint16_t  encode16_half( float fl );
/*/ 16-bit half to float
/*/ static float     decode16_half( uint16_t half );
/*/ float [0..1] to 8-bit byte
/*/ static uint8_t   encode8_unorm( float    x );
/*/ 8-bit byte to float [0..1]
/*/ static float     decode8_unorm( uint8_t  x );
/*/ float [-1..1] to 8-bit byte
/*/ static uint8_t   encode8_snorm( float    x );
/*/ 8-bit byte to float [-1..1]
/*/ static float     decode8_snorm( uint8_t  x );
/*/ float [-1..1] to 8-bit byte (OpenGL version)
/*/ static uint8_t   encode8_snorm_gl2( float   x );
/*/ 8-bit byte to float [-1..1] (OpenGL version)
/*/ static float     decode8_snorm_gl2( uint8_t x );

/*/ generic N-bits[1..16] half encoder
/*/ template<unsigned N> static uint16_t  encode_half( float fl ) { return encode16_half(fl) >> (16-N); }
/*/ generic N-bits[1..16] half decoder
/*/ template<unsigned N> static float     decode_half( uint16_t fl ) { return decode16_half(fl << (16-N)); }

/*/ generic N-bit unorm encoder (based on D3D10_UNORM)
/*/ template<unsigned N> static uint16_t  encode_unorm( float    x ) { return uint16_t( int (x * ((1<<(N))-1) + 0.5f) ); }
/*/ generic N-bit unorm decoder (based on D3D10_UNORM)
/*/ template<unsigned N> static float     decode_unorm( uint16_t x ) { return x / float((1<<(N))-1); }
/*/ generic N-bit snorm encoder (based on D3D10_UNORM)
/*/ template<unsigned N> static uint16_t  encode_snorm( float    x ) { return (x < 0) | (encode_unorm<N-1>(x < 0 ? -x : x) << 1); }
/*/ generic N-bit snorm decoder (based on D3D10_UNORM)
/*/ template<unsigned N> static float     decode_snorm( uint16_t x ) { return decode_unorm<N-1>(x>>1) * (x & 1 ? -1:1); }
}

// [api: specialized, per type]
namespace quant {
/*/ quaternion to 32-bit integer
/*/ static void      encode101010_quat( uint32_t &out, float x, float y, float z, float w );
/*/ 32-bit integer to quaternion
/*/ static void      decode101010_quat( float &x, float &y, float &z, float &w, uint32_t in );
/*/ position or scale to 48-bit integer
/*/ static void      encode161616_vec( uint64_t &out, float x, float y, float z );
/*/ 48-bit integer to position or scale
/*/ static void      decode161616_vec( float &x, float &y, float &z, uint64_t in );
/*/ position or scale to 32-bit integer
/*/ static void      encode8814_vec( uint32_t &out, float x, float y, float z );
/*/ 32-bit integer to position or scale
/*/ static void      decode8814_vec( float &x, float &y, float &z, uint32_t in );
/*/ position or scale to 16-bit integer
/*/ static void      encode555_vec( uint16_t &out, float x, float y, float z );
/*/ 16-bit integer to position or scale
/*/ static void      decode555_vec( float &x, float &y, float &z, uint16_t in );

/*/ quaternion to 32-bit integer (struct version)
/*/ template<typename T> static void encode101010_quat( uint32_t &out, const T &q );
/*/ 32-bit integer to quaternion (struct version)
/*/ template<typename T> static void decode101010_quat( T &q, uint32_t in );

/*/ position or scale to 48-bit integer (struct version) (good compromise)
/*/ template<typename vec3> static void encode161616_vec( uint64_t &out, const vec3 &v );
/*/ 48-bit integer to position or scale (struct version) (good compromise)
/*/ template<typename vec3> static void decode161616_vec( vec3 &v, uint64_t in );
/*/ position or scale to 16-bit integer (struct version)
/*/ template<typename vec3> static void encode555_vec( uint16_t &out, const vec3 &v );
/*/ 16-bit integer to position or scale (struct version)
/*/ template<typename vec3> static void decode555_vec( vec3 &v, uint16_t in );
}

// [api: specialized, scale/positions variants]
namespace quant {
/*/ position or scale to 32-bit integer (struct version)(best rotation, poor distance)
/*/ template<typename vec3> static void encode11118_vec( uint32_t &out, const vec3 &v );
/*/ 32-bit integer to position or scale (struct version)(best rotation, poor distance)
/*/ template<typename vec3> static void decode11118_vec( vec3 &v, uint32_t in );
/*/ position or scale to 32-bit integer (struct version)(large rotation, small distance)
/*/ template<typename vec3> static void encode101010_vec( uint32_t &out, const vec3 &v );
/*/ 32-bit integer to position or scale (struct version)(large rotation, small distance)
/*/ template<typename vec3> static void decode101010_vec( vec3 &v, uint32_t in );
/*/ position or scale to 32-bit integer (struct version)(medium orientation, medium distance)
/*/ template<typename vec3> static void encode9912_vec( uint32_t &out, const vec3 &v );
/*/ 32-bit integer to position or scale (struct version)(medium orientation, medium distance)
/*/ template<typename vec3> static void decode9912_vec( vec3 &v, uint32_t in );
/*/ position or scale to 32-bit integer (struct version)(small orientation, large distance) (good compromise)
/*/ template<typename vec3> static void encode8814_vec( uint32_t &out, const vec3 &v );
/*/ 32-bit integer to position or scale (struct version)(small orientation, large distance) (good compromise)
/*/ template<typename vec3> static void decode8814_vec( vec3 &v, uint32_t in );
/*/ position or scale to 32-bit integer (struct version)(poor orientation, best distance)
/*/ template<typename vec3> static void encode7716_vec( uint32_t &out, const vec3 &v );
/*/ 32-bit integer to position or scale (struct version)(poor orientation, best distance)
/*/ template<typename vec3> static void decode7716_vec( vec3 &v, uint32_t in );
}

// [api: misc]
namespace quant {
/*/ helper function. reverse/inverse square root (-DQUANT_USE_STD_SQRT to use standard sqrt() instead)
/*/ static float     rsqrt( float number );
/*/ helper function. remap floating number in range [min1..max1] to range [min2..max2]
/*/ static float     remap( float x, float min1, float max1, float min2, float max2 );
}




// implementation, beware of dog

namespace quant {

//
// Helper: inverse square root (default: fast, less precision, predictable results across platforms)
static float rsqrt( float number ) {
#if defined(QUANT_USE_STD_SQRT)
    // a more precise inverse square root.
    // also, this provides less predictable floating results across platforms.
    return float(1 / sqrt(number));
#else
    // fast inverse square root.
    // [ref] http://en.wikipedia.org/wiki/Fast_inverse_square_root
    long i;
    float x2, y;
    const float threehalfs = 1.5F;
    x2 = number * 0.5F;
    y  = number;
    i  = * ( long * ) &y;                       // evil floating point bit level hacking
    i  = 0x5f3759df - ( i >> 1 );               // what the fuck? 
    y  = * ( float * ) &i;
    y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
    y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed
    y  = y * ( threehalfs - ( x2 * y * y ) );   // 3rd iteration, this can be removed
    return y;
#endif
}

//
// Helper: remap floating number in range [min1..max1] to range [min2..max2]
static float remap( float x, float min1, float max1, float min2, float max2 ) {
    return ( ( ( x - min1 ) / ( max1 - min1 ) ) * ( max2 - min2 ) ) + min2;
}

// 
// [src]: https://gist.github.com/rygorous/2156668
union FP32 {
    uint32_t u;
    float f;
    struct {
        uint32_t Mantissa : 23;
        uint32_t Exponent : 8;
        uint32_t Sign : 1;
    };
};
 
union FP16 {
    uint16_t u;
    struct {
        uint32_t Mantissa : 10;
        uint32_t Exponent : 5;
        uint32_t Sign : 1;
    };
};
 
static uint16_t encode16_half( float fl ) {
    FP16 o = { 0 };
    FP32 f; f.f = fl;
    // Based on ISPC reference code (with minor modifications)
    if (f.Exponent == 0) { 
        // if Signed zero/denormal (which will underflow)
        o.Exponent = 0;
    } else if (f.Exponent == 255) { 
        // if Inf or NaN (all exponent bits set)
        o.Exponent = 31;
        o.Mantissa = f.Mantissa ? 0x200 : 0; // NaN->qNaN and Inf->Inf
    }
    else { 
        // if Normalized number
        // Exponent unbias the single, then bias the halfp
        int newexp = f.Exponent - 127 + 15;
        if (newexp >= 31) { 
            // if Overflow, return signed infinity
            o.Exponent = 31;
        } else if (newexp <= 0) { 
            // if Underflow
            if ((14 - newexp) <= 24) { 
                // Mantissa might be non-zero
                uint32_t mant = f.Mantissa | 0x800000; // Hidden 1 bit
                o.Mantissa = mant >> (14 - newexp);
                // Check for rounding
                if ((mant >> (13 - newexp)) & 1) { 
                    // Round, might overflow into exp bit, but this is OK
                    o.u++; 
            }   }
        } else {
            o.Exponent = newexp;
            o.Mantissa = f.Mantissa >> 13;
            // Check for rounding
            if (f.Mantissa & 0x1000) { 
                // Round, might overflow to inf, this is OK
                o.u++; 
        }   }
    }
    o.Sign = f.Sign;
    return o.u;
}
 
static float decode16_half( uint16_t half ) {
    static const FP32 magic = { 113 << 23 };
    static const uint32_t shifted_exp = 0x7c00 << 13; // exponent mask after shift
    FP32 o;
    FP16 h; h.u = half;
 
    // exponent/mantissa bits
    o.u = (h.u & 0x7fff) << 13;     
    // just the exponent
    uint32_t exp = shifted_exp & o.u;   
    // exponent adjust
    o.u += (127 - 15) << 23;        
    // handle exponent special cases
    if (exp == shifted_exp) { 
        // if Inf/NaN
        // extra exp adjust
        o.u += (128 - 16) << 23;    
    } else if (exp == 0) { 
        // if Zero/Denormal
        // extra exp adjust
        o.u += 1 << 23;             
        // renormalize
        o.f -= magic.f;             
    }
    // sign bit
    o.u |= (h.u & 0x8000) << 16;    
    return o.f;
}

//
// [ref] http://zeuxcg.org/2010/12/14/quantizing-floats/
// D3D10, GL2: for *_UNORM formats of n-bit length, the decode function is decode(x) = x / (2^n - 1)
// Unsigned quantization: input: [0..1] float; output: [0..255] integer
// D3D10: for *_SNORM formats of n-bit length the decode function is decode(x) = clamp(x / (2^(n-1) - 1), -1, 1)
// Signed quantization for D3D10 rules: input: [-1..1] float; output: [-127..127] integer
static uint8_t   encode8_unorm(float    x) { return uint8_t( int (x * 255.f + 0.5f) ); }
static float     decode8_unorm(uint8_t  x) { return x / 255.f; }
static uint8_t   encode8_snorm(float    x) { return uint8_t( int (x * 127.f + (x > 0 ? 0.5f : -0.5f)) ); }
static float     decode8_snorm(uint8_t  x) { float f = x / 127.f; return f <= -1 ? -1.f : (f >= 1 ? 1.f : f); }
//
// [ref] http://zeuxcg.org/2010/12/14/quantizing-floats/
// OpenGL2: same decoding function for unsigned numbers, but a different one for signed: decode(x) = (2x + 1) / (2^n - 1)
// Signed quantization for OpenGL rules: input: [-1..1] float; output: [-128..127] integer
// Warning: This has slightly better precision (all numbers encode distinct values), but can't represent 0 exactly. 
static uint8_t   encode8_snorm_gl2(float   x) { return uint8_t( int (x * 127.5f + (x >= 0.f ? 0.f : -1.f)) ); }
static float     decode8_snorm_gl2(uint8_t x) { return (2*x + 1) / 255.f; }

//
// Quantize rotation as 10+10+10 bits for quaternion components

static void encode101010_quat( uint32_t &out, float x, float y, float z, float w ) {
    // [ref] http://bitsquid.blogspot.com.es/2009/11/bitsquid-low-level-animation-system.html 
    // "For quaternions we use 2 bits to store the index of the largest component,
    // then 10 bits each to store the value of the remaining three components. We use the knowledge
    // that 1 = x^2 + y^2 + z^2 + w^2 to restore the largest component, so we don't actually have to store its value.
    // Since we don't store the largest component we know that the remaining ones must be in the range (-1/sqrt(2), 1/sqrt(2))
    // (otherwise, one of them would be largest). So we use the 10 bits to quantize a value in that range, giving us a precision of 0.0014.
    // The quaternions (x, y, z, w) and (-x, -y, -z, -w) represent the same rotation, so I flip the signs so that the largest component
    // is always positive." - Niklas Frykholm / bitsquid.se
    static const float rmin = -rsqrt(2), rmax = rsqrt(2);
    float xx = x*x, yy = y*y, zz = z*z, ww = w*w;
    /**/ if( xx >= yy && xx >= zz && xx >= ww ) {
        y = remap( y, rmin, rmax, -1, 1 );
        z = remap( z, rmin, rmax, -1, 1 );
        w = remap( w, rmin, rmax, -1, 1 );
        out = x >= 0 ? uint32_t( (0<<30) | (encode_snorm<10>( y)<<20) | (encode_snorm<10>( z)<<10) | (encode_snorm<10>( w)<<0) )
                     : uint32_t( (0<<30) | (encode_snorm<10>(-y)<<20) | (encode_snorm<10>(-z)<<10) | (encode_snorm<10>(-w)<<0) );
    }
    else if( yy >= zz && yy >= ww ) {
        x = remap( x, rmin, rmax, -1, 1 );
        z = remap( z, rmin, rmax, -1, 1 );
        w = remap( w, rmin, rmax, -1, 1 );
        out = y >= 0 ? uint32_t( (1<<30) | (encode_snorm<10>( x)<<20) | (encode_snorm<10>( z)<<10) | (encode_snorm<10>( w)<<0) )
                     : uint32_t( (1<<30) | (encode_snorm<10>(-x)<<20) | (encode_snorm<10>(-z)<<10) | (encode_snorm<10>(-w)<<0) );
    }
    else if( zz >= ww ) {
        x = remap( x, rmin, rmax, -1, 1 );
        y = remap( y, rmin, rmax, -1, 1 );
        w = remap( w, rmin, rmax, -1, 1 );
        out = z >= 0 ? uint32_t( (2<<30) | (encode_snorm<10>( x)<<20) | (encode_snorm<10>( y)<<10) | (encode_snorm<10>( w)<<0) )
                     : uint32_t( (2<<30) | (encode_snorm<10>(-x)<<20) | (encode_snorm<10>(-y)<<10) | (encode_snorm<10>(-w)<<0) );
    }
    else {
        x = remap( x, rmin, rmax, -1, 1 );
        y = remap( y, rmin, rmax, -1, 1 );
        z = remap( z, rmin, rmax, -1, 1 );
        out = w >= 0 ? uint32_t( (3<<30) | (encode_snorm<10>( x)<<20) | (encode_snorm<10>( y)<<10) | (encode_snorm<10>( z)<<0) )
                     : uint32_t( (3<<30) | (encode_snorm<10>(-x)<<20) | (encode_snorm<10>(-y)<<10) | (encode_snorm<10>(-z)<<0) );
    }
}

//
// DeQuantize rotation as 10+10+10 bits for quaternion components

static void decode101010_quat( float &x, float &y, float &z, float &w, uint32_t in ) {
    // [ref] http://bitsquid.blogspot.com.es/2009/11/bitsquid-low-level-animation-system.html 
    // See encode101010_quat() function above.
    static const float rmin = -rsqrt(2), rmax = rsqrt(2);
    switch( in >> 30 ) {
        default: case 0:
            y = decode_snorm<10>( ( in >> 20 ) & 0x3FF );
            z = decode_snorm<10>( ( in >> 10 ) & 0x3FF );
            w = decode_snorm<10>( ( in >>  0 ) & 0x3FF );
            y = remap( y, -1, 1, rmin, rmax );
            z = remap( z, -1, 1, rmin, rmax );
            w = remap( w, -1, 1, rmin, rmax );
            x = 1 / rsqrt( 1 - y * y - z * z - w * w );
        break; case 1:
            x = decode_snorm<10>( ( in >> 20 ) & 0x3FF );
            z = decode_snorm<10>( ( in >> 10 ) & 0x3FF );
            w = decode_snorm<10>( ( in >>  0 ) & 0x3FF );
            x = remap( x, -1, 1, rmin, rmax );
            z = remap( z, -1, 1, rmin, rmax );
            w = remap( w, -1, 1, rmin, rmax );
            y = 1 / rsqrt( 1 - x * x - z * z - w * w );
        break; case 2:
            x = decode_snorm<10>( ( in >> 20 ) & 0x3FF );
            y = decode_snorm<10>( ( in >> 10 ) & 0x3FF );
            w = decode_snorm<10>( ( in >>  0 ) & 0x3FF );
            x = remap( x, -1, 1, rmin, rmax );
            y = remap( y, -1, 1, rmin, rmax );
            w = remap( w, -1, 1, rmin, rmax );
            z = 1 / rsqrt( 1 - x * x - y * y - w * w );
        break; case 3:
            x = decode_snorm<10>( ( in >> 20 ) & 0x3FF );
            y = decode_snorm<10>( ( in >> 10 ) & 0x3FF );
            z = decode_snorm<10>( ( in >>  0 ) & 0x3FF );
            x = remap( x, -1, 1, rmin, rmax );
            y = remap( y, -1, 1, rmin, rmax );
            z = remap( z, -1, 1, rmin, rmax );
            w = 1 / rsqrt( 1 - x * x - y * y - z * z );
    }
}

//
// Quantize scale/position as X+Y bits for unit vector rotation + Z bits for vector length (note: 2 bits reserved)
// Requires: each X, Y, Z in [7..16] range && X+Y+Z <= 30

template<unsigned X, unsigned Y, unsigned Z>
static void encode_vec( uint32_t &out, float x, float y, float z ) { 
    // Somehow similar to encode101010_quat() function above.
    // We decompose given vector into unit vector and length (magnitude). Then we discard the largest component, as unit vectors
    // follow 1 = x^2 + y^2 + z^2 expression (similar to encode101010_quat() function above). The unit vectors (x, y, z) and
    // (-x, -y, -z) represent the same direction, so I flip the magnitude so that the largest component is always positive.
    static const float rmin = -rsqrt(2), rmax = rsqrt(2);
    float xx = x*x, yy = y*y, zz = z*z;
    float len = rsqrt( xx + yy + zz ); // float len = sqrt( xx + yy + zz );
    x *= len; y *= len; z *= len;      // x /= len; y /= len; z /= len;
    /****/ if( xx >= yy && xx >= zz ) {
        // y = remap( y, rmin, rmax, -1, 1 );
        // z = remap( z, rmin, rmax, -1, 1 );
        out = ( x >= 0 ? uint32_t( (0<<30)|(encode_snorm<X>( y)<<(Z+Y))|(encode_snorm<Y>( z)<<Z)|(encode_half<Z>( 1/len)) ) :  // ( len)
                         uint32_t( (0<<30)|(encode_snorm<X>(-y)<<(Z+Y))|(encode_snorm<Y>(-z)<<Z)|(encode_half<Z>(-1/len)) ) ); // (-len)
    } else if( yy >= zz ) {
        // x = remap( x, rmin, rmax, -1, 1 );
        // z = remap( z, rmin, rmax, -1, 1 );
        out = ( y >= 0 ? uint32_t( (1<<30)|(encode_snorm<X>( x)<<(Z+Y))|(encode_snorm<Y>( z)<<Z)|(encode_half<Z>( 1/len)) ) :  // ( len)
                         uint32_t( (1<<30)|(encode_snorm<X>(-x)<<(Z+Y))|(encode_snorm<Y>(-z)<<Z)|(encode_half<Z>(-1/len)) ) ); // (-len)
    } else {
        // x = remap( x, rmin, rmax, -1, 1 );
        // y = remap( y, rmin, rmax, -1, 1 );
        out = ( z >= 0 ? uint32_t( (2<<30)|(encode_snorm<X>( x)<<(Z+Y))|(encode_snorm<Y>( y)<<Z)|(encode_half<Z>( 1/len)) ) :  // ( len)
                         uint32_t( (2<<30)|(encode_snorm<X>(-x)<<(Z+Y))|(encode_snorm<Y>(-y)<<Z)|(encode_half<Z>(-1/len)) ) ); // (-len)
    }
}

//
// DeQuantize scale/position as X+Y bits for unit vector rotation + Z bits for vector length (note: 2 bits reserved)
// Requires: each X, Y, Z in [7..16] range && X+Y+Z <= 30

template<unsigned X, unsigned Y, unsigned Z>
static void decode_vec( float &x, float &y, float &z, uint32_t in ) { 
    // See encode_vec() function above.
    static const float rmin = -rsqrt(2), rmax = rsqrt(2);
    switch( in >> 30 ) {
        default: case 0:
            y = decode_snorm<X>( ( in >> (Z+Y) ) & ((1<<X)-1) );
            z = decode_snorm<Y>( ( in >>    Z  ) & ((1<<Y)-1) );
            // y = remap( y, -1, 1, rmin, rmax );
            // z = remap( z, -1, 1, rmin, rmax );
            x = 1 / rsqrt( 1 - y * y - z * z ); // x = sqrt( 1 - y * y - z * z );
        break; case 1:
            x = decode_snorm<X>( ( in >> (Z+Y) ) & ((1<<X)-1) );
            z = decode_snorm<Y>( ( in >>    Z  ) & ((1<<Y)-1) );
            // x = remap( x, -1, 1, rmin, rmax );
            // z = remap( z, -1, 1, rmin, rmax );
            y = 1 / rsqrt( 1 - x * x - z * z ); // y = sqrt( 1 - x * x - z * z );
        break; case 2:
            x = decode_snorm<X>( ( in >> (Z+Y) ) & ((1<<X)-1) );
            y = decode_snorm<Y>( ( in >>    Z  ) & ((1<<Y)-1) );
            // x = remap( x, -1, 1, rmin, rmax );
            // y = remap( y, -1, 1, rmin, rmax );
            z = 1 / rsqrt( 1 - x * x - y * y ); // z = sqrt( 1 - x * x - y * y );
    }
    float len = decode_half<Z>( uint16_t(in & ((1<<Z)-1)) ); 
    x *= len, y *= len, z *= len;
}

// position or scale to 48-bit integer
static void      encode161616_vec( uint64_t &out, float x, float y, float z ) {
    out  = ((uint64_t)encode_half<16>(x)) << 32;
    out |= ((uint64_t)encode_half<16>(y)) << 16;
    out |= ((uint64_t)encode_half<16>(z)) <<  0;
}
// 48-bit integer to position or scale
static void      decode161616_vec( float &x, float &y, float &z, uint64_t in ) {
    x = decode_half<16>((in >> 32) & 0xffff);
    y = decode_half<16>((in >> 16) & 0xffff);
    z = decode_half<16>((in >>  0) & 0xffff);
}
// position or scale to 48-bit integer (struct version)
template<typename vec3> static void      encode161616_vec( uint64_t &out, const vec3 &v ) { encode161616_vec( out, v.x, v.y, v.z ); }
// 48-bit integer to position or scale (struct version)
template<typename vec3> static void      decode161616_vec( vec3 &v, uint64_t in ) { decode161616_vec( v.x, v.y, v.z, in ); }


// in-place class sugars

template<typename quat>
static void encode101010_quat( uint32_t &out, const quat &q ) {
    encode101010_quat( out, q.x, q.y, q.z, q.w );
}
template<typename quat>
static void decode101010_quat( quat &q, uint32_t in ) {
    decode101010_quat( q.x, q.y, q.z, q.w, in );
}


// position or scale to 16-bit integer (struct version)
static void      encode555_vec( uint16_t &out, float x, float y, float z ) {
    uint16_t x5 = encode_snorm<5>( x );
    uint16_t y5 = encode_snorm<5>( y );
    uint16_t z5 = encode_snorm<5>( z );
    out = (x5 << 10) | (y5 << 5) | (z5 << 0);
}
// 16-bit integer to position or scale (struct version)
static void      decode555_vec( float &x, float &y, float &z, uint16_t in ) {
    x = decode_snorm<5>( (in>>10) & 0x1f );
    y = decode_snorm<5>( (in>> 5) & 0x1f );
    z = decode_snorm<5>( (in>> 0) & 0x1f );
}
// position or scale to 16-bit integer (struct version)
template<typename vec3> static void encode555_vec( uint16_t &out, const vec3 &v ) {
    encode555_vec( out, v.x, v.y, v.z );
}
// 16-bit integer to position or scale (struct version)
template<typename vec3> static void decode555_vec( vec3 &v, uint16_t in ) {
    decode555_vec( v.x, v.y, v.z, in );
}




template<unsigned X, unsigned Y, unsigned Z, typename vec3>
static void encode_vec( uint32_t &out, const vec3 &v ) { 
    encode_vec<X,Y,Z>( out, v.x, v.y, v.z );
}
template<unsigned X, unsigned Y, unsigned Z, typename vec3>
static void decode_vec( vec3 &v, uint32_t in ) { 
    decode_vec<X,Y,Z>( v.x, v.y, v.z, in );
}

template<typename vec3>
static void encode7716_vec( uint32_t &out, const vec3 &v ) { 
    encode_vec<7,7,16>( out, v.x, v.y, v.z );
}
template<typename vec3>
static void decode7716_vec( vec3 &v, uint32_t in ) { 
    decode_vec<7,7,16>( v.x, v.y, v.z, in );
}

template<typename vec3>
static void encode8814_vec( uint32_t &out, const vec3 &v ) { 
    encode_vec<8,8,14>( out, v.x, v.y, v.z );
}
template<typename vec3>
static void decode8814_vec( vec3 &v, uint32_t in ) { 
    decode_vec<8,8,14>( v.x, v.y, v.z, in );
}

template<typename vec3>
static void encode9912_vec( uint32_t &out, const vec3 &v ) { 
    encode_vec<9,9,12>( out, v.x, v.y, v.z );
}
template<typename vec3>
static void decode9912_vec( vec3 &v, uint32_t in ) { 
    decode_vec<9,9,12>( v.x, v.y, v.z, in );
}

template<typename vec3>
static void encode101010_vec( uint32_t &out, const vec3 &v ) { 
    encode_vec<10,10,10>( out, v.x, v.y, v.z );
}
template<typename vec3>
static void decode101010_vec( vec3 &v, uint32_t in ) { 
    decode_vec<10,10,10>( v.x, v.y, v.z, in );
}

template<typename vec3>
static void encode11118_vec( uint32_t &out, const vec3 &v ) { 
    encode_vec<11,11,8>( out, v.x, v.y, v.z );
}
template<typename vec3>
static void decode11118_vec( vec3 &v, uint32_t in ) { 
    decode_vec<11,11,8>( v.x, v.y, v.z, in );
}

}

// } quant.hpp
 
#ifdef QUANT_TESTS
 
#include <cmath>
#include <iostream>
 
float acc[3] = {};
float hit[3] = {};
 
void verify( float x ) { 
    using namespace quant;
    const int y1 = encode8_unorm(x);     const float x1 = decode8_unorm(y1);
    const int y2 = encode8_snorm(x);     const float x2 = decode8_snorm(y2);
    const int y3 = encode8_snorm_gl2(x); const float x3 = decode8_snorm_gl2(y3);
    
    float err;
    err = std::abs(x - x1); acc[0] += err; hit[0]++; std::cout << "v1 " << x1 << " vs " << x << " (error: " << err << ")" << std::endl;
    err = std::abs(x - x2); acc[1] += err; hit[1]++; std::cout << "v2 " << x2 << " vs " << x << " (error: " << err << ")" << std::endl;
    err = std::abs(x - x3); acc[2] += err; hit[2]++; std::cout << "v3 " << x3 << " vs " << x << " (error: " << err << ")" << std::endl;
}
 
// run 'sample | sort'
int main() {
    using namespace quant;
    for( short x = 0, e = 1000; x <= e; ++x ) {
        verify( x * 0.001f );
    }
    std::cout << "average error (lower is better)" << std::endl;
    std::cout << "avg-1:" << (acc[0] / hit[0]) << std::endl; 
    std::cout << "avg-2:" << (acc[1] / hit[1]) << std::endl; 
    std::cout << "avg-3:" << (acc[2] / hit[2]) << std::endl;

    std::cout << decode16_half( encode16_half(3.14159f) ) << std::endl;
}
 
#endif
