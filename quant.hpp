// Quantization suite: floats <-> 16-bit halfs, s/unorm floats <-> 8/10-bit shorts, quaternions <-> 32-bit integers, positions <-> 32-bit integers
// r-lyeh, public domain

#pragma once
#include <math.h>
#include <stdint.h>

#define QUANT_VERSION "1.0.1" /* (2015/06/02) improved quant precision; fixed rotation bug
#define QUANT_VERSION "1.0.0" // (2015/05/29) initial revision */

// [usage]
// - Quantize animations from standard 64-bytes (matrix4x4f) to 12-bytes (uint32_t pos,rot,sca) per bone.
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

// [api]
namespace quant {
/*/ float to 16-bit half
/*/ static uint16_t  encode16_float( float fl );
/*/ 16-bit half to float
/*/ static float     decode16_float( uint16_t half );
/*/ float [0..1] to 8-bit byte
/*/ static uint8_t   encode8_unorm(float    x);
/*/ 8-bit byte to float [0..1]
/*/ static float     decode8_unorm(uint8_t  x);
/*/ float [-1..1] to 8-bit byte
/*/ static uint8_t   encode8_snorm(float    x);
/*/ 8-bit byte to float [-1..1]
/*/ static float     decode8_snorm(uint8_t  x);
/*/ float [-1..1] to 8-bit byte (OpenGL version)
/*/ static uint8_t   encode8_snorm_gl2(float   x);
/*/ 8-bit byte to float [-1..1] (OpenGL version)
/*/ static float     decode8_snorm_gl2(uint8_t x);
/*/ float [0..1] to 10-bit short
/*/ static uint16_t  encode10_unorm(float    x);
/*/ 10-bit short to float [0..1]
/*/ static float     decode10_unorm(uint16_t x);
/*/ float [-1..1] to 10-bit short
/*/ static uint16_t  encode10_snorm(float    x);
/*/ 10-bit short to float [-1..1]
/*/ static float     decode10_snorm(uint16_t x);
/*/ float [0..1] to 7-bit byte
/*/ static uint8_t   encode7_unorm(float    x);
/*/ 7-bit byte to float [0..1]
/*/ static float     decode7_unorm(uint8_t  x);
/*/ float [-1..1] to 7-bit byte
/*/ static uint8_t   encode7_snorm(float    x);
/*/ 7-bit byte to float [-1..1]
/*/ static float     decode7_snorm(uint8_t  x);
/*/ float [0..1] to 8-bit byte (alt version)
/*/ static uint8_t   encode8_unorm_alt(float    x);
/*/ 8-bit byte to float [0..1] (alt version)
/*/ static float     decode8_unorm_alt(uint8_t  x);
/*/ float [-1..1] to 8-bit byte (alt version)
/*/ static uint8_t   encode8_snorm_alt(float    x);
/*/ 8-bit byte to float [-1..1] (alt version)
/*/ static float     decode8_snorm_alt(uint8_t  x);
/*/ quaternion to 32-bit integer
/*/ static void      encode101010_quat( uint32_t &out, float a, float b, float c, float d );
/*/ 32-bit integer to quanternion
/*/ static void      decode101010_quat( float &a, float &b, float &c, float &d, uint32_t in );
/*/ vector3 to 32-bit integer (larger distance)
/*/ static void      encode7716_vec( uint32_t &out, float x, float y, float z );
/*/ 32-bit integer to vector3 (larger distance)
/*/ static void      decode7716_vec( float &x, float &y, float &z, uint32_t in );
/*/ vector3 to 32-bit integer (better orientation)
/*/ static void      encode8814_vec( uint32_t &out, float x, float y, float z );
/*/ 32-bit integer to vector3 (better orientation)
/*/ static void      decode8814_vec( float &x, float &y, float &z, uint32_t in );
/*/ quaternion to 32-bit integer (struct version)
/*/ template<typename T> static void encode101010_quat( uint32_t &out, const T &q );
/*/ 32-bit integer to quanternion (struct version)
/*/ template<typename T> static void decode101010_quat( T &q, uint32_t in );
/*/ vector3 to 32-bit integer (struct version)(larger distance)
/*/ template<typename T> static void encode7716_vec( uint32_t &out, const T &v );
/*/ 32-bit integer to vector3 (struct version)(larger distance)
/*/ template<typename T> static void decode7716_vec( T &v, uint32_t in );
/*/ vector3 to 32-bit integer (struct version)(better orientation)
/*/ template<typename T> static void encode8814_vec( uint32_t &out, const T &v );
/*/ 32-bit integer to vector3 (struct version)(better orientation)
/*/ template<typename T> static void decode8814_vec( T &v, uint32_t in );
/*/ helper function. reverse/inverse square root (-DQUANT_USE_STD_SQRT to use standard sqrt() instead)
/*/ static float     rsqrt( float number );
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
 
static uint16_t encode16_float( float fl ) {
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
 
static float decode16_float( uint16_t half ) {
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
// Custom: 10-bit based on D3D10_UNORM
// Used to de/quantize quat rotations (see de/encode101010_quat functions below)
static uint16_t  encode10_unorm(float    x) { return uint16_t( int (x * 1023.f + 0.5f) ); }
static float     decode10_unorm(uint16_t x) { return x / 1023.f; }
static uint16_t  encode10_snorm(float    x) { return encode10_unorm( ( x + 1 ) * 0.5f ); }
static float     decode10_snorm(uint16_t x) { return ( decode10_unorm( x ) * 2 ) - 1 ; }

//
// Helper: 7-bit, 8-bit based on D3D10_UNORM
// Used to de/quantize vec3 positions (see de/encode7716_vec && de/encode8814_vec functions below)
static uint8_t   encode7_unorm(float    x) { return uint8_t( int (x * 127.f + 0.5f) ); }
static float     decode7_unorm(uint8_t  x) { return x / 127.f; }
static uint8_t   encode7_snorm(float    x) { return encode7_unorm( ( x + 1 ) * 0.5f ); }
static float     decode7_snorm(uint8_t  x) { return ( decode7_unorm( x ) * 2 ) - 1 ; }

static uint8_t   encode8_unorm_alt(float    x) { return uint8_t( int (x * 255.f + 0.5f) ); }
static float     decode8_unorm_alt(uint8_t  x) { return x / 255.f; }
static uint8_t   encode8_snorm_alt(float    x) { return encode8_unorm_alt( ( x + 1 ) * 0.5f ); }
static float     decode8_snorm_alt(uint8_t  x) { return ( decode8_unorm_alt( x ) * 2 ) - 1 ; }

//
// Quantize rotation as 10+10+10 bits for quaternion components

static void encode101010_quat( uint32_t &out, float a, float b, float c, float d ) {
    // [ref] http://bitsquid.blogspot.com.es/2009/11/bitsquid-low-level-animation-system.html 
    // "For quaternions we use 2 bits to store the index of the largest component,
    // then 10 bits each to store the value of the remaining three components. We use the knowledge
    // that 1 = x^2 + y^2 + z^2 + w^2 to restore the largest component, so we don't actually have to store its value.
    // Since we don't store the largest component we know that the remaining ones must be in the range (-1/sqrt(2), 1/sqrt(2))
    // (otherwise, one of them would be largest). So we use the 10 bits to quantize a value in that range, giving us a precision of 0.0014.
    // The quaternions (x, y, z, w) and (-x, -y, -z, -w) represent the same rotation,
    // so I flip the signs so that the largest component is always positive." - Niklas Frykholm / bitsquid.se
    static const float rmin = -rsqrt(2), rmax = rsqrt(2);
    float aa = a*a, bb = b*b, cc = c*c, dd = d*d;
    /**/ if( aa >= bb && aa >= cc && aa >= dd ) {
		b = remap( b, rmin, rmax, -1, 1 );
		c = remap( c, rmin, rmax, -1, 1 );
		d = remap( d, rmin, rmax, -1, 1 );
        out = a >= 0 ? uint32_t( (0<<30) | (encode10_snorm( b)<<20) | (encode10_snorm( c)<<10) | (encode10_snorm( d)<<0) )
                     : uint32_t( (0<<30) | (encode10_snorm(-b)<<20) | (encode10_snorm(-c)<<10) | (encode10_snorm(-d)<<0) );
    }
    else if( bb >= cc && bb >= dd ) {
		a = remap( a, rmin, rmax, -1, 1 );
		c = remap( c, rmin, rmax, -1, 1 );
		d = remap( d, rmin, rmax, -1, 1 );
        out = b >= 0 ? uint32_t( (1<<30) | (encode10_snorm( a)<<20) | (encode10_snorm( c)<<10) | (encode10_snorm( d)<<0) )
                     : uint32_t( (1<<30) | (encode10_snorm(-a)<<20) | (encode10_snorm(-c)<<10) | (encode10_snorm(-d)<<0) );
    }
    else if( cc >= dd ) {
		a = remap( a, rmin, rmax, -1, 1 );
		b = remap( b, rmin, rmax, -1, 1 );
		d = remap( d, rmin, rmax, -1, 1 );
        out = c >= 0 ? uint32_t( (2<<30) | (encode10_snorm( a)<<20) | (encode10_snorm( b)<<10) | (encode10_snorm( d)<<0) )
                     : uint32_t( (2<<30) | (encode10_snorm(-a)<<20) | (encode10_snorm(-b)<<10) | (encode10_snorm(-d)<<0) );
    }
    else {
		a = remap( a, rmin, rmax, -1, 1 );
		b = remap( b, rmin, rmax, -1, 1 );
		c = remap( c, rmin, rmax, -1, 1 );
        out = d >= 0 ? uint32_t( (3<<30) | (encode10_snorm( a)<<20) | (encode10_snorm( b)<<10) | (encode10_snorm( c)<<0) )
                     : uint32_t( (3<<30) | (encode10_snorm(-a)<<20) | (encode10_snorm(-b)<<10) | (encode10_snorm(-c)<<0) );
    }
}

//
// DeQuantize rotation as 10+10+10 bits for quaternion components

static void decode101010_quat( float &a, float &b, float &c, float &d, uint32_t in ) {
    // [ref] http://bitsquid.blogspot.com.es/2009/11/bitsquid-low-level-animation-system.html 
    // See encode101010_quat() function above.
    static const float rmin = -rsqrt(2), rmax = rsqrt(2);
    switch( in >> 30 ) {
        default: case 0:
            b = decode10_snorm( ( in >> 20 ) & 0x3FF );
            c = decode10_snorm( ( in >> 10 ) & 0x3FF );
            d = decode10_snorm( ( in >>  0 ) & 0x3FF );
			b = remap( b, -1, 1, rmin, rmax );
			c = remap( c, -1, 1, rmin, rmax );
			d = remap( d, -1, 1, rmin, rmax );
			a = 1 / rsqrt( 1 - b * b - c * c - d * d );
        break; case 1:
            a = decode10_snorm( ( in >> 20 ) & 0x3FF );
            c = decode10_snorm( ( in >> 10 ) & 0x3FF );
            d = decode10_snorm( ( in >>  0 ) & 0x3FF );
			a = remap( a, -1, 1, rmin, rmax );
			c = remap( c, -1, 1, rmin, rmax );
			d = remap( d, -1, 1, rmin, rmax );
            b = 1 / rsqrt( 1 - a * a - c * c - d * d );
        break; case 2:
            a = decode10_snorm( ( in >> 20 ) & 0x3FF );
            b = decode10_snorm( ( in >> 10 ) & 0x3FF );
            d = decode10_snorm( ( in >>  0 ) & 0x3FF );
			a = remap( a, -1, 1, rmin, rmax );
			b = remap( b, -1, 1, rmin, rmax );
			d = remap( d, -1, 1, rmin, rmax );
            c = 1 / rsqrt( 1 - a * a - b * b - d * d );
        break; case 3:
            a = decode10_snorm( ( in >> 20 ) & 0x3FF );
            b = decode10_snorm( ( in >> 10 ) & 0x3FF );
            c = decode10_snorm( ( in >>  0 ) & 0x3FF );
			a = remap( a, -1, 1, rmin, rmax );
			b = remap( b, -1, 1, rmin, rmax );
			c = remap( c, -1, 1, rmin, rmax );
            d = 1 / rsqrt( 1 - a * a - b * b - c * c );
    }
}

//
// Quantize scale/position as 7+7 bits for unit vector rotation + 16 bits for vector length

static void encode7716_vec( uint32_t &out, float x, float y, float z ) { 
    // Somehow similar to encode101010_quat() function above.
    // We decompose given vector into unit vector and length (magnitude). Then we discard the largest component, as unit vectors
    // follow 1 = x^2 + y^2 + z^2 expression (similar to encode101010_quat() function above). The unit vectors (x, y, z, w) and
    // (-x, -y, -z) represent the same direction, so I flip the magnitude so that the largest component is always positive.
    float xx = x*x;
    float yy = y*y;
    float zz = z*z;
    float len = rsqrt( xx + yy + zz ); // float len = sqrt( xx + yy + zz );
    x *= len; y *= len; z *= len;      // x /= len; y /= len; z /= len;
    /****/ if( xx >= yy && xx >= zz ) {
        out = ( x >= 0 ? uint32_t( (0<<30)|(encode7_snorm( y)<<23)|(encode7_snorm( z)<<16)|(encode16_float( 1/len)) ) :  // ( len)
                         uint32_t( (0<<30)|(encode7_snorm(-y)<<23)|(encode7_snorm(-z)<<16)|(encode16_float(-1/len)) ) ); // (-len)
    } else if( yy >= zz ) {
        out = ( y >= 0 ? uint32_t( (1<<30)|(encode7_snorm( x)<<23)|(encode7_snorm( z)<<16)|(encode16_float( 1/len)) ) :  // ( len)
                         uint32_t( (1<<30)|(encode7_snorm(-x)<<23)|(encode7_snorm(-z)<<16)|(encode16_float(-1/len)) ) ); // (-len)
    } else {
        out = ( z >= 0 ? uint32_t( (2<<30)|(encode7_snorm( x)<<23)|(encode7_snorm( y)<<16)|(encode16_float( 1/len)) ) :  // ( len)
                         uint32_t( (2<<30)|(encode7_snorm(-x)<<23)|(encode7_snorm(-y)<<16)|(encode16_float(-1/len)) ) ); // (-len)
    }
}

//
// DeQuantize scale/position as 7+7 bits for unit vector rotation + 16 bits for vector length

static void decode7716_vec( float &x, float &y, float &z, uint32_t in ) { 
    // See encode7716_vec() function above.
    switch( in >> 30 ) {
        default: case 0:
            y = decode7_snorm( ( in >> 23 ) & 0x7F );
            z = decode7_snorm( ( in >> 16 ) & 0x7F );
            x = 1 / rsqrt( 1 - y * y - z * z ); // x = sqrt( 1 - y * y - z * z );
        break; case 1:
            x = decode7_snorm( ( in >> 23 ) & 0x7F );
            z = decode7_snorm( ( in >> 16 ) & 0x7F );
            y = 1 / rsqrt( 1 - x * x - z * z ); // y = sqrt( 1 - x * x - z * z );
        break; case 2:
            x = decode7_snorm( ( in >> 23 ) & 0x7F );
            y = decode7_snorm( ( in >> 16 ) & 0x7F );
            z = 1 / rsqrt( 1 - x * x - y * y ); // z = sqrt( 1 - x * x - y * y );
    }
    float len = decode16_float( uint16_t(in & 0xFFFF) ); 
    x *= len, y *= len, z *= len;
}

//
// Quantize scale/position as 8+8 bits for unit vector rotation + 14 bits for vector length

static void encode8814_vec( uint32_t &out, float x, float y, float z ) { 
    // Similar to encode7716_vec() function above.
    // We degrade precision while representing distance from origin in favor of a more precise direction vector.
    // We move the 2-bit selector into the less significant bits of the mantissa in the half-float (that holds the vector length).
    float xx = x*x;
    float yy = y*y;
    float zz = z*z;
    float len = rsqrt( xx + yy + zz );  // float len = sqrt( xx + yy + zz );
    x *= len; y *= len; z *= len;       // x /= len; y /= len; z /= len;
    /****/ if( xx >= yy && xx >= zz ) {
        out = ( x >= 0 ? uint32_t( (encode8_snorm_alt( y)<<24)|(encode8_snorm_alt( z)<<16)|(encode16_float( 1/len)&0xFFFC)|(0) ) :  // ( len)
                         uint32_t( (encode8_snorm_alt(-y)<<24)|(encode8_snorm_alt(-z)<<16)|(encode16_float(-1/len)&0xFFFC)|(0) ) ); // (-len)
    } else if( yy >= zz ) {
        out = ( y >= 0 ? uint32_t( (encode8_snorm_alt( x)<<24)|(encode8_snorm_alt( z)<<16)|(encode16_float( 1/len)&0xFFFC)|(1) ) :  // ( len)
                         uint32_t( (encode8_snorm_alt(-x)<<24)|(encode8_snorm_alt(-z)<<16)|(encode16_float(-1/len)&0xFFFC)|(1) ) ); // (-len)
    } else {
        out = ( z >= 0 ? uint32_t( (encode8_snorm_alt( x)<<24)|(encode8_snorm_alt( y)<<16)|(encode16_float( 1/len)&0xFFFC)|(2) ) :  // ( len)
                         uint32_t( (encode8_snorm_alt(-x)<<24)|(encode8_snorm_alt(-y)<<16)|(encode16_float(-1/len)&0xFFFC)|(2) ) ); // (-len)
    }
}

//
// DeQuantize scale/position as 8+8 bits for unit vector rotation + 14 bits for vector length

static void decode8814_vec( float &x, float &y, float &z, uint32_t in ) { 
    // See encode8814_vec() function above.
    switch( in & 3 ) {
        default: case 0:
            y = decode8_snorm_alt( ( in >> 24 ) & 0xFF );
            z = decode8_snorm_alt( ( in >> 16 ) & 0xFF );
            x = 1 / rsqrt( 1 - y * y - z * z ); // x = sqrt( 1 - y * y - z * z );
        break; case 1:
            x = decode8_snorm_alt( ( in >> 24 ) & 0xFF );
            z = decode8_snorm_alt( ( in >> 16 ) & 0xFF );
            y = 1 / rsqrt( 1 - x * x - z * z ); // y = sqrt( 1 - x * x - z * z );
        break; case 2:
            x = decode8_snorm_alt( ( in >> 24 ) & 0xFF );
            y = decode8_snorm_alt( ( in >> 16 ) & 0xFF );
            z = 1 / rsqrt( 1 - x * x - y * y ); // z = sqrt( 1 - x * x - y * y );
    }
    float len = decode16_float( uint16_t(in & 0xFFFC) ); 
    x *= len, y *= len, z *= len;
}

// in-place class sugars

template<typename quat>
static void encode101010_quat( uint32_t &out, const quat &q ) {
    encode101010_quat( out, q.x, q.y, q.z, q.w );
}
template<typename quat>
static void decode101010_quat( quat &q, uint32_t in ) {
    decode101010_quat( q.x, q.y, q.z, q.w, in );
}

template<typename vec3>
static void encode7716_vec( uint32_t &out, const vec3 &v ) { 
    encode7716_vec( out, v.x, v.y, v.z );
}
template<typename vec3>
static void decode7716_vec( vec3 &v, uint32_t in ) { 
    decode7716_vec( v.x, v.y, v.z, in );
}

template<typename vec3>
static void encode8814_vec( uint32_t &out, const vec3 &v ) { 
    encode8814_vec( out, v.x, v.y, v.z );
}
template<typename vec3>
static void decode8814_vec( vec3 &v, uint32_t in ) { 
    decode8814_vec( v.x, v.y, v.z, in );
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

    std::cout << decode16_float( encode16_float(3.14159f) ) << std::endl;
}
 
#endif
