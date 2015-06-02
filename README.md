# Quant :fish_cake: <a href="https://travis-ci.org/r-lyeh/quant"><img src="https://api.travis-ci.org/r-lyeh/quant.svg?branch=master" align="right" /></a>
- Quant is a quantization suite supporting conversion to/from halves, s/unorm bytes, quats and vec3s.
- Quant is tiny, header-only, cross-platform.
- Quant is public domain.

## Motivation
- Support a few de-facto standards: s/unorm (both d3d10/opengl), half-floats.
- Support quaternions and vectors for both position and scales.
- Keep conversions as cross-platform and architecture-friendly as possible.
- Keep good visual quality while aiming to smallest types.
- Pack a standard 3D transform matrix (matrix4x4f 64 bytes) into 12 bytes (3x uint32_t).

## Usage
- Quantize animations.
- Quantize colors
- Quantize sounds.
- Quantize user input.
- Quantize network packets.
- Etc

## Todos
- Integrate (and finish) curve simplification and hermite splines that I have somewhere around.
- Lossless/lossy animation format proposal:
  - demultiplex vertex streams to a single giant (mono) stream
  - quantize whole stream `(iter - min) / ( max - min ) -> [0..1]`
  - apply lossless audio (FLAC) or lossy audio codec (mp3/ogg)
  - decode, upscale stream `iter * ( max - min ) + min`, and multiplex vertex
  - profit (?)

## API
```c++
namespace quant {
// float to 16-bit half
uint16_t  encode16_float( float fl );
// 16-bit half to float
float     decode16_float( uint16_t half );
// float [0..1] to 8-bit byte
uint8_t   encode8_unorm(float    x);
// 8-bit byte to float [0..1]
float     decode8_unorm(uint8_t  x);
// float [-1..1] to 8-bit byte
uint8_t   encode8_snorm(float    x);
// 8-bit byte to float [-1..1]
float     decode8_snorm(uint8_t  x);
// float [-1..1] to 8-bit byte (OpenGL version)
uint8_t   encode8_snorm_gl2(float   x);
// 8-bit byte to float [-1..1] (OpenGL version)
float     decode8_snorm_gl2(uint8_t x);
// float [0..1] to 10-bit short
uint16_t  encode10_unorm(float    x);
// 10-bit short to float [0..1]
float     decode10_unorm(uint16_t x);
// float [-1..1] to 10-bit short
uint16_t  encode10_snorm(float    x);
// 10-bit short to float [-1..1]
float     decode10_snorm(uint16_t x);
// float [0..1] to 7-bit byte
uint8_t   encode7_unorm(float    x);
// 7-bit byte to float [0..1]
float     decode7_unorm(uint8_t  x);
// float [-1..1] to 7-bit byte
uint8_t   encode7_snorm(float    x);
// 7-bit byte to float [-1..1]
float     decode7_snorm(uint8_t  x);
// float [0..1] to 8-bit byte (alt version)
uint8_t   encode8_unorm_alt(float    x);
// 8-bit byte to float [0..1] (alt version)
float     decode8_unorm_alt(uint8_t  x);
// float [-1..1] to 8-bit byte (alt version)
uint8_t   encode8_snorm_alt(float    x);
// 8-bit byte to float [-1..1] (alt version)
float     decode8_snorm_alt(uint8_t  x);
// quaternion to 32-bit integer
void      encode101010_quat( uint32_t &out, float a, float b, float c, float d );
// 32-bit integer to quanternion
void      decode101010_quat( float &a, float &b, float &c, float &d, uint32_t in );
// vector3 to 32-bit integer (larger distance)
void      encode7716_vec( uint32_t &out, float x, float y, float z );
// 32-bit integer to vector3 (larger distance)
void      decode7716_vec( float &x, float &y, float &z, uint32_t in );
// vector3 to 32-bit integer (better orientation)
void      encode8814_vec( uint32_t &out, float x, float y, float z );
// 32-bit integer to vector3 (better orientation)
void      decode8814_vec( float &x, float &y, float &z, uint32_t in );
// quaternion to 32-bit integer (struct version)
template <typename T> void encode101010_quat( uint32_t &out, const T &q );
// 32-bit integer to quanternion (struct version)
template <typename T> void decode101010_quat( T &q, uint32_t in );
// vector3 to 32-bit integer (struct version)(larger distance)
template <typename T> void encode7716_vec( uint32_t &out, const T &v );
// 32-bit integer to vector3 (struct version)(larger distance)
template <typename T> void decode7716_vec( T &v, uint32_t in );
// vector3 to 32-bit integer (struct version)(better orientation)
template <typename T> void encode8814_vec( uint32_t &out, const T &v );
// 32-bit integer to vector3 (struct version)(better orientation)
template <typename T> void decode8814_vec( T &v, uint32_t in );
// helper function. reverse/inverse square root (-DQUANT_USE_STD_SQRT to use standard sqrt() instead)
float     rsqrt( float number );
// helper function. remap floating number in range [min1..max1] to range [min2..max2]
float     remap( float x, float min1, float max1, float min2, float max2 );
}
```

## Showcase (demo not provided)
![image](https://raw.github.com/r-lyeh/depot/master/skull-quant.png)

## Changelog
- v1.0.1 (2015/06/02)
  - Improved quant precision
  - Fixed rotation bug
- v1.0.0 (2015/05/29)
  - Initial revision

## References and links
- https://gist.github.com/rygorous/2156668
- http://zeuxcg.org/2010/12/14/quantizing-floats/
- http://en.wikipedia.org/wiki/Fast_inverse_square_root
- http://bitsquid.blogspot.com.es/2009/11/bitsquid-low-level-animation-system.html 

## Licenses
- [Quant](https://github.com/r-lyeh/quant) (Public Domain).
- [float from/to half variants](https://gist.github.com/rygorous/2156668) by Fabian "ryg" Giesen (Public Domain).
