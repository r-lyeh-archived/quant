# Quant :fish_cake: <a href="https://travis-ci.org/r-lyeh/quant"><img src="https://api.travis-ci.org/r-lyeh/quant.svg?branch=master" align="right" /></a>
- Quant is a quantization suite supporting many targets and unit types.
- Quant is tiny, header-only, cross-platform.
- Quant is public domain.

## Features
- Support conversion from/to signed normalized floats to/from N-bits shorts (d3d10 variant).
- Support conversion from/to signed normalized floats to/from N-bits shorts (opengl variant).
- Support conversion from/to unsigned normalized floats to/from N-bits shorts.
- Support conversion from/to quaternions to/from 32-bits integers.
- Support conversion from/to positions to/from 48/32/16-bits integers.
- Support conversion from/to scales to/from 48/32/16-bits integers.
- Pack standard 3D transform matrix (matrix4x4f 64 bytes) into 10 or 12 bytes.
- Conversions done as cross-platform and architecture-friendly as possible.
- Good visual quality while aiming to smallest types.

## Usages
- To de/quantize animations.
- To de/quantize colors.
- To de/quantize sounds.
- To de/quantize user input.
- To de/quantize network packets.
- Etc

## API
```c++
namespace quant {
  // For generic floats
  uint16_t encode16_half(float);    // 16-bit
  float    decode16_half(uint16_t); // 16-bit

  // For signed normalized [-1..1] floats
  uint8_t encode8_snorm(float);     // 8-bit
  float   decode8_snorm(uint8_t);   // 8-bit

  // For unsigned normalized [0..1] floats
  uint8_t encode8_unorm(float);    // 8-bit
  float   decode8_unorm(uint8_t);  // 8-bit

  // For rotation quaternions
  encode101010_quant(uint32_t &q, float x, y, z, w);    // 32-bit
  decode101010_quant(float &x, &y, &z, &w, uint32_t q); // 32-bit

  // For position vectors
  encode161616_vec(uint64_t &q, float x, y, z);   // 48-bit version
  decode161616_vec(float &x, &y, &z, uint64_t q); // 48-bit version
  encode8814_vec(uint32_t &q, float x, y, z);     // 32-bit version
  decode8814_vec(float &x, &y, &z, uint32_t q);   // 32-bit version
  encode555_vec(uint16_t q, float x, y, z);       // 16-bit version
  decode555_vec(float &x, &y, &z, uint16_t q);    // 16-bit version

  // For scale vectors
  // Scale tip:
  //- Try to de/quantize scale vectors as `fn(q,1-x,1-y,1-z)` if possible, rather than `fn(q,x,y,z)`
  //- So, scales close to one will be numerically stabler (~less visual glitches at larger scales)
  //- So, scales close to zero will be numerically unstabler (~more visual glitches at smaller scales)
  encode8814_vec(uint32_t q, float x, y, z);      // 32-bit version
  decode8814_vec(float &x, &y, &z, uint32_t q);   // 32-bit version
  encode555_vec(uint16_t q, float x, y, z);       // 16-bit version
  decode555_vec(float &x, &y, &z, uint16_t q);    // 16-bit version
}
```
For more variants, please check [quant.hpp header](quant.hpp).

## Todos
- Finish (and integrate) curve simplification and hermite splines that I have somewhere lying around.
- Lossless/lossy animation format proposal:
  - demultiplex vertex streams to a single giant (mono) stream
  - quantize whole stream `(iter - min) / ( max - min ) -> [0..1]`
  - apply lossless audio (FLAC) or lossy audio codec (mp3/ogg)
  - decode, upscale stream `iter * ( max - min ) + min`, and multiplex vertex
  - profit (?)

## Showcase (demo not provided)
![image](https://raw.github.com/r-lyeh/depot/master/skull-quant.png)

## Changelog
- v1.0.2 (2015/06/05)
  - Expanded api: 48-bit and 16-bit vector support.
  - Improved numerical stability
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
