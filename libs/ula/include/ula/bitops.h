#ifndef ULA_BITOPS_H
#define ULA_BITOPS_H

//#define WE_HAVE_64_BITS

#ifdef cray

#include <intrinsics.h>
#define popcnt _popcnt
#define gbit _gbit
#define gbits _gbits
#define maskr _maskr

#else

#include <ula/bittypes.h>

/** extract a bit from a word.
   @param x the word
   @param n position of the bit to be extracted
   @return the n-th bit of the word x */
inline uint32 gbit(int32 x, int32 n) { return (x>>n)&(uint32(1)); }
inline uint64 gbit(uint64 x,int32 n) {return (x>>n)&(uint64(1));}

/** extract bits from a word.
   @param x the word
   @param m the number of bits to be extracted
   @param n position of the bit to be extracted
   @return the m bits starting at bit n  */
inline uint32 gbits(int32 x, int32 m, int32 n) { return (x>>n)&((uint32(1)<<m)-1); }
inline uint64 gbits(uint64 x, int32 m, int32 n) { return (x>>n)&((uint64(1)<<m)-1); }

/** create a right-justified mask
    creates a bitpattern with the rightmost i bits set
    @param i the number of bits to be set to 1
    @return a word with the rightmost bits set to 1 */
inline uint32 maskr(int32 i) {return (int32(1)<<i)-1;}
inline uint64 maskr(uint64 i) {return (uint64(1)<<i)-1;}

inline int32 BX_(int32 x) { 
  return ((x) 
          - (((x)>>1)&0x77777777)
          - (((x)>>2)&0x33333333)
          - (((x)>>3)&0x11111111));
}

inline uint64 BX_(uint64 x) {
#ifdef WE_HAVE_64_BITS
  return (
	  (x)
	  -(((x)>>1)&uint64(0x7777777777777777))
	  -(((x)>>2)&uint64(0x3333333333333333))
	  -(((x)>>3)&uint64(0x1111111111111111)));
#else
  return ((x) 
          - (((x)>>1)&0x77777777)
          - (((x)>>2)&0x33333333)
          - (((x)>>3)&0x11111111));
#endif
}

inline uint32 popcnt(int32 x) 
{ return    (((BX_(x)+(BX_(x)>>4)) & 0x0F0F0F0F) % 255);}

inline uint32 popcnt(uint64 x)
{ 
#ifdef WE_HAVE_64_BITS
  return    (((BX_(x)+(BX_(x)>>4)) & uint64(0x0F0F0F0F0F0F0F0F)) % 255);
#else
  return    (((BX_(x)+(BX_(x)>>4)) & 0x0F0F0F0F) % 255);
#endif
}

#endif

#endif
