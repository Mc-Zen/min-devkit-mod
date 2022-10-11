
// Optimized mathematical functions

#pragma once
#include <concepts>


namespace Butterfly {

template<std::integral T>
constexpr bool isPowerOf2(T x) {
	return !(x & (x - 1)) && x > 0;
}


//
// template<std::integral T> requires(sizeof(T) == 4)
// T log2OfPowerOf2(T x) {
//	assert((x & (x - 1)) == 0 && x > 0 && "Size N has to be a power of 2 and must not be 0");
//	int p = 0;
//	while (x > 1) {
//		x >>= 1;
//		++p;
//	}
//	return p;
//}

namespace detail {
struct coeffs
{
	static constexpr int deBruijnBitPosition2[32] = {
		0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
		31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
	};
};
}

template<std::integral T>
requires(sizeof(T) == 4) constexpr T log2OfPowerOf2(T x) {
	assert(isPowerOf2(x) && "This function is only designed for numbers which are a power of 2");
	// static const int deBruijnBitPosition2[32] = {
	//	0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
	//	31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
	// };
	return detail::coeffs::deBruijnBitPosition2[(uint32_t)(x * 0x077CB531U) >> 27];
}


template<std::integral T>
requires(sizeof(T) == 4) constexpr T bitReverse(T x, int nb) {
	assert(nb > 0 && 32 > nb && "invalid bit count");
	x = (x << 16) | (x >> 16);
	x = ((x & 0x00FF00FF) << 8) | ((x & 0xFF00FF00) >> 8);
	x = ((x & 0x0F0F0F0F) << 4) | ((x & 0xF0F0F0F0) >> 4);
	x = ((x & 0x33333333) << 2) | ((x & 0xCCCCCCCC) >> 2);
	x = ((x & 0x55555555) << 1) | ((x & 0xAAAAAAAA) >> 1);

	return ((x >> (32 - nb)) & (0xFFFFFFFF >> (32 - nb)));
}


}
