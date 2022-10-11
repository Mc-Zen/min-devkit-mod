
// EXPERIMENTAL:
//
// Tiny library for fixed point arithmetic and modular arithmetic for
// arbitrary ranges of fixed floats


#pragma once

#include <cstdint>


namespace Butterfly {


template<int size>
struct FixedPointSizeType
{
};
template<>
struct FixedPointSizeType<8>
{
	using type = uint8_t;
};
template<>
struct FixedPointSizeType<16>
{
	using type = uint16_t;
};
template<>
struct FixedPointSizeType<32>
{
	using type = uint32_t;
};
template<>
struct FixedPointSizeType<64>
{
	using type = uint64_t;
};


/// @brief Wrapper for fixed point types
///
/// @tparam size        size of underlying unsigned integer type (needs to be a power of 2)
/// @tparam integerBits number of bits to assign to the integer part of the fixed float
template<int size, int integerBits>
class Fixed
{
	static_assert((size & (size - 1)) == 0 && size > 0, "Size needs to be a power of 2");
	static_assert(integerBits <= size, "The number of integer bits needs to be less then the size");

public:
	using int_type = typename FixedPointSizeType<size>::type;

	constexpr Fixed() = default;

	template<std::floating_point T>
	constexpr Fixed(T value) : value(static_cast<int_type>(std::round(value * fixedPointMultiplicator))) {}


	template<class FP>
	constexpr FP to() const { return static_cast<FP>(value * fixedPointMultiplicatorInv); }

	constexpr int_type integer() const { return value >> fixedPointFractionalBits; }
	constexpr double fractional() const { return (value & fractionalMask) * fixedPointMultiplicatorInv; }


	constexpr Fixed operator+(Fixed a) const { return Fixed{ value + a.value }; }
	constexpr Fixed operator-(Fixed a) const { return Fixed{ value - a.value }; }
	constexpr Fixed operator<<(Fixed a) const { return Fixed{ value << a.value }; }
	constexpr Fixed operator>>(Fixed a) const { return Fixed{ value >> a.value }; }
	constexpr Fixed operator%(Fixed a) const { return Fixed{ value % a.value }; }
	constexpr bool operator==(Fixed a) const { return value == a.value; }
	constexpr Fixed& operator+=(Fixed a) {
		value += a.value;
		return *this;
	}
	constexpr Fixed& operator-=(Fixed a) {
		value -= a.value;
		return *this;
	}

	template<std::integral T>
	constexpr Fixed& operator*=(T a) { value *= a; return *this; }
	template<std::integral T>
	constexpr Fixed& operator/=(T a) { value /= a; return *this; }
	template<std::integral T>
	constexpr Fixed operator*(T a) const { auto tmp = *this; return tmp *= a; }
	template<std::integral T>
	constexpr Fixed operator/(T a) const { auto tmp = *this; return tmp /= a; }
	template<std::floating_point T>
	constexpr Fixed& operator*=(T a) { value *= a; return *this; }
	template<std::floating_point T>
	constexpr Fixed& operator/=(T a) { value /= a; return *this; }
	template<std::floating_point T>
	constexpr Fixed operator*(T a) const { auto tmp = *this; return tmp *= a; }
	template<std::floating_point T>
	constexpr Fixed operator/(T a) const { auto tmp = *this; return tmp /= a; }

private:
	constexpr Fixed(int_type value) : value(value) {}
	int_type value{};
	static constexpr int_type fixedPointFractionalBits = 8 * sizeof(int_type) - integerBits;
	static constexpr double fixedPointMultiplicator = integerBits != 0 ? (int_type(1) << fixedPointFractionalBits) : (int_type(1) << (fixedPointFractionalBits - 1)) * 2.;
	static constexpr double fixedPointMultiplicatorInv = 1.0 / fixedPointMultiplicator;
	static constexpr int_type fractionalMask = (int_type(1) << fixedPointFractionalBits) - 1;
};


// Wrapper for real numbers (fixed point) with modular arithmetic.
// An interval [0,max] or [0,max) is defined which maps the range of
// an unsigned

/// @brief Wrapper for real numbers (fixed point) with modular arithmetic.
///        An interval [0, max] or [0, max) is defined which maps the range of
///        an unsigned integer to the specified range.
///
/// @tparam size        size of the underlying type, may be 8, 16, 32 or 64
/// @tparam max         max of the range
/// @tparam maxExcluded determines if the maximum value may be reached
template<int size, double max, bool maxExcluded = true>
class Wrapping_Fixed
{
	static_assert((size & (size - 1)) == 0 && size > 0, "size needs to be a power of 2");


public:
	using int_type = typename FixedPointSizeType<size>::type;


	template<class FP>
	constexpr Wrapping_Fixed(FP value)
		: value(static_cast<int_type>(std::round(value * scale))) {}

	template<class FP>
	constexpr FP to() const { return static_cast<FP>(value * scaleInv); }

	Wrapping_Fixed& operator+=(Wrapping_Fixed a) {
		value += a.value;
		return *this;
	}
	Wrapping_Fixed& operator-=(Wrapping_Fixed a) {
		value -= a.value;
		return *this;
	}
	Wrapping_Fixed operator+(Wrapping_Fixed a) { return { value + a.value }; }
	Wrapping_Fixed operator-(Wrapping_Fixed a) { return { value - a.value }; }
	Wrapping_Fixed operator<<(Wrapping_Fixed a) { return { value << a.value }; }
	Wrapping_Fixed operator>>(Wrapping_Fixed a) { return { value >> a.value }; }
	Wrapping_Fixed operator%(Wrapping_Fixed a) { return { value % a.value }; }
	Wrapping_Fixed operator==(Wrapping_Fixed a) const { return value == a.value; }

	template<std::floating_point T>
	Wrapping_Fixed operator*(T a) { return { static_cast<int_type>(value * a) }; }

	template<std::floating_point T>
	Wrapping_Fixed operator/(T a) { return { static_cast<int_type>(value / a) }; }

	template<std::integral T>
	Wrapping_Fixed operator*(T a) { return { static_cast<int_type>(value * a) }; }

	template<std::integral T>
	Wrapping_Fixed operator/(T a) { return { static_cast<int_type>(value / a) }; }

private:
	Wrapping_Fixed(int_type value) : value(value) {}
	int_type value{};

	static constexpr double scale = static_cast<double>(std::numeric_limits<int_type>::max()) / max + (maxExcluded ? (1. / max) : 0.);
	static constexpr double scaleInv = 1.0 / scale;
	static_assert(!maxExcluded || static_cast<double>(std::numeric_limits<int_type>::max()) * scaleInv < max, "Internal error: range max is not excluded");
	static_assert(maxExcluded || static_cast<double>(std::numeric_limits<int_type>::max()) * scaleInv == max, "Internal error: range max is not included");
};







}
