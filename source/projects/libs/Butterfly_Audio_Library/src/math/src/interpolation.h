
// Interpolation algorithms and interpolation structs with static members providing
// access to interpolation method meta-information and the algorithm itself.


#pragma once

#include <concepts>

namespace Butterfly {


/// @brief Linear interpolation between values y0 and y1
///
/// @tparam T Value type
/// @param t Interpolation parameter in the interval [0, 1]
/// @param y0 First value
/// @param y1 Second value
/// @return interpolated value
template<class T>
constexpr T linear_interpolation(const T& t, const T& y0, const T& y1) {
	return y0 + t * (y1 - y0);
}


/// @brief 3rd-order hermite interpolation between four values y-1, y0, y1 and y2.
///        The interpolation is performed between y0 and y1.
///
/// @tparam T Value type
/// @param t Interpolation parameter in the interval [0, 1]
/// @param ym1 First value
/// @param y0 Second value
/// @param y1 Third value
/// @param y2 Fourth value
/// @return interpolated value
template<class T>
constexpr T hermite_interpolation(const T& t, const T& ym1, const T& y0, const T& y1, const T& y2) {
	const auto c0 = y0;
	const auto c1 = T(0.5) * (y1 - ym1);
	const auto c2 = ym1 - T(2.5) * y0 + T(2.) * y1 - T(0.5) * y2;
	const auto c3 = T(1.5) * (y0 - y1) + T(0.5) * (y2 - ym1);

	return ((c3 * t + c2) * t + c1) * t + c0;
}


/// @brief Cubic interpolation between four values y-1, y0, y1 and y2.
///        The interpolation is performed between y0 and y1.
///
/// @tparam T Value type
/// @param t Interpolation parameter in the interval [0, 1]
/// @param ym1 First value
/// @param y0 Second value
/// @param y1 Third value
/// @param y2 Fourth value
/// @return interpolated value
template<class T>
constexpr T cubic_interpolation(const T& t, const T& ym1, const T& y0, const T& y1, const T& y2) {
	const auto c3 = y2 - y1 + y0 - ym1;
	const auto c2 = ym1 - y0 - c3;
	const auto c1 = y1 - ym1;
	const auto c0 = y0;
	const auto t_sq = t * t;

	return ((c3 * t + c2) * t + c1) * t + c0;
}

/// @brief Evaluation of cubic a Bézier function with knots x0 through x3. 
///	      Also, some sort of point class can be passed in that supports basic arithmetic. 
/// 
/// @tparam T Value type
/// @param t Curve parameter in the interval [0, 1]
/// @param x0 Start point
/// @param x1 First handle
/// @param x2 Second handle
/// @param x3 End point
/// @return evaluated value
template<class T>
constexpr T cubic_bezier(const T& t, const T& x0, const T& x1, const T& x2, const T& x3) {
	const auto x0i = linear_interpolation(t, x0, x1);
	const auto x1i = linear_interpolation(t, x1, x2);
	const auto x2i = linear_interpolation(t, x2, x3);
	const auto x0ii = linear_interpolation(t, x0i, x1i);
	const auto x1ii = linear_interpolation(t, x1i, x2i);
	return linear_interpolation(t, x0ii, x1ii);
}

//template<class T>
//constexpr T bezier1(T t, T x0, T x1, T x2, T x3) {
//	// 12 Additions, 6 Multiplications
//	const auto x0i = x0 + t * (x1 - x0);
//	const auto x1i = x1 + t * (x2 - x1);
//	const auto x2i = x2 + t * (x3 - x2);
//	const auto x0ii = x0i + t * (x1i - x0i);
//	const auto x1ii = x1i + t * (x2i - x1i);
//	return x0ii + t * (x1ii - x0ii);
//}
//
// actually slower!
//template<class T>
//constexpr T bezier2(T t, T x0, T x1, T x2, T x3) {
//	// 12 Additions, 6 Multiplications
//	const auto c0 = x1 - x0;
//	const auto c1 = c0 + t * (x2 - x1 - c0);
//	const auto c2 = x2 - x1 + t * (x3 - 2 * x2 + x1); // x2 - 2 * x1 + x0;
//	return x0 + t * (c0 + 2 * c1 + t * (c2 - c1));
//}


template<class T, class U>
concept interpolator = requires(U* data, size_t index, U offset) {
	{ T::interpolate(data, index, offset) } -> std::convertible_to<U>;
	{ T::getLookbehindLength() } -> std::convertible_to<size_t>;
	{ T::getLookaheadLength() } -> std::convertible_to<size_t>;
};


struct LinearInterpolator
{
	template<class RandomAccessContainer, class T>
	static constexpr T interpolate(const RandomAccessContainer& data, size_t index, T offset) {
		static_assert(std::is_same_v<typename RandomAccessContainer::value_type, T>);
		return linear_interpolation(offset, data[index], data[index + 1]);
	}
	static constexpr size_t getLookbehindLength() { return 0; }
	static constexpr size_t getLookaheadLength() { return 1; }
};

struct HermiteInterpolator
{
	template<class RandomAccessContainer, class T>
	static constexpr T interpolate(const RandomAccessContainer& data, size_t index, T offset) {
		static_assert(std::is_same_v<typename RandomAccessContainer::value_type, T>);
		return hermite_interpolation(offset, data[index - 1], data[index], data[index + 1], data[index + 2]);
	}
	static constexpr size_t getLookbehindLength() { return 1; }
	static constexpr size_t getLookaheadLength() { return 2; }
};

struct CubicInterpolator
{
	template<class RandomAccessContainer, class T>
	static constexpr T interpolate(const RandomAccessContainer& data, size_t index, T offset) {
		static_assert(std::is_same_v<typename RandomAccessContainer::value_type, T>);
		return cubic_interpolation(offset, data[index - 1], data[index], data[index + 1], data[index + 2]);
	}
	static constexpr size_t getLookbehindLength() { return 1; }
	static constexpr size_t getLookaheadLength() { return 2; }
};


}
