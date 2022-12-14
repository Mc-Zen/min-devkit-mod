
// Routines for computing the fast fourier transform and its inverse.

#pragma once

#include <cmath>
#include <cstdint>
#include <cassert>
#include <array>
#include <vector>
#include <complex>
#include <numbers>
#include "optimized_math.h"

namespace Butterfly {

namespace detail {

enum class FFTDirection {
	Forward = 1,
	Backward = -1
};

template<std::random_access_iterator InIt, std::random_access_iterator OutIt>
constexpr void fft(InIt in_first, InIt in_last, OutIt out, FFTDirection dir) {
	using complex = typename OutIt::value_type;
	using T = typename complex::value_type;

	const auto size = std::distance(in_first, in_last);
	assert(isPowerOf2(size) && "Input size needs to be a power of 2");

	const int log_n = log2OfPowerOf2(static_cast<int>(size));
	const T nrm = T(1) / std::sqrt(T(size));

	for (unsigned int i = 0; i < size; ++i) {
		out[i] = nrm * in_first[bitReverse(i, log_n)];
	}

	for (int s = 0; s < log_n; ++s) {
		const int m2 = 1 << s; // 2^(s-1)
		const int m = m2 << 1; // 2^s
		std::complex<T> w{ 1, 0 };

		const auto wm = std::polar(T(1), T(dir) * std::numbers::pi_v<T> / T(m2));
		for (int j = 0; j < m2; ++j) {
			for (int k = j; k < size; k += m) {
				const auto t = w * out[k + m2];
				const auto u = out[k];
				out[k] = u + t;
				out[k + m2] = u - t;
			}
			w *= wm;
		}
	}
}

template<std::random_access_iterator InIt, std::random_access_iterator OutIt>
void fft2(InIt in_first, InIt in_last, OutIt out, FFTDirection dir) {
	using complex = typename OutIt::value_type;
	using T = typename complex::value_type;

	const ptrdiff_t size = std::distance(in_first, in_last);
	assert(isPowerOf2(size) && "Input size needs to be a power of 2");

	const auto size_half = size / 2;
	const int log_n = log2OfPowerOf2(static_cast<int>(size));
	const T nrm = T(1) / std::sqrt(T(size));

	for (int j = 0; j < size; ++j)
		out[j] = nrm * in_first[bitReverse(j, log_n)];

	for (int i = 0; i < log_n; ++i) {
		const int bm = 1 << i;
		const int bw = 2 << i;
		const T ang = T((int)(dir)) * std::numbers::pi_v<T> / T(bm);

		for (int j = 0; j < size_half; ++j) {
			const int i1 = ((j >> i) << (i + 1)) + j % bm;
			const int i2 = i1 ^ bm;
			const std::complex<T> z1 = std::polar(T(1), ang * T(i1 ^ bw));
			const std::complex<T> z2 = std::polar(T(1), ang * T(i2 ^ bw));
			const std::complex<T> tmp = out[i1];

			out[i1] += z1 * out[i2];
			out[i2] = tmp + z2 * out[i2];
		}
	}
}

}



/// @brief Computes the fast fourier transform of a signal that is given by the range
/// [in_first, in_last) whose length needs to be a power of 2. The result is written via
/// the given output iterator.
///
/// @tparam InIt      Input iterator to real or complex values. Must meet the requirements of LegacyRandomAccessIterator
/// @tparam OutIt     Output iterator to complex values. Must meet the requirements of LegacyRandomAccessIterator.
/// @param  in        iterator to input data
/// @param  in_first  iterator to begin of input data
/// @param  in_last   iterator to end of input data
/// @param  out       iterator to output data
template<std::random_access_iterator InIt, std::random_access_iterator OutIt>
void fft(InIt in_first, InIt in_last, OutIt out) {
	detail::fft(in_first, in_last, out, detail::FFTDirection::Forward);
}

/// @brief Computes the inverse fast fourier transform of a signal that is given by the range
/// [in_first, in_last) whose length needs to be a power of 2. The result is written via
/// the given output iterator.
///
/// @tparam InIt     Input iterator to complex values. Must meet the requirements of LegacyRandomAccessIterator
/// @tparam OutIt    Output iterator to complex values. Must meet the requirements of LegacyRandomAccessIterator.
/// @param  in       iterator to input data
/// @param  in_first iterator to begin of input data
/// @param  in_last  iterator to end of input data
/// @param  out      iterator to output data
template<std::random_access_iterator InIt, std::random_access_iterator OutIt>
void ifft(InIt in_first, InIt in_last, OutIt out) {
	detail::fft(in_first, in_last, out, detail::FFTDirection::Backward);
}

template<std::floating_point T>
std::vector<std::complex<T>> fft(std::vector<std::complex<T>>& data) {
	std::vector<std::complex<T>> result(data.size());
	fft(data.begin(), data.end(), result.begin());
	return result;
}

template<std::floating_point T>
std::vector<std::complex<T>> ifft(std::vector<std::complex<T>>& data) {
	std::vector<std::complex<T>> result(data.size());
	ifft(data.begin(), data.end(), result.begin());
	return result;
}



template<std::random_access_iterator InIt, std::random_access_iterator OutIt>
void fft2(InIt in_first, InIt in_last, OutIt out) {
	detail::fft(in_first, in_last, out, detail::FFTDirection::Forward);
}

template<std::random_access_iterator InIt, std::random_access_iterator OutIt>
void ifft2(InIt in_first, InIt in_last, OutIt out) {
	detail::fft(in_first, in_last, out, detail::FFTDirection::Backward);
}

template<std::floating_point T>
std::vector<std::complex<T>> fft2(std::vector<std::complex<T>>& data) {
	std::vector<std::complex<T>> result(data.size());
	fft2(data.begin(), data.end(), result.begin());
	return result;
}

template<std::floating_point T>
std::vector<std::complex<T>> ifft2(std::vector<std::complex<T>>& data) {
	std::vector<std::complex<T>> result(data.size());
	ifft2(data.begin(), data.end(), result.begin());
	return result;
}


/// @brief FFT calculator for one specific size (which needs to be a power of two).
///
/// Internal factors and numbers are precalculated during construction so
/// that when calling the transform functions only a few additions and multiplications
/// need to be performed improving the performance if many ffts with the same
/// size need to be computed.
///
/// @tparam T Data type
/// @tparam N Size of data
template<std::floating_point T, int N>
class FFTCalculator
{
public:
	using complex = std::complex<T>;

	static_assert(isPowerOf2(N), "Size N has to be a power of 2");

	constexpr FFTCalculator() {
		for (int i = 0; i < N; ++i)
			butterfly_indices[i] = bitReverse(i, log_n);

		for (int s = 0; s < log_n; ++s) {
			const int m2 = 1 << s; // 2^(s-1)
			z[s] = std::polar(T(1), std::numbers::pi_v<T> / m2);
		}
	}


	/// @brief
	///  Fourier transform for range based input. The memory that input and output point to
	///  must not overlap.
	///
	/// @tparam InIt   Input iterator to real or complex values. Must meet the requirements of LegacyRandomAccessIterator
	/// @tparam OutIt  Output iterator to complex values. Must meet the requirements of LegacyRandomAccessIterator.
	/// @param  in     iterator to input data
	/// @param  out    iterator to output data
	template<std::random_access_iterator InIt, std::random_access_iterator OutIt>
	constexpr void fft(InIt in, OutIt out) const {
		for (int i = 0; i < N; ++i) {
			out[i] = nrm * in[butterfly_indices[i]];
		}
		for (int s = 0; s < log_n; ++s) {
			const int m2 = 1 << s;
			const int m = m2 << 1;
			const auto wm = z[s];
			complex w(1, 0);

			for (int j = 0; j < m2; ++j) {
				for (int k = j; k < N; k += m) {
					const complex t = w * out[k + m2];
					const complex u = out[k];
					out[k] = u + t;
					out[k + m2] = u - t;
				}
				w *= wm;
			}
		}
	}


	/// @brief
	///  Inverse fourier transform for range based input. The memory that input and output point to
	///  must not overlap.
	///
	/// @tparam InIt   Input iterator to complex values. Must meet the requirements of LegacyRandomAccessIterator
	/// @tparam OutIt  Output iterator to complex values. Must meet the requirements of LegacyRandomAccessIterator.
	/// @param  in     iterator to input data
	/// @param  out    iterator to output data
	template<std::random_access_iterator InIt, std::random_access_iterator OutIt>
	constexpr void ifft(InIt in, OutIt out) const {
		for (int i = 0; i < N; ++i) {
			out[i] = nrm * in[butterfly_indices[i]];
		}
		for (int s = 0; s < log_n; ++s) {
			const int m2 = 1 << s;
			const int m = m2 << 1;
			const auto tmp = z[s];
			const auto wm = complex{ tmp.real(), -tmp.imag() };
			complex w(1, 0);

			for (int j = 0; j < m2; ++j) {
				for (int k = j; k < N; k += m) {
					const complex t = w * out[k + m2];
					const complex u = out[k];
					out[k] = u + t;
					out[k + m2] = u - t;
				}
				w *= wm;
			}
		}
	}


	/// @brief
	///  Inverse fourier transform for range based input. The memory that input and output point to
	///	 must not overlap. Only the real part is output and the imaginary part discarded. Therefore
	///  the input should be hermitian symmetric.
	///
	///  Note: This function allocates memory in order to convert to real values because the
	///  intermediate values need to be complex
	///
	/// @tparam InIt   Input iterator to complex values. Must meet the requirements of LegacyRandomAccessIterator
	/// @tparam OutIt  Output iterator to real values. Must meet the requirements of LegacyRandomAccessIterator.
	/// @param  in     iterator to input data
	/// @param  out    iterator to output data
	template<std::random_access_iterator InIt, std::random_access_iterator OutIt>
	constexpr void ifft_real(InIt in, OutIt out) const {
		std::vector<complex> ifft(N);
		for (int i = 0; i < N; ++i) {
			ifft[i] = nrm * in[butterfly_indices[i]];
		}
		for (int s = 0; s < log_n; ++s) {
			const int m2 = 1 << s;
			const int m = m2 << 1;
			const auto tmp = z[s];
			const auto wm = complex{ tmp.real(), -tmp.imag() };
			complex w(1, 0);

			for (int j = 0; j < m2; ++j) {
				for (int k = j; k < N; k += m) {
					const complex t = w * ifft[k + m2];
					const complex u = ifft[k];
					ifft[k] = u + t;
					ifft[k + m2] = u - t;
				}
				w *= wm;
			}
		}
		for (int i = 0; i < N; i++) {
			out[i] = ifft[i].real();
		}
	}

private:
	static constexpr int log_n{ log2OfPowerOf2(N) };
	T nrm{ T(1) / std::sqrt(T(N)) };

	std::array<int, N> butterfly_indices;
	std::array<complex, log_n> z;
};




/// @brief FFT calculator (slower version) for one specific size (which needs to be a power of two).
///
/// Internal factors and numbers are precalculated during construction so
/// that when calling the transform functions only a few additions and multiplications
/// need to be performed improving the performance if many ffts with the same
/// size need to be computed.
///
/// @tparam T Data type
/// @tparam N Size of data
template<class T, int N>
class FFTCalculator2
{
public:
	using complex = std::complex<T>;

	static_assert(isPowerOf2(N), "Size N has to be a power of 2");

	constexpr FFTCalculator2() {
		for (int i = 0; i < N; ++i)
			butterfly_indices[i] = bitReverse(i, log_n);

		for (int i = 0; i < log_n; ++i) {
			const int bm = 1 << i;
			const int bw = 2 << i;
			const T ang = std::numbers::pi_v<T> / T(bm);

			for (int j = 0; j < n_half; ++j) {
				const auto index = i * n_half + j;
				const int i1 = ((j >> i) << (i + 1)) + j % bm;
				const int i2 = i1 ^ bm;
				indices1[index] = i1;
				indices2[index] = i2;
				zs1[index] = std::polar(T(1), ang * T(i1 ^ bw));
				zs2[index] = std::polar(T(1), ang * T(i2 ^ bw));
			}
		}
	}


	/// @brief
	///  Fourier transform for range based input. The memory that input and output point to
	///  must not overlap.
	///
	/// @tparam InIt   Input iterator to real or complex values. Must meet the requirements of LegacyRandomAccessIterator
	/// @tparam OutIt  Output iterator to complex values. Must meet the requirements of LegacyRandomAccessIterator.
	/// @param  in     iterator to input data
	/// @param  out    iterator to output data
	template<std::random_access_iterator InIt, std::random_access_iterator OutIt>
	constexpr void fft(InIt in, OutIt out) const {
		for (int i = 0; i < N; ++i) {
			out[i] = nrm * in[butterfly_indices[i]];
		}

		for (int i = 0; i < log_n * n_half; ++i) {
			const auto i1 = indices1[i];
			const auto i2 = indices2[i];
			const auto tmp = out[i1];
			out[i1] += zs1[i] * out[i2];
			out[i2] = tmp + zs2[i] * out[i2];
		}
	}


	/// @brief
	///  Inverse fourier transform for range based input. The memory that input and output point to
	///  must not overlap.
	///
	/// @tparam InIt   Input iterator to complex values. Must meet the requirements of LegacyRandomAccessIterator
	/// @tparam OutIt  Output iterator to complex values. Must meet the requirements of LegacyRandomAccessIterator.
	/// @param  in     iterator to input data
	/// @param  out    iterator to output data
	template<std::random_access_iterator InIt, std::random_access_iterator OutIt>
	constexpr void ifft(InIt in, OutIt out) const {
		for (int i = 0; i < N; ++i) {
			out[i] = nrm * in[butterfly_indices[i]];
		}

		for (int i = 0; i < log_n * n_half; ++i) {
			const auto i1 = indices1[i];
			const auto i2 = indices2[i];
			const auto tmp = out[i1];
			complex z1_conj{ zs1[i].real(), -zs1[i].imag() };
			complex z2_conj{ zs2[i].real(), -zs2[i].imag() };
			out[i1] += z1_conj * out[i2];
			out[i2] = tmp + z2_conj * out[i2];
		}
	}


	/// @brief
	///  Inverse fourier transform for range based input. The memory that input and output point to
	///	 must not overlap. Only the real part is output and the imaginary part discarded. Therefore
	///  the input should be hermitian symmetric.
	///
	///  Note: This function allocates memory in order to convert to real values because the
	///  intermediate values need to be complex
	///
	/// @tparam InIt   Input iterator to complex values. Must meet the requirements of LegacyRandomAccessIterator
	/// @tparam OutIt  Output iterator to real values. Must meet the requirements of LegacyRandomAccessIterator.
	/// @param  in     iterator to input data
	/// @param  out    iterator to output data
	template<std::random_access_iterator InIt, std::random_access_iterator OutIt>
	constexpr void ifft_real(InIt in, OutIt out) const {
		std::vector<complex> ifft(N);
		for (int i = 0; i < N; ++i) {
			ifft[i] = nrm * in[butterfly_indices[i]];
		}

		for (int i = 0; i < log_n * n_half; ++i) {
			const auto i1 = indices1[i];
			const auto i2 = indices2[i];
			const auto tmp = ifft[i1];
			complex z1_conj{ zs1[i].real(), -zs1[i].imag() };
			complex z2_conj{ zs2[i].real(), -zs2[i].imag() };
			ifft[i1] += z1_conj * ifft[i2];
			ifft[i2] = tmp + z2_conj * ifft[i2];
		}
		for (int i = 0; i < N; i++) {
			out[i] = ifft[i].real();
		}
	}

	/// @brief       Fourier transform for complex-valued array input.
	///
	/// @param data  input data array
	/// @return      Fourier transform of data
	constexpr std::array<complex, N> fft(const std::array<complex, N>& data) const {
		std::array<complex, N> result;
		fft(data.begin(), result.begin());
		return result;
	}

	/// @brief       Fourier transform for real-valued array input.
	///
	/// @param data  input data array
	/// @return      Fourier transform of data
	constexpr std::array<complex, N> fft(const std::array<T, N>& data) const {
		std::array<complex, N> result;
		fft(data.begin(), result.begin());
		return result;
	}

	/// @brief       Inverse fourier transform for array input.
	///
	/// @param data  input data array
	/// @return      Inverse fourier transform.
	constexpr std::array<complex, N> ifft(const std::array<complex, N>& data) const {
		std::array<complex, N> result;
		ifft(data.begin(), result.begin());
		return result;
	}

	/// @brief       Inverse fourier transform for array input which discards the imaginary part.
	///
	/// @param data  hermitian symmetric input data array
	/// @return      Inverse fourier transform.
	constexpr std::array<T, N> ifft_real(const std::array<complex, N>& data) const {
		std::array<T, N> result;
		ifft_real(data.begin(), result.begin());
		return result;
	}



private:
	static constexpr int log_n{ log2OfPowerOf2(N) };
	static constexpr int n_half{ N / 2 };
	T nrm{ T(1) / std::sqrt(T(N)) };

	std::array<int, N> butterfly_indices;
	std::array<int, n_half * log_n> indices1, indices2;
	std::array<complex, n_half * log_n> zs1, zs2;
};

}
