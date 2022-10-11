#pragma once


#include "fft.h"

namespace Butterfly {



/// @brief Remove spectral components in given DFT so the signal will not aliase in time domain when
/// played back with given maximum playback frequency at given samplerate.
///
/// @tparam It                     Input/output iterator, must meet the requirements of LegacyRandomAccessIterator
/// @param  first                  Signal range start
/// @param  last                   Signal range end
/// @param  samplerate             Sampling rate
/// @param  max_playback_frequency Maximum frequency at which the buffer may be played back periodically without aliasing at given samplerate
template<std::random_access_iterator It, std::floating_point T>
void antialiase_dft(
	It first, It last,
	T samplerate,
	T max_playback_frequency) {

	const auto size = static_cast<size_t>(std::distance(first, last));

	const auto nyquist = samplerate * 0.5;
	const auto nyquist_index = nyquist / max_playback_frequency;
	const auto cutoff_index = static_cast<size_t>(std::floor(nyquist_index)) + 1;

	if (cutoff_index > size / 2) return;

	(*first).imag(0);

	for (auto it = first + cutoff_index; it != last - cutoff_index + 1; ++it) {
		*it = {};
	}
}


/// @brief Antialiase given signal for a number of maximum frequencies in [freq_first, freq_last) using fourier bandlimiting.
/// The signal length needs to be a power of 2 and must match the size of the `FFTCalculator`. The latter defines type and size of the signal.
/// In order to write the antialiased signals in `std::distance(freq_first, freq_last)` outputs, the iterator
/// @p out_table_first needs to be dereferencable to a type which has a random access iterator
/// that can be fetched by `std::begin(*out_table_first + n)`.
///
/// @tparam T                      Sample type (i.e. float or double)
/// @tparam size                   Size of signal (needs to be a power of 2)
/// @tparam SignalIt               Input data iterator, needs to fulfill the requirements of ForwardIterator
/// @tparam FrequencyIt            Frequency input iterator, needs to fulfill the requirements of ForwardIterator
/// @tparam OutputIterator         Output data iterator, needs to fulfill the requirements of LegacyForwardIterator and must point to something that <code>std::begin(T&)</code> can be called upon and is itself a LegacyIt
/// @param  signal_first           Iterator to first input data
/// @param  freq_first             Iterator to first frequency input element
/// @param  freq_last              Iterator to last frequency input element
/// @param  out_table_first        Iterator to first range
/// @param  samplerate             Sampling rate
/// @param  fft_calculator         FFT calculator
template<std::floating_point T, int size, std::forward_iterator SignalIt, std::forward_iterator FrequencyIt, std::forward_iterator OutputIterator>
requires requires(OutputIterator it) { requires std::forward_iterator<decltype(std::begin(*it))>; }
void antialiase(
	SignalIt signal_first,
	FrequencyIt freq_first, FrequencyIt freq_last,
	OutputIterator out_table_first,
	T samplerate,
	const Butterfly::FFTCalculator<T, size>& fft_calculator) {

	std::vector<T> data(size);
	std::vector<std::complex<T>> fft(size);
	std::copy(signal_first, signal_first + size, data.begin());

	fft_calculator.fft(signal_first, fft.begin());

	auto freq_it = freq_first;
	auto table_it = out_table_first;
	for (; freq_it != freq_last; ++freq_it, ++table_it) {
		auto copy = fft;
		antialiase_dft(copy.begin(), copy.end(), samplerate, *freq_it);
		const auto iter = std::begin(*table_it);
		fft_calculator.ifft_real(copy.begin(), iter);
	}
}


}