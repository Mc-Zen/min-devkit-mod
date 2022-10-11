#pragma once
#include <algorithm>
#include <iterator>
#include <vector>
#include <cmath>


namespace Butterfly {

/// @brief  Get absolute maximum of signal (peak value).
/// @tparam It    Input iterator, must meet the requirements of LegacyForwardIterator
/// @param  first Signal range start
/// @param  last  Signal range end
/// @return Peak value
template<std::forward_iterator It>
typename std::iterator_traits<It>::value_type peak(It first, It last) {
	typename std::iterator_traits<It>::value_type abs_max{};
	std::for_each(first, last, [&abs_max](const auto& v) { abs_max = std::max(abs_max, std::abs(v)); });
	return abs_max;
}


/// @brief  Get RMS (root mean square) of signal,
/// @tparam It    Input iterator, must meet the requirements of LegacyForwardIterator
/// @param  first Signal range start
/// @param  last  Signal range end
/// @return RMS
template<std::forward_iterator It>
typename std::iterator_traits<It>::value_type rms(It first, It last) {
	typename std::iterator_traits<It>::value_type rms{};
	std::for_each(first, last, [&rms](const auto& v) { rms += v * v; });
	return static_cast<typename std::iterator_traits<It>::value_type>(std::sqrt(rms / static_cast<double>(std::distance(first, last))));
}


/// @brief  Normalize signal by peak to range [-\value,\value].
/// @tparam It    Input-output iterator, must meet the requirements of LegacyForwardIterator
/// @param  first start of range
/// @param  last  end of range
/// @param  value normalization value (default is 1.0)
template<std::forward_iterator It>
void peak_normalize(It first, It last, typename std::iterator_traits<It>::value_type value = 1.0) {
	const auto inv = value / peak(first, last);
	for (auto it = first; it != last; ++it)
		*it *= inv;
}


/// @brief  Normalize signal by RMS to range [-\value,\value].
/// @tparam It    Input-output iterator, must meet the requirements of LegacyForwardIterator
/// @param  first start of range
/// @param  last  end of range
/// @param  value normalization value (default is 1.0)
template<std::forward_iterator It>
void rms_normalize(It first, It last, typename std::iterator_traits<It>::value_type value = 1.0) {
	const auto inv = value / rms(first, last);
	for (auto it = first; it != last; ++it)
		*it *= inv;
}



/// @brief Get points where given data crosses the given value. 
/// 
/// @tparam InIt           Input iterator
/// @tparam T              output data type (defaults to double)
/// @param first           start of range
/// @param last            end of range
/// @param value           crossing value
/// @param maxNumberToFind stop searching once this number of crossings have been found
/// @return vector of zero crossings
template<std::forward_iterator InIt, std::floating_point T = double>
std::vector<T> getCrossings(InIt first, InIt last, T value = 0, size_t maxNumberToFind = std::numeric_limits<size_t>::max()) {
	if (first == last) return {};

	std::vector<T> crossings;

	auto previousValue = *first;
	bool isAbove = previousValue > value;
	auto it = first;
	++it;
	auto count = 1;

	while (crossings.size() < maxNumberToFind && it != last) {
		const auto val = *it;
		if ((val > value) != isAbove) {
			isAbove = !isAbove;
			const T dy = val - previousValue;
			crossings.push_back(-(val - dy * count - value) / dy);
		}
		previousValue = val;
		++it;
		++count;
	}
	return crossings;
}


/// @brief Compute the discrete difference given by out[i] = in[i+1] - in[i]
///
/// @param firstIn  Iterator to first input element
/// @param lastIn   Iterator to last input element
/// @param firstOut Iterator to first output element. This must point to some storage that is at least as big as the input range minus 1.
template<std::forward_iterator InIt, std::forward_iterator OutIt>
void differentiate(InIt firstIn, InIt lastIn, OutIt firstOut) {
	auto inIt = firstIn;
	auto previousValue = *inIt;
	++inIt;
	auto outIt = firstOut;
	for (; inIt != lastIn; ++inIt, ++outIt) {
		auto value = *inIt;
		*outIt = value - previousValue;
		previousValue = value;
	}
}



/// @brief Average magnitude difference function
///
/// @param firstIn Iterator to first input element
/// @param lastIn Iterator to last input element
/// @param firstOut Iterator to first output element. This must point to some storage that is at least as big as the input range.
template<std::random_access_iterator InIt, std::random_access_iterator OutIt>
void amdf(InIt firstIn, InIt lastIn, OutIt firstOut) {
	const auto size = static_cast<size_t>(std::distance(firstIn, lastIn));
	for (size_t i = 0; i < size; ++i) {
		firstOut[i] = 0;
		for (size_t j = i; j < size; ++j) {
			firstOut[i] += std::abs(firstIn[j - i] - firstIn[j]);
		}
	}
}

}