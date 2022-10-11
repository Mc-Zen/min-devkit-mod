
// Pitch detection algorithms


#pragma once

#include "waveform_processing.h"
//#include <coroutine>
#include <numeric>
#include <iostream>
#include <optional>
#include <algorithm>


namespace Butterfly {

namespace detail {

template<std::forward_iterator It, std::floating_point T = double>
std::pair<T, T> meanAndStandardDeviation(It first, It last) {
	using data_type = typename std::iterator_traits<It>::value_type;
	data_type m{}, sd{};
	size_t count{};
	for (auto it = first; it != last; ++it, ++count) {
		m += *it;
	}
	const T mean = m / static_cast<data_type>(count);
	for (auto it = first; it != last; ++it) {
		const auto tmp = *it - mean;
		sd += tmp * tmp;
	}
	return { mean, std::sqrt(sd / static_cast<data_type>(count)) };
}

} // namespace detail



struct PitchInfo
{
	double frequency;
	double standardDeviation;
	double maxDeviation;
};

struct PitchFindingParameters
{
	// Period finder tolerance. Increasing the tolerance helps to find pitches os noisy
	// data but can decrease the accuracy and lead to false positives.
	double tolerance{ 0.3 };

	// Filters outliers.
	// The algorithm finds a number of periods. Some values might be off by a large factor
	// which can be filtered out using a filter. All periods which do not lie in the
	// filter range [avg+deviationFilter*avg,avg-deviationFilter*avg] of the intermediate
	// average are discarded, increasing the accuracy.
	double deviationFilter{ 0.3 };

	// When a long sample is analyzed, the pitch might vary. This parameter constrains
	// the search to the first n periods.
	size_t maxPeriodsToAverage = std::numeric_limits<size_t>::max();
};

template<std::random_access_iterator InIt>
std::optional<PitchInfo> getPitch(InIt first, InIt last, const PitchFindingParameters& parameters = {}) {
	using data_type = typename std::iterator_traits<InIt>::value_type;
	const auto size = std::distance(first, last);
	const auto devFilter = parameters.deviationFilter;
	const auto tol = parameters.tolerance;

	if (size < 10) return std::nullopt;

	std::vector<data_type> data(first, last);

	std::vector<data_type> amdf(size);
	Butterfly::amdf(data.begin(), data.end(), amdf.begin());
	peak_normalize(amdf.begin(), amdf.end());

	// get extrema of amdf by getting zero crossings of diff(amdf)
	std::vector<data_type> diff(size - 1);
	differentiate(amdf.begin(), amdf.end(), diff.begin());
	peak_normalize(diff.begin(), diff.end());

	const auto crossings = getCrossings(diff.begin(), diff.end(), data_type{ 0 });

	// only get extrema that are close to 0 (determined by tolerance)
	std::vector<data_type> filteredCrossings;
	std::copy_if(crossings.begin(), crossings.end(), std::back_inserter(filteredCrossings),
		[&amdf, tol](auto crossing) {
			crossing += 0.5;
			auto z = std::abs(amdf[static_cast<int>(crossing)]);
			const auto corrected_tolerance = (1 - crossing / amdf.size()) * tol;
			return std::abs(amdf[static_cast<int>(crossing)]) < corrected_tolerance && crossing > 3.;
		});
	if (filteredCrossings.size() < 2) return std::nullopt;

	// get period lengths as differences between periods
	std::vector<data_type> values(filteredCrossings.size() - 1);
	differentiate(filteredCrossings.begin(), filteredCrossings.end(), values.begin());
	// transform to frequency
	std::for_each(values.begin(), values.end(), [](auto& a) { a = 1 / a; });

	const auto f_and_sdv = detail::meanAndStandardDeviation(values.begin(), values.end());
	const auto f0 = f_and_sdv.first;
	std::erase_if(values, [f0, devFilter](auto a) { return std::abs(a - f0) > f0 * devFilter; });
	const auto f_and_sdv1 = detail::meanAndStandardDeviation(values.begin(), values.end());
	const auto f1 = f_and_sdv1.first;
	const auto sdv1 = f_and_sdv1.second;
	if (values.empty()) return std::nullopt;

	//std::cout << "sdv " << sdv1 << "\n";
	const auto maxDeviationIt = std::max_element(values.begin(), values.end(),
		[f1](auto a, auto b) { return std::abs(a - f1) < std::abs(b - f1); });
	const auto maxDeviation = std::abs(*maxDeviationIt - f1);

	return PitchInfo{ f1, sdv1, maxDeviation };
}

}
