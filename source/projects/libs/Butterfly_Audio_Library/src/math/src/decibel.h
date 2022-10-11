
// Decibel conversion functions

#pragma once
#include <cmath>


namespace Butterfly {



/// @brief Convert decibels to normalized volume.
/// @tparam T Sample type
/// @param dB decibel value
/// @return normalized value
template<class T>
T dBToNormalized(T dB) {
	return static_cast<T>(std::pow(10.0, 0.05f * static_cast<double>(dB)));
}

/// @brief Convert normalized volume to decibels.
/// @tparam T Sample tpe
/// @param volume normalized value
/// @return decibel
template<class T>
T normalizedTodB(T volume) {
	return static_cast<T>(20.0 * std::log10(static_cast<double>(volume)));
}

}
