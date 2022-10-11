
#pragma once

#include <cmath>
#include <cassert>


namespace Butterfly {


enum class RampingType {
	Linear,
	Exponential
};


template<class T, RampingType rampingType = RampingType::Linear>
class RampedValue
{
public:

	constexpr RampedValue() = default;
	constexpr RampedValue(T value, int steps = 100) : value(value), target(value) {
		if constexpr (rampingType == RampingType::Exponential) {
			assert(value > T{ 0 } && "Value needs to be positive for exponential ramping.");
		}
		setSteps(steps);
	};

	constexpr T operator()() const { return value; }

	constexpr T operator++() {
		if (countDown <= 0) return value = target; // set to target because the incremental approach may lead to an imprecise result
		countDown--;
		if constexpr (rampingType == RampingType::Linear) {
			return value += inc;
		} else {
			return value *= inc;
		}
	}

	constexpr T operator++(int) {
		const auto tmp = value;
		++(*this);
		return tmp;
	}

	/// @brief Set the new target value.
	/// @param v New value
	/// @return Whether ramping is needed.
	constexpr bool set(T v) {
		target = v;
		if (steps == 0 || value == target) {
			value = target;
			return false;
		}
		if constexpr (rampingType == RampingType::Linear) {
			inc = (target - value) / static_cast<T>(steps);
		} else {
			assert(value > T{ 0 } && "Value needs to be positive for exponential ramping.");
			inc = std::pow(target / value, T{ 1 } / static_cast<T>(steps));
		}
		countDown = steps;
		return true;
	}

	constexpr void setImmediately(T v) {
		value = target = v;
		countDown = 0;
	}

	constexpr void setSteps(int steps) {
		this->steps = steps;
	}

	constexpr void setTime(T milliseconds, T sampleRate) {
		steps = static_cast<int>(milliseconds / T{ 1000 } * sampleRate);
	}

	constexpr T getTarget() const { return target; }
	constexpr int getSteps() const { return steps; }
	constexpr bool isRamping() const { return countDown > 0; }
	constexpr RampingType getRampingType() const { return rampingType; }

private:
	T value{ rampingType == RampingType::Linear ? 0 : 1 };
	T target{ value };
	int steps{ 100 };
	int countDown{};
	T inc{};
};


}
