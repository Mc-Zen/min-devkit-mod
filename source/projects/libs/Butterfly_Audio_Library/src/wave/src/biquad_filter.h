#pragma once

#include <numbers>
#include <concepts>

namespace Butterfly {


template<std::floating_point T>
class BiquadBase
{
public:
	using value_type = T;

	constexpr BiquadBase() = default;
	constexpr BiquadBase(double samplerate) {}

	// y[n] = b0·x[n] + b1·x[n-1] + b2·x[n-2] - a1·y[n-1] - a2·y[n-2]
	constexpr T operator()(T x) {
		const auto y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
		x2 = x1;
		x1 = x;
		y2 = y1;
		y1 = y;
		return y;
	}

	void reset() {
		x1 = x2 = y1 = y2 = T(0);
		a1 = a2 = b1 = b2 = T(0);
		b0 = T(1);
	}

protected:
	T x1{}, x2{}, y1{}, y2{};
	T a1{}, a2{};		// poles
	T b0{}, b1{}, b2{}; // zeros
};


template<std::floating_point T>
class BiquadFilter : BiquadBase<T>
{
public:
	enum class Type {
		Lowpass,
		Highpass,
		Bandpass,
		Notch,
		Peak,
		Lowshelf,
		Highshelf,
		Allpass
	};
	
	using value_type = T;

	constexpr BiquadFilter() = default;

	constexpr BiquadFilter(double samplerate, double frequency)
		: samplerate(samplerate),
		  invSamplerate(T(1) / samplerate),
		  frequency(frequency) {
	}

	void setFrequency(T frequency) {
		this->frequency = frequency;
		update();
	}

	void setQ(T q) {
		this->q = q;
		update();
	}

	void setGain(T gain) {
		this->gain = gain;
		update();
	}

	void setType(Type type) {
		this->type = type;
		update();
	}

	void getSamplerate() const { return samplerate; }
	void getFrequency() const { return frequency; }
	void getGain() const { return gain; }
	void getQ() const { return q; }
	Type getType() const { return type; }

protected:
	void update() {
		switch (type) {
		case Type::Lowpass: updateLowpass(); return;
		case Type::Highpass: updateHighpass(); return;
		case Type::Bandpass: updateBPF_constantPeakGain0(); return;
		case Type::Notch: updateNotch(); return;
		case Type::Peak: updatePeak(); return;
		case Type::Lowshelf: updateLowshelf(); return;
		case Type::Highshelf: updateHighshelf(); return;
		case Type::Allpass: updateAllpass(); return;
		}
	}

	void updateHighpass() { updateHPF_or_LPF(1.); }
	void updateLowpass() { updateHPF_or_LPF(-1.); }
	void updateHighshelf() { updateLow_or_HighShelf(-1.); }
	void updateLowshelf() { updateLow_or_HighShelf(1.); }

	void updateNotch() {
		const auto w0 = T(2) * std::numbers::pi_v<T> * frequency * invSamplerate;
		const auto a = std::sin(w0) / (2 * q);
		const auto a0_inv = 1. / (1. + a);

		this->b0 = b2 = 1. * a0_inv;
		b1 = a1 = -2. * std::cos(w0) * a0_inv;
		a0 = (1. + a) * a0_inv;
		a2 = (1. - a) * a0_inv;
	}

	void updateBPF_constantSkirtGainQ() {
		const auto w0 = T(2) * std::numbers::pi_v<T> * frequency * invSamplerate;
		const auto a = std::sin(w0) / (2 * q);
		const auto a0_inv = 1. / (1. + a);

		this->b0 = (a * q) * a0_inv;
		b1 = 0;
		b2 = -this->b0;
		a1 = -2 * std::cos(w0) * a0_inv;
		a2 = (1 - a) * a0_inv;
	}

	void updateBPF_constantPeakGain0() {
		const auto w0 = T(2) * std::numbers::pi_v<T> * frequency * invSamplerate;
		const auto a = std::sin(w0) / (2 * q);
		const auto a0_inv = 1. / (1. + a);

		this->b0 = a * a0_inv;
		b1 = 0;
		b2 = -this->b0;
		a1 = -2 * std::cos(w0) * a0_inv;
		a2 = (1 - a) * a0_inv;
	}

	void updateAllpass() {
		const auto w0 = T(2) * std::numbers::pi_v<T> * frequency * invSamplerate;
		const auto a = std::sin(w0) / (2 * q);
		const auto a0_inv = 1. / (1. + a);

		b1 = a1 = -2. * std::cos(w0) * a0_inv;
		b2 = a0 = (1. + a) * a0_inv;
		this->b0 = a2 = (1. - a) * a0_inv;
	}
	void updateLow_or_HighShelf(T lowOrHigh) { // low: 1, high: -1
		const auto A = std::pow(10, gain * 0.025);
		const auto w0 = T(2) * std::numbers::pi_v<T> * frequency * invSamplerate;
		const auto cosw = std::cos(w0);
		const auto a = std::sin(w0) / (2 * q);
		const auto a0_inv = 1. / (1. + a);

		const auto Ap1 = A + 1.;
		const auto Am1 = A - 1.;
		const auto sqA = std::sqrt(A);
		const auto sqAa2 = sqA * a * 2.;
		const auto e = Ap1 - lowOrHigh * Am1 * cosw;
		const auto f = Ap1 + lowOrHigh * Am1 * cosw;

		const auto a0_inv = 1. / (f + sqAa2);

		this->b0 = A * (e + sqAa2) * a0_inv;
		b1 = lowOrHigh * 2 * A * (Am1 - Ap1 * cosw) * a0_inv;
		b2 = A * (e - sqAa2) * a0_inv;

		a1 = lowOrHigh * -2. * (Am1 + Ap1 * cosw) * a0_inv;
		a2 = (f - sqAa2) * a0_inv;
	}

	void updatePeak() {
		const auto A = std::pow(10, gain * 0.025);
		const auto w0 = T(2) * std::numbers::pi_v<T> * frequency * invSamplerate;
		const auto a = std::sin(w0) / (2 * q);
		const auto cosw = std::cos(w0);
		const auto a0_inv = 1. / (1. + a / A);

		this->b0 = (1 + a * A) * a0_inv;
		b2 = (1 - a * A) * a0_inv;
		b1 = -2 * cosw * a0_inv;
		a1 = b1;
		a2 = (1 - a / A) * a0_inv;
	}

	void updateHPF_or_LPF(T hpf_lpf) { // hpf: 1, lpf: -1
		const auto w0 = T(2) * std::numbers::pi_v<T> * frequency * invSamplerate;
		const auto a = std::sin(w0) / (2 * q);
		const auto cosw = std::cos(w0);
		const auto a0_inv = 1. / (1. + a);
		b1 = (1. + hpf_lpf * cosw) * a0_inv;
		this->b0 = b2 = -hpf_lpf * 0.5 * b1;
		a1 = -2 * cosw * a0_inv;
		a2 = (1. - a) * a0_inv;
	}

	T samplerate{}, invSamplerate{};
	T frequency{};
	T q{};
	T gain{};

	Type type{ Type::Lowpass };
};


}
