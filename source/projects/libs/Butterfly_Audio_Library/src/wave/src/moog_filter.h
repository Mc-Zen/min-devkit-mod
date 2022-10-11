#pragma once

#include <numbers>


//findCrossings(begin, end, value=0) -> vector<double>
// 
// // optinoal or throw exception on failure and catch it
//
// std::optional<std::pair<float, float>> getPeriod(data) {
//	pitch = getPitch(data.begin(), data.begin() + 1000);
//	if (!pitch) fail;
//	period = 1 / pitch;
//	crossings = findCrossings(begin(), begin() + 2 * period);
//	if (crossings.size() < 2) fail;
//	c1 = crossings[0];
// 
//	c1 = findFirstCrossing(data.rbegin(), data.rend())
//
//	c2 = 0;
//	for (int i = 1; i < crossings.size(); i++) {
//		// next to pitch
//	}
//
//}


namespace Butterfly {


template<class T, class ParamType = double>
class MoogFilter
{
public:
	MoogFilter(ParamType sampleRate) : sampleRate(sampleRate), sampleRateInv(1. / sampleRate) {
		init();
	}

	void setSampleRate(ParamType fs) {
		sampleRate = fs;
		sampleRateInv = ParamType(1) / fs;
		calc();
	}

	void setResonance(ParamType filterRezo) {
		resonance = filterRezo;
		calc();
	}

	void setFrequency(ParamType filterCutoff) {
		cutoff = 2. * filterCutoff * sampleRateInv;
		calc();
	}

	constexpr T operator()(T input) {
		// process input
		x = input - r * y4;

		// four cascaded one-pole filters (bilinear transform)
		y1 = x * p + oldx * p - k * y1;
		y2 = y1 * p + oldy1 * p - k * y2;
		y3 = y2 * p + oldy2 * p - k * y3;
		y4 = y3 * p + oldy3 * p - k * y4;

		// clipper band limited sigmoid
		//constexpr T inv_6 = T(1.0 / 6.0);
		//	y4 -= (y4 * y4 * y4) * inv_6;

		oldx = x;
		oldy1 = y1;
		oldy2 = y2;
		oldy3 = y3;

		return y4;
	}


	constexpr void reset() {
		init();
	}
	//// filter an input sample using normalized params
	//T filter(T input, T cutoff, T resonance) {
	//	// set params first
	//	cutoff = filterCutoff;
	//	resonance = filterRezo;
	//	calc();

	//	return process(input);
	//}

	constexpr T getSampleRate() const { return sampleRate; }
	constexpr T getResonance() const { return resonance; }
	constexpr T getCutoff() const { return cutoff; }
	constexpr T getCutoffHz() const { return cutoff * sampleRate * 0.5; }

protected:
	void init() {
		y1 = y2 = y3 = y4 = oldx = oldy1 = oldy2 = oldy3 = T(0.0);
		calc();
	}

	void calc() {
		// empirical tuning
		p = cutoff * (T(1.8) - T(0.8) * cutoff);
		// k = p + p - T(1.0);
		// A much better tuning seems to be:
		k = T(2.0) * std::sin(cutoff * std::numbers::pi_v<T> * T(0.5)) - T(1.0);

		T t1 = (T(1.0) - p) * T(1.386249);
		T t2 = T(12.0) + t1 * t1;
		r = resonance * (t2 + T(6.0) * t1) / (t2 - T(6.0) * t1);
	}

private:
	// cutoff and resonance [0 - 1]
	T cutoff{ 1 };
	T resonance{};
	T sampleRate{}, sampleRateInv{};
	T fs{};
	T y1{}, y2{}, y3{}, y4{};
	T oldx{};
	T oldy1{}, oldy2{}, oldy3{};
	T x{}, r{}, p{}, k{};
};






}