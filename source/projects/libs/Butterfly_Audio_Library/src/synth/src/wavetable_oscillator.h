
// Implementation of a wavetable oscillator.
// More overloads for antialiase to simplify creation of antialiased wavetables.


#pragma once
#include <cmath>
#include <iostream>
#include "wavetable.h"
#include "antialiase.h"
#include "fixed_point.h"

namespace Butterfly {


/// @brief Table selector for a wavetable oscillator which iterates forward through
///        the given range until a table is found which has a playback frequency above the
///        given frequency.
struct ForwardSearchTableSelector
{
	template<class Iterator, class T>
	static constexpr Iterator selectTable(Iterator begin, Iterator end, T frequency) {
		size_t i = 0;
		for (auto it = begin; it != end; ++it, ++i) {
			if (it->getMaximumPlaybackFrequency() >= frequency) {
				return it;
			}
		}
		return end;
	}
};


/// @brief Wavetable oscillator wrapping access to multiple wavetables depending on the
///        frequency used (i.e. in order to prevent aliasing).
///
///        It is assumed that the tables are sorted in an ascending order, so that the table with
///        the lowest frequency is the first one.
///
///        The oscillator is only in a valid and usable state when the tables and the frequency is
///        set (as done in the first constructor).
///
///        The frequency may technically exceed the frequency of the last table which may result
///        in aliasing. However, the frequency shall not exceed the sample rate (which is also
///        asserted internally).
///
/// Usage example using a std::vector of wavetables (which is the default storage type):
/// \code
///	    std::vector<Wavetable<double>> tables(4);
///	    // fill tables
///	    WavetableOscillator<Wavetable<double>> osc{ &tables, 40100., 200. };
///	    auto sample = ++osc;
/// \endcode
///
/// You can also use a std::array (or any other random access storage type) for your tables.
/// \code
///     using Table = Wavetable<double>;
///     using MultiTable = std::array<Table, 5>;
///     using Osc = WavetableOscillator<Table, MultiTable>;
///
///     MultiTable tables;
///     // fill tables
///	    Osc{ tables, 40100., 200 };
/// \endcode
///
/// @tparam Wavetable Wavetable class (needs to feature size_t size(), T get(T) and getMaximumPlaybackFrequency()).
/// @tparam TableSelector Functor which is used to select the ideal table for the current playback frequency.
template<class Wavetable, class StorageType = std::vector<Wavetable>, class TableSelector = ForwardSearchTableSelector, class ParamType = double>
class WavetableOscillator
{
public:
	using wavetable_type = Wavetable;
	using multiwavetable_type = StorageType;
	using value_type = typename Wavetable::value_type;
	using param_type = ParamType;
	using T = value_type;

	constexpr WavetableOscillator() = default;
	constexpr WavetableOscillator(ParamType sampleRate) : sampleRateInv(ParamType{ 1 } / sampleRate) {}

	constexpr WavetableOscillator(StorageType* wavetables, ParamType sampleRate, ParamType frequency)
		: sampleRateInv(1.0 / sampleRate),
		  frequency(frequency) {
		setTable(wavetables);
	}

	constexpr void setTable(StorageType* wavetables) {
		this->wavetables = wavetables;
		topFreq = bottomFreq = 0.0;
		setFrequency(frequency);
	}

	constexpr void setSampleRate(ParamType sampleRate) {
		this->sampleRateInv = 1.0 / sampleRate;
		setFrequency(frequency);
	}

	constexpr void setFrequency(ParamType frequency) {
		this->frequency = frequency;
		assert(frequency * sampleRateInv < 1.0 && "The frequency needs to be lower that the sample rate");
		selectTable();
		updateDelta();
	}
	

	/// @brief Increment the oscillator by one step and get the current value.
	/// @return current value
	constexpr T operator++() {
		currentSamplePosition += delta;
		if (currentSamplePosition >= currentTableSize) {
			currentSamplePosition -= currentTableSize;
		}
		value = (*currentTable)(static_cast<T>(currentSamplePosition));
		return value;
	}

	/// @brief Increment the oscillator by one step and get the former current value.
	/// @return former value
	constexpr T operator++(int) {
		const auto tmp = value;
		currentSamplePosition += delta;
		if (currentSamplePosition >= currentTableSize) {
			currentSamplePosition -= currentTableSize;
		}
		value = (*currentTable)(static_cast<T>(currentSamplePosition));
		return tmp;
	}

	/// @brief Get current value of the oscillator without changing its state.
	/// @return current value
	constexpr T operator()() const { return value; }

	/// @brief Reset the position/phase to 0. Also updates the current value.
	constexpr void retrigger() {
		currentSamplePosition = 0;
		value = (*currentTable)(static_cast<T>(currentSamplePosition));
	}

	constexpr void reset() { retrigger(); }
	constexpr Wavetable* getSelectedTable() const { return currentTable; }
	constexpr ParamType getFrequency() const { return frequency; }
	constexpr ParamType getSampleRate() const { return ParamType{ 1 } / sampleRateInv; }

private:
	constexpr void selectTable() {
		// chances are good that the frequency has changed just a little bit, so
		// only select table if antialiasing condition really does not hold.
		if (frequency <= topFreq && frequency > bottomFreq) return;

		assert(wavetables);
		assert(wavetables->size() > 0);
		auto it = TableSelector::selectTable(wavetables->begin(), wavetables->end(), frequency);

		if (it == wavetables->end()) { // No table can be selected without aliasing -> then we just get aliasing
			it--;
		}

		auto newTable = &(*it);

		if (currentTableSize != 0) {
			const auto oldTableSize = currentTableSize;
			const auto newTableSize = newTable->size();
			currentSamplePosition = currentSamplePosition * static_cast<double>(newTableSize) / static_cast<double>(oldTableSize);
			currentSamplePosition = std::max(0., std::min((newTableSize - 1.e-7), currentSamplePosition));
		}

		currentTable = newTable;
		currentTableSize = currentTable->size();
		assert(currentTableSize > 0 && "Size of wavetables may not be zero");
		value = (*currentTable)(static_cast<T>(currentSamplePosition));


		topFreq = currentTable->getMaximumPlaybackFrequency();
		if (it == wavetables->begin()) {
			bottomFreq = T(0);
		} else {
			bottomFreq = (it - 1)->getMaximumPlaybackFrequency();
		}
	}

	constexpr void updateDelta() {
		delta = static_cast<double>(frequency * currentTableSize * sampleRateInv);
	}


	ParamType sampleRateInv{};
	ParamType frequency{};
	double delta{};
	double currentSamplePosition{};

	T value{};

	StorageType* wavetables{};
	Wavetable* currentTable{};
	size_t currentTableSize{};

	double topFreq{}, bottomFreq{}; // Frequency interval of current table
};
//
// Oscillator with phase in [0,1)
//
// template<class Wavetable, class StorageType = std::vector<Wavetable>, class TableSelector = ForwardSearchTableSelector>
// class WavetableOscillator1
//{
// public:
//	using wavetable_type = Wavetable;
//	using value_type = typename Wavetable::value_type;
//	using T = typename Wavetable::value_type;
//	using multiwavetable_type = StorageType;
//
//	constexpr WavetableOscillator1(multiwavetable_type* wavetables, T sampleRate, T frequency)
//		: sampleRateInv(T(1.) / sampleRate),
//		  frequency(frequency) {
//		setTable(wavetables);
//	}
//
//	constexpr WavetableOscillator1(T sampleRate) : sampleRateInv(T(1.) / sampleRate) {}
//
//	constexpr void setTable(multiwavetable_type* wavetables) {
//		this->wavetables = wavetables;
//		topFreq = bottomFreq = T(0);
//		setFrequency(frequency);
//	}
//
//	constexpr void setFrequency(T frequency) {
//		this->frequency = frequency;
//		assert(frequency * sampleRateInv < 1 && "The frequency needs to be lower that the sample rate");
//		selectTable();
//		updateDelta();
//	}
//
//	/// @brief Increment the oscillator by one step and get the current value.
//	/// @return current value
//	constexpr T operator++() {
//		phase += phaseInc;
//		if (phase >= 1.0) {
//			phase -= 1.0;
//		}
//		value = (*currentTable)(phase * static_cast<T>(currentTableSize));
//		return value;
//	}
//
//	/// @brief Increment the oscillator by one step and get the former current value.
//	/// @return former value
//	constexpr T operator++(int) {
//		const auto tmp = value;
//		phase += phaseInc;
//		if (phase >= 1.0) {
//			phase -= 1.0;
//		}
//		value = (*currentTable)(phase * static_cast<T>(currentTableSize));
//		return tmp;
//	}
//
//	/// @brief Get current value of the oscillator without changing its state.
//	/// @return current value
//	constexpr T operator()() const { return value; }
//
//	/// @brief Reset the position/phase to 0. Also updates the current value.
//	constexpr void retrigger() {
//		phase = 0;
//		value = (*currentTable)(phase * static_cast<T>(currentTableSize));
//	}
//
//	constexpr void reset() { retrigger(); }
//	constexpr Wavetable* getSelectedTable() const { return currentTable; }
//	constexpr T getFrequency() const { return frequency; }
//	constexpr T getSampleRate() const { return 1.0 / sampleRateInv; }
//
// private:
//	constexpr void selectTable() {
//		// chances are good that the frequency has changed just a little bit, so
//		// only select table if antialiasing condition really does not hold.
//		if (frequency <= topFreq && frequency > bottomFreq) return;
//
//		assert(wavetables);
//		assert(wavetables->size() > 0);
//		auto it = TableSelector::selectTable(wavetables->begin(), wavetables->end(), frequency);
//
//		if (it == wavetables->end()) { // No table can be selected without aliasing -> then we just get aliasing
//			it--;
//		}
//
//		auto newTable = &(*it);
//
//		currentTable = newTable;
//		currentTableSize = currentTable->size();
//		assert(currentTableSize > 0 && "Size of wavetables may not be zero");
//		value = (*currentTable)(phase * static_cast<T>(currentTableSize));
//
//
//		topFreq = currentTable->getMaximumPlaybackFrequency();
//		if (it == wavetables->begin()) {
//			bottomFreq = T(0);
//		} else {
//			bottomFreq = (it - 1)->getMaximumPlaybackFrequency();
//		}
//	}
//
//	constexpr void updateDelta() {
//		phaseInc = frequency * sampleRateInv;
//	}
//
//
//	T sampleRateInv{};
//	T frequency{};
//	T phaseInc{};
//
//	T phase{};
//	T value{};
//
//	multiwavetable_type* wavetables{};
//	Wavetable* currentTable{};
//	size_t currentTableSize{};
//
//	T topFreq{}, bottomFreq{}; // Frequency interval of current table
//};
//
//
//
// Oscillator with fixed point phase counting
// template<class Wavetable, class StorageType = std::vector<Wavetable>, class TableSelector = ForwardSearchTableSelector>
// class WavetableOscillator2
//{
// public:
//	using wavetable_type = Wavetable;
//	using value_type = typename Wavetable::value_type;
//	using T = typename Wavetable::value_type;
//	using multiwavetable_type = StorageType;
//
//	constexpr WavetableOscillator2(multiwavetable_type* wavetables, T sampleRate, T frequency)
//		: sampleRateInv(T(1.) / sampleRate),
//		  frequency(frequency) {
//		setTable(wavetables);
//	}
//
//	constexpr WavetableOscillator2(T sampleRate) : sampleRateInv(T(1.) / sampleRate) {}
//
//	constexpr void setTable(multiwavetable_type* wavetables) {
//		this->wavetables = wavetables;
//		topFreq = bottomFreq = T(0);
//		setFrequency(frequency);
//	}
//
//	constexpr void setFrequency(T frequency) {
//		this->frequency = frequency;
//		assert(frequency * sampleRateInv < 1 && "The frequency needs to be lower that the sample rate");
//		selectTable();
//		updateDelta();
//	}
//
//	/// @brief Increment the oscillator by one step and get the current value.
//	/// @return current value
//	constexpr T operator++() {
//		phase += phaseInc;
//		// value = (*currentTable)(phase.to<T>() * (currentTableSize));
//		value = (*currentTable)(static_cast<T>(phase * inv_multiplicator));
//		return value;
//	}
//
//	/// @brief Increment the oscillator by one step and get the former current value.
//	/// @return former value
//	constexpr T operator++(int) {
//		const auto tmp = value;
//		phase += phaseInc;
//		// value = (*currentTable)(phase.to<T>() * (currentTableSize));
//		value = (*currentTable)(static_cast<T>(phase * inv_multiplicator));
//		return tmp;
//	}
//
//	/// @brief Get current value of the oscillator without changing its state.
//	/// @return current value
//	constexpr T operator()() const { return value; }
//
//	/// @brief Reset the position/phase to 0. Also updates the current value.
//	constexpr void retrigger() {
//		phase = 0;
//		// value = (*currentTable)(phase.to<T>() * (currentTableSize));
//		value = (*currentTable)(static_cast<T>(phase * inv_multiplicator));
//	}
//
//	constexpr void reset() { retrigger(); }
//	constexpr Wavetable* getSelectedTable() const { return currentTable; }
//	constexpr T getFrequency() const { return frequency; }
//	constexpr T getSampleRate() const { return 1.0 / sampleRateInv; }
//
// private:
//	constexpr void selectTable() {
//		// chances are good that the frequency has changed just a little bit, so
//		// only select table if antialiasing condition really does not hold.
//		if (frequency <= topFreq && frequency > bottomFreq) return;
//
//		assert(wavetables);
//		assert(wavetables->size() > 0);
//		auto it = TableSelector::selectTable(wavetables->begin(), wavetables->end(), frequency);
//
//		if (it == wavetables->end()) { // No table can be selected without aliasing -> then we just get aliasing
//			it--;
//		}
//
//		auto newTable = &(*it);
//
//		currentTable = newTable;
//		currentTableSize = static_cast<double>(currentTable->size());
//		assert(currentTableSize > 0 && "Size of wavetables may not be zero");
//		// value = (*currentTable)(phase.to<T>() * (currentTableSize));
//		value = (*currentTable)(static_cast<T>(phase * inv_multiplicator));
//
//		inv_multiplicator = currentTableSize / multiplicator;
//
//		if (inv_multiplicator * static_cast<T>(std::numeric_limits<uint64_t>::max()) >= currentTableSize) {
//			std::cout << std::setprecision(22) << "Oh no " << inv_multiplicator * static_cast<double>(std::numeric_limits<uint64_t>::max()) << " " << (double)currentTableSize << "\n";
//		}
//		assert(inv_multiplicator * static_cast<double>(std::numeric_limits<uint64_t>::max()) < currentTableSize);
//
//		topFreq = currentTable->getMaximumPlaybackFrequency();
//		if (it == wavetables->begin()) {
//			bottomFreq = T(0);
//		} else {
//			bottomFreq = (it - 1)->getMaximumPlaybackFrequency();
//		}
//	}
//
//	constexpr void updateDelta() {
//		// phaseInc = Fixed(frequency * sampleRateInv);
//		phaseInc = static_cast<uint64_t>(multiplicator * static_cast<double>(frequency) * sampleRateInv);
//	}
//
//	using Fixed = Fixed<64, 0>;
//
//	T sampleRateInv{};
//	T frequency{};
//	// Fixed phaseInc, phase;
//	T value{};
//
//	const double multiplicator = (1ull << 63) * 2.0; // = 1 << 64
//	double inv_multiplicator;
//	uint64_t phaseInc{}, phase{};
//
//	multiwavetable_type* wavetables{};
//	Wavetable* currentTable{};
//	double currentTableSize{};
//
//	T topFreq{}, bottomFreq{}; // Frequency interval of current table
//};



/// @brief Wavetable oscillator for morphing between two wavetables. A parameter in the interval
///        [0, 1] is used to blend between the first and the second table.
///
/// @tparam WavetableOscillator Wavetable oscillator class
template<class WavetableOscillator>
class MorpingWavetableOscillator
{

public:
	using oscillator_type = WavetableOscillator;
	using value_type = typename WavetableOscillator::value_type;
	using param_type = typename WavetableOscillator::param_type;
	using wavetable_type = typename WavetableOscillator::wavetable_type;
	using multiwavetable_type = typename WavetableOscillator::multiwavetable_type;
	using T = value_type;


	constexpr MorpingWavetableOscillator() = default;
	constexpr MorpingWavetableOscillator(param_type sampleRate) : osc1(sampleRate), osc2(sampleRate) {}

	constexpr MorpingWavetableOscillator(multiwavetable_type* firstTable, multiwavetable_type* secondTable, T sampleRate, T frequency)
		: osc1(firstTable, sampleRate, frequency),
		  osc2(secondTable, sampleRate, frequency) {}

	constexpr void setTable(multiwavetable_type* firstTable, multiwavetable_type* secondTable) {
		osc1.setTable(firstTable);
		osc2.setTable(secondTable);
	}

	constexpr void setSampleRate(param_type sampleRate) {
		osc1.setSampleRate(sampleRate);
		osc2.setSampleRate(sampleRate);
	}

	constexpr void setFrequency(param_type frequency) {
		osc1.setFrequency(frequency);
		osc2.setFrequency(frequency);
	}

	constexpr void setParam(param_type param) {
		this->param = static_cast<T>(param);
	}

	constexpr T operator++() {
		return (T(1) - param) * ++osc1 + param * ++osc2;
	}

	constexpr T operator++(int) {
		return (T(1) - param) * osc1++ + param * osc2++;
	}

	constexpr T operator()() const {
		return (T(1) - param) * osc1() + param * osc2();
	}

	constexpr void retrigger() {
		osc1.retrigger();
		osc2.retrigger();
	}

	constexpr void reset() {
		osc1.reset();
		osc2.reset();
	}

	constexpr param_type getFrequency() const { return osc1.getFrequency(); }
	constexpr param_type getSampleRate() const { return osc1.getSampleRate(); }
	constexpr param_type getParam() const { return static_cast<param_type>(param); }

private:
	T param{};
	WavetableOscillator osc1, osc2;
};






template<class InputIterator, class FreqIterator, class T, int size, class multiwavetable_type>
void antialiase(InputIterator signalBegin,
	FreqIterator freqsBegin, FreqIterator freqsEnd,
	T samplerate,
	const FFTCalculator<T, size>& fftCalculator,
	multiwavetable_type& tables) {

	const size_t freqs_size = std::distance(freqsBegin, freqsEnd);
	assert(tables.size() == freqs_size);

	std::vector<std::array<T, size>> data(freqs_size);
	antialiase(signalBegin, freqsBegin, freqsEnd, data.begin(), samplerate, fftCalculator);

	for (size_t i = 0; i < tables.size(); ++i) {
		tables[i].setData(data[i].begin(), data[i].end(), freqsBegin[i]);
	}
}


template<class T, int size>
class Antialiaser
{
public:
	Antialiaser(T samplerate, const FFTCalculator<T, size>& fftCalculator) : samplerate(samplerate), fftCalculator(fftCalculator) {}

	template<class Container1, class Container2, class multiwavetable_type>
	void antialiase(const Container1& signal, const Container2& freqs, multiwavetable_type& tables) {
		std::vector<std::array<T, size>> data(freqs.size());
		Butterfly::antialiase(signal.begin(), freqs.begin(), freqs.end(), data.begin(), samplerate, fftCalculator);

		for (size_t i = 0; i < tables.size(); ++i) {
			tables[i].setData(data[i].begin(), data[i].end(), freqs[i]);
		}
	}

	template<class InputIterator, class FreqIterator, class multiwavetable_type>
	void antialiase(InputIterator signalBegin, FreqIterator freqsBegin, FreqIterator freqsEnd, multiwavetable_type& tables) {
		const size_t freqs_size = std::distance(freqsBegin, freqsEnd);
		assert(tables.size() == freqs_size);

		std::vector<std::array<T, size>> data(freqs_size);
		Butterfly::antialiase(signalBegin, freqsBegin, freqsEnd, data.begin(), samplerate, fftCalculator);

		for (size_t i = 0; i < tables.size(); ++i) {
			tables[i].setData(data[i].begin(), data[i].end(), freqsBegin[i]);
		}
	}


private:
	const T samplerate;
	const FFTCalculator<T, size>& fftCalculator;
};


}
