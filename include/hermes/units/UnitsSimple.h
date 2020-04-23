#ifndef HERMES_UNITSSIMPLE_H
#define HERMES_UNITSSIMPLE_H

#define HERMES_UNITSDISABLE

#include <array>
#include <cmath>

#include "UnitsSIPrefixes.h"

namespace hermes { namespace units {

typedef double QNumber;
typedef double QAngle;
typedef double QSolidAngle;
typedef double QLength;
typedef double QTime;
typedef double QEnergy;
typedef double QMField;
typedef double QMass;
typedef double QTemperature;
typedef double QIntensity;
typedef double QFrequency;
typedef double QECurrent;
typedef double QSubstance;
typedef double QLIntensity;
typedef double QArea;
typedef double QVolume;
typedef double QForce;
typedef double QMomentum;
typedef double QAMomentum;
typedef double QPressure;
typedef double QECharge;
typedef double QPower;
typedef double QEResistance;
typedef double QECapacitance;
typedef double QEPotential;
typedef double QPDensityPerEnergy;
typedef double QAcceleration;
typedef double QSpeed;
typedef double QDiffFlux;
typedef double QDiffIntensity;
typedef double QEmissivity;
typedef double QPDensity;
typedef double QInverseLength;
typedef double QColumnDensity;
typedef double QDispersionMeasure;
typedef double QRotationMeasure;
typedef double QDiffCrossSection;
typedef double QEnergyDensity;

typedef double QRMIntegral;
typedef double QSynchroConstant;
typedef double QPiZeroIntegral;
typedef double QRingCOIntensity;
typedef double QRingX0Unit;
typedef double QICInnerIntegral;
typedef double QGREmissivity;

template <int power>
constexpr double pow(const QNumber &num) {
	return std::pow(static_cast<double>(num), power);
}
constexpr double squared(const QNumber &num) { return num * num; }

}}  // namespace hermes::units

#endif  // HERMES_UNITSSIMPLE_H
