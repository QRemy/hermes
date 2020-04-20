#ifndef HERMES_DIFFERENTIALCROSSSECTION_H
#define HERMES_DIFFERENTIALCROSSSECTION_H

#include "hermes/Units.h"
#include "hermes/ParticleID.h"

namespace hermes { namespace interactions {

using namespace units;

class DifferentialCrossSection {
protected:
	bool cachingEnabled;
public:
	DifferentialCrossSection(bool cachingEnabled = false);
	~DifferentialCrossSection();
	
	void enableCaching();
	void disableCaching();
	
	virtual QDifferentialCrossSection getDiffCrossSection(
			const QEnergy &E_proton,
			const QEnergy &E_gamma) const;
        virtual QDifferentialCrossSection getDiffCrossSection(
                        const QEnergy &E_electron,
                        const QEnergy &E_photon,
                        const QEnergy &E_gamma) const;
	virtual QNumber getSigma(
			const PID &projectile,
			const PID &target) const;
};

} // namespace interactions
} // namespace hermes

#endif // HERMES_DIFFERENTIALCROSSSECTION_H
