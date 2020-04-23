#ifndef HERMES_BREMSSTRAHLUNGSIMPLE_H
#define HERMES_BREMSSTRAHLUNGSIMPLE_H

#include <memory>

#include "hermes/CacheTools.h"
#include "hermes/interactions/DiffCrossSection.h"

namespace hermes { namespace interactions {

class BremsstrahlungSimple : public DifferentialCrossSection {
  private:
	bool cachingEnabled;
	std::unique_ptr<CacheStorageCrossSection> cache;

  public:
	BremsstrahlungSimple();

	void enableCaching();
	void disableCaching();

	QDiffCrossSection getDiffCrossSection(
	    const QEnergy &E_photon, const QEnergy &E_gamma) const override;
	QDiffCrossSection getDiffCrossSectionDirectly(const QEnergy &E_photon,
	                                              const QEnergy &E_gamma) const;

	QNumber ElwertFactor(const QNumber &beta_i, const QNumber &beta_f,
	                     int Z) const;
	QNumber xiFunc(const QNumber &T_electron_i, const QNumber &k, int Z,
	               int N) const;
	QNumber Phi_u(const QNumber &gamma_i, const QNumber &gamma_f,
	              const QNumber &k) const;

	inline double R_1(double q, double Z) const;
	inline double R_2(double q, double Z) const;
	QNumber I_Phi_1(const QNumber &delta, int Z, int N) const;
	QNumber I_Phi_2(const QNumber &delta, int Z, int N) const;
	QNumber Phi_1(const QNumber &gamma_i, const QNumber &gamma_f,
	              const QNumber &k, const QNumber &delta, int Z, int N) const;
	QNumber Phi_2(const QNumber &gamma_i, const QNumber &gamma_f,
	              const QNumber &k, const QNumber &delta, int Z, int N) const;

	QArea dsdk_LowEnergy(const QNumber &p_i, const QNumber &p_f, QNumber k,
	                     int Z) const;
	QArea dsdk_IntermediateEnergy(const QNumber &gamma_i,
	                              const QNumber &gamma_f, const QNumber &p_i,
	                              const QNumber &p_f, const QNumber &k,
	                              int Z) const;
	QArea dsdk_HighEnergy(const QNumber &gamma_i, const QNumber &gamma_f,
	                      const QNumber &k, int Z, int N) const;
};

}}  // namespace hermes::interactions

#endif  // HERMES_BREMSSTRAHLUNGSIMPLE_H
