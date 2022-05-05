#ifdef HERMES_HAVE_CFITSIO

#include "hermes/neutralgas/RingData.h"

#include <cmath>
#include <iostream>

#include "hermes/Common.h"

namespace hermes { namespace neutralgas {

RingData::RingData(GasType gas) : type(gas) {
	if (gas == GasType::HI) {
		readDataFile("NHtot_rings.fits");
	}
	if (gas == GasType::H2) {
		readDataFile("WCOrings_COGAL.fits.gz");
	}
}

void RingData::readDataFile(const std::string &filename) {
	ffile = std::make_unique<FITSFile>(FITSFile(getDataPath("GasDensity/Remy18/" + filename)));
	ffile->openFile(FITS::READ);

	n_lon = ffile->readKeyValueAsInt("NAXIS1");
	n_lat = ffile->readKeyValueAsInt("NAXIS2");
	n_rings = ffile->readKeyValueAsInt("NAXIS3");

	crpix_lon = ffile->readKeyValueAsDouble("CRPIX1");
	min_lon = ffile->readKeyValueAsDouble("CRVAL1");
	delta_lon = ffile->readKeyValueAsDouble("CDELT1");
	crpix_lat = ffile->readKeyValueAsDouble("CRPIX2") ;
	min_lat = ffile->readKeyValueAsDouble("CRVAL2");
	delta_lat = ffile->readKeyValueAsDouble("CDELT2");

	std::cerr << "Number of rings: " << n_rings << std::endl;

	int firstElement = 1;
	int nElements = n_lon * n_lat * n_rings;
	dataVector = ffile->readImageAsFloat(firstElement, nElements);
}

GasType RingData::getGasType() const { return type; }

int RingData::getRingNumber() const { return n_rings; }

double RingData::getRawValue(int ring, const QDirection &dir) const {
//	QAngle lat = 180_deg - dir[0];
//	QAngle lon = 180_deg - dir[1];  // because the galactic centre of
//	                                // the ring model is in the middle of the map

//	int pxl_lat = static_cast<int>(round(lat / 180_deg * n_lat));
//	int pxl_lon = static_cast<int>(round(lon / 360_deg * n_lon));

	auto gal_dir = toGalCoord(dir);
	QAngle lat = gal_dir[0];
	QAngle lon = gal_dir[1];
	//if ((lon > 180_deg) && (min_lon<0))
	//	lon = lon - 360_deg;

	int pxl_lat = static_cast<int>(round((lat - min_lat * 1_deg) / (delta_lat * 1_deg)));
	int pxl_lon = static_cast<int>(round((lon - min_lon * 1_deg) / (delta_lon * 1_deg)));
	pxl_lat = pxl_lat + crpix_lat - 1;
	pxl_lon = pxl_lon + crpix_lon - 1;

     if ((pxl_lat < 0) || (pxl_lat > n_lat)) return 0;
     //if ((pxl_lon < 0) || (pxl_lon > n_lon)) return 0;
	// NAXIS1 x NAXIS2 x NAXIS3 => lon x lat x ring
	return dataVector[(ring * n_lat + pxl_lat) * n_lon + pxl_lon];
}

QColumnDensity RingData::getHIColumnDensityInRing(int ring, const QDirection &dir) const {
	// the data is given in cm^-2
	if (ring == 11) return QColumnDensity(0);
	return getRawValue(ring, dir) / 1_cm2;
}

QRingCOIntensity RingData::getCOIntensityInRing(int ring, const QDirection &dir) const {
	// the data is given in K km s^-2
	if (ring == 0)
		return (getRawValue(10, dir) + getRawValue(11, dir)) * 1_K * 1_km / 1_s;
	else if (ring == 11)
		return QRingCOIntensity(0);
	else
		return getRawValue(ring - 1, dir) * 1_K * 1_km / 1_s;
}

}}  // namespace hermes::neutralgas

#endif  // HERMES_HAVE_CFITSIO