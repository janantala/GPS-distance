###	Run: python distance.py

###===================================================================================================
### get distance in metres between 2 points:
### Vincenty Formula http://www.movable-type.co.uk/scripts/latlong-vincenty.html

import math

def getDistance(lat1, lon1, lat2, lon2):

	a = 6378137.0 
	b = 6356752.314245 
	f = 1 / 298.257223563

	L = math.radians(lon2 - lon1)

	U1 = math.atan((1 - f) * math.tan(math.radians(lat1)));
	U2 = math.atan((1 - f) * math.tan(math.radians(lat2)));
	sinU1 = math.sin(U1)
	cosU1 = math.cos(U1)
	sinU2 = math.sin(U2)
	cosU2 = math.cos(U2)
	cosSqAlpha = float()
	sinSigma = float()
	cos2SigmaM = float()
	cosSigma = float()
	sigma = float()



	# l == lambda
	l = L
	lambdaP  = float() 
	iterLimit = 100


	while True: 

		sinLambda = math.sin(l)
		cosLambda = math.cos(l)
		sinSigma = math.sqrt((cosU2 * sinLambda) * (cosU2 * sinLambda) + (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) * (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda))
		if (sinSigma == 0):
			return 0;

		cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda
		sigma = math.atan2(sinSigma, cosSigma)
		sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
		cosSqAlpha = 1 - sinAlpha * sinAlpha
		cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha

		C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha))
		lambdaP = l
		l = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma	* (-1 + 2 * cos2SigmaM * cos2SigmaM)));
	
	
		if (iterLimit == 0) or ((math.fabs(l - lambdaP) > 1e-12) and (iterLimit > 0)):
			break

		iterLimit = iterLimit - 1

	# end while

	if (iterLimit == 0):
		return 0

	uSq = cosSqAlpha * (a * a - b * b) / (b * b)
	A = 1 + uSq / 16384	* (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
	B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
	deltaSigma = B * sinSigma * (cos2SigmaM + B / 4	* (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) - B / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)))
	
	s = b * A * (sigma - deltaSigma)
	return s

###===================================================================================================

print getDistance(48.154563, 17.072561, 48.154564, 17.072562) #     0.13378944235175813
print getDistance(48.154563, 17.072561, 48.158800, 17.064064) #   788.4148295245418
print getDistance(48.148636, 17.107558, 48.208810, 16.372477) # 55073.682463659854
