#ifndef _VETOR3F_H_
#define _VETOR3F_H_

#include <cmath>

// vetor 3D representado como 3 floats x, y, z
class Vetor3f
{
public:
	Vetor3f(float x, float y, float z);

	float x, y, z;

	// soma de vetores
	Vetor3f operator+(const Vetor3f& outro);

	static float prod_interno(const Vetor3f& v1, const Vetor3f& v2);
	static Vetor3f prod_vetorial(const Vetor3f& v1, const Vetor3f& v2);

	// normalizar este vetor
	void normalizar();
};

#endif
