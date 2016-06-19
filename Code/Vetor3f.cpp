#include "Vetor3f.h"

Vetor3f::Vetor3f(float x, float y, float z) : x(x), y(y), z(z){};

Vetor3f Vetor3f::operator+(const Vetor3f& v2)
{
	return Vetor3f(x + v2.x, y + v2.y, z + v2.z);
}

static float prod_interno(const Vetor3f& v1, const Vetor3f& v2)
{
	return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

static Vetor3f prod_vetorial(const Vetor3f& v1, const Vetor3f& v2)
{
	return Vetor3f(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

void Vetor3f::normalizar()
{
	float norma = (float)sqrt(x*x + y*y + z*z);
	x /= norma;
	y /= norma;
	z /= norma;
}
