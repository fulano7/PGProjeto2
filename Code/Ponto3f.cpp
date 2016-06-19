#include "Ponto3f.h"

Ponto3f::Ponto3f(float x, float y, float z) : x(x), y(y), z(z) {};

Vetor3f Ponto3f::operator-(const Ponto3f& p2)
{
	return Vetor3f(this->x - p2.x, this->y - p2.y, this->z - p2.z);
}
