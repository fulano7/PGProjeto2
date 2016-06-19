#ifndef _PONTO3F_H_
#define _PONTO3F_H_

#include "Vetor3f.h"

// ponto 3D representado como 3 floats x, y, z
class Ponto3f
{
public:
	Ponto3f(float x, float y, float z);

	float x, y, z;

	// subtracao de 2 pontos
	Vetor3f operator-(const Ponto3f& p2);
};

#endif // _PONTO3F_H_
