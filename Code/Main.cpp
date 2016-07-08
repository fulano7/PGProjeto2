/*
-----------------------------------------------------------------------------
OpenGL Tutorial
VOXAR Labs
Computer Science Center - CIn
Federal University of Pernambuco - UFPE
http://www.cin.ufpe.br/~voxarlabs

-----------------------------------------------------------------------------
*/

#include "openGL_tutorial.h"
#define PI  3.14159265

struct Triangulo
{
	int v1, v2, v3;
};

struct Ponto3D
{
	float x, y, z;
};

struct Ponto2D
{
	float x, y;
};

struct Camera
{
	Ponto3D C, N, V, U;
	float d, hx, hy;
};

struct Cor
{
	float r, g, b;
};

struct Luz
{
	Ponto3D Pl;
	float Od[3];
	Cor Ia, Il;
	float ka, kd, ks, n;
};

Ponto3D* pontos_objeto_vista;
Ponto2D* pontos_objeto_tela;
Ponto3D* normais_vertices;
Triangulo* triangulos;
Camera c;
Luz luz;
int num_pontos, num_triangulos;
float z_buffer[WIDTH][HEIGHT];
float theta_x, theta_y, theta_z_camera;
float y_auxiliar;//pontos para realizar as rotações
float x_auxiliar;
float z_auxiliar;
double angulo;
Ponto3D centroideGlobal;
GLfloat window_width = 600.0;
GLfloat window_height = 600.0;

static void produto_vetorial(Ponto3D* resultado, Ponto3D* p1, Ponto3D* p2)
{
	resultado->x = p1->y * p2->z - p1->z * p2->y;
	resultado->y = p1->z * p2->x - p1->x * p2->z;
	resultado->z = p1->x * p2->y - p1->y * p2->x;
}

static void normalizar(Ponto3D* v)
{
	float l = (float)(sqrt(v->x * v->x + v->y * v->y + v->z * v->z));
	v->x /= l;
	v->y /= l;
	v->z /= l;
}

static void normal_triangulo(Ponto3D* resultado, Ponto3D* v1, Ponto3D* v2, Ponto3D* v3)
{
	Ponto3D v21, v31;

	v21.x = v2->x - v1->x;
	v21.y = v2->y - v1->y;
	v21.z = v2->z - v1->z;

	v31.x = v3->x - v1->x;
	v31.y = v3->y - v1->y;
	v31.z = v3->z - v1->z;

	produto_vetorial(resultado, &v21, &v31);

	normalizar(resultado);

}

static float produto_escalar(Ponto3D* p1, Ponto3D* p2)
{
	return (p1->x*p2->x + p1->y*p2->y + p1->z*p2->z);
}

static void ler_camera(char * path)
{
	FILE * p_arq = fopen(path, "r");

	fscanf(p_arq, "%f %f %f", &c.C.x, &c.C.y, &c.C.z);

	fscanf(p_arq, "%f %f %f", &c.N.x, &c.N.y, &c.N.z);

	fscanf(p_arq, "%f %f %f", &c.V.x, &c.V.y, &c.V.z);

	fscanf(p_arq, "%f %f %f", &c.d, &c.hx, &c.hy);

	fclose(p_arq);
}

static void inicializar_camera()
{
	normalizar(&c.N);

	// ortogonalizar V
	float coeficiente = produto_escalar(&c.V, &c.N) / produto_escalar(&c.N, &c.N);
	c.V.x -= (coeficiente * c.N.x);
	c.V.y -= (coeficiente * c.N.y);
	c.V.z -= (coeficiente * c.N.z);

	normalizar(&c.V);

	// U = N x V
	produto_vetorial(&c.U, &c.N, &c.V);
}

static void ler_luz(char * path)
{
	FILE * p_arq = fopen(path, "r");

	fscanf(p_arq, "%f %f %f", &luz.Pl.x, &luz.Pl.y, &luz.Pl.z);

	fscanf(p_arq, "%f", &luz.ka);

	fscanf(p_arq, "%f %f %f", &luz.Ia.r, &luz.Ia.g, &luz.Ia.b);

	fscanf(p_arq, "%f", &luz.kd);

	fscanf(p_arq, "%f %f %f", &luz.Od[0], &luz.Od[1], &luz.Od[2]);

	fscanf(p_arq, "%f", &luz.ks);

	fscanf(p_arq, "%f %f %f", &luz.Il.r, &luz.Il.g, &luz.Il.b);

	fscanf(p_arq, "%f", &luz.n);

	fclose(p_arq);
}

static void ler_objeto(char * path)
{
	FILE * p_arq = fopen(path, "r");

	fscanf(p_arq, "%d %d", &num_pontos, &num_triangulos);

	pontos_objeto_vista = (Ponto3D *)malloc((num_pontos + 1) * sizeof(Ponto3D)); // indice 1

	pontos_objeto_tela = (Ponto2D *)malloc((num_pontos + 1) * sizeof(Ponto2D)); // indice 1

	normais_vertices = (Ponto3D *)malloc((num_pontos + 1) * sizeof(Ponto3D)); // indice 1
	for (int i = 1; i <= num_pontos; i++) normais_vertices[i].x = normais_vertices[i].y = normais_vertices[i].z = 0.0;

	triangulos = (Triangulo *)malloc((num_triangulos + 1)*sizeof(Triangulo)); // indice 1

	float x, y, z;
	for (int i = 1; i <= num_pontos; i++)
	{
		fscanf(p_arq, "%f %f %f", &x, &y, &z);

		pontos_objeto_vista[i].x = x;
		pontos_objeto_vista[i].y = y;
		pontos_objeto_vista[i].z = z;
	}

	int v1, v2, v3;

	for (int i = 1; i <= num_triangulos; i++)
	{
		fscanf(p_arq, "%d %d %d", &v1, &v2, &v3);

		triangulos[i].v1 = v1;
		triangulos[i].v2 = v2;
		triangulos[i].v3 = v3;
	}

	fclose(p_arq);
}

static void calculaCentroide()
{
	Ponto3D centroide;
	centroide.x = centroide.y = centroide.z = 0.0f;
	for (int i = 1; i <= num_pontos; i++) {
		centroide.x += pontos_objeto_vista[i].x;
		centroide.y += pontos_objeto_vista[i].y;
		centroide.z += pontos_objeto_vista[i].z;
	}
	centroide.x = (centroide.x / num_pontos);
	centroide.y = (centroide.y / num_pontos);
	centroide.z = (centroide.z / num_pontos);

	centroideGlobal = centroide;
}

static void rotacaoZ()
{
	double cosseno = cos(angulo * 180 / PI);
	double seno = sin(angulo * 180 / PI);
	for (int i = 1; i <= num_pontos; i++)
	{
		x_auxiliar = pontos_objeto_vista[i].x;
		y_auxiliar = pontos_objeto_vista[i].y;
		pontos_objeto_vista[i].x = x_auxiliar*cosseno - y_auxiliar*seno;
		pontos_objeto_vista[i].y = x_auxiliar*seno + y_auxiliar*cosseno;
	}
}

static void rotacaoY()
{
	calculaCentroide();
	double cosseno = cos(angulo * 180 / PI);
	double seno = sin(angulo * 180 / PI);

	for (int i = 1; i <= num_pontos; i++) {
		x_auxiliar = pontos_objeto_vista[i].x;
		z_auxiliar = pontos_objeto_vista[i].z;
		pontos_objeto_vista[i].x = x_auxiliar*cosseno + z_auxiliar*seno - centroideGlobal.x*cosseno - centroideGlobal.z*seno + centroideGlobal.x;
		pontos_objeto_vista[i].z = z_auxiliar*cosseno - x_auxiliar*seno + centroideGlobal.x*seno - centroideGlobal.z*cosseno + centroideGlobal.z;
	}
}

static void rotacaoX()
{
	calculaCentroide();
	double cosseno = cos(angulo * 180 / PI);
	double seno = sin(angulo * 180 / PI);

	for (int i = 1; i <= num_pontos; i++) {
		y_auxiliar = pontos_objeto_vista[i].y;
		z_auxiliar = pontos_objeto_vista[i].z;
		pontos_objeto_vista[i].y = y_auxiliar*cosseno - z_auxiliar*seno - centroideGlobal.y*cosseno + centroideGlobal.z*seno + centroideGlobal.y;
		pontos_objeto_vista[i].z = z_auxiliar*cosseno + y_auxiliar*seno - centroideGlobal.y*seno - centroideGlobal.z*cosseno + centroideGlobal.z;
	}
}

static void mudanca_base_luz()
{
	// Luz de coordenada mundial para coordenada de vista

	float x_vista, y_vista, z_vista;

	// P <- P-C
	luz.Pl.x -= c.C.x;
	luz.Pl.y -= c.C.y;
	luz.Pl.z -= c.C.z;

	// multiplicando pela matriz de mudanca de base
	x_vista = produto_escalar(&luz.Pl, &c.U);
	y_vista = produto_escalar(&luz.Pl, &c.V);
	z_vista = produto_escalar(&luz.Pl, &c.N);

	luz.Pl.x = x_vista;
	luz.Pl.y = y_vista;
	luz.Pl.z = z_vista;
}

static void mudanca_base_objeto()
{
	// Objeto de coordenadas mundiais (nas variaveis pontos_objeto_vista e centroideGlobal) para coordenada de vista
	float x_vista, y_vista, z_vista;

	for (int i = 1; i <= num_pontos; i++)
	{
		// P <- P-C
		pontos_objeto_vista[i].x -= c.C.x;
		pontos_objeto_vista[i].y -= c.C.y;
		pontos_objeto_vista[i].z -= c.C.z;

		// multiplicando pela matriz de mudanca de base
		x_vista = produto_escalar(&pontos_objeto_vista[i], &c.U);
		y_vista = produto_escalar(&pontos_objeto_vista[i], &c.V);
		z_vista = produto_escalar(&pontos_objeto_vista[i], &c.N);

		pontos_objeto_vista[i].x = x_vista;
		pontos_objeto_vista[i].y = y_vista;
		pontos_objeto_vista[i].z = z_vista;
	}

	centroideGlobal.x -= c.C.x;
	centroideGlobal.y -= c.C.y;
	centroideGlobal.z -= c.C.z;

	x_vista = produto_escalar(&centroideGlobal, &c.U);
	y_vista = produto_escalar(&centroideGlobal, &c.V);
	z_vista = produto_escalar(&centroideGlobal, &c.N);

	centroideGlobal.x = x_vista;
	centroideGlobal.y = y_vista;
	centroideGlobal.z = z_vista; 
}

static void calcula_normais()
{
	Ponto3D normal;
	int v1, v2, v3;
	for (int i = 1; i <= num_triangulos; i++)
	{
		v1 = triangulos[i].v1;
		v2 = triangulos[i].v2;
		v3 = triangulos[i].v3;

		// ja vem normalizado.
		normal_triangulo(&normal, &pontos_objeto_vista[v1], &pontos_objeto_vista[v2], &pontos_objeto_vista[v3]);

		normais_vertices[v1].x += normal.x;
		normais_vertices[v2].x += normal.x;
		normais_vertices[v3].x += normal.x;

		normais_vertices[v1].y += normal.y;
		normais_vertices[v2].y += normal.y;
		normais_vertices[v3].y += normal.y;

		normais_vertices[v1].z += normal.z;
		normais_vertices[v2].z += normal.z;
		normais_vertices[v3].z += normal.z;
	}

	// normaliza as normais dos vertices
	for (int i = 1; i <= num_pontos; i++) normalizar(&normais_vertices[i]);
}

static void calcula_pontos_tela()
{
	for (int i = 1; i <= num_pontos; i++)
	{
		pontos_objeto_tela[i].x = (c.d / c.hx) * (pontos_objeto_vista[i].x / pontos_objeto_vista[i].z);
		pontos_objeto_tela[i].y = (c.d / c.hy) * (pontos_objeto_vista[i].y / pontos_objeto_vista[i].z);

		pontos_objeto_tela[i].x = (int)((pontos_objeto_tela[i].x + 1) * WIDTH / 2);
		pontos_objeto_tela[i].y = (int)((1 - pontos_objeto_tela[i].y) * HEIGHT / 2);
	}
}

// inicializa z-buffer com valor maximo de float em todas as posicoes
static void inicializa_z_buffer()
{
	for (int i = 0; i < WIDTH; i++)
	{
		for (int j = 0; j < HEIGHT; j++) z_buffer[i][j] = std::numeric_limits<float>::infinity();
	}
}

// testar
static void ordenar_pontos(Ponto2D* desordenado[])
{
	Ponto2D* aux;
	if (desordenado[1]->y < desordenado[0]->y)
	{
		aux = desordenado[0];
		desordenado[0] = desordenado[1];
		desordenado[1] = aux;
	}
	if (desordenado[2]->y < desordenado[0]->y)
	{
		aux = desordenado[0];
		desordenado[0] = desordenado[2];
		desordenado[2] = aux;
	}
	if (desordenado[1]->y > desordenado[2]->y)
	{
		aux = desordenado[2];
		desordenado[2] = desordenado[1];
		desordenado[1] = aux;
	}
}

//subtraction of two points
static Ponto2D sub_points(Ponto2D *a, Ponto2D *b) {
	Ponto2D result;
	result.x = a->x - b->x;
	result.y = a->y - b->y;
	return result;
}

static float produto_escalar_2d(Ponto2D *a, Ponto2D *b) {
	return (a->x * b->x) + (a->y * b->y);
}

//barycentric coordinates
static void coord_baricentricas(Ponto2D *p, Ponto2D *ponto1, Ponto2D *ponto2, Ponto2D *ponto3, Ponto3D *baricentrica) {
	Ponto2D v0 = sub_points(ponto2, ponto1), v1 = sub_points(ponto3, ponto1), v2 = sub_points(p, ponto1);
	float d00 = produto_escalar_2d(&v0, &v0);
	float d01 = produto_escalar_2d(&v0, &v1);
	float d11 = produto_escalar_2d(&v1, &v1);
	float d20 = produto_escalar_2d(&v2, &v0);
	float d21 = produto_escalar_2d(&v2, &v1);

	float denom = d00 * d11 - d01 * d01;
	baricentrica->y = (d11 * d20 - d01 * d21) / denom;
	baricentrica->z = (d00 * d21 - d01 * d20) / denom;
	baricentrica->x = 1.0f - baricentrica->y - baricentrica->z;
}

static void scanline(float xMin, float xMax, int yScan, Triangulo* t)// Ka * Ia + Kd*(N.L) * (Od * IL) + Ks *(R.V)^n * IL
{
	Ponto3D p3d, vetorNormal, V, L, R, baricentrica;
	Ponto2D p;
	float r, g, b;
	float aux;
	for (int x = xMin; x <= xMax; x++)
	{
		p.x = x;
		p.y = yScan;
		coord_baricentricas(&p, &pontos_objeto_tela[t->v1], &pontos_objeto_tela[t->v2], &pontos_objeto_tela[t->v3], &baricentrica);

		p3d.x = baricentrica.x * pontos_objeto_vista[t->v1].x + baricentrica.y * pontos_objeto_vista[t->v2].x + baricentrica.z * pontos_objeto_vista[t->v3].x;
		p3d.y = baricentrica.x * pontos_objeto_vista[t->v1].y + baricentrica.y * pontos_objeto_vista[t->v2].y + baricentrica.z * pontos_objeto_vista[t->v3].y;
		p3d.z = baricentrica.x * pontos_objeto_vista[t->v1].z + baricentrica.y * pontos_objeto_vista[t->v2].z + baricentrica.z * pontos_objeto_vista[t->v3].z;
		
		// consulta ao z-buffer
		if (x >= 0 && yScan >= 0 && x < WIDTH && yScan < HEIGHT && p3d.z < z_buffer[x][yScan])
		{
			z_buffer[x][yScan] = p3d.z;
			vetorNormal.x = baricentrica.x * normais_vertices[t->v1].x + baricentrica.y * normais_vertices[t->v2].x + baricentrica.z * normais_vertices[t->v3].x;
			vetorNormal.y = baricentrica.x * normais_vertices[t->v1].y + baricentrica.y * normais_vertices[t->v2].y + baricentrica.z * normais_vertices[t->v3].y;
			vetorNormal.z = baricentrica.x * normais_vertices[t->v1].z + baricentrica.y * normais_vertices[t->v2].z + baricentrica.z * normais_vertices[t->v3].z;

			V.x = -p3d.x;
			V.y = -p3d.y;
			V.z = -p3d.z;
			L.x = luz.Pl.x - p3d.x;
			L.y = luz.Pl.y - p3d.y;
			L.z = luz.Pl.z - p3d.z;
			normalizar(&vetorNormal);
			normalizar(&V);
			normalizar(&L);

			if (produto_escalar(&vetorNormal, &V) < 0) {
				vetorNormal.x = -vetorNormal.x;
				vetorNormal.y = -vetorNormal.y;
				vetorNormal.z = -vetorNormal.z;
			}

			//colorir apenas com a ambiental KA
			r = luz.ka * luz.Ia.r;
			g = luz.ka * luz.Ia.g;
			b = luz.ka * luz.Ia.b;

			float nL = produto_escalar(&vetorNormal, &L);
			// colorir tambem com a difusa e a especular
			if (nL >= 0) {
				R.x = 2 * nL * (vetorNormal.x - L.x);
				R.y = 2 * nL * (vetorNormal.y - L.y);
				R.z = 2 * nL * (vetorNormal.z - L.z);
				normalizar(&R);
				aux = produto_escalar(&R, &V);
				r += (luz.kd * nL * luz.Od[0] * luz.Il.r) + (luz.ks * pow(aux, luz.n) * luz.Il.r);
				g += (luz.kd * nL * luz.Od[1] * luz.Il.g) + (luz.ks * pow(aux, luz.n) * luz.Il.g);
				b += (luz.kd * nL * luz.Od[2] * luz.Il.b) + (luz.ks * pow(aux, luz.n) * luz.Il.b);
			}

			if (r > 255.0f) r = 255.0f;
			if (g > 255.0f) g = 255.0f;
			if (b > 255.0f) b = 255.0f;

			glColor3f(r / 255.0f, g / 255.0f, b / 255.0f);
			glBegin(GL_POINTS);
			glVertex2i(x, yScan);
			glEnd();
		}



		// calcular coordenadas baricentricas
		// calcular ponto 3D
		// consultar z-buffer
		// iluminacao - pintar com phong
	}
}

// testar
// tem que passar o triangulo 2d para o scanline.
static void scan_triangulo_flat_top(Ponto2D* p1, Ponto2D* p2, Ponto2D* p3, Triangulo* triangulo_2d)
{
	float dx_esquerda, dx_direita, altura;
	float x1 = p1->x, x2 = p2->x, x3 = p3->x, y1 = p1->y, y3 = p3->y;
	float xMin, xMax;
	float aux;

	altura = y3 - y1;

	if (x2 < x1)
	{
		aux = x2;
		x2 = x1;
		x1 = aux;
	}

	dx_esquerda = (x3 - x1) / altura;
	dx_direita = (x3 - x2) / altura;

	xMin = (float)x1;
	xMax = (float)x2 + (float)0.5;

	for (int yScan = y1; yScan <= y3; yScan++)
	{
		scanline(xMin, xMax, yScan, triangulo_2d);
		xMin += dx_esquerda;
		xMax += dx_direita;
	}
}

// testar
// tem que passar o triangulo 2d para o scanline.
static void scan_triangulo_flat_bottom(Ponto2D* p1, Ponto2D* p2, Ponto2D* p3, Triangulo* triangulo_2d)
{
	float dx_esquerda, dx_direita, altura;
	float x1 = p1->x, x2 = p2->x, x3 = p3->x, y1 = p1->y, y3 = p3->y;
	float xMin, xMax;
	float aux;

	altura = y3 - y1;

	if (x3 < x2)
	{
		aux = x3;
		x3 = x2;
		x2 = aux;
	}

	dx_esquerda = (x2 - x1) / altura;
	dx_direita = (x3 - x1) / altura;

	xMin = (float)x1;
	xMax = (float)x1; // vertice de cima

	for (int yScan = y1; yScan <= y3; yScan++)
	{
		scanline(xMin, xMax, yScan, triangulo_2d);
		xMin += dx_esquerda;
		xMax += dx_direita;
	}
}

// testar
static void scan_conversion()
{
	int v1, v2, v3;
	for (int i = 1; i <= num_triangulos; i++)
	{
		v1 = triangulos[i].v1;
		v2 = triangulos[i].v2;
		v3 = triangulos[i].v3;

		Ponto2D* p1 = &pontos_objeto_tela[v1];
		Ponto2D* p2 = &pontos_objeto_tela[v2];
		Ponto2D* p3 = &pontos_objeto_tela[v3];

		Ponto2D* desordenado[3] = { p1, p2, p3 };

		ordenar_pontos(desordenado);

		Ponto2D* menor = desordenado[0];
		Ponto2D* medio = desordenado[1];
		Ponto2D* maior = desordenado[2];
		/*

		FLAT TOP

		------
		\    /
		\  /
		V

		*/
		if (menor->y == medio->y) scan_triangulo_flat_top(menor, medio, maior, &triangulos[i]);

		/*

		FLAT BOTTOM

		^
		/   \
		/     \
		-------

		*/


		else if (medio->y == maior->y) scan_triangulo_flat_bottom(menor, medio, maior, &triangulos[i]);


		// caso geral - divide o triangulo em 2 triangulos faceis de preencher
		else
		{
			Ponto2D novo;
			novo.x = menor->x + (int)(0.5 + (float)(medio->y - menor->y)*(float)(maior->x - menor->x) / (float)(maior->y - menor->y));
			novo.y = medio->y;
			scan_triangulo_flat_bottom(menor, &novo, medio, &triangulos[i]);
			scan_triangulo_flat_top(medio, &novo, maior, &triangulos[i]);
		}
	}
}

// Função callback chamada para fazer o desenho
void Desenha()
{
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();

	glPointSize(1.f);

	calcula_pontos_tela();
	inicializa_z_buffer();
	scan_conversion();
	printf("yay, pontos\n");
	// glBegin(GL_POINTS);
	
	//glColor3f(0.0f, 0.0f, 0.0f);
	/*for (int i = 1; i <= num_pontos; i++) {
		glVertex2f(pontos_objeto_tela[i].x, pontos_objeto_tela[i].y);
	}
	glEnd();*/

	// libera recursos
	glFlush();
	printf("acabou\n");
}

// Inicializa parâmetros de rendering
/*void Inicializa()
{
	// Define a cor de fundo da janela de visualização como preta
	glClearColor(1.0f, 0.0f, 1.0f, 1.0f);
	//para ver os parametros da função (e de qualquer outra) usar ctrl+shift+spacebar
	//dentro dos parênteses 
}*/


void myKeyboard(unsigned char key, int x, int y)
{
	if (key == GLUT_KEY_F1) {
		rotacaoX();
		glutPostRedisplay();
	}
	if (key == GLUT_KEY_F2) {
		rotacaoY();
		glutPostRedisplay();
	}
	if (key == GLUT_KEY_F3) {
		rotacaoZ();
		glutPostRedisplay();
	}
	if(key == ESC) {
		free(pontos_objeto_vista);
		pontos_objeto_vista = 0;
		free(pontos_objeto_tela);
		pontos_objeto_tela = 0;
		free(normais_vertices);
		normais_vertices = 0;
		free(triangulos);
		triangulos = 0;
		exit(0);
	}
}

void myreshape(GLsizei w, GLsizei h) {
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	window_width = (GLfloat)w;
	window_height = (GLfloat)h;
	glOrtho(0, window_width, window_height, 0.f, -5.0, 5.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

// Programa Principal 
int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	angulo = 30;
	
	ler_camera("Cameras/01_Camera.cfg");
	inicializar_camera();
	ler_luz("iluminacao.txt");
	ler_objeto("Objetos/01_Objeto.byu");
	printf("leu objeto\n");
	mudanca_base_luz();
	printf("mudou base luz\n");
	mudanca_base_objeto();
	printf("mudou base objeto\n");
	calcula_normais();
	printf("calculou normais\n");
	//calcula_pontos_tela();
	printf("calculou pontos de tela\n");
	//inicializa_z_buffer();
	printf("inicializou z-buffer\n");
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	//setar modo de exibição, nesse caso buffer simples e modelo de cor RGB
	glutInitWindowSize(window_width, window_height);
	//tamanho da janela
	glutInitWindowPosition(10, 10);
	//onde a janela vai aparecer na tela do PC
	glutCreateWindow("janelinha");
	//não lembro como criar 2 janelas, vejam ai isso se quiserem
	//criar um janela
	glutDisplayFunc(Desenha);
	glutReshapeFunc(myreshape);
	//callback da função que desenha na tela
	glutKeyboardFunc(myKeyboard);
	// Inicializa();
	//inicializar alguns parametros do glut (nessa caso a cor do fundo da tela).
	//cor que vai limpar o buffer
	glMatrixMode(GL_MODELVIEW); // estou alterando a matrix do modelo da cena
	glLoadIdentity();
	glutMainLoop();
	//começa a execução da maquina de estados do glut/opengl que controla as funções
	//de callback (controlador de mouse, teclado, callback de controle de tela, etc).
}
