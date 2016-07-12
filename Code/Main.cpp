/*
-----------------------------------------------------------------------------
OpenGL Tutorial
VOXAR Labs
Computer Science Center - CIn
Federal University of Pernambuco - UFPE
http://www.cin.ufpe.br/~voxarlabs

-----------------------------------------------------------------------------
*/

#pragma comment( linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"" )

#include "openGL_tutorial.h"
#include <cstring>
#include <string>
#define PI  3.14159265

using namespace std;

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
Triangulo* triangulos; // armazena os IDs dos vertices de cada triangulo.
Camera c;
Luz luz;
int num_pontos, num_triangulos;
int auxObj = 0;
float z_buffer[WIDTH][HEIGHT];
float y_auxiliar;//pontos para realizar as rotações
float x_auxiliar;
float z_auxiliar;
double angulo;
char caminhoArquivo[256] = "Cameras/";
char caminhoArquivoObj[256] = "Objetos/";
char NomeCamera[256] = "01_Camera.cfg";
char NomeObjeto[256] = "01_Objeto.byu";

Ponto3D centroideGlobal;
GLfloat window_width = 600.0;
GLfloat window_height = 600.0;
GLdouble orthoMatrix[16];
GLfloat identidade[16];



void ortoMatriz(GLdouble left, GLdouble right, GLdouble bot, GLdouble top, GLdouble zNear, GLdouble zFar) {

	orthoMatrix[0] = (2 / (right - left));
	orthoMatrix[1] = 0;
	orthoMatrix[2] = 0;
	orthoMatrix[3] = 0;
	orthoMatrix[4] = 0;
	orthoMatrix[5] = (2 / (top - bot));
	orthoMatrix[6] = 0;
	orthoMatrix[7] = 0;
	orthoMatrix[8] = 0;
	orthoMatrix[9] = 0;
	orthoMatrix[10] = ((-2) / (zFar - zNear));
	orthoMatrix[11] = 0;
	orthoMatrix[12] = -((right + left) / (right - left));
	orthoMatrix[13] = -((top + bot) / (top - bot));
	orthoMatrix[14] = -((zFar + zNear) / (zFar - zNear));
	orthoMatrix[15] = 1;
	
	
		
}

// calcula o produto vetorial (a.k.a. produto cruzado) entre 2 vetores e armazena o resultado em "resultado".
static void produto_vetorial(Ponto3D* resultado, Ponto3D* p1, Ponto3D* p2)
{
	resultado->x = p1->y * p2->z - p1->z * p2->y;
	resultado->y = p1->z * p2->x - p1->x * p2->z;
	resultado->z = p1->x * p2->y - p1->y * p2->x;
}

// normaliza o vetor "v".
static void normalizar(Ponto3D* v)
{
	float l = (float)(sqrt(v->x * v->x + v->y * v->y + v->z * v->z));
	v->x /= l;
	v->y /= l;
	v->z /= l;
}

// calcula a normal de um triangulo, normaliza-a e armazena o resultado em "resultado".
// para este projeto, nao nos preocupamos com o sentido horario/anti-horario dos vertices do triangulo.
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

// calcula o produto escalar entre 2 vetores
static float produto_escalar(Ponto3D* p1, Ponto3D* p2)
{
	return (p1->x*p2->x + p1->y*p2->y + p1->z*p2->z);
}

// lê um arquivo .cfg com parametros de camera e atribui os valores lidos às variáveis correspondentes.
// a camera eh devidamente inicializada na funcao "inicializar_camera".
static void ler_camera(char * path)
{
	FILE * p_arq = fopen(path, "r");

	fscanf(p_arq, "%f %f %f", &c.C.x, &c.C.y, &c.C.z);

	fscanf(p_arq, "%f %f %f", &c.N.x, &c.N.y, &c.N.z);

	fscanf(p_arq, "%f %f %f", &c.V.x, &c.V.y, &c.V.z);

	fscanf(p_arq, "%f %f %f", &c.d, &c.hx, &c.hy);

	fclose(p_arq);
}

// inicializa os parametros UVN da camera:
//     N: o vetor "look at" da camera, que parte da camera ao seu alvo, correspondendo ao eixo Oz.
//     V: o vetor "up" da camera, que corresponde ao eixo Oy.
//     U: o vetor "right side" da camera, que corresponde ao eixo Ox.
static void inicializar_camera()
{
	// normaliza N.
	normalizar(&c.N);

	// projetamos V ortogonalmente em N com gram-schmidt.
	float coeficiente = produto_escalar(&c.V, &c.N) / produto_escalar(&c.N, &c.N);
	c.V.x -= (coeficiente * c.N.x);
	c.V.y -= (coeficiente * c.N.y);
	c.V.z -= (coeficiente * c.N.z);

	// normaliza V.
	normalizar(&c.V);

	// U = N x V
	produto_vetorial(&c.U, &c.N, &c.V);
	// agora U, V e N sao uma base ortonormal do R3.
}

// lê um arquivo com parametros de iluminaçao e inicializa as variáveis correspondentes.
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

// lê um arquivo .byu que representa um objeto e inicializa as variáveis desse objeto.
static void ler_objeto(char * path)
{
	FILE * p_arq = fopen(path, "r");

	// leitura da quantidade de vertices e de triangulos.
	fscanf(p_arq, "%d %d", &num_pontos, &num_triangulos);

	pontos_objeto_vista = (Ponto3D *)malloc((num_pontos + 1) * sizeof(Ponto3D)); // indice 1

	pontos_objeto_tela = (Ponto2D *)malloc((num_pontos + 1) * sizeof(Ponto2D)); // indice 1

	normais_vertices = (Ponto3D *)malloc((num_pontos + 1) * sizeof(Ponto3D)); // indice 1
	// inicializando as normais dos vertices com o vetor nulo.
	for (int i = 1; i <= num_pontos; i++) normais_vertices[i].x = normais_vertices[i].y = normais_vertices[i].z = 0.0;

	triangulos = (Triangulo *)malloc((num_triangulos + 1)*sizeof(Triangulo)); // indice 1

	// leitura dos vertices.
	float x, y, z;
	for (int i = 1; i <= num_pontos; i++)
	{
		fscanf(p_arq, "%f %f %f", &x, &y, &z);

		pontos_objeto_vista[i].x = x;
		pontos_objeto_vista[i].y = y;
		pontos_objeto_vista[i].z = z;
	}

	// leitura dos triangulos.
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

static void calculaCentroide() //faz o cálculo do centroide, que seria atribuir a média 
							  //de cada coordenada dos pontos de vista à coordenada correspondente
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

static void rotacaoZ() //realiza a rotação no eixo OZ, calculando previamente o seno e o cosseno do ângulo (em graus)
					//para depois aplicar ao x e ao y dos pontos de vista a rotação 
{
	double cosseno = cos(angulo * PI / 180);
	double seno = sin(angulo * PI / 180);
	for (int i = 1; i <= num_pontos; i++)
	{
		x_auxiliar = pontos_objeto_vista[i].x;
		y_auxiliar = pontos_objeto_vista[i].y;
		pontos_objeto_vista[i].x = x_auxiliar*cosseno - y_auxiliar*seno;
		pontos_objeto_vista[i].y = x_auxiliar*seno + y_auxiliar*cosseno;
	}
}

static void rotacaoY() //calcula a rotação no eixo paralelo a OY, levando em conta translações necessárias 
						//usando as coordenadas do centroide
{
	calculaCentroide();
	double cosseno = cos(angulo * PI / 180);
	double seno = sin(angulo * PI / 180);

	for (int i = 1; i <= num_pontos; i++) {
		x_auxiliar = pontos_objeto_vista[i].x;
		z_auxiliar = pontos_objeto_vista[i].z;
		pontos_objeto_vista[i].x = x_auxiliar*cosseno + z_auxiliar*seno - centroideGlobal.x*cosseno - centroideGlobal.z*seno + centroideGlobal.x;
		pontos_objeto_vista[i].z = z_auxiliar*cosseno - x_auxiliar*seno + centroideGlobal.x*seno - centroideGlobal.z*cosseno + centroideGlobal.z;
	}
}

static void rotacaoX() //calcula a rotação no eixo paralelo a OX, levando em conta translações necessárias
						//usando as coordenadas do centroide
{
	calculaCentroide();
	double cosseno = cos(angulo * PI / 180);
	double seno = sin(angulo * PI / 180);

	for (int i = 1; i <= num_pontos; i++) {
		y_auxiliar = pontos_objeto_vista[i].y;
		z_auxiliar = pontos_objeto_vista[i].z;
		pontos_objeto_vista[i].y = y_auxiliar*cosseno - z_auxiliar*seno - centroideGlobal.y*cosseno + centroideGlobal.z*seno + centroideGlobal.y;
		pontos_objeto_vista[i].z = z_auxiliar*cosseno + y_auxiliar*seno - centroideGlobal.y*seno - centroideGlobal.z*cosseno + centroideGlobal.z;
	}
}

// Luz de coordenada mundial para coordenada de vista
static void mudanca_base_luz()
{

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

// Objeto de coordenadas mundiais (nas variaveis pontos_objeto_vista e centroideGlobal) para coordenada de vista
static void mudanca_base_objeto()
{
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
	//aplicando a mudança de base também no centroide
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

// calcula as normais dos vertices do objeto.
static void calcula_normais()
{
	Ponto3D normal;
	int v1, v2, v3;
	for (int i = 1; i <= num_triangulos; i++)
	{
		v1 = triangulos[i].v1;
		v2 = triangulos[i].v2;
		v3 = triangulos[i].v3;

		// obtem a normal normalizada do triangulo.
		normal_triangulo(&normal, &pontos_objeto_vista[v1], &pontos_objeto_vista[v2], &pontos_objeto_vista[v3]);

		// soma essa normal à normal de cada vértice.

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

// passa os pontos de coordenadas de vista (3d) para coordenadas de tela (2d)
static void calcula_pontos_tela()
{
	for (int i = 1; i <= num_pontos; i++)
	{
		// como visto em sala, projetamos os pontos para 2d de acordo com os parametros d, hx e hy da camera.
		// x' = d/hx * x/z
		pontos_objeto_tela[i].x = (c.d / c.hx) * (pontos_objeto_vista[i].x / pontos_objeto_vista[i].z);
		// y' = d/hy * y/z
		pontos_objeto_tela[i].y = (c.d / c.hy) * (pontos_objeto_vista[i].y / pontos_objeto_vista[i].z);
		// neste ponto estamos com o ponto em coordenadas 2d que independem da resoluçao da tela (intervalo [-1, 1])

		// agora parametrizamos para a resolucao da tela. (intervalos [0, window_width] para coordenada x e [0, window_height] para coordenada y):
		pontos_objeto_tela[i].x = (int)((pontos_objeto_tela[i].x + 1) * window_width / 2);
		pontos_objeto_tela[i].y = (int)((1 - pontos_objeto_tela[i].y) * window_height / 2);
	}
}

// inicializa z-buffer com +infinito em todas as posicoes
static void inicializa_z_buffer()
{
	for (int i = 0; i < WIDTH; i++)
	{
		for (int j = 0; j < HEIGHT; j++) z_buffer[i][j] = std::numeric_limits<float>::infinity();
	}
}

// ordena 3 pontos em ordem crescente da coordenada y.
static void ordenar_pontos(Ponto2D* desordenado[])
{
	Ponto2D* aux;

	// compara p1 e p2, ordenando-os entre si.
	if (desordenado[1]->y < desordenado[0]->y)
	{
		aux = desordenado[0];
		desordenado[0] = desordenado[1];
		desordenado[1] = aux;
	}

	// compara p1 e p3, ordenando-os entre si.
	if (desordenado[2]->y < desordenado[0]->y)
	{
		aux = desordenado[0];
		desordenado[0] = desordenado[2];
		desordenado[2] = aux;
	}
	// neste ponto, o ponto com menor y é mesmo p1.

	// por ultimo, compara p2 e p3, ordenando-os entre si, para que o maior y seja p3.
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

// preenche uma linha de um triangulo com pixels,
// e define a cor de cada pixel usando o modelo de iluminaçao de phong:
// I = Ka * Ia + Kd*(N.L) * (Od * IL) + Ks *(R.V)^n * Il
static void scanline(float xMin, float xMax, int yScan, Triangulo* t)
{
	Ponto3D p3d, vetorNormal, V, L, R, baricentrica;
	Ponto2D p;
	float r, g, b;
	float aux;
	// para cada x de xMin a xMax
	for (int x = xMin; x <= xMax; x++)
	{
		// nosso ponto 2d é (x, yScan)
		p.x = x;
		p.y = yScan;

		// calculamos as coordenadas baricentricas do ponto 2d em relaçao aos vertices 2d.
		coord_baricentricas(&p, &pontos_objeto_tela[t->v1], &pontos_objeto_tela[t->v2], &pontos_objeto_tela[t->v3], &baricentrica);

		// multiplicamos as coordenadas baricentricas do ponto 2d pelos vertices 3d originais,
		// obtendo uma aproximação para o ponto 3d.
		p3d.x = baricentrica.x * pontos_objeto_vista[t->v1].x + baricentrica.y * pontos_objeto_vista[t->v2].x + baricentrica.z * pontos_objeto_vista[t->v3].x;
		p3d.y = baricentrica.x * pontos_objeto_vista[t->v1].y + baricentrica.y * pontos_objeto_vista[t->v2].y + baricentrica.z * pontos_objeto_vista[t->v3].y;
		p3d.z = baricentrica.x * pontos_objeto_vista[t->v1].z + baricentrica.y * pontos_objeto_vista[t->v2].z + baricentrica.z * pontos_objeto_vista[t->v3].z;
		
		// consulta ao z-buffer.
		// primeiro, checamos os limites do array.
		if (x >= 0 && yScan >= 0 && x < WIDTH && yScan < HEIGHT 
			// depois, a consulta tradicional pra saber se o ponto deve ser visto por estar mais à frente.
			&& p3d.z < z_buffer[x][yScan] 
			// por ultimo, se o z é negativo, o ponto está atrás da câmera e não deve ser visto.
			&& p3d.z >= 0)
		{
			// atualizando o z-buffer
			z_buffer[x][yScan] = p3d.z;

			// multiplicando as coordenadas baricentricas do ponto 2d pelas normais dos vertices 3d originais,
			// para obter uma aproximação para a normal do ponto 3d.
			vetorNormal.x = baricentrica.x * normais_vertices[t->v1].x + baricentrica.y * normais_vertices[t->v2].x + baricentrica.z * normais_vertices[t->v3].x;
			vetorNormal.y = baricentrica.x * normais_vertices[t->v1].y + baricentrica.y * normais_vertices[t->v2].y + baricentrica.z * normais_vertices[t->v3].y;
			vetorNormal.z = baricentrica.x * normais_vertices[t->v1].z + baricentrica.y * normais_vertices[t->v2].z + baricentrica.z * normais_vertices[t->v3].z;

			// V = -P			
			V.x = -p3d.x;
			V.y = -p3d.y;
			V.z = -p3d.z;

			// L = Pl - P
			L.x = luz.Pl.x - p3d.x;
			L.y = luz.Pl.y - p3d.y;
			L.z = luz.Pl.z - p3d.z;

			// normaliza todos.
			normalizar(&vetorNormal);
			normalizar(&V);
			normalizar(&L);

			// se N.V < 0, a normal esta "pra dentro". invertemos ela fazendo N <- -N.
			if (produto_escalar(&vetorNormal, &V) < 0) {
				vetorNormal.x = -vetorNormal.x;
				vetorNormal.y = -vetorNormal.y;
				vetorNormal.z = -vetorNormal.z;
			}

			// colorir com a ambiental Ka, sempre!
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

			// se alguma coordenada de cor der acima de 255, truncamos em 255.
			if (r > 255.0f) r = 255.0f;
			if (g > 255.0f) g = 255.0f;
			if (b > 255.0f) b = 255.0f;

			// os resultados r, g e b estao no intervalo [0, 255]. para passar para o intervalo [0, 1] divide-se por 255.
			// essa sera a cor do ponto, representada como uma tripla de pontos flutuantes
			// no intervalo [0, 1].
			glColor3f(r / 255.0f, g / 255.0f, b / 255.0f);
			
			// dizendo ao opengl que sera desenhado um ou mais pontos.
			glBegin(GL_POINTS);
			// desenhamos o ponto, cujo tamanho eh 1 pixel:  (x, yScan) 
			glVertex2i(x, yScan);
			glEnd(); // dizendo ao opengl que terminou.
		}
	}
}

// percorre os pixels de um triangulo flat top. recebe como parametros:
//    p1, p2, p3 os vertices desse triangulo.
//    triangulo_2d, com o objetivo de passar ele como parametro para o scanline, que precisara recuperar os vertices 3d originais.
static void scan_triangulo_flat_top(Ponto2D* p1, Ponto2D* p2, Ponto2D* p3, Triangulo* triangulo_2d)
{
	// dx_esquerda eh o inverso do coeficiente angular da reta a esquerda e
	// dx_direita eh o inverso do coeficiente angular da reta a direita.
	float dx_esquerda, dx_direita, altura;
	float x1 = p1->x, x2 = p2->x, x3 = p3->x, y1 = p1->y, y3 = p3->y;
	float xMin, xMax;
	float aux;

	// para garantir que p1 é o ponto à esquerda:
	if (x2 < x1)
	{
		aux = x2;
		x2 = x1;
		x1 = aux;
	}

	// p1 e p2 tem o mesmo y entao tanto faz.
	altura = y3 - y1;

	/*
                            _
 p1 *............* p2       |
      .        .            |
   \   .      .   /         |  altura
    v   .    .   v          |
         .  .               |
          *                 _ 
         p3                 
	
	para percorrer cada reta, de cima para baixo, incrementamos o y em 1 e o x em (inverso do coeficiente angular).

	*/

	// inverso do coeficiente angular da reta da esquerda (reta p1p3):
	dx_esquerda = (x3 - x1) / altura;
	// inverso do coeficiente angular da reta da direita (reta p2p3):
	dx_direita = (x3 - x2) / altura;

	// ponto de partida para percorrer a reta da esquerda: p1
	xMin = (float)x1;
	// ponto de partida para percorrer a reta da direita: p2
	xMax = (float)x2 + (float)0.5;

	for (int yScan = y1; yScan <= y3; yScan++)
	{
		// preenchemos a linha atual
		scanline(xMin, xMax, yScan, triangulo_2d);
		// incrementamos o x à esquerda e à direita, e incrementamos o y de 1 em 1 para ir à linha seguinte.
		xMin += dx_esquerda;
		xMax += dx_direita;
	}
}

// percorre os pixels de um triangulo flat bottom. recebe como parametros:
//    p1, p2, p3 os vertices desse triangulo.
//    triangulo_2d, com o objetivo de passar ele como parametro para o scanline, que precisara recuperar os vertices 3d originais.
static void scan_triangulo_flat_bottom(Ponto2D* p1, Ponto2D* p2, Ponto2D* p3, Triangulo* triangulo_2d)
{
	// dx_esquerda eh o inverso do coeficiente angular da reta a esquerda e
	// dx_direita eh o inverso do coeficiente angular da reta a direita.
	float dx_esquerda, dx_direita, altura;
	float x1 = p1->x, x2 = p2->x, x3 = p3->x, y1 = p1->y, y3 = p3->y;
	float xMin, xMax;
	float aux;

	// para garantir que p2 é o ponto à esquerda:
	if (x3 < x2)
	{
		aux = x3;
		x3 = x2;
		x2 = aux;
	}

	// p2 e p3 tem o mesmo y entao tanto faz.
	altura = y3 - y1;

	/*

           * p1                 _
    /    .    .    \            |
   v    .      .    v           |
       .        .               |  altura
      .          .              |
  p2 *............* p3          _

   para percorrer cada reta, de cima para baixo, incrementamos o y em 1 e o x em (inverso do coeficiente angular).
   
	*/

	// inverso do coeficiente angular da reta da esquerda (reta p1p2):
	dx_esquerda = (x2 - x1) / altura;
	// inverso do coeficiente angular da reta da direita (reta p1p3):
	dx_direita = (x3 - x1) / altura;

	// ponto de partida para percorrer a reta da esquerda: p1
	xMin = (float)x1;
	// ponto de partida para percorrer a reta da direita: p1
	xMax = (float)x1;

	for (int yScan = y1; yScan <= y3; yScan++)
	{
		// preenchemos a linha atual
		scanline(xMin, xMax, yScan, triangulo_2d);
		// incrementamos o x à esquerda e à direita, e incrementamos o y de 1 em 1 para ir à linha seguinte.
		xMin += dx_esquerda;
		xMax += dx_direita;
	}
}

// conversao por varredura: percorre os pixels de um triangulo.
static void scan_conversion()
{
	// armazenam os vertices de cada triangulo
	int v1, v2, v3;

	// para cada triangulo,
	for (int i = 1; i <= num_triangulos; i++)
	{
		v1 = triangulos[i].v1;
		v2 = triangulos[i].v2;
		v3 = triangulos[i].v3;

		Ponto2D* p1 = &pontos_objeto_tela[v1];
		Ponto2D* p2 = &pontos_objeto_tela[v2];
		Ponto2D* p3 = &pontos_objeto_tela[v3];

		Ponto2D* desordenado[3] = { p1, p2, p3 };

		// ordena os vertices do menor y (mais em cima) ao maior y (mais embaixo)
		ordenar_pontos(desordenado);

		Ponto2D* menor = desordenado[0];
		Ponto2D* medio = desordenado[1];
		Ponto2D* maior = desordenado[2];
		/*

		FLAT TOP

		trata o caso simples em que existem 2 vertices mais em cima com o mesmo y

       -----
       \   /
        \ /
         v
		
		*/
		if (menor->y == medio->y) scan_triangulo_flat_top(menor, medio, maior, &triangulos[i]);

		/*

		FLAT BOTTOM

		trata o caso simples em que existem 2 vertices mais embaixo com o mesmo y
		
           ^
         /   \
        /     \
        -------

		*/
		else if (medio->y == maior->y) scan_triangulo_flat_bottom(menor, medio, maior, &triangulos[i]);


		// caso geral - divide o triangulo em dois: um flat top e um flat bottom.
		else
		{
			/*

menor->          *
               .   .
              .      .
             .         .
            .            .
medio->    * ------------  .   <- novo
                 .           .
                         .     
maior->                        *

			*/

			Ponto2D novo;
			
			// mesmo y de "medio".
			novo.y = medio->y;

			// para encontrar o x do ponto novo,
			// calculamos ele na reta "menor -> maior"
			// eq da reta:
			// x(t) = (maior.x - menor.x)t + menor.x
			// y(t) = (maior.y - menor.y)t + menor.y
			// sabendo que y(t) = medio.y:
			// medio.y = (maior.y - menor.y)t + menor.y
			// t = (medio.y - menor.y) / (maior.y - menor.y)
			// portanto x(t) = a linha abaixo:
			novo.x = menor->x + (int)(0.5 + (float)(medio->y - menor->y)*(float)(maior->x - menor->x) / (float)(maior->y - menor->y));
			// (int)(0.5 + f) arredonda o ponto flutuante f para o inteiro mais proximo.

			// percorremos um triangulo flat bottom e depois um flat top.
			scan_triangulo_flat_bottom(menor, &novo, medio, &triangulos[i]);
			scan_triangulo_flat_top(medio, &novo, maior, &triangulos[i]);
		}
	}
}

void mudar_objeto(char *cam_camera, char *cam_objeto) {
	ler_camera(cam_camera);
	inicializar_camera();
	ler_luz("iluminacao.txt");
	ler_objeto(cam_objeto);
	mudanca_base_luz();

	mudanca_base_objeto();

	glutPostRedisplay();
}

void integracao_java() {
	STARTUPINFOW si;
	PROCESS_INFORMATION pi;
	char comando[] = "javaw -jar pg.jar";
	wchar_t w_comando[20];
	mbstowcs(w_comando, comando, 18);
	LPWSTR lpw_comando = w_comando;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));
	if (CreateProcessW(NULL, lpw_comando, NULL, NULL, FALSE, CREATE_NO_WINDOW, NULL, NULL, &si, &pi))
	{
		WaitForSingleObject(pi.hProcess, INFINITE);
		CloseHandle(pi.hProcess);
		CloseHandle(pi.hThread);
	}

	FILE * p_arq = fopen("config.txt", "r");
	
	char acao;
	fscanf(p_arq, "%c", &acao);
	if (acao == 'r') {
		char eixo[3];
		fscanf(p_arq, " %lf", &angulo);
		fscanf(p_arq, " %s", eixo);

		if (eixo[1] == 'x') {
			rotacaoX();
		} else if(eixo[1] == 'y') {
			rotacaoY();
		} else if (eixo[1] == 'z') {
			rotacaoZ();
		}
		glutPostRedisplay();
	}
	else if (acao == 'm') {
		char cam_camera_f[1000], cam_objeto_f[1000];
		fscanf(p_arq, " %s", cam_camera_f);
		fscanf(p_arq, " %s", cam_objeto_f);
		
		mudar_objeto(cam_camera_f, cam_objeto_f);
	}
	else if (acao == 'q')
	{
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

// Função callback chamada para fazer o desenho
void Desenha()
{
	// cor de fundo: preta.
	glClearColor(0.0, 0.0, 0.0, 1.0);
	
	// limpa buffer de cor
	glClear(GL_COLOR_BUFFER_BIT);
	
	// carrega matriz identidade na matriz do modelo da cena
	glLoadIdentity();

	// tamanho de um ponto a ser desenhado eh 1 pixel!
	glPointSize(1.f);

	calcula_normais();
	calcula_pontos_tela();
	inicializa_z_buffer();
	scan_conversion();

	// libera recursos
	glFlush();

	// exibe a janela de interacao com o usuario.
	integracao_java();
}

void myKeyboardAscii(unsigned char key, int x, int y)
{
	if (key == ESC) 
	{
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
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	window_width = (GLfloat)w;
	window_height = (GLfloat)h;
	ortoMatriz(0, window_width, window_height, 0.f, -5.0, 5.0);
	glLoadMatrixd(orthoMatrix);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

// Programa Principal 
int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	strcat(caminhoArquivo, NomeCamera);
	strcat(caminhoArquivoObj, NomeObjeto);
	
	ler_camera(caminhoArquivo);
	inicializar_camera();
	ler_luz("iluminacao.txt");
	ler_objeto(caminhoArquivoObj);
	mudanca_base_luz();
	mudanca_base_objeto();

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	//setar modo de exibição, nesse caso buffer simples e modelo de cor RGB
	glutInitWindowSize(window_width, window_height);
	//tamanho da janela
	glutInitWindowPosition(10, 10);
	//onde a janela vai aparecer na tela do PC
	glutCreateWindow("janelinha");
	//criar um janela
	glutDisplayFunc(Desenha);
	glutReshapeFunc(myreshape);
	//callback da função que desenha na tela
	glutKeyboardFunc(myKeyboardAscii);
	glMatrixMode(GL_MODELVIEW); // estou alterando a matriz do modelo da cena
	glLoadIdentity(); // carrega identidade na matriz do modelo da cena
	glutMainLoop();
	//começa a execução da maquina de estados do glut/opengl que controla as funções
	//de callback (controlador de mouse, teclado, callback de controle de tela, etc).
}
