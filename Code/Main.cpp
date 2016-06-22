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

struct Triangulo
{
	int v1, v2, v3;
};

struct Ponto3D
{
	float x, y, z;
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
	Ponto3D Pl, Od;
	Cor Ia, Il;
	float ka, kd, ks, n;
};

Ponto3D* pontos_objeto_mundo;
Ponto3D* pontos_objeto_vista;
Ponto3D* pontos_objeto_tela;
Ponto3D* normais_vertices;
Triangulo* triangulos;
Camera c;
Luz luz;
int num_pontos, num_triangulos;

static Ponto3D* produto_vetorial(Ponto3D* p1, Ponto3D* p2)
{
	Ponto3D *v = (Ponto3D *)malloc(sizeof(Ponto3D));
	v->x = p1->y * p2->z - p1->z * p2->y;
	v->y = p1->z * p2->x - p1->x * p2->z;
	v->z = p1->x * p2->y - p1->y * p2->x;
	return v;
}

static void normalizar(Ponto3D* v)
{
	float l = (float)(sqrt(v->x * v->x + v->y * v->y + v->z * v->z));
	v->x /= l;
	v->y /= l;
	v->z /= l;
}

static Ponto3D* normal_triangulo(Ponto3D* v1, Ponto3D* v2, Ponto3D* v3)
{
	Ponto3D v21, v31;
	
	v21.x = v2->x - v1->x;
	v21.y = v2->y - v1->y;
	v21.z = v2->z - v1->z;

	v31.x = v3->x - v1->x;
	v31.y = v3->y - v1->y;
	v31.z = v3->z - v1->z;
	
	Ponto3D* resp = produto_vetorial(&v21, &v31);
	
	normalizar(resp);
	
	return resp;
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
	c.U = *(produto_vetorial(&c.N, &c.V));
}

static void ler_luz(char * path)
{
	FILE * p_arq = fopen(path, "r");

	fscanf(p_arq, "%f %f %f", &luz.Pl.x, &luz.Pl.y, &luz.Pl.z);

	fscanf(p_arq, "%f", &luz.ka);

	fscanf(p_arq, "%f %f %f", &luz.Ia.r, &luz.Ia.g, &luz.Ia.b);

	fscanf(p_arq, "%f", &luz.kd);

	fscanf(p_arq, "%f %f %f", &luz.Od.x, &luz.Od.y, &luz.Od.z);

	fscanf(p_arq, "%f", &luz.ks);

	fscanf(p_arq, "%f %f %f", &luz.Il.r, &luz.Il.g, &luz.Il.b);

	fscanf(p_arq, "%f", &luz.n);

	fclose(p_arq);
}

static void ler_objeto(char * path)
{
	FILE * p_arq = fopen(path, "r");
	fscanf(p_arq, "%d %d", &num_pontos, &num_triangulos);
	
	pontos_objeto_mundo = (Ponto3D *)malloc((num_pontos + 1) * sizeof(Ponto3D)); // indice 1

	pontos_objeto_vista = (Ponto3D *) malloc((num_pontos + 1) * sizeof(Ponto3D)); // indice 1
	
	normais_vertices = (Ponto3D *) malloc((num_pontos + 1) * sizeof(Ponto3D)); // indice 1
	for (int i = 1; i <= num_pontos; i++) normais_vertices[i].x = normais_vertices[i].y = normais_vertices[i].z = 0.0;
	
	triangulos = (Triangulo *) malloc((num_triangulos + 1)*sizeof(Triangulo)); // indice 1
	
	float x, y, z;
	for (int i = 1; i <= num_pontos; i++)
	{
		fscanf(p_arq, "%f %f %f", &x, &y, &z);
		
		pontos_objeto_mundo[i].x = x;
		pontos_objeto_mundo[i].y = y;
		pontos_objeto_mundo[i].z = z;
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

static void rotacao()
{
	// TODO
}

static void mudanca_base()
{
	// TODO
	// Luz de coordenada mundial para coordenada de vista

	// Objeto para coordenada de vista
}

static void calcula_normais()
{
	Ponto3D* normal;
	int v1, v2, v3;
	for (int i = 1; i <= num_triangulos; i++)
	{
		v1 = triangulos[i].v1;
		v2 = triangulos[i].v2;
		v3 = triangulos[i].v3;

		// ja vem normalizado.
		normal = normal_triangulo(&pontos_objeto_vista[v1], &pontos_objeto_vista[v2], &pontos_objeto_vista[v3]);

		normais_vertices[v1].x += normal->x;
		normais_vertices[v2].x += normal->x;
		normais_vertices[v3].x += normal->x;

		normais_vertices[v1].y += normal->y;
		normais_vertices[v2].y += normal->y;
		normais_vertices[v3].y += normal->y;

		normais_vertices[v1].z += normal->z;
		normais_vertices[v2].z += normal->z;
		normais_vertices[v3].z += normal->z;
	}

	// normaliza as normais dos vertices
	for (int i = 1; i <= num_pontos; i++) normalizar(&normais_vertices[i]);
}

// Função callback chamada para fazer o desenho
void Desenha()
{
/*glMatrixMode(GL_MODELVIEW);
//definir que todas as tranformações vão ser em cena (no desenho)
     glLoadIdentity();*/
                   
     glFlush();
	//não sei exatamente o que faz, umas das coisas que não funciona sem.
}

// Inicializa parâmetros de rendering
void Inicializa (void)
{   
    // Define a cor de fundo da janela de visualização como preta
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	//para ver os parametros da função (e de qualquer outra) usar ctrl+shift+spacebar
	//dentro dos parênteses 
}

// Programa Principal 
int main()
{
	/*ler_camera("Cameras/01_Camera.cfg");
	inicializar_camera();
	ler_luz("iluminacao.txt");
	ler_objeto("Objetos/01_Objeto.byu");*/
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	//setar modo de exibição, nesse caso buffer duplo e modelo de cor RGB
	glutInitWindowSize(WIDTH,HEIGHT);
	//tamanho da janela
	glutInitWindowPosition(10,10);
	//onde a janela vai aparecer na tela do PC (acho isso inutil)
	glutCreateWindow("janelinha");
	//não lembro como criar 2 janelas, vejam ai isso se quiserem
	//criar um janela
	glutDisplayFunc(Desenha);
	//callback da função que dezenha na tela
	Inicializa();
	//inicializar alguns parametros do glut (nessa caso a cor do fundo da tela).
	//cor que vai limpar o buffer
	glutMainLoop();
	//começa a execução da maquina de estados do glut/opengl que controla as funções
	//de callback (controlador de mouse, teclado, callback de controle de tela, etc).
}
