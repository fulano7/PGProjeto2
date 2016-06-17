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

/////////////////////////////////////////////////////////////////////////////
/*
http://www.inf.pucrs.br/~manssour/OpenGL/Tutorial.html
esse link!!!!
n�oo vou parar aquui pra mostrar mais coisas de C/C++ hj..
*/
/////////////////////////////////////////////////////////////////////////////



// Fun��o callback chamada para fazer o desenho
void Desenha(void)
{
glMatrixMode(GL_MODELVIEW);
//definir que todas as tranforma��es v�o ser em cena (no desenho)
     glLoadIdentity();
                   
     // Limpa a janela de visualiza��o com a cor de fundo especificada
     glClear(GL_COLOR_BUFFER_BIT);

     // Especifica que a cor corrente � vermelha
     //         R     G     B
     glColor3f(1.0f, 0.0f, 0.0f);
	 //cor que vai usar para pintar

     // Desenha um quadrado preenchido com a cor corrente
     glBegin(GL_QUADS);
	 //inicializar um denzenho, nesse caso um quadrado preenchido, existem muitas outras
	 //primitivas co opengl (GL_TRIANGLES,GL_LINE,GL_LINE_STIPPLE,GL_POLYGON,etc)
               glVertex2i(100,150);
               glVertex2i(100,100);
               // Especifica que a cor corrente � azul
               glColor3f(0.0f, 0.0f, 1.0f);
               //cor que vai usar para pintar (mudar a cor no meio do desenho
			   //vai ficar legal
			   glVertex2i(150,100);
               glVertex2i(150,150);               
     glEnd();

     glColor3f(1, 1, 1);
	 glPointSize(15);
	 glBegin(GL_POINTS);
               glVertex2i(10,10);
               glVertex2i(200,200);
			   glColor3f(1, 0, 1);
               glVertex2i(200,200);
	 glEnd();
	 //termina um desenho
     
	 // Executa os comandos OpenGL
     glFlush();
	//n�o sei exatamente o que faz, umas das coisas que n�o funciona sem.
}

// Inicializa par�metros de rendering
void Inicializa (void)
{   
    // Define a cor de fundo da janela de visualiza��o como preta
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	//para ver os parametros da fun��o (e de qualquer outra) usar ctrl+shift+spacebar
	//dentro dos par�nteses 
}

void AlteraTamanhoJanela(GLsizei w, GLsizei h)
{
                   // Evita a divisao por zero
                   if(h == 0) h = 1;
                           
                   // Especifica as dimens�es da Viewport
                   glViewport(0, 0, w, h);
				   //define a �rea dentro janela, em coordenadas de tela,
				   //que OpenGL pode usar para fazer o desenho. O volume de
				   //visualiza��o �, ent�o, mapeado para a nova viewport
                   
				   // Inicializa o sistema de coordenadas
                   glMatrixMode(GL_PROJECTION);
				   //definir que todas as tranforma��es v�o ser em camera (usuario)
                   glLoadIdentity();

                   // Estabelece a janela de sele��o (left, right, bottom, top)
                   if (w <= h) 
                           gluOrtho2D (0.0f, 250.0f, 0.0f, 250.0f*h/w);
                   else 
                           gluOrtho2D (0.0f, 250.0f*w/h, 0.0f, 250.0f);
				   //� usada para determinar que a proje��o ortogr�fica (2D) 
				   //ser� utilizada para exibir na tela a imagem
}

// Programa Principal 
int main(void)
{
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	//setar modo de exibi��o, nesse caso um unico buffer e modelo de cor RGB
	////recomendo usar buffer duplo (GLUT_DOUBLE).
	glutInitWindowSize(640,480);
	//tamanho da janela
	glutInitWindowPosition(10,10);
	//onde a janela vai aparecer na tela do PC (acho isso inutil)
	glutCreateWindow("janelinha");
	//n�o lembro como criar 2 janelas, vejam ai isso se quiserem
	//criar um janela
	glutDisplayFunc(Desenha);
	//callback da fun��o que dezenha na tela
	glutReshapeFunc(AlteraTamanhoJanela);
	//callback que controla o tamanho da janela (maximizada e minimizada com e 
	// sem essa fun��o e vc vai entender)
	Inicializa();
	//inicializar alguns parametros do glut (nessa caso a cor do fundo da tela).
	//cor que vai limpar o buffer
	glutMainLoop();
	//come�a a execu��o da maquina de estados do glut/opengl que controla as fun��es
	//de callback (controlador de mouse, teclado, callback de controle de tela, etc).
}