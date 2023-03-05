#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <GL/glut.h>
#include <iostream>
#define VECTOR_LIMIT 6
#include <math.h>
#define K 100
#define VECTOR_LENGTH_LIMIT 5
#define POS_SCALE 20
#define POS_OFFSET_X 0
#define POS_OFFSET_Z 0
#include <GL/glut.h>
#define MSPF 16
#define TIME_CONSTANT 0.005
#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286

struct Vector {
	float x, y, z;
};

float dotProduct(Vector a, Vector b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

float magnitude(Vector a) {
	return sqrt(dotProduct(a, a));
}

int magnitude(int a) {
	return a > 0 ? 1 : -1;
}

Vector crossProduct(Vector a, Vector b) {
	return Vector{ a.y * b.z - b.y * a.z, a.z * b.x - a.x * b.z, a.x * b.y - b.x * a.y };
}

class Charge;

Charge* scene[5];

class Charge {
public:
	float x, y, z; // pozisyon
	float vx = 0, vy = 0, vz = 0; // hız
	float ax = 0, ay = 0, az = 0; // ivme vektörü
	int coulomb;
	Charge(float x, float y, float z, int coulomb) {
		std::cout << vx << vy << vz << std::endl;
		this->x = x;
		this->y = y;
		this->z = z;
		this->coulomb = coulomb;
	}
	void accelerate(float time) { // ivme fonksiyonu
		int a = magnitude(this->coulomb);
		this->vx += this->ax * time * a;
		this->vy += this->ay * time * a;
		this->vz += this->az * time * a;
	}
	void translate(float time) {
		this->x += this->vx * time;
		this->y += this->vy * time;
		this->z += this->vz * time;
	}
	void computeAcceleration() {
		Vector temp;
		this->ax = 0, this->ay = 0, this->az = 0;
		for (int l = 0; l < 5; l++) {
			if (scene[l] == NULL) break;
			if (scene[l] == this) continue;
			temp = scene[l]->fieldAt(this->x, this->y, this->z);
			this->ax += temp.x, this->ay += temp.y, this->az += temp.z;
		}
	}
	void advance(float time) {
		translate(time);
		accelerate(time);
	}
	Vector fieldAt(float x, float y, float z) { // alan hesabı
		Vector r{ x - this->x, y - this->y, z - this->z }; // r vektörü
		float distanceSquared = dotProduct(r, r);
		float factor = K * this->coulomb / distanceSquared;
		r.x *= factor, r.y *= factor, r.z *= factor;
		return r;
	}
};

void drawline(int from_x, int from_y, int from_z, int to_x, int to_y, int to_z)
{
	// From coordinate position
	glVertex3i(from_x, from_y, from_z);

	// To coordinate position
	glVertex3i(to_x, to_y, to_z);
}
void drawArrow(int x1, int y1, int z1, int x2, int y2, int z2, Vector c)
{

	glLineWidth(4.0); // Set line width to 4.0
	glBegin(GL_LINES);
	glColor3f(c.x, c.y, c.z); //çizginin renk ayarları
	drawline(x1, y1, z1, x2, y2, z2);
	glEnd();
	Vector a{ x1,y1,z1 };
	Vector b{ x2, y2,z2 };
	Vector v{ b.x - a.x, b.y - a.y, b.z - a.z };

	double baseScale = magnitude(v) * -5;
	Vector base{ v.x / baseScale, v.y / baseScale, v.z / baseScale };

	Vector basePoint{ base.x + b.x, base.y + b.y, base.z + b.z };

	double theta = 0;
	Vector rotationVector = crossProduct(v, { 0,0, +1 });
	double crossMag = magnitude(rotationVector);
	glTranslated(basePoint.x, basePoint.y, basePoint.z);
	if (crossMag != 0) {
		theta = acos(dotProduct(v, { 0, 0, 1 }) / magnitude(v)) * 180 / PI;
		glRotated(-theta, rotationVector.x, rotationVector.y, rotationVector.z);
	}
	glutSolidCone(1, 2, 20, 10);
	if (crossMag != 0)
		glRotated(+theta, rotationVector.x, rotationVector.y, rotationVector.z);
	glTranslated(-basePoint.x, -basePoint.y, -basePoint.z);
	glLineWidth(2.0); // set line width back to 2.0
}
void drawSphere(int x, int y, int z, int size, Vector c)
{
	if (size < 0) size = -size;
	glTranslatef(x, y, z); // küre yeri belirleme
	//küreyi çizme
	glColor4f(c.x, c.y, c.z, 0.05);
	glutSolidSphere(size, 200, 200);
	glTranslatef(-x, -y, -z);
}


void renderFrame() {
	for (int i = 0; i < VECTOR_LIMIT; i++) {
		for (int j = -VECTOR_LIMIT / 2; j < VECTOR_LIMIT / 2; j++) {
			for (int k = 0; k < VECTOR_LIMIT; k++) {
				Vector temp, total = Vector{ 0,0,0 };
				for (int l = 0; l < 5; l++) {
					if (scene[l] == NULL) break;
					temp = scene[l]->fieldAt(i * POS_SCALE, j * POS_SCALE, k * POS_SCALE);
					total.x += temp.x, total.y += temp.y, total.z += temp.z;
				}
				float mag = magnitude(total);
				Vector color = Vector{ 1, 1, 0 };
				if (mag > VECTOR_LENGTH_LIMIT) {
					color.x = VECTOR_LENGTH_LIMIT / mag;
				}
				total = { total.x * VECTOR_LENGTH_LIMIT / mag, total.y * VECTOR_LENGTH_LIMIT / mag, total.z * VECTOR_LENGTH_LIMIT / mag };
				drawArrow(i * POS_SCALE, j * POS_SCALE, k * POS_SCALE, total.x + i * POS_SCALE, total.y + j * POS_SCALE, total.z + k * POS_SCALE, color);
			}
		}
	}
	for (int i = 0; i < 5; i++) {
		if (scene[i] == NULL) break;
		scene[i]->advance(TIME_CONSTANT);
	}
	for (int i = 0; i < 5; i++) {
		if (scene[i] == NULL) break;
		scene[i]->computeAcceleration();
		drawSphere(scene[i]->x, scene[i]->y, scene[i]->z, scene[i]->coulomb, scene[i]->coulomb > 0 ? Vector{ 0, 0, 1 } : Vector{ 1, 0, 0 });
	}
}

using namespace std;

#define WHEEL_UP 3
#define WHEEL_DOWN 4
#define MAX_ROWS 100
#define MAX_COLS 100
#define CAMERA_DISTANCE_MIN 1.0
#define CAMERA_DISTANCE_MAX 500.0
#define SBSIZE 256
#define M_PI 3.14159265358979323846264338327950288

class wcPt3d {
public:
	GLfloat x, y, z;
};

enum {
	DL_BOX = 1,
	DL_GRID,
};

GLsizei pencereEni = 700, pencereBoyu = 500;
static GLint MouseX = 0;
static GLint MouseY = 0;
static GLint MouseButton;
static GLint MouseState;
static double CameraLatitude = 45.0;
static double CameraLongitude = 25.0;
static double CameraDistance = 50.0;
static bool CameraOrthographic = false;
static double EyeX = -100.0;
static double EyeY = 50.0;
static double EyeZ = 100.0;
static double DirX = 0.0;
static double DirY = 0.0;
static double DirZ = 0.0;
static GLuint SelectBuffer[SBSIZE];
static int FrameRate = 60;
static wcPt3d CameraFocus = { MAX_ROWS / 2, 0.0, MAX_COLS / 2 };
static wcPt3d SelectedBox = { MAX_ROWS / 2, 0.0, MAX_COLS / 2 };
static bool ShowGrid = true;
static vector<wcPt3d> Controls;
static GLint ControlIndex = -1;

// Update camera
void kameraKonumunuGuncelle() {
	double L = CameraDistance * cos(M_PI * CameraLongitude / 180.0);
	double X = MAX_ROWS / 2;
	double Z = MAX_COLS / 2;
	EyeX = X + L * -sin(M_PI * CameraLatitude / 180.0);
	EyeY = CameraDistance * sin(M_PI * CameraLongitude / 180.0);
	EyeZ = Z + L * cos(M_PI * CameraLatitude / 180.0);
	DirX = X;
	DirY = 0.0;
	DirZ = Z;
}
// Initialize scene
void init_scene() {
	glEnable(GL_DEPTH_TEST);
	GLfloat diffuse0[] = { 0.8, 0.8, 0.8, 1.0 };
	GLfloat ambient0[] = { 0.2, 0.2, 0.2, 1.0 };
	GLfloat specular0[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat position0[] = { 100.0, 100.0, 100.0, 1.0 };
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular0);
	glLightfv(GL_LIGHT0, GL_POSITION, position0);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_NORMALIZE);
	glEnable(GL_SMOOTH);
	glShadeModel(GL_SMOOTH);
	kameraKonumunuGuncelle();
}
// To set the projection
void projeksiyonuAyarla(bool reset = true) {
	glMatrixMode(GL_PROJECTION);
	if (reset)
		glLoadIdentity();
	if (CameraOrthographic)
		glOrtho(-MAX_ROWS, MAX_ROWS, -MAX_ROWS, MAX_ROWS, -MAX_COLS,
			MAX_COLS);
	else
		gluPerspective(45.0, (GLdouble)(pencereEni) / (GLdouble)(pencereBoyu), 0.1, 1000.0);
}

void pencereYenidenSekillendir(int yeniEn, int yeniBoy) {
	glViewport(0, 0, yeniEn, yeniBoy);
	glMatrixMode(GL_PROJECTION); // projeksiyon parametrelerini ayarla
	glLoadIdentity();
	gluOrtho2D(0, (GLdouble)yeniEn, 0, (GLdouble)yeniBoy);
	// Mevcut pencere eni ve boyunu g�ncelle
	pencereEni = yeniEn;
	pencereBoyu = yeniBoy;
	projeksiyonuAyarla();
}

// Box �izmek i�in 
void boxCizmek() {
	glBegin(GL_QUADS); {
		glVertex3f(-0.5, 0.0, -0.5);
		glVertex3f(0.5, 0.0, -0.5);
		glVertex3f(0.5, 0.0, 0.5);
		glVertex3f(-0.5, 0.0, 0.5);
	} glEnd();
}

// Grid �izmek i�in 
void gridCizmek() {
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glDisable(GL_LIGHTING);
	for (int r = 0; r < MAX_ROWS; r++) {
		glPushName(r);
		for (int c = 0; c < MAX_COLS; c++) {
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glPushName(c);
			glPushMatrix(); {
				glTranslatef(r, 0, c);
				glCallList(DL_BOX);
			} glPopMatrix();
			glPopName();
		}
		glPopName();
	}
	glPopAttrib();
}

// G�r�nt�leme listeleri yapmak i�in
void make_display_lists() {
	glNewList(DL_BOX, GL_COMPILE); {
		boxCizmek();
	}; glEndList();
	glNewList(DL_GRID, GL_COMPILE); {
		gridCizmek();
	}; glEndList();
}

// Draw scene 
void draw_scene() {
	if (CameraOrthographic) {
		double d = CAMERA_DISTANCE_MAX - CameraDistance;
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glRotatef(90, 1.0, 0.0, 0.0);
		glScalef(d / 10.0, d / 10.0, d / 10.0);
		glTranslatef(-CameraFocus.x, 0.0, -CameraFocus.z);
	}
	else {
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(EyeX, EyeY, EyeZ, DirX, DirY, DirZ, 0.0, 1.0, 0.0);
		renderFrame();
	}
	if (ShowGrid) {
		glColor3d(0.3, 0.3, 0.3);
		glCallList(DL_GRID);
	}
	glFlush();
}

void timer(int val) {
	glutPostRedisplay();
	glutTimerFunc(MSPF, timer, 0);
}

// display function
void display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	projeksiyonuAyarla();
	draw_scene();
	glutSwapBuffers();
}

// Keyboard fonksiyonu 
void keyboard(GLubyte egriCizmeTusu, GLint xFare, GLint yFare) {
	GLint x = xFare;
	GLint y = pencereBoyu - yFare;
	switch (egriCizmeTusu) {
	case 'e':
	case 'E':
		exit(EXIT_SUCCESS);
		break;
	case 'p':
	case 'P':
		CameraOrthographic = !CameraOrthographic;
		CameraFocus.x = SelectedBox.x;
		CameraFocus.z = SelectedBox.z;
		if (!CameraOrthographic)
			ControlIndex = -1;
		break;
	case 'g':
	case 'G':
		ShowGrid = !ShowGrid;
		break;
	}
	kameraKonumunuGuncelle();
}

// Mouse fonksiyonu 
void mouse(int button, int state, int x, int y) {
	MouseX = x;
	MouseY = y;
	MouseButton = button;
	MouseState = state;
	switch (button) {
	case WHEEL_UP:
		if (ControlIndex >= 0 && !CameraOrthographic)
			Controls[ControlIndex].y += 0.1;
		else
			CameraDistance = (CameraDistance > CAMERA_DISTANCE_MIN ? CameraDistance -
				1.0 : CAMERA_DISTANCE_MIN);
		break;
	case WHEEL_DOWN:
		if (ControlIndex >= 0 && !CameraOrthographic)
			Controls[ControlIndex].y -= 0.1;
		else
			CameraDistance = (CameraDistance < CAMERA_DISTANCE_MAX ? CameraDistance
				+ 1.0 : CAMERA_DISTANCE_MAX);
		break;
	}
	kameraKonumunuGuncelle();
}

// Motion callback 
void motion(int x, int y) {
	double dx = (double)(x - MouseX) / (double)(pencereEni);
	double dy = (double)(y - MouseY) / (double)(pencereBoyu);
	if (MouseButton == GLUT_RIGHT_BUTTON) {
		CameraLatitude += 180.0 * dx;
		CameraLongitude += 180.0 * dy;
		if (CameraLongitude < -90.0)
			CameraLongitude = -90.0 + 0.0001;
		if (CameraLongitude > 90.0)
			CameraLongitude = 90.0 - 0.0001;
		kameraKonumunuGuncelle();
	}
	MouseX = x;
	MouseY = y;
}

// Main execution 
int main(int argc, char* argv[]) {
	scene[0] = new Charge(2.5 * POS_SCALE, -2.5 * POS_SCALE, 2.5 * POS_SCALE, 3);
	//scene[0]->vz = 5;
	scene[1] = new Charge(2.5 * POS_SCALE, +2.5 * POS_SCALE, 2.5 * POS_SCALE, -3);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(pencereEni, pencereBoyu);
	glutCreateWindow("proje");
	glutDisplayFunc(display);
	glutReshapeFunc(pencereYenidenSekillendir);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	init_scene();
	make_display_lists();
	timer(0);
	glutMainLoop();
	return (EXIT_SUCCESS);
}