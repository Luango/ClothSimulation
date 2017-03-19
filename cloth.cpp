//
//  main.cpp
//  MyCloth
//
//  Created by Lu Ang on 15/3/15.
//  Copyright (c) 2015 Lu Ang. All rights reserved.
//

#include <iostream>
#include <GLUT/GLUT.h>
#include <math.h>
#include <vector>
using namespace std;

#define TIME_STEP 0.25;
#define DAMP 0.9;
class Vec3;

class Mat3
{
public:
    float f[3][3];
    
    Mat3(float e00, float e01, float e02, float e10, float e11, float e12, float e20, float e21, float e22)
    {
        f[0][0] = e00;
        f[0][1] = e01;
        f[0][2] = e02;
        f[1][0] = e10;
        f[1][1] = e11;
        f[1][2] = e12;
        f[2][0] = e20;
        f[2][1] = e21;
        f[2][2] = e22;
    }
    Mat3() {}
    Mat3 operator/ (const float &a)
    {
        return Mat3(f[0][0]/a,f[0][1]/a,f[0][2]/a,
                    f[1][0]/a,f[1][1]/a,f[1][2]/a,
                    f[2][0]/a,f[2][1]/a,f[2][2]/a);
    }
    
    Mat3 operator- (const Mat3 &v)
    {
        return Mat3(f[0][0]-v.f[0][0],f[0][1]-v.f[0][1],f[0][2]-v.f[0][2],
                    f[1][0]-v.f[1][0],f[1][1]-v.f[1][1],f[1][2]-v.f[1][2],
                    f[2][0]-v.f[2][0],f[2][1]-v.f[2][1],f[2][2]-v.f[2][2]);
    }
    
    Mat3 operator+ (const Mat3 &v)
    {
        return Mat3(f[0][0]+v.f[0][0],f[0][1]+v.f[0][1],f[0][2]+v.f[0][2],
                    f[1][0]+v.f[1][0],f[1][1]+v.f[1][1],f[1][2]+v.f[1][2],
                    f[2][0]+v.f[2][0],f[2][1]+v.f[2][1],f[2][2]+v.f[2][2]);
    }
    
    Mat3 operator* (const float &a)
    {
        return Mat3(f[0][0]*a,f[0][1]*a,f[0][2]*a,
                    f[1][0]*a,f[1][1]*a,f[1][2]*a,
                    f[2][0]*a,f[2][1]*a,f[2][2]*a);
    }
    friend Mat3 operator *(Vec3 &a);
    
    void show_mat()
    {
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                cout<<f[i][j]<<" ";
            }
            cout<<endl;
        }
    }
};

class Vec3
{
public:
    float f[3];
    
    Vec3(float x, float y, float z) { f[0] =x; f[1] =y; f[2] =z; }
    Vec3() {}
    float length() { return sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]); }
    Vec3 normalized(){ float l = length(); return Vec3(f[0]/l,f[1]/l,f[2]/l); }
    void operator+= (const Vec3 &v) {
        f[0]+=v.f[0];
        f[1]+=v.f[1];
        f[2]+=v.f[2];
    }
    Vec3 operator/ (const float &a){ return Vec3(f[0]/a,f[1]/a,f[2]/a); }
    Vec3 operator- (const Vec3 &v){ return Vec3(f[0]-v.f[0],f[1]-v.f[1],f[2]-v.f[2]); }
    Vec3 operator+ (const Vec3 &v){ return Vec3(f[0]+v.f[0],f[1]+v.f[1],f[2]+v.f[2]); }
    Vec3 operator* (const float &a){ return Vec3(f[0]*a,f[1]*a,f[2]*a); }
    Vec3 operator-(){ return Vec3(-f[0],-f[1],-f[2]); }
    float dot(const Vec3 &v) { return f[0]*v.f[0] + f[1]*v.f[1] + f[2]*v.f[2]; }
    
    Mat3 outerProduct( const Vec3 &v )
    {
        return Mat3(f[0]*v.f[0],f[0]*v.f[1],f[0]*v.f[2],
                    f[1]*v.f[0],f[1]*v.f[1],f[1]*v.f[2],
                    f[2]*v.f[0],f[2]*v.f[1],f[2]*v.f[2]);
    }
    Vec3 multiply_Mat3(const Mat3 &v)
    {
        return Vec3(v.f[0][0]*f[0]+v.f[0][1]*f[1]+v.f[0][2]*f[2],
                    v.f[1][0]*f[0]+v.f[1][1]*f[1]+v.f[1][2]*f[2],
                    v.f[2][0]*f[0]+v.f[2][1]*f[1]+v.f[2][2]*f[2]);
    }
    void show()
    {
        cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<endl<<endl;
    }
};

class Particle
{
private:
    bool movable;   // determine if particle can move
    float mass;    // the mass of each particle
    Vec3 pos;   // position
    Vec3 old_pos;   // previous position
    
public:
    Vec3 acceleration;
    Vec3 force;
    Vec3 velocity;
    
    Particle(Vec3 pos): pos(pos),old_pos(pos),mass(1),movable(true),acceleration(Vec3 (0,0,0)),velocity(Vec3 (0,0,0)),force(Vec3 (0,0,0)) {}
    Particle(){}
    
    void addForce(Vec3 force)
    {
        acceleration += force/mass;
    }
    void timeStep()
    {
        if (movable) {
            Vec3 temp = pos;
            pos += acceleration * TIME_STEP;
            old_pos = temp;
            acceleration = Vec3(0,0,0);
        }
    }
    Vec3& getPos()  // return position
    {
        return pos;
    }
    void makeUnmovable()
    {
        movable = false;
    }
    void moveParticle(Vec3 vec)
    {
        if(movable)
        {
            pos += vec*(1);
        }
    }
};

class Constraint
{
private:
    float rest_distance;  // the distance between particle p1 and p2
public:
    Particle *p1, *p2;
    Constraint(Particle *p1, Particle *p2): p1(p1),p2(p2)
    {
        Vec3 vec = p1->getPos() - p2->getPos();
        rest_distance = vec.length();
    }
    // verlet method to satisfy constraint
    void verletSatisfyConstraint()
    {
        Vec3 p1_to_p2 = p2->getPos() - p1->getPos();
        float current_distance = p1_to_p2.length();
        // stretch vector
        if(rest_distance<current_distance){
            Vec3 stretch_vector = p1_to_p2 * (1 - rest_distance/current_distance);
            p1->moveParticle(stretch_vector * 0.5);
            p2->moveParticle(-stretch_vector * 0.5);
        }
    }
    // explicit method to satisfy constraint
    void explicitSatisfyConstraint()
    {
        Vec3 p1_to_p2 = p2->getPos() - p1->getPos();
        float current_distance = p1_to_p2.length();
        Vec3 force = p1_to_p2 * (1 - rest_distance/current_distance);
        p1->velocity = p1->velocity + force * TIME_STEP;
        p2->velocity = p2->velocity - force * TIME_STEP;
        p1->velocity = p1->velocity * DAMP;
        p2->velocity = p2->velocity * DAMP;
        
        p1->moveParticle(p1->velocity);
        p2->moveParticle(p2->velocity);
    }
    // implicit method to satisfy stretch constraint
    void implicitSatisfyStretchConstraint()
    {
        Vec3 p1_to_p2 = p2->getPos() - p1->getPos();
        float current_distance = p1_to_p2.length();
        Vec3 force = p1_to_p2 * (1 - rest_distance/current_distance);
        Mat3 multiply_result_matrix;
        Mat3 identity_matrix(1, 0, 0, 0, 1, 0, 0, 0, 1);
        Mat3 div_result_matrix;
        Mat3 result_matrix;
        float dot_result;
        Vec3 b;
        Vec3 delta_x;
        if(rest_distance<current_distance){
            // First calculate A matrix
            dot_result = p1_to_p2.dot(p1_to_p2);
            multiply_result_matrix = p1_to_p2.outerProduct(p1_to_p2);
            multiply_result_matrix = multiply_result_matrix/dot_result;
            div_result_matrix = multiply_result_matrix+(identity_matrix-multiply_result_matrix)*(1-(rest_distance/current_distance));
            div_result_matrix = div_result_matrix*TIME_STEP;
            result_matrix = identity_matrix - div_result_matrix;
            // Second calculate b vector
            p1->velocity = p1->velocity + force * TIME_STEP;
            p2->velocity = p2->velocity - force * TIME_STEP;
            b = (p1->velocity+force)*TIME_STEP;
            
            conjugate_gradient(&delta_x, result_matrix, b, Vec3(0,0,0));
            p1->moveParticle(delta_x);
            p2->moveParticle(-delta_x);
            
            p1->velocity = p1->velocity * 0.9;
            p2->velocity = p2->velocity * 0.9;
        }
    }
    // explicit method to satisfy bending constraint
    void explicitSatisfyBendingConstraint()
    {
        Vec3 p1_to_p2 = p2->getPos() - p1->getPos();
        float current_distance = p1_to_p2.length();
        float k = (2*(current_distance/rest_distance)/(sin(current_distance/rest_distance)))/(current_distance/rest_distance);
        float fb = (1.0*k*k)/(cos(k*rest_distance/2)-(k*rest_distance/2)/sin(k*rest_distance/2));
        Vec3 force = (p1_to_p2/current_distance)*fb*0.001;
        Vec3 force2 = (p1_to_p2/current_distance)*0.01;
        if(force.length()<force2.length())
        {
            force = force2;
        }
        if (rest_distance>current_distance)
        {
            p1->velocity = p1->velocity*0;
            p2->velocity = p2->velocity*0;
            p1->velocity = p1->velocity + force * TIME_STEP;
            p2->velocity = p2->velocity - force * TIME_STEP;
            p1->velocity = p1->velocity * DAMP;
            p2->velocity = p2->velocity * DAMP;
            
            p1->moveParticle(p1->velocity);
            p2->moveParticle(p2->velocity);
        }
    }
    
    // conjugate gradient function algorithm
    // Start with some x(0).
    // Set p(0) = r(0) = b - Ax(0)
    // For k = 0, 1, 2, ...
    // x(k+1) = x(k) + a(k)p(k),
    // a(k) = r(k)(T)r(k)/p(k)(T)Ap(k)
    // r(k+1) = b - Ax(k+1) = r(k) - a(k)Ap(k)
    // p(k+1) = r(k+1) + beta(k)p(k),
    // beta(k) = r(k+1)(T)r(k+1)/r(k)(T)r(k)
    void conjugate_gradient(Vec3* delta_x, Mat3 A, Vec3 b, Vec3 x0)
    {
        Vec3 x_next;
        Vec3 x_previous;
        float a;
        Vec3 r_next;
        Vec3 r_previous;
        Vec3 p_next;
        Vec3 p_previous;
        float beta;
        
        int k=0;
        p_previous = b - x0.multiply_Mat3(A);
        r_previous = p_previous;
        x_previous = x0;
        while(k<10)
        {
            a = (r_previous.dot(r_previous))/(p_previous.multiply_Mat3(A).dot(p_previous));
            x_next = x_previous+p_previous*a;
            r_next = p_previous.multiply_Mat3(A);
            r_next = r_previous - r_next * a;
            
            beta = (r_next.dot(r_next))/(r_previous.dot(r_previous));
            p_next = r_next + p_previous * beta;
            x_previous = x_next;
            r_previous = r_next;
            p_previous = p_next;
            
            k++;
            if(r_next.length()<1.0)
            {
                *delta_x = x_next;
                break;
            }
        }
    }
};

class Cloth
{
private:
    int cloth_width_num;
    int cloth_height_num;
    vector<Particle> particles;
    vector<Constraint> stretch_constraints;
    vector<Constraint> bending_constraints;
    
    void drawPoints(Particle *p1)
    {
        glVertex3fv((GLfloat*) &(p1->getPos() ));
    }
    Particle* getParticle(int x, int y)
    {
        return &particles[y*cloth_width_num+x];
    }
    void makeStretchConstraint(Particle *p1, Particle *p2)
    {
        stretch_constraints.push_back(Constraint(p1, p2));
    }
    void makeBendingConstraint(Particle *p1, Particle *p2)
    {
        bending_constraints.push_back(Constraint(p1, p2));
    }
public:
    Cloth(float width, float height, int width_num, int height_num) : cloth_width_num(width_num),cloth_height_num(height_num)
    {
        particles.resize(cloth_height_num*cloth_width_num);
        for (int x=0; x<cloth_width_num; x++) {
            for (int y=0; y<cloth_height_num; y++) {
                Vec3 pos = Vec3(width*(x/(float)cloth_width_num), -height * (y/(float)cloth_height_num), 0);
                particles[y*cloth_width_num + x] = Particle(pos);
            }
        }
        for (int x=0; x<cloth_width_num; x++)
        {
            for (int y=0; y<cloth_height_num; y++)
            {
                if (x<cloth_width_num-1){
                    makeStretchConstraint(getParticle(x, y), getParticle(x+1, y));
                }
                if (y<cloth_height_num-1) {
                    makeStretchConstraint(getParticle(x, y), getParticle(x, y+1));
                }
                if (x<cloth_width_num-1 && y<cloth_height_num-1) {
                    makeStretchConstraint(getParticle(x, y), getParticle(x+1, y+1));
                }
                if (x<cloth_width_num-1 && y<cloth_height_num-1) {
                    makeStretchConstraint(getParticle(x+1, y), getParticle(x, y+1));
                }
            }
        }
        for (int x=0; x<cloth_width_num; x++)
        {
            for (int y=0; y<cloth_height_num; y++)
            {
                if (x<cloth_width_num-2){
                    makeBendingConstraint(getParticle(x, y), getParticle(x+2, y));
                }
                if (y<cloth_height_num-2) {
                    makeBendingConstraint(getParticle(x, y), getParticle(x, y+2));
                }
                if (x<cloth_width_num-2 && y<cloth_height_num-2) {
                    makeBendingConstraint(getParticle(x, y), getParticle(x+2, y+2));
                }
                if (x<cloth_width_num-2 && y<cloth_height_num-2) {
                    makeBendingConstraint(getParticle(x+2, y), getParticle(x, y+2));
                }
            }
        }
        for (int i =0; i<4; i++) {
            getParticle(i, 0)->makeUnmovable();
            getParticle(cloth_width_num-i, 0)->makeUnmovable();
        }
    }
    void drawCloth()
    {
        glBegin(GL_LINES);
        for (int x = 0; x<cloth_width_num-1; x++) {
            for (int y = 0; y<cloth_height_num-1; y++) {
                drawPoints(getParticle(x, y));
                drawPoints(getParticle(x+1, y+1));
            }
        }
        glEnd();
    }
    void addForce(Vec3 force)
    {
        vector<Particle>::iterator particle;
        for (particle = particles.begin(); particle != particles.end(); particle++) {
            (*particle).addForce(force);
        }
    }
    // cloth timeStep
    void timeStep()
    {
        vector<Constraint>::iterator stretch_constraint;
        vector<Constraint>::iterator bending_constraint;
        vector<Particle>::iterator particle;
        // the larger i the stronger cloth
        for (int i = 0; i<30; i++)
        {
            for (stretch_constraint = stretch_constraints.begin(); stretch_constraint != stretch_constraints.end(); stretch_constraint++)
            {
                //(*constraint).verletSatisfyConstraint();
                //(*constraint).explicitSatisfyConstraint();
                (*stretch_constraint).implicitSatisfyStretchConstraint();
            }
            for (bending_constraint = bending_constraints.begin(); bending_constraint != bending_constraints.end(); bending_constraint++)
            {
                //(*bending_constraint).verletSatisfyConstraint();
                //(*bending_constraint).explicitSatisfyConstraint();
                (*bending_constraint).explicitSatisfyBendingConstraint();
            }
            for (particle = particles.begin(); particle != particles.end(); particle++)
            {
                (*particle).timeStep();
            }
        }
    }
    //collide with a ball
    void ballCollision(Vec3 center, float radius)
    {
        vector<Particle>::iterator particle_iterator;
        for(particle_iterator = particles.begin(); particle_iterator!=particles.end();particle_iterator++)
        {
            Vec3 v=(*particle_iterator).getPos()-center;
            float distance = v.length();
            if (v.length() < radius) {
                (*particle_iterator).moveParticle(v.normalized()*(radius-distance));
            }
        }
    }
};
Cloth cloth(10,10,30,30);
Vec3 ball_pos(5,-11,0);
float ball_radius = 4;
float cube_length =8;
float ball_time = 0;

void display(void)
{
    ball_time++;
    ball_pos.f[2] = -cos(ball_time/30.0)*10;
    cloth.addForce(Vec3(0, -0.511, 0));
    cloth.timeStep();
    cloth.ballCollision(ball_pos,ball_radius);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    
    glDisable(GL_LIGHTING);
    glBegin(GL_POLYGON);
    glColor3f(1.0f,1.0f,1.0f);
    glVertex3f(-200.0f,-100.0f,-100.0f);
    glVertex3f(200.0f,-100.0f,-100.0f);
    glColor3f(0.8f,0.8f,0.8f);
    glVertex3f(200.0f,100.0f,-100.0f);
    glVertex3f(-200.0f,100.0f,-100.0f);
    glEnd();
    glEnable(GL_LIGHTING);
    
    glTranslatef(-6.5,6,-12.0f);
    glRotated(45, 0, 1, 0);
    cloth.drawCloth();
    
    glPushMatrix();
    glTranslatef(ball_pos.f[0],ball_pos.f[1],ball_pos.f[2]);
    glColor3f(1.0f,1.0f,1.0f);
    glutSolidSphere(ball_radius-0.01,50,50);
    glPopMatrix();
    
    glutSwapBuffers();
    glutPostRedisplay();
}
void init(void)
{
    glShadeModel(GL_SMOOTH);
    glClearColor(0.2f, 0.2f, 0.4f, 0.5f);
    glClearDepth(1.0f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_COLOR_MATERIAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    GLfloat lightPos[4] = {1.0,1.0,0.5,0.0};
    glLightfv(GL_LIGHT0,GL_POSITION,(GLfloat *) &lightPos);
    
    glEnable(GL_LIGHT1);
    
    GLfloat lightAmbient1[4] = {0.0,0.0,0.0,0.0};
    GLfloat lightPos1[4] = {-1.0,0.0,-0.2,0.0};
    GLfloat lightDiffuse1[4] = {0.3,0.3,0.3,0.0};
    
    glLightfv(GL_LIGHT1,GL_POSITION,(GLfloat *) &lightPos1);
    glLightfv(GL_LIGHT1,GL_AMBIENT,(GLfloat *) &lightAmbient1);
    glLightfv(GL_LIGHT1,GL_DIFFUSE,(GLfloat *) &lightDiffuse1);
    
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
}
void reshape(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(80, (float)w/(float)h, 1.0, 5000.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void timer(int p)
{
    glutPostRedisplay();
    glutTimerFunc(10, timer, 0);
}

int main(int argc, char ** argv) {
    cout << "Hello, World!\n";
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB|GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutCreateWindow("Assignmen 3 Cloth");
    init();
    glutTimerFunc(400, timer, 0);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMainLoop();
    
    return 0;
}
