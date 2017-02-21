
#include "loop_subdivision.h"

PFNGLGENQUERIESARBPROC        glGenQueriesARB        = NULL;
PFNGLDELETEQUERIESARBPROC     glDeleteQueriesARB     = NULL;
PFNGLISQUERYARBPROC           glIsQueryARB           = NULL;
PFNGLBEGINQUERYARBPROC        glBeginQueryARB        = NULL;
PFNGLENDQUERYARBPROC          glEndQueryARB          = NULL;
PFNGLGETQUERYIVARBPROC        glGetQueryivARB        = NULL;
PFNGLGETQUERYOBJECTIVARBPROC  glGetQueryObjectivARB  = NULL;
PFNGLGETQUERYOBJECTUIVARBPROC glGetQueryObjectuivARB = NULL;

GLuint BoundingBoxQuery=-1,numberOfFragmentsRendered=0;

void reshape ( int w, int h);
void keyboard( unsigned char key, int x, int y);
void idle(void);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void display();

int last_x, last_y, spin_x, spin_y;
GLboolean pauseFlag = FALSE;

Model *model1, *model2, *model3, *model4, *model5, *model6;
Scene *scene;
int zoomCounter=0;
int STEP_SIZE = 4;
GLboolean first_time=TRUE;


int current_mode = BASIC_MODE;
int initial_subdivision_level  = 0;


int main(int argc, char* argv[])
{
    printf("\n----------------------------------------------------------------------------------------------\n");
    printf("\nThe are three modes supported:\n\t1) Basic Mode : You start with the base-mesh. You can subdivide to any level by pressing L as many times as you want to subdivide.");
    printf("\n\t2) Zoom Mode : You start with a base mesh far away so that it looks reasonably smooth from this distance. Then you can press Z to zoom into the scene.\n\tWhile you zoom in, the model will be automatically subdivided so as to look reasonably smooth for that distance.");
    printf("\n\t You can also press L anytime to explicitly subdivide the model.");
    printf("\n\t3) Level Mode : Here, you directly start with the model subdivided to the specified level, instead of the base mesh.");
    printf("\n\t Again, you can also press L anytime to explicitly subdivide the model.");
    printf("\n----------------------------------------------------------------------------------------------\n");
    printf("\nYou can also ROTATE the model using mouse drags.");
    fflush(stdout);
    


    if ( argc == 3 )
    {
        if ( strcmp(argv[1],"level") == 0 )
        {
            initial_subdivision_level = atoi(argv[2]);
            current_mode = LEVEL_MODE;
            printf("\nProgram started in Level Mode. You would be seeing the model subdivided to level %d. ",initial_subdivision_level ); fflush(stdout);
            printf("\nYou may wish to subdivide further by pressing L.");
        }
        else 
        {
            printf("\nPossible Usages :\n\t <program_name> level <number_of_subdivision_levels>. \n\t <program_name> zoomMode \n\t <program_name> basicMode\n");
            fflush(stdout);
            return 0;
        }

    }
    else if ( argc == 2 )
    {
        if ( strcmp(argv[1],"zoomMode") == 0 )
        {
            current_mode = ZOOM_MODE;
        }
        else if ( strcmp(argv[1],"basicMode") == 0 )
        {
            current_mode = BASIC_MODE;
        }
        else 
        {
            printf("\nPossible Usages :\n\t <program_name> level <number_of_subdivision_levels>. \n\t <program_name> zoomMode \n\t <program_name> basicMode\n");
            fflush(stdout);
            return 0;
        }
    }



    
    
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(512, 512);
    glutInitWindowPosition(50, 50);
    glutInit(&argc, argv);

    glutCreateWindow("OpenGL Viewer for subdivision");


    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutIdleFunc(idle);
    glClearDepth(20.0);
    
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_CULL_FACE);
    glEnable(GL_COLOR_MATERIAL);

    GLfloat light_position[] = {0.0, 0.0, -6.0, 0.0};
    glLightfv(GL_LIGHT0, GL_POSITION,light_position ); 
    GLfloat lightColor[] = {0.4f, 0.2f, 0.0f, 1.0f};
    GLfloat spec_lightColor[] = {1.0f, 1.0f, 1.0f, 1.0f};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);
    glLightfv(GL_LIGHT0, GL_SPECULAR, spec_lightColor);

    
    glShadeModel(GL_SMOOTH);
    model1 = new Model;
    model1->readFromOBJ("obj_data//1.obj");
    
    model2 = new Model;
    model2->readFromOBJ("obj_data//2.obj");

    model3 = new Model;
    model3->readFromOBJ("obj_data//3.obj");

    model4 = new Model;
    model4->readFromOBJ("obj_data//4.obj");

    model5 = new Model;
    model5->readFromOBJ("obj_data//5.obj");

    model6 = new Model;
    model6->readFromOBJ("obj_data//6.obj");

    scene = new Scene();
    scene->addModel(model1);
    scene->addModel(model2);
    scene->addModel(model3);
    scene->addModel(model4);
    scene->addModel(model5);
    scene->addModel(model6);

    //model->printModel();
 
    //model->assertFullNeighbourhoodInfo();
    //getch();


    
    glutMainLoop();
    
    return 0;

}




void printError(char *errorString)
{
    printf("\r\nError : %s",errorString );
    fflush(stdout);
    exit(1);
}

void idle(void)
{

    if ( pauseFlag == FALSE ) glutPostRedisplay();
}

void keyboard( unsigned char key, int x, int y)
{
    switch(key)
    {
    case 'L':
    case 'l':
        scene->subdivide(1);
        glutPostRedisplay();
        break;
    case 'T':
    case't':
        scene->refineTopology();
        glutPostRedisplay();
        break;
    case 'G':
    case'g':
        scene->refineGeometry();
        glutPostRedisplay();
        break;
    case 'w':
    case 'W':
        //rotateVector[0] = 1.0; rotateVector[1] = 0.0; rotateVector[2] = 0.0;
        glutPostRedisplay();
        break;
    case 'p':
        pauseFlag = ~(pauseFlag);
        glutPostRedisplay();
        break;
    case 'Z':
    case 'z':
        if ( current_mode == ZOOM_MODE )
        {
            zoomCounter++;
            glTranslatef(-1,0.0,0.0);
        }
        glutPostRedisplay();
        break;
    case 'q':
    case 'Q':
        exit(0);
        break;

    

    default: 
        printf("\n----------------------------------------------------------------------------------------------\n");
        printf("\nThe are three modes supported:\n\t1) Basic Mode : You start with the base-mesh. You can subdivide to any level by pressing L as many times as you want to subdivide.");
        printf("\n\t2) Zoom Mode : You start with a base mesh far away so that it looks reasonably smooth from this distance. Then you can press Z to zoom into the scene.\n\tWhile you zoom in, the model will be automatically subdivided so as to look reasonably smooth for that distance.");
        printf("\n\t You can also press L anytime to explicitly subdivide the model.");
        printf("\n\t3) Level Mode : Here, you directly start with the model subdivided to the specified level, instead of the base mesh.");
        printf("\n\t Again, you can also press L anytime to explicitly subdivide the model.");
        printf("\n----------------------------------------------------------------------------------------------\n");
        printf("\nYou can also ROTATE the model using mouse drags.");
        fflush(stdout);
        break;
    }
}

void mouse(int button, int state, int x, int y)
{
    if (state & GLUT_UP)
    {
        if (x != last_x || y != last_y)
        {
            //animate = 1;
        }
    }
    else
    {
        //animate = 0;
        last_x = x;
        last_y = y;
    }
    glutPostRedisplay();
}

void motion(int x, int y)
{
    spin_x += x-last_x;
    spin_y += y-last_y;
    last_x = x;
    last_y = y;
    glutPostRedisplay();
}

void reshape ( int w, int h)
{
    glViewport(0,0,w,h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40.0,1, 1.0, 400.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0.0, 0.0, 6.0, 0.0,0.0,0.0,0.0,1.0,0.0);
    glRotatef(90,0.0,1.0,0.0);
    glRotatef(90,1.0,0.0,0.0);
    
    if ( current_mode == ZOOM_MODE ) glTranslatef(20.0,0.0,0.0);

}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if ( current_mode == ZOOM_MODE  && zoomCounter != 0 && (first_time || zoomCounter%STEP_SIZE == 0 )) 
    {
        scene->subdivide(1);
        zoomCounter = 0;
        first_time = FALSE;
    }



    if ( current_mode == LEVEL_MODE )
    {
        glPushMatrix();
        glTranslatef(0.0, 1.0,0.0);
        scene->drawScene();
        glutSwapBuffers();
        glPopMatrix();
        scene->subdivide(initial_subdivision_level);
        current_mode = BASIC_MODE;
    }

    glPushMatrix();
    glRotatef((float)-spin_x,0.0,0.0,1.0);
    glRotatef((float)spin_y,0.0,1.0,0.0);
    scene->drawScene();
    glPopMatrix();



    glutSwapBuffers();
}


