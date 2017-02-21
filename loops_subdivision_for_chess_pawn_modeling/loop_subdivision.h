#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include <GLUT/glut.h>
#include <GLUT/glext.h> /*for occlusion querries*/


#define FALSE 0
#define TRUE 1

#ifndef M_PI
#define M_PI 3.14159265f
#endif

#define BASIC_MODE 0
#define LEVEL_MODE 1
#define ZOOM_MODE 2

#define ALPHA_MAX 20
#define ALPHA_LIMIT 0.469 /* converges to ~ 0.469 */

#define BETA_MAX 20
#define BETA_LIMIT 0.469 /* converges to ~ 0.469 */


class Edge;
class Face;
class Triangle;

#define SMOOTH_VERTEX 0
#define DART_VERTEX 1
#define REG_CREASE_VERTEX 2
#define NONREG_CREASE_VERTEX  3
#define CORNER_VERTEX 4

#define SMOOTH_EDGE_SUBDIVISION_MASK 1
#define REG_CREASE_EDGE_SUBDIVISION_MASK 2
#define NONREG_CREASE_EDGE_SUBDIVISION_MASK 3

class Vertex
{
    public:
    float x,y,z; // position
    Vertex *refinedV; // when we refine its position (in refineGeometry()), we use this as accumulator for 
    Edge* edge;
    int sharpEdgeIncidence;// number of sharp edges incident on this vertex
    Vertex* LLnext;

    Vertex()
    {
        x = y = z = 0.0;
        refinedV = NULL; // will be allocated on need basis. 
        edge = NULL ;
        sharpEdgeIncidence = 0;
        LLnext = NULL;

    }
    int getCategory();
};

class Edge
{
    public:
    GLboolean divided;
    GLboolean isSharp;
    Vertex *vertex;
    Edge* twin;
    Face *face;
    Edge *prev,*next;
    Edge *LLnext;

    
    Edge()
    {
        divided = FALSE;
        isSharp = FALSE;
        vertex = NULL ;  
        twin =  prev = next = LLnext = NULL; 
        face = NULL;
    }

    int getSubdivisionMaskType();
};

class Face
{
    public:
    Edge *edge;
    Face *LLnext;
    Vertex *normal;
    Face()
    {
        normal = NULL;
        edge = NULL; 
        LLnext = NULL;
    }
    void addNormal();
};

class Model
{
    public:
    int old_n_vertices;
    int n_vertices;
    int n_faces;
    int n_sharp_edges;

    Vertex *HVertex;
    Edge *HEdge;
    Face *HFace;
    Model *LLnext;
    
    Model()
    {
        old_n_vertices = n_vertices = n_faces = n_sharp_edges = 0;
        HVertex = NULL;
        HEdge = NULL ;
        HFace = NULL;
        LLnext = NULL;
    }

    void addFaceNormals();
    void printModel();
    void drawModel();
    void readFromOBJ(char filename[]);
    void buildModel(char* filename);
    void assertFullNeighbourhoodInfo();
    void subdivide();
    void refineTopology();
    void deleteEdgeFromLL(Edge *e);
    void deleteFaceFromLL(Face *f);
    Vertex* addVertexToLL();
    Edge* addEdgeToLL();
    Face* addFaceToLL();

    void refineGeometry();

};


class Scene
{
public:
    Model* HModel;
    int n_models;

    Scene()
    {
        HModel = NULL;
        n_models = 0;
    }

    void drawScene();
    void subdivide(int k);
    void addModel(Model *model);
    void refineTopology();
    void refineGeometry();
};


class Triangle
{
public:
    int verts[3];
};



void printError(char *errorString);
void printVLinkedList ( Vertex *HVertex);
void printFLinkedList ( Face *HFace);
void printELinkedList ( Edge *HEdge);
void printTwins(Edge *HEdge);


