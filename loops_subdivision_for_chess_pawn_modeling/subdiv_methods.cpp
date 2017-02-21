#include "loop_subdivision.h"

void Scene::refineTopology()
{
    Model *tempM = HModel;
    while(tempM != NULL)
    {
        tempM->refineTopology();
        tempM->addFaceNormals();
        tempM = tempM->LLnext;
    }
}
void Scene::refineGeometry()
{
    Model *tempM = HModel;
    while(tempM != NULL)
    {
        tempM->refineGeometry();
        tempM = tempM->LLnext;
    }
}

void Scene::drawScene()
{
    Model *tempM = HModel;
    glPushMatrix();
    glScalef(0.4,0.4,0.4);
    
    //1
    
    glPushMatrix();
    glTranslatef(0.0,0.0, 0.2);
    tempM->drawModel();
    glPopMatrix();
    tempM = tempM->LLnext;

    //2
    glPushMatrix();
    glTranslatef(0.0,0.0,-0.05);
    tempM->drawModel();
    glPopMatrix();
    tempM = tempM->LLnext;
    
    //3
    glPushMatrix();
    glTranslatef(0.0,0.0,0.05);
    tempM->drawModel();
    glPopMatrix();
    tempM = tempM->LLnext;

    //4
    glPushMatrix();
    glTranslatef(0.0,0.0,2.00);
    tempM->drawModel();
    glPopMatrix();
    tempM = tempM->LLnext;

    //5
    glPushMatrix();
    glTranslatef(0.0,0.0,2.15);
    tempM->drawModel();
    glPopMatrix();
    tempM = tempM->LLnext;

    //6
    glPushMatrix();
    glTranslatef(0.0,0.0,3.10);
    tempM->drawModel();
    glPopMatrix();
    tempM = tempM->LLnext;

    glPopMatrix();
}

void Scene::addModel(Model *model)
{
    if ( HModel == NULL)
    {
        HModel = model;
        n_models++;
        return;
    }
    else
    {
        Model* tempM = HModel;
        while ( tempM->LLnext != NULL)
            tempM = tempM->LLnext;
        tempM->LLnext = model;
        n_models++;
        return;
    }
}

void Scene::subdivide(int k)
{
    for ( int level = 0 ; level < k ; level++)
    {
        Model *tempM = HModel;
        for ( int  i=0; i < n_models ; i++)
        {
            tempM->subdivide();
            tempM->addFaceNormals();
            tempM = tempM->LLnext;
        }
    }
}


int Edge::getSubdivisionMaskType()
{
    if (    vertex->getCategory() == DART_VERTEX ||  twin->vertex->getCategory() == DART_VERTEX 
        ||  vertex->getCategory() == SMOOTH_VERTEX ||  twin->vertex->getCategory() == SMOOTH_VERTEX )
        return SMOOTH_EDGE_SUBDIVISION_MASK;

    if ( vertex->getCategory() == twin->vertex->getCategory() )
        return REG_CREASE_EDGE_SUBDIVISION_MASK;

    if (    (vertex->getCategory() == CORNER_VERTEX && twin->vertex->getCategory() == NONREG_CREASE_VERTEX)
        ||  (vertex->getCategory() == NONREG_CREASE_VERTEX && twin->vertex->getCategory() == CORNER_VERTEX)  )
        return REG_CREASE_EDGE_SUBDIVISION_MASK;

    else return NONREG_CREASE_EDGE_SUBDIVISION_MASK;
    
}


int Vertex::getCategory()
{
    
    if( sharpEdgeIncidence == 0 )
        return SMOOTH_VERTEX;

    if( sharpEdgeIncidence == 1)
        return DART_VERTEX;
    if( sharpEdgeIncidence > 2 ) 
        return CORNER_VERTEX;

    assert(sharpEdgeIncidence == 2);

    /*determine whether currentV is extraordinary */
    int valency =0;
    GLboolean isExtraordinary = FALSE;
    Edge *tempE = edge;
    do
    {
        if ( tempE->twin == NULL)
        {
            isExtraordinary = TRUE;
            break;
        }
        assert(tempE->twin->next != NULL);
        valency++;
        tempE = tempE->twin->next;
    }while ( tempE != edge );
    if(!isExtraordinary)
    {
        if ( valency != 6 ) return NONREG_CREASE_VERTEX;
        tempE = edge;
        do
        {
            if (tempE->isSharp)
                if (tempE->twin->next->isSharp == FALSE &&
                    tempE->twin->next->twin->next->isSharp == FALSE &&
                    tempE->twin->next->twin->next->twin->next->isSharp == TRUE )
                    return REG_CREASE_VERTEX;
                else 
                    return NONREG_CREASE_VERTEX;
                tempE = tempE->twin->next;
        }while ( tempE != edge);
    }
        valency = 0;
    if(isExtraordinary)
    {
        //First, we go clock-wise round the link
        tempE = edge;
        do
        {               
            //assert(tempE->twin->next != NULL);
            valency++;
            tempE = tempE->twin->next;
            
        }while ( tempE->twin != NULL );
        
        //Then, we go counter-clock-wise round the link
        tempE = edge;
        do
        {
            //assert(tempE->prev->twin != NULL);
            valency++;
            tempE = tempE->prev->twin;
            
        }while ( tempE->twin != NULL );
            if ( valency == 4 ) return REG_CREASE_VERTEX;
        else return NONREG_CREASE_VERTEX;
    }
    
    assert(0); // should never come here
}

void Model::subdivide()
{
    refineTopology();
    refineGeometry();

}


void Model::refineTopology()
{
    //add new vertices ( one for each edge )
    //for each edge in the linked list
    Edge *tempE;//for temporary use
    Vertex *tempV;//for temporary use
    Face *tempF;//for temporary use
    Vertex *v1,*v2,*v3, *v4; //for temporary use
    Edge *e1,*e2,*e3,*e4,*oe1,*oe2,*oe3,*oe4,*oe5,*oe6, *ne1, *ne2, *ne3, *ne4 ,*ne5, *ne6, *parentEdge;
    Face *QF1, *QF2, *QF3, *QF4, *parentFace;
    
    /* Part (I) : create ONE vertex for EVERY TWO half-edges : each edge and its twin*/
    tempE = HEdge;
    while(tempE != NULL)
    {
        /* add new vertex for only those edges for which adding has not already been done.*/
        /* we add a single new vertex for a couple of edges : a half edge plus its twin*/
        if ( tempE->divided == FALSE)
        {
            //for normal edges
            if ( tempE->twin != NULL)
            {
                //add vertex
                Vertex *newVertex = addVertexToLL();
                //horizontal vertices
                v1 = tempE->vertex;
                v2 = tempE->twin->vertex;
                //vertical vertices
                v3 = tempE->next->vertex;
                v4 = tempE->twin->next->vertex;

                switch( tempE->getSubdivisionMaskType())
                {
                case SMOOTH_EDGE_SUBDIVISION_MASK:
                    newVertex->x = ( (float)3/8 * ( v1->x + v2->x) )+ ((float)1/8 * (v3->x + v4->x ));
                    newVertex->y = ( (float)3/8 * ( v1->y + v2->y) )+ ((float)1/8 * (v3->y + v4->y ));
                    newVertex->z = ( (float)3/8 * ( v1->z + v2->z) )+ ((float)1/8 * (v3->z + v4->z ));
                    break;
                case REG_CREASE_EDGE_SUBDIVISION_MASK:
                    newVertex->x = ( (float)1/2 * ( v1->x + v2->x) );
                    newVertex->y = ( (float)1/2 * ( v1->y + v2->y) );
                    newVertex->z = ( (float)1/2 * ( v1->z + v2->z) );
                    break;
                case NONREG_CREASE_EDGE_SUBDIVISION_MASK:
                    newVertex->x = ((float)5/8 * (v1->x )) + ((float)3/8 * (v2->x));
                    newVertex->y = ((float)5/8 * (v1->y )) + ((float)3/8 * (v2->y));
                    newVertex->z = ((float)5/8 * (v1->z )) + ((float)3/8 * (v2->z));
                    break;
                default: assert(0);// should not happen
                    
                }

                newVertex->edge = tempE; // just to keep reference of the old topology; it will be used at **THIS** place ( find "**THIS**")
                //DO NOT update n_vertices
            }
            else
            {
                //for boundary edges
                //add vertex
                Vertex *newVertex = addVertexToLL();
                //horizontal vertices
                v1 = tempE->vertex;
                v2 = tempE->prev->vertex;

                //no vertical vertices for boundary

                newVertex->x = ( (float)1/2 * ( v1->x + v2->x) );
                newVertex->y = ( (float)1/2 * ( v1->y + v2->y) );
                newVertex->z = ( (float)1/2 * ( v1->z + v2->z) );

                newVertex->edge = tempE; // just to keep reference of the old topology
                //DO NOT update n_vertices
            }
            tempE->divided = TRUE;
            if ( tempE->twin != NULL ) tempE->twin->divided = TRUE;
            
        }

        tempE = tempE->LLnext;
    }
    /*end of Part I */

    /* Part (II) : for EACH NEW vertex, create 4 half edges, and delete two old edges. */

    //while loop to reach the last old vertex ( n_vertices)
    tempV = HVertex;
    for ( int i = 0 ; i < n_vertices; i++)
    {
        tempV = tempV->LLnext;
    }//so we now point to the FIRST new vertex in the linked list

    old_n_vertices = n_vertices;

    //goal : add 4 half edges and delete old two
    //for each new vertex ( follow the linked list  till null )
    while( tempV != NULL)
    {
        parentEdge = tempV->edge; //**THIS**
        
        //e1 = FROM newvertex TO parentEdge->vertex 
        e1= addEdgeToLL();
        //e2 = FROM parentEdge->prev->vertex TO newvertex 
        //e2= addEdgeToLL();
        
        e1->prev = parentEdge;
        e1->next = parentEdge->next;
        parentEdge->next->prev = e1;
        e1->face =  parentEdge->face;
        e1->vertex = parentEdge->vertex;

        e2 = parentEdge;
        e2->next = e1;
        e2->vertex = tempV; //tempV means the new vertex
        /* e2's face is already parentEdge's face and e2's prev is already parentEdge's prev */
        /* e2's twin is still parentEdge->twin*/
        
        if ( parentEdge->isSharp == TRUE ) 
        {
            e1->isSharp = e2->isSharp = TRUE;
            parentEdge->vertex->sharpEdgeIncidence = 2;
            /*both e1 and e2 are incident on this vertex, both of which are sharp. And no other edge incident on it are sharp.*/
        }


        //adding e3 and e4 now
        if ( parentEdge->twin != NULL )
        {
            //e3 = FROM newvertex TO parentEdge->prev->vertex
            e3 = addEdgeToLL();
            //e4 = FROM parentEdge->vertex TO newvertex
            //e4 = addEdgeToLL();
            e3->prev = parentEdge->twin;
            e3->vertex = parentEdge->twin->vertex;
            e3->next = parentEdge->twin->next;
            parentEdge->twin->next->prev = e3;
            e3->face = parentEdge->twin->face;

            e4 = parentEdge->twin;
            e4->next = e3;
            e4->vertex = tempV; //tempV means the new vertex
            
            e3->twin = e2;
            e2->twin = e3;

            e1->twin = e4;
            e4->twin = e1;
            if ( parentEdge->twin->isSharp == TRUE ) 
            {
                e3->isSharp = e4->isSharp = TRUE;
            }
        }

        tempV->edge = e1; // we have to set ANY ONE edge emerging from the vertex
        
        //face to edge pointers updation is not required as the old edge are retained
        
        n_vertices++; /* vertices have been added in pass I only. Here we are doing a deferred update of n_vertices*/
        tempV = tempV->LLnext;

    }
    /* end of Part (II) */

    /* Part (III) : for EACH face(will be called parentFace) : add beechwali 6 ( 3*2 ) edges; make parentFace= QF1 & add 3 faces */

    int old_n_faces = n_faces;
    tempF = HFace;
    for ( int i=0 ;i < old_n_faces ; i++)
    {
        parentFace = tempF;

        oe1 = parentFace->edge;
        oe2 = oe1->next;
        oe3 = oe2->next;
        oe4 = oe3->next;
        oe5 = oe4->next;
        oe6 = oe5->next;


        //join following 3 vertices
        v1 = oe1->vertex;
        v2 = oe3->vertex;
        v3 = oe5->vertex;

        //create three edges by joining v1,v2,v3
        //ne1 = v1-v2
        ne1 = addEdgeToLL();
        //ne2 = v2-v3
        ne2 = addEdgeToLL();
        //ne3 = v3-v1  ( we have kept the direction of central new face ( QF4 ) CCW ( the same as the original triangle )
        ne3 = addEdgeToLL();
        
        //ne4 = v2-v1;
        ne4 = addEdgeToLL();
        //ne5 = v3-v2;
        ne5 = addEdgeToLL();
        //ne6 = v1- v3;
        ne6 = addEdgeToLL();

        /*now we have 6 new edges, and 6 old edges, togetherly making 4 different triangles. So we have to define prev-next and twin
        relationships for all these 12 edges.*/
        
        ne1->twin = ne4;
        ne4->twin = ne1;

        ne2->twin = ne6;
        ne6->twin = ne2;

        ne3->twin = ne5;
        ne5->twin = ne3;

        //QF1: prev next relationships

        oe1->next = ne5;
        ne5->prev = oe1;

        ne5->next = oe6;
        oe6->prev = ne5;

        oe6->next = oe1;
        oe1->prev = oe6;

        //QF2 : prev next relationships

        oe2->next = oe3;
        oe3->prev = oe2;

        oe3->next = ne4;
        ne4->prev = oe3;

        ne4->next = oe2;
        oe2->prev = ne4;

        //QF3 : prev next

        ne6->next = oe4;
        oe4->prev = ne6;

        oe4->next = oe5;
        oe5->prev = oe4;

        oe5->next = ne6;
        ne6->prev = oe5;

        //QF4 : prev next
        
        ne1->next = ne2;
        ne2->prev = ne1;

        ne2->next = ne3;
        ne3->prev = ne2;

        ne3->next = ne1;
        ne1->prev = ne3;

        // setting vertex

        ne1->vertex = ne6->vertex = v2;
        ne2->vertex = ne5->vertex = v3;
        ne3->vertex = ne4->vertex = v1;

        //NOW THE ONLY PART WHICH MUST REMAINS IS EDGE->FACE POINTERS OF ALL 12 EDGES AS NO QUARTER FACES HAVE BEEN CREATED YET
        
        
        //-create 4 new quarter faces : QF1, QF2, QF3, QF4 ; and add them to the linked list of faces
        QF1 = parentFace;
        QF2 = addFaceToLL();
        QF3 = addFaceToLL();
        QF4 = addFaceToLL();


        QF1->edge = oe1;
        QF2->edge = oe2;
        QF3->edge = oe4;
        QF4->edge = ne1;

        oe1->face = ne5->face = oe6->face = QF1;
        oe2->face = oe3->face = ne4->face = QF2;
        oe4->face = oe5->face = ne6->face = QF3;
        ne1->face = ne2->face = ne3->face = QF4;

        
        
        n_faces += 3;
        tempF = tempF->LLnext;
    }

    tempE = HEdge;
    int i =0;
    while (tempE != NULL)
    {
        tempE->divided = FALSE;
        i++;
        tempE = tempE->LLnext;
    }
}

Vertex* Model::addVertexToLL()
{
    Vertex *newVertex;

    newVertex = new Vertex;
    newVertex->LLnext = NULL;

    if ( HVertex == NULL )
    {
        HVertex = newVertex;
        return newVertex;
    }

    else
    {
        Vertex *temp = HVertex;
        while ( temp->LLnext != NULL )
        {
            //do timepass
            temp = temp->LLnext;
        }
        temp->LLnext = newVertex;
        return newVertex;

    }
}


Edge* Model::addEdgeToLL()
{
    Edge *newEdge;

    newEdge = new Edge;
    newEdge->LLnext = NULL;

    if ( HEdge == NULL )
    {
        HEdge = newEdge;
        return newEdge;
    }

    else
    {
        Edge *temp = HEdge;
        while ( temp->LLnext != NULL )
        {
            //do timepass
            temp = temp->LLnext;
        }
        temp->LLnext = newEdge;
        return newEdge;

    }
}

Face* Model::addFaceToLL()
{
    Face *newFace;

    newFace = new Face;
    newFace->LLnext = NULL;

    if ( HFace == NULL )
    {
        HFace = newFace;
        return newFace;
    }

    else
    {
        Face *temp = HFace;
        while ( temp->LLnext != NULL )
        {
            //do timepass
            temp = temp->LLnext;
        }
        temp->LLnext = newFace;
        return newFace;

    }
}


void Model::deleteEdgeFromLL(Edge *e)
{
    if ( e == HEdge)
    {
        assert( e->LLnext != NULL );
        HEdge = e->LLnext;
        delete e;
    }
    else
    {
        Edge *temp = HEdge;
        while ( temp->LLnext != e )
        {
            //do timepass
            temp = temp->LLnext;
        }
        if ( e->LLnext != NULL ) temp ->LLnext = e->LLnext;
        else temp->LLnext = NULL;
        
        //if ( e->twin != NULL && e->twin->twin == e ) e->twin->twin = NULL;
        
        if (e != NULL ) delete e;
        //e = NULL;
    }
}

void Model::deleteFaceFromLL(Face *f)
{
    if ( f == HFace)
    {
        assert( f->LLnext != NULL );
        HFace = f->LLnext;
        delete f;
    }
    else
    {
        Face *temp = HFace;
        while ( temp->LLnext != f )
        {
            //do timepass
            temp = temp->LLnext;
        }
        if ( f->LLnext != NULL ) temp ->LLnext = f->LLnext;
        else temp->LLnext = NULL;

        delete f;
    }
}


//--------------------------------------------------------------------------------------------------------------------------

/* return some pre-generated betas */
GLfloat beta(GLuint n)
{
    return (5.0/4.0 - (3+2*cos(2*M_PI / n)) * (3+2*cos(2*M_PI / n)) / 32);
}


GLfloat alpha(GLuint n)
{
    GLfloat b;
 
    b = beta(n);
    
    return n * (1 - b) / b; 
}



void Model::refineGeometry()
{
    Vertex *currentV, *tempV, *extremeRightV, *extremeLeftV;
    Edge *tempE;
    int noOfNeighbourVertices=0;
    float result_alpha=0.0;

    GLboolean isExtraordinary=FALSE;
    
    /* goal : for each vertex, find its new position*/
    currentV = HVertex;
    for ( int i=0 ; i< old_n_vertices; i++ )
    {
        if ( currentV->refinedV == NULL ) currentV->refinedV = new Vertex; // allocate it only on the first time

        /*determine whether currentV is extraordinary */
        noOfNeighbourVertices =0;
        isExtraordinary = FALSE;
        currentV->refinedV->x = currentV->refinedV->y = currentV->refinedV->z = 0.0;
        tempE = currentV->edge;
        do
        {
            if ( tempE->twin == NULL)
            {
                isExtraordinary = TRUE;
                break;
            }
            assert(tempE->twin->next != NULL);
            noOfNeighbourVertices++;
            tempE = tempE->twin->next;
        }while ( tempE != currentV->edge );

        if(isExtraordinary)
        {
            //printf("\n(%f,%f,%f) is Extraordinary. noOfNeighbourVertices = %d ",currentV->x, currentV->y, currentV->z, noOfNeighbourVertices);
            //fflush(stdout);
            //First, we go clock-wise round the link
            tempE = currentV->edge;
            while ( tempE->twin != NULL )
            {               
                //assert(tempE->twin->next != NULL);
                tempE = tempE->twin->next;
                
            }
            extremeRightV = tempE->vertex;
            
            //Then, we go counter-clock-wise round the link
            tempE = currentV->edge;
            while ( tempE->twin != NULL )
            {
                //assert(tempE->prev->twin- != NULL);
                if ( tempE->prev->twin == NULL ) break;
                tempE = tempE->prev->twin;
                
            }
            extremeLeftV = tempE->prev->prev->vertex;

            currentV->refinedV->x = (float)1/8 * ( extremeRightV->x + extremeLeftV->x );
            currentV->refinedV->y = (float)1/8 * ( extremeRightV->y + extremeLeftV->y );
            currentV->refinedV->z = (float)1/8 * ( extremeRightV->z + extremeLeftV->z );

            currentV->refinedV->x = ( (float)6/8 * currentV->x ) + currentV->refinedV->x;
            currentV->refinedV->y = ( (float)6/8 * currentV->y ) + currentV->refinedV->y;
            currentV->refinedV->z = ( (float)6/8 * currentV->z ) + currentV->refinedV->z;

        }

        else // NOT-extraordinary case
        {
            //printf("\n(%f,%f,%f) is NOT Extraordinary. noOfNeighbourVertices = %d",currentV->x, currentV->y, currentV->z, noOfNeighbourVertices);
            //fflush(stdout);

            if ( currentV->getCategory() == SMOOTH_VERTEX  || currentV->getCategory() == DART_VERTEX )
            {

                currentV->refinedV->x = currentV->refinedV->y = currentV->refinedV->z = 0.0;
                
                /* first, take summation of all the vertices in the neighborhood of currentV*/
                tempE = currentV->edge;
                do
                {
                    currentV->refinedV->x += tempE->vertex->x; 
                    currentV->refinedV->y += tempE->vertex->y;
                    currentV->refinedV->z += tempE->vertex->z;

                tempE = tempE->twin->next;
                }while ( tempE != currentV->edge);

                result_alpha = alpha(noOfNeighbourVertices);

                currentV->refinedV->x = ( (result_alpha) * currentV->x ) + currentV->refinedV->x;
                currentV->refinedV->y = ( (result_alpha) * currentV->y ) + currentV->refinedV->y;
                currentV->refinedV->z = ( (result_alpha) * currentV->z ) + currentV->refinedV->z;
            
                currentV->refinedV->x =  currentV->refinedV->x / (result_alpha+noOfNeighbourVertices);
                currentV->refinedV->y =  currentV->refinedV->y / (result_alpha+noOfNeighbourVertices);
                currentV->refinedV->z =  currentV->refinedV->z / (result_alpha+noOfNeighbourVertices);
            }
            else if ( currentV->getCategory() == REG_CREASE_VERTEX )
            {
                
                tempE = currentV->edge;
                while (tempE->isSharp != TRUE )
                {
                    tempE = tempE->twin->next;
                }
                
                currentV->refinedV->x = ((float)4/6 * currentV->x) + ( (float)1/6 * (tempE->vertex->x + tempE->twin->next->twin->next->twin->next->vertex->x)) ;
                currentV->refinedV->y = ((float)4/6 * currentV->y) + ( (float)1/6 * (tempE->vertex->y + tempE->twin->next->twin->next->twin->next->vertex->y)) ;
                currentV->refinedV->z = ((float)4/6 * currentV->z) + ( (float)1/6 * (tempE->vertex->z + tempE->twin->next->twin->next->twin->next->vertex->z)) ;
            }

            else if ( currentV->getCategory() == NONREG_CREASE_VERTEX )
            {
                tempE = currentV->edge;
                while (tempE->isSharp != TRUE )
                {
                    tempE = tempE->twin->next;
                }
                //first sharp edge detected 
                
                currentV->refinedV->x = ((float)3/5 * currentV->x) + ( (float)1/5 * (tempE->vertex->x) ) ;
                currentV->refinedV->y = ((float)3/5 * currentV->y) + ( (float)1/5 * (tempE->vertex->y) ) ;
                currentV->refinedV->z = ((float)3/5 * currentV->z) + ( (float)1/5 * (tempE->vertex->z) ) ;
                
                //continued....
                tempE = tempE->twin->next;
                while (tempE->isSharp != TRUE )
                {
                    tempE = tempE->twin->next;
                }
                //the next sharp edge detected 
                
                currentV->refinedV->x = ( currentV->refinedV->x) + ( (float)1/5 * (tempE->vertex->x) ) ;
                currentV->refinedV->y = ( currentV->refinedV->y) + ( (float)1/5 * (tempE->vertex->y) ) ;
                currentV->refinedV->z = ( currentV->refinedV->z) + ( (float)1/5 * (tempE->vertex->z) ) ;

            }

        }

        currentV = currentV->LLnext;
    }//while loop ends

    tempV = HVertex;
    for ( int i=0 ; i< old_n_vertices ; i++ )
    {
        tempV->x = tempV->refinedV->x;
        tempV->y = tempV->refinedV->y;
        tempV->z = tempV->refinedV->z;
        tempV = tempV->LLnext;
    } 
}



//--------------------------------------------------------------------------------------------------------------------------

GLfloat dotProduct(Vertex* u, Vertex* v)
{
    assert(u); assert(v);
    return (u->x * v->x + u->y * v->y + u->z * v->z);
}

Vertex* crossProduct(Vertex* u, Vertex* v)
{
    Vertex* result;
    assert(u); assert(v);

    result = new Vertex ;
    result->x = u->y * v->z - u->z * v->y;
    result->y = u->z * v->x - u->x * v->z;
    result->z = u->x * v->y - u->y * v->x;

    return result;
}

void normalizeVector(Vertex* v)
{
    GLfloat l;

    assert(v);
    l = (GLfloat)sqrt(dotProduct(v,v));
    v->x /= l;
    v->y /= l;
    v->z /= l;
}


void Face::addNormal()
{
    Vertex *edge1 = new Vertex ;
    Vertex *edge2 = new Vertex ;

    edge1->x = edge->prev->vertex->x - edge->vertex->x;
    edge1->y = edge->prev->vertex->y - edge->vertex->y;
    edge1->z = edge->prev->vertex->z - edge->vertex->z;
    
    edge2->x = edge->prev->prev->vertex->x - edge->prev->vertex->x;
    edge2->y = edge->prev->prev->vertex->y - edge->prev->vertex->y;
    edge2->z = edge->prev->prev->vertex->z - edge->prev->vertex->z;

    normal = crossProduct(edge1, edge2);
    normalizeVector(normal);
         
    free(edge1);
    free(edge2);
}

void Model::addFaceNormals()
{
    Face *tempF = HFace;
    while ( tempF != NULL )
    {
        tempF->addNormal();
        tempF = tempF->LLnext;
    }
}

void Model::drawModel()
{
    assert(HVertex != NULL);
    assert(HEdge!= NULL);
    assert(HFace!= NULL);
    
    
    Face *tempF;
    Edge *tempE;
    tempF = HFace;
    while ( tempF != NULL)
    {
        glBegin(GL_TRIANGLES);
        glNormal3f(tempF->normal->x, tempF->normal->y, tempF->normal->z );
        tempE = tempF->edge;
        do
        {
            //if ( tempE->isSharp == TRUE  ) glColor3f(1.0, 0.0,0.0 );
            //if ( tempE->twin != NULL ) 
            //  if ( tempE->twin->isSharp == TRUE ) glColor3f(1.0, 0.0,0.0 );
            //else glColor3f(1.0,1.0,1.0);
            
            glVertex3f(tempE->vertex->x,tempE->vertex->y,tempE->vertex->z);
            tempE = tempE->next;

        }while ( tempE != tempF->edge);

        glEnd();
        tempF = tempF->LLnext;
    }


}


void Model::printModel()
{
    printf("\nVerts : %d\nfaces : %d\nSharp Edges : %d", n_vertices, n_faces , n_sharp_edges);
    printf("following are sharp edges: \n------------------------------");
    Edge *tempE = HEdge;
    for ( int i=0 ; i < n_faces *3 ; i++)
    {
        if ( tempE->isSharp == TRUE )
            printf("\n( %f, %f, %f ) TO ( %f, %f, %f )", tempE->vertex->x, tempE->vertex->y, tempE->vertex->z, tempE->twin->vertex->x, tempE->twin->vertex->y, tempE->twin->vertex->z);

        tempE = tempE->LLnext;
    }
    

}

//-------------------------------------------------------------------------------------------------------------------


void Model::readFromOBJ(char filename[])
{
    FILE* fpOBJ = fopen(filename,"r");
    char buffer[128];
    if ( fpOBJ == NULL )
    {
        printError("fpOBJ NULL");
    }

    
    while( fgets(buffer,sizeof(buffer),fpOBJ) != NULL )
    {
        switch(buffer[0])
        {
        case 'v':
            n_vertices++;
            break;

        case 'f':
            n_faces++;
            break;
        case 's':
            n_sharp_edges++;
            break;

        default:
            break;

        }
    }

    fclose(fpOBJ);
    //now we build the model
    buildModel(filename);

}


void Model::buildModel(char* filename)
/* buildModel follows after readFromOBJ(). In readFromOBJ, we count the number of vertices and faces specified in the .obj file.
    The, subsequently buildModel() is automatically called which extracts the vertex and face information from the .obj file into
    the Vertices and Triangles array.

    The vertices[] and faces[] collections in the model are basically linked lists. Although, now that we have the n_vertices and 
    n_faces ( w.r.t. the base mesh ) ready with us ( readFromOBJ() has alredy done it ), we allocate the known number of the vertices
    and faces as arrays. ( We make sure that the link list nature is stil preserved, as we *do* maintains the links between the array
    nodes, so that, later when we subdivide and add more vertices and faces, we can easily extend the current link-lists.
    Allocating the base-mesh vertices and faces as arrays makes the neighborhood detection and initialisation part to be simple to 
    implement.
*/
{
    FILE* fpOBJ = fopen(filename,"r");
    char buffer[128];
    if ( fpOBJ == NULL )
    {
        printError("fpOBJ NULL");
    }
    
    /*linklist allocations as arrays*/
    Vertex *arrVertices = new Vertex[n_vertices];

    Edge *arrEdges = new Edge[ 3* n_faces ];
    
    Face *arrFaces = new Face[n_faces];

    Triangle *arrTriangles = new Triangle[n_faces];

    /*get the OBJ info into arrays*/
    n_vertices = n_faces = 0;

    int sharp_edge_vertex1=-1,sharp_edge_vertex2=-1, SEcounter=0; //sharp-Edge counter

    int **sharpEdges = new int*[n_sharp_edges];
    for ( int i=0 ; i < n_sharp_edges ; i++)
        sharpEdges[i] = new int[2];
    

    while( fscanf(fpOBJ,"%s",buffer) != EOF )
    {
        switch(buffer[0])
        {
        case 'v':
            fscanf(fpOBJ, "%f %f %f", &(arrVertices[n_vertices].x), &(arrVertices[n_vertices].y), &(arrVertices[n_vertices].z) );
            //printf( "\n%f %f %f", arrVertices[n_vertices].x, arrVertices[n_vertices].y, arrVertices[n_vertices].z );
            n_vertices++;
            break;

        case 'f':
            fscanf(fpOBJ, "%d %d %d", &(arrTriangles[n_faces].verts[0]), &(arrTriangles[n_faces].verts[1]), &(arrTriangles[n_faces].verts[2]) );
            arrTriangles[n_faces].verts[0]--; //as indices in OBJ start from 1; and indices of our array start from 0
            arrTriangles[n_faces].verts[1]--;
            arrTriangles[n_faces].verts[2]--;
            //printf( "\n%d %d %d", arrTriangles[n_faces].verts[0], arrTriangles[n_faces].verts[1], arrTriangles[n_faces].verts[2] );
            n_faces++;          
            break;
        
        case 's':
            fscanf(fpOBJ, "%d %d", &sharp_edge_vertex1, &sharp_edge_vertex2 );
            sharp_edge_vertex1--; sharp_edge_vertex2--; //as indices in OBJ start from 1; and indices of our array start from 0
            arrVertices[sharp_edge_vertex1].sharpEdgeIncidence++;
            arrVertices[sharp_edge_vertex2].sharpEdgeIncidence++;
            sharpEdges[SEcounter][0] = sharp_edge_vertex1;
            sharpEdges[SEcounter][1] = sharp_edge_vertex2;
            SEcounter++;
            break;


        default:
            fgets(buffer,sizeof(buffer),fpOBJ);
            break;

        }
    }

    /*forming the edges linked list and relating edges to faces and faces to edges*/
    for (int i=0; i <= (n_faces-1) ; i++)
    {
        for ( int j=0 ; j < 3 ; j++)
        {
            /*relating every edge to the vertex-that-it-points-to and the-vertex-it-starts-from to the edge.*/
            arrEdges[i*3+j].vertex = &(arrVertices[ arrTriangles[i].verts[j] ]);
            //printf("\n%dth face : %dth edge : making %dth VERTEX of %ldth edge",i,j,arrTriangles[i].verts[j],&(arrEdges[i*3+j]) );
            if ( arrVertices[ arrTriangles[i].verts[(j+2)%3]].edge == NULL )
            {
                arrVertices[ arrTriangles[i].verts[(j+2)%3]].edge  =  &(arrEdges[i*3 + j]);
                //printf("\n%dth face : %dth edge : making %ld EDGE of %dth vertex",i,j, &(arrEdges[i*3 + j]), arrTriangles[i].verts[(j+2)%3] );
            }

            arrEdges[i*3+j].prev = &(arrEdges[i*3+ ((j+2)%3)]);
            arrEdges[i*3+j].next = &(arrEdges[i*3+ ((j+1)%3)]);

            arrEdges[i*3+j].face = &(arrFaces[i]);
            arrFaces[i].edge = &(arrEdges[i*3]);
        }
    }

    /*link the nodes */
    HVertex = arrVertices;
    HEdge = arrEdges;
    HFace = arrFaces;
    for (int i=0; i <= (n_vertices-2) ; i++)
        arrVertices[i].LLnext = &(arrVertices[i+1]);
    arrVertices[n_vertices-1].LLnext = NULL;

    for (int i=0; i <= (n_faces-2) ; i++)
        arrFaces[i].LLnext = &(arrFaces[i+1]);
    arrFaces[n_faces-1].LLnext = NULL;
    
    for (int i=0; i <= (n_faces*3 -2) ; i++)
        arrEdges[i].LLnext = &(arrEdges[i+1]);
    arrEdges[n_faces*3-1].LLnext = NULL;

    /*setting up twin edges*/
    for ( int i=0 ; i < n_faces*3 ; i++)
    {
        if ( arrEdges[i].twin == NULL )
        {
            for ( int j=i+1; j < n_faces*3 ; j++)
            {
                if ( arrEdges[i].prev->vertex == arrEdges[j].vertex &&  arrEdges[j].prev->vertex == arrEdges[i].vertex )
                {
                    arrEdges[i].twin = &(arrEdges[j]);
                    arrEdges[j].twin = &(arrEdges[i]);
                }
            }
        }
        //assert(arrEdges[i].twin != NULL );
    }

    /*based on the sharp-edges data obtained from OBJ, mark sharp edges in the linked list as sharp*/
    for (int i=0; i< n_faces*3 ; i++)
    {
        if ( arrEdges[i].isSharp == FALSE )
        {
            for ( int j=0; j < n_sharp_edges; j++)
            {
                if (    ( arrEdges[i].vertex == &(arrVertices[sharpEdges[j][0]]) && arrEdges[i].twin->vertex == &(arrVertices[sharpEdges[j][1]]) )
                    ||  (arrEdges[i].vertex == &(arrVertices[sharpEdges[j][1]]) && arrEdges[i].twin->vertex == &(arrVertices[sharpEdges[j][0]]) ) )
                {
                    arrEdges[i].isSharp = TRUE;
                    arrEdges[i].twin->isSharp = TRUE;
                }
            }
        }
    }




    addFaceNormals();


    fclose(fpOBJ);



}

//------------------------------------------ debugging functions ----------------------------------------------------------------


void Model::assertFullNeighbourhoodInfo()
{
    Edge *EDGE = HFace->edge;

    Face *tempF= HFace;
    int i =0, j=0; // only for debugging
    Edge *tempE = HFace->edge ;
    while ( tempF != NULL )
    {
        tempE = tempF->edge ;
        j=0;
        do
        {
            //printf("\r\n # %dth face , %dth edgge",i,j);
            assert( tempE->prev->next == tempE );
            assert(tempE->face == tempF );
            //assert( tempE->twin != NULL );

            assert( tempE->vertex != NULL );
            assert( tempE->vertex->edge != NULL );
            assert( tempE->next != NULL );
            
            //assert( tempE->twin->vertex == tempE->prev->vertex );

            tempE = tempE->next ; 
            j++;
        } while ( tempE != tempF->edge );
        
        tempF = tempF->LLnext ;
        i++;
    }
    Vertex *tempV = HVertex;
    for ( int i =0 ; i < n_vertices ; i++ )
    {
        //printf("\ncategory of (%f, %f, %f) is %d ", tempV->x, tempV->y, tempV->z, tempV->getCategory() );
        tempV = tempV->LLnext;
    }

}


void printTwins(Edge *HEdge)
{
    Edge *temp= HEdge;
    while ( temp != NULL)
    {
        //printf("\r\n%ld %ld", temp, temp->twin );
        temp = temp->LLnext ;
    }
}







void printVLinkedList ( Vertex *HVertex)
{
    Vertex *temp= HVertex;
    while ( temp != NULL)
    {
        //printf( "\n%f %f %f", temp->x,  temp->y,  temp->z );
        temp = temp->LLnext;
    }
}

void printFLinkedList ( Face *HFace)
{
    Face *temp= HFace;
    while ( temp != NULL)
    {
        //printf( "\n%f %f %f", temp->edge->vertex->x,  temp->edge->vertex->y,  temp->edge->vertex->z );
        temp = temp->LLnext;
    }
}
void printELinkedList ( Edge *HEdge)
{
    Edge *temp= HEdge;
    while ( temp != NULL)
    {
        //printf("\n%ld %ld",temp->vertex->edge, temp->next);
        //printf( "\n%f %f %f ", temp->vertex->x,  temp->vertex->y,  temp->vertex->z );
        //printf( "\tNext : %f %f %f ", temp->next->vertex->x,  temp->next->vertex->y,  temp->next->vertex->z );
        temp = temp->LLnext;
    }
}
