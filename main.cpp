//============================================================================
// Name        : main.cpp
// Author      : David Briemann
// Version     :
// Copyright   :
// Description :
//============================================================================


#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <climits>
using namespace std;

#include <tclap/CmdLine.h>
#include <tclap/SwitchArg.h>
#include <tclap/ValueArg.h>
#include <tclap/ArgException.h>
#include <GL/glut.h>
#include <boost/thread.hpp>

#include "Conformation.hpp"
#include "Protein.hpp"
#include "Population.hpp"



//FUNCTION HEADERS
void drawAminoAcid(double, double, double, char);
void render();
void calculation(Population *);
void drawConnector(double, double, double, double);
void renderBitmapString(float, float, void *, char *);
void usekeys(unsigned char, int, int);


//GLOBAL VARIABLES
Conformation *globalFittestPtr;
int isTerminated = 0;
int xWindowSize = 1000;
int yWindowSize = 700;
float xView = xWindowSize*2/3;
float yView = yWindowSize*2/3;
int xScaleFactor = xWindowSize / 10;
int yScaleFactor = yWindowSize / 10;
float scaleFactor = 1.0;

//command line switches
bool switch_enable_graphics = false;
unsigned int switch_max_evaluations = 100000;
float switch_crossover_prob = 0.6;
float switch_mutation_prob = 0.08;
unsigned int switch_popsize = 100;
string switch_protein = "";
int switch_minen = INT_MIN;



int main(int argc, char* argv[]) {

    //process command line arguments
    try {
        TCLAP::CmdLine cmd("Folding Visualizer - HP 2d Model", ' ', "0.5");
        TCLAP::SwitchArg enableGraphics("g", "graphics", "Enables OpenGL Window", cmd, false);
        TCLAP::ValueArg<float> setCros("c", "crossprob", "Sets the crossover probability", false, 1.0, "float 0..1", cmd);
        TCLAP::ValueArg<float> setMut("m", "mutateprob", "Sets the mutation probability", false, 0.08, "float 0..1", cmd);
        TCLAP::ValueArg<unsigned int> setPopsize("s", "popsize", "Sets the population size", false, 100, "pos. int", cmd);
        TCLAP::ValueArg<string> setProtein("p", "protein", "Sets the protein", true, "", "string", cmd);
        TCLAP::ValueArg<int> setMinen("y", "minenergy", "Sets the minimum energy", false, INT_MIN, "neg. int", cmd);
        TCLAP::ValueArg<unsigned int> setMaxeval("e", "maxeval", "Sets the maximum number of energy evaluations",
                                                false, 100000, "pos. int", cmd);

        // Parse the argv array.
        cmd.parse( argc, argv );
        switch_enable_graphics = enableGraphics.getValue();
        switch_crossover_prob = setCros.getValue();
        switch_mutation_prob = setMut.getValue();
        switch_popsize = setPopsize.getValue();
        switch_protein = setProtein.getValue();
        switch_max_evaluations = setMaxeval.getValue();

    } catch (TCLAP::ArgException &e) { // catch any exceptions
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
    }

    //random seed
    srand(time(NULL));

    //Protein p("WBWwB");
    //Protein p("WBWWBWWBBWWB");
    //Protein p("BWBWWBBWBWWBWBBWWBWB"); //MAXEN = -9
    //Protein p("WWBWWBBWWBBWWWWWBBBBBBBBBBWWWWWWBBWWBBWWBWWBBBBB"); //MAXEN = -23
    //Protein p("BBBBBBBBBBBBWBWBWWBBWWBBWWBWWBBWWBBWWBWWBBWWBBWWBWBWBBBBBBBBBBBB"); //MAXEN = -40
    //Protein p("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB");
    //Protein p("BBBBWWWWBBBBBBBBBBBBWWWWWWBBBBBBBBBBBBWWWBBBBBBBBBBBBWWWBBBBBBBBBBBBWWWBWWBBWWBBWWBWB"); //MAXEN = -52

    Protein p( switch_protein );

    Population pop( switch_popsize, p, switch_mutation_prob, switch_crossover_prob);

    //create and start thread for calculation
    boost::thread calcThread(&calculation, &pop);


    if( switch_enable_graphics ) {
        //openGl part
        glutInit(&argc, argv);									//init glut
        glutInitWindowSize(xWindowSize, yWindowSize);			//define the window size
        glutInitWindowPosition(0,0);							//Position the window
        glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB|GLUT_DOUBLE);	//Define the drawing mode
        glutCreateWindow("Folding Visualizer - HP 2d Model");
        gluOrtho2D(-xView, xView, -yView, yView); 				//define view area
        glClearColor(0.7,0.7,0.7,0);							//Define our background color (grey)

        //register render function
        glutDisplayFunc(render);

        glutIdleFunc(render);

        //register keyboard function
        glutKeyboardFunc(usekeys);

        glutMainLoop();
    }
    //if opengl is closed wait for calculation
    calcThread.join();

    return 0;
}


void drawAminoAcid(double x, double y, double size, char kind) {
    if( kind == 'B' ) {
        glColor3f(0,0,0);//Change the color to black
    } else if ( kind == 'W' ) {
        glColor3f(1,1,1);//Change the color to white
    }

    glBegin(GL_QUADS);//Start drawing quads
        glVertex2f(x-(size/2),y-(size/2));//first coordinate
        glVertex2f(x-(size/2),y+(size/2));//second coordinate
        glVertex2f(x+(size/2),y+(size/2));//third coordinate
        glVertex2f(x+(size/2),y-(size/2));//last coordinate*/
    glEnd();//Stop drawing quads
}

void drawConnector(double x1, double y1, double x2, double y2) {
    glBegin(GL_LINES);
        glColor3f(0,0,0);//Change the color to black
        glVertex2f(x1, y1);
        glVertex2f(x2, y2);
    glEnd();
}


void render() {
    int size = 16;
    double x = 0;
    double y = 0;
    double x2 = 0;
    double y2 = 0;

    string tt = globalFittestPtr->getStatusString();

    glClear(GL_COLOR_BUFFER_BIT);//Clear the screen

    //draw info text
    if(!isTerminated) {
        glColor3f(1,1,1); //white
    }
    /*TODO FONT VIEW*/
    renderBitmapString( (-xView + xScaleFactor) * (2.0 - scaleFactor),
                        (yView - yScaleFactor) * (2.0 - scaleFactor),
                        GLUT_BITMAP_HELVETICA_18, const_cast<char *>(tt.c_str()));

    if(isTerminated) {
        tt = "TERMINATED";
        glColor3f(0,0,1);//blue
        /*TODO FONT VIEW*/
          renderBitmapString( (-xView + xScaleFactor) * (2.0 - scaleFactor),
                                (yView - yScaleFactor*1.5) * (2.0 - scaleFactor),
                                GLUT_BITMAP_HELVETICA_18, const_cast<char *>(tt.c_str()));
    }

    //draw connector to first amino acid
    x = Conformation::extractX(globalFittestPtr->getAbsAt(0)) * size * 2;
    y = Conformation::extractY(globalFittestPtr->getAbsAt(0)) * size * 2;
    drawConnector(x,y,xView,yView);

    //for each amino acid connector
    for( int i = 1; i < globalFittestPtr->getProtein()->getLength() - 1; i++ ) {
        x2 = Conformation::extractX(globalFittestPtr->getAbsAt(i)) * size * 2;
        y2 = Conformation::extractY(globalFittestPtr->getAbsAt(i)) * size * 2;

        //draw connector to previous acid
        drawConnector(x,y,x2,y2);

        x = x2;
        y = y2;
    }

    x = Conformation::extractX(globalFittestPtr->getAbsAt( globalFittestPtr->getProtein()->getLength() - 1 )) * size * 2;
    y = Conformation::extractY(globalFittestPtr->getAbsAt( globalFittestPtr->getProtein()->getLength() - 1 )) * size * 2;

    //draw connector to the last acid
    drawConnector(x,y,x2,y2);

    //draw first amino acid
    x = Conformation::extractX(globalFittestPtr->getAbsAt(0)) * size * 2;
    y = Conformation::extractY(globalFittestPtr->getAbsAt(0)) * size * 2;
    drawAminoAcid(x , y, size, globalFittestPtr->getProtein()->getNth(0));

    //for each amino acid from second to second last
    for( int i = 1; i < globalFittestPtr->getProtein()->getLength() - 1; i++ ) {
        x2 = Conformation::extractX(globalFittestPtr->getAbsAt(i)) * size * 2;
        y2 = Conformation::extractY(globalFittestPtr->getAbsAt(i)) * size * 2;

        //draw acid
        drawAminoAcid(x2, y2, size, globalFittestPtr->getProtein()->getNth(i));

        x = x2;
        y = y2;
    }

    //draw last amino acid
    x = Conformation::extractX(globalFittestPtr->getAbsAt( globalFittestPtr->getProtein()->getLength() - 1 )) * size * 2;
    y = Conformation::extractY(globalFittestPtr->getAbsAt( globalFittestPtr->getProtein()->getLength() - 1 )) * size * 2;
    drawAminoAcid(x , y, size, globalFittestPtr->getProtein()->getNth(globalFittestPtr->getProtein()->getLength() - 1));

    glutSwapBuffers();
}

void calculation( Population *pop) {
    globalFittestPtr = pop->getFittest();

    while(pop->getFittest()->getFitness() > switch_minen && Conformation::energyEvalSteps < switch_max_evaluations ) {
        pop->crossover();

        if( pop->getFittest()->getFitness() < globalFittestPtr->getFitness() ) {
            globalFittestPtr = pop->getFittest();

            //if graphics is disabled: output ascii status and picture to console
            if( !switch_enable_graphics ) {
                cout << globalFittestPtr->getStatusString() << endl;
                globalFittestPtr->printAsciiPicture();
            }
        }
    }

    isTerminated = true;
}


void renderBitmapString(float x, float y, void *font, char *str) {
    char *c;
    glRasterPos2f(x, y);
    for (c=str; *c != '\0'; c++) {
        glutBitmapCharacter(font, *c);
    }
}

// keyboard function
void usekeys(unsigned char key, int x, int y) {
    switch(key) {
    case '+':
        scaleFactor = 1.1;
        glScalef(scaleFactor, scaleFactor, 1.0);
        break;

    case '-':case'r':
        scaleFactor = 0.9;
        glScalef(scaleFactor, scaleFactor, 1.0);
        break;

    default:
        break;
    }
}
