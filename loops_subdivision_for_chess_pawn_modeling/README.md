Application of Loop's Subdivision Surface to Model a Chess Pawn
==========

* Author:  Neeraj Kulkarni <neerajkulkarni.iitd@gmail.com>
* License: GNU General Public License, version 3 (GPL-3.0)
* Development date: September 2010

A very simple implementation for a class project in C++

INSTALL NOTES
---------

Tested on Ubuntu 12.04
       
       cd into directory
       make
       ./loop.out


IMPLEMENTATION NOTES
---------

The program is written in C++ and uses the openGL/GLUT/GLU framework for rendering.

The pawn consists of six parts as shown in the image -- each designed separately by manually designing the base shapes (programmatically generated base shapes) and then joined them together to make the final model. The code is written from scratch to generate k levels of Loop subdivision, where k (an integer) is a command line argument. The base mesh is stored in obj format with markers for the boundaries and edges in an input file. An OpenGL viewer is provided for viewing subdivided model and on zooming the model up on screen, it continues to look smooth (with appropriate k).

