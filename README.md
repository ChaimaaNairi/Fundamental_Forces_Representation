﻿# Fundamental_Forces_Representation
The Fundamental_Forces_Representation project is a C++ and OpenGL application that aims to represent fundamental forces, such as electric or gravitational fields, as three-dimensional vector fields by giving different parameters. The project implements a vector field representation of fundamental forces using numerical methods and allows users to manipulate the camera view during the simulation phase. </br>

 <h2>Project Scope </h2>
The project has the following features: </br>

<h3>Interaction</h3>
Users can select the force to be displayed as input to the application. For the electric field, users can select the charges, while for the gravitational field, users can select the masses and determine their positions. During the simulation phase, users can manipulate the camera view using their mouse. </br>

<h3>Animation</h3>
Force interactions are time-varying, and the application detects these changes at certain time intervals using numerical methods. The animation system then plans the frames and sends this frame information to the drawing algorithm. If force interaction calculations are expensive, the animation system does the necessary interpolation to draw more frames. </br>

<h3>Drawing</h3>
The drawing system statically draws vectors and objects, such as charges and masses, in vector fields using OpenGL functions for the given environment information. </br>
