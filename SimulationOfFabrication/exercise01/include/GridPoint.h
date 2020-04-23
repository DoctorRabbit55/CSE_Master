#pragma once

#include <stdint.h>
#include <vector>

#include "geometry.h"

#include <cmath>
#include <algorithm>

#include <iostream>

enum Direction{x, y};
enum BC {periodic, reflective};
enum Derivative {forwards, backwards, central};


struct Vector2d {
  double x;
  double y;
};



class GridPoint {

  public:
    GridPoint(float x, float y);
    GridPoint(){};
    
    void calculateDistanceToRectangle(Rectangle rec, uint32_t grid_size_x, uint32_t grid_size_y, float spacing, BC bc);
    void calculateDistanceToSphere(Sphere sphere, uint32_t grid_size_x, uint32_t grid_size_y, float spacing, BC bc);
    double getDistance();


    float x_ = 0;
    float y_ = 0;
    
    double distance_ = 0;
    
};
