#pragma once

#include <iostream>
#include <vector>

#include "geometry.h"
#include "GridPoint.h"


class Grid {

  public:
    Grid(int size_x, int size_y, float spacing, BC bc);    
    
    std::vector<std::vector<GridPoint> > grid_points_;
    
    void calculateDistancesToRectangle(Rectangle rec);
    void calculateDistancesToSphere(Sphere sphere);
    std::vector<std::vector<double> > getDistances();   
  
    double getDerivative(int x, int y, Direction direction, Derivative typ);
    Vector2d getNormalVector(int x, int y);
    double getCurvature(int x, int y);
  
  
    float spacing_ = 0;
    unsigned int size_x_;
    unsigned int size_y_;
    BC bc_;
    
    int checkXandUpdate(int x);
    int checkYandUpdate(int y);
    
};
