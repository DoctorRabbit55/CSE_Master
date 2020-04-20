#pragma once

#include <iostream>
#include <vector>

#include "geometry.h"
#include "GridPoint.h"

#define FORWARD_DER           0
#define BACKWARD_DER          1
#define CENTRAL_DER           2

#define X_DIR                 0
#define Y_DIR                 1

#define PERIODIC_BC           1
#define REFLECTIVE_BC         0

class Grid {

  public:
    Grid(int size_x, int size_y, float spacing, bool BC);    
    
    void calculateDistancesToRectangle(Rectangle rec);
    void calculateDistancesToSphere(Sphere sphere);
    std::vector<std::vector<double> > getDistances();   
  
    double getDerivative(int x, int y, bool direction, uint8_t kind);
  
  private:  
    std::vector<std::vector<GridPoint> > grid_points_;
    float spacing_ = 0;
    bool bc_ = 0;
};
