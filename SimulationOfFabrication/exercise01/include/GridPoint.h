#pragma once

#include <stdint.h>
#include <vector>

#include "geometry.h"

#include <cmath>
#include <algorithm>


class GridPoint {

  public:
    GridPoint(float x, float y);
    GridPoint(){};
    
    void calculateDistanceToRectangle(Rectangle rec, uint32_t grid_size_x, uint32_t grid_size_y, bool bc_periodic);
    void calculateDistanceToSphere(Sphere sphere, uint32_t grid_size_x, uint32_t grid_size_y, bool bc_periodic);
    double getDistance();

  private:
    float x_ = 0;
    float y_ = 0;
    
    double distance_ = 0;
    
};
