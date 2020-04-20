#include "GridPoint.h"

#include <iostream>

GridPoint::GridPoint(float x, float y) {

  x_ = x;
  y_ = y;

}

void GridPoint::calculateDistanceToRectangle(Rectangle rec, uint32_t grid_size_x, uint32_t grid_size_y, bool bc_periodic) {

  double delta_x;
  double delta_y;

  if (x_ > rec.x_min && x_ < rec.x_max && y_ > rec.y_min && y_ < rec.y_max) {

    delta_x = std::min( abs(rec.x_min - x_), abs(rec.x_max - x_) );
    delta_y = std::min( abs(rec.y_min - y_), abs(rec.y_max - y_) );
    
    distance_ = - std::min(delta_x, delta_y);   
  
  }
  else {
  
    delta_x = (x_ > rec.x_min && x_ < rec.x_max) ? 0 : std::min( abs(rec.x_min - x_), abs(rec.x_max - x_) );    
    delta_y = (y_ > rec.y_min && y_ < rec.y_max) ? 0 : std::min( abs(rec.y_min - y_), abs(rec.y_max - y_) );
  
    distance_ = sqrt( pow(delta_x, 2) + pow(delta_y, 2) );
  }

  std::cout << "x: " << x_ << " y: " << y_ << " distance: " << distance_ << std::endl;
}

/*
|x x x x|x x x x
|x.x x x|x.x x x
|x x x x|x x x x
|x x x x|x x x x
*/

void GridPoint::calculateDistanceToSphere(Sphere sphere, uint32_t grid_size_x, uint32_t grid_size_y, bool bc_periodic) {

  
  double radius = sqrt( pow(x_ - sphere.center_x, 2) + pow(y_ - sphere.center_y, 2)); 

  std::vector<double> radi = std::vector<double>(5, 1000000);
  radi[0] = radius;
  
  if (bc_periodic) {
  
    if (x_ < grid_size_x/2) {
      radi[1] = sqrt( pow(x_ - sphere.center_x-grid_size_x, 2) + pow(y_ - sphere.center_y, 2));
    }
    if (x_ > grid_size_x/2) {
      radi[2] = sqrt( pow(x_ - sphere.center_x+grid_size_x, 2) + pow(y_ - sphere.center_y, 2));
    }
    if (y_ < grid_size_y/2) {
      radi[3] = sqrt( pow(x_ - sphere.center_x, 2) + pow(y_ - sphere.center_y-grid_size_y, 2));
    }
    if (y_ > grid_size_y/2) {
      radi[4] = sqrt( pow(x_ - sphere.center_x, 2) + pow(y_ - sphere.center_y+grid_size_y, 2));
    }
  
  
  }

  distance_ = *std::min_element(radi.begin(), radi.end()) - sphere.radius;

}

double GridPoint::getDistance() {
  return distance_;
}
