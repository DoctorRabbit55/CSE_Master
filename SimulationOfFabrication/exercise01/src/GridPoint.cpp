#include "GridPoint.h"


GridPoint::GridPoint(float x, float y) {

  x_ = x;
  y_ = y;
  distance_ = 0;

}

void GridPoint::calculateDistanceToRectangle(Rectangle rec, uint32_t grid_size_x, uint32_t grid_size_y, float spacing, BC bc) {

  double delta_x;
  double delta_y;
  
  // point is inside of rectangle
  if (x_ >= rec.x_min && x_ <= rec.x_max && y_ >= rec.y_min && y_ <= rec.y_max) {

    delta_x = std::min( abs(rec.x_min - x_), abs(rec.x_max - x_) );
    delta_y = std::min( abs(rec.y_min - y_), abs(rec.y_max - y_) );
    
    distance_ = - std::min(delta_x, delta_y);   
  
  }
  // point is outside of rectangle
  else {
  
    delta_x = (x_ > rec.x_min && x_ < rec.x_max) ? 0 : std::min( abs(rec.x_min - x_), abs(rec.x_max - x_) );    
    delta_y = (y_ > rec.y_min && y_ < rec.y_max) ? 0 : std::min( abs(rec.y_min - y_), abs(rec.y_max - y_) );
  
    distance_ = sqrt( pow(delta_x, 2) + pow(delta_y, 2) );
  }
  
  if (bc == BC::periodic) {
  
    // lists of other rectangles, we have to use
    float per_x_min[] = { rec.x_min - grid_size_x, rec.x_min, rec.x_min + grid_size_x };
    float per_x_max[] = { rec.x_max - grid_size_x, rec.x_max, rec.x_max + grid_size_x };  
    float per_y_min[] = { rec.y_min - grid_size_y, rec.y_min, rec.y_min + grid_size_y };  
    float per_y_max[] = { rec.y_max - grid_size_y, rec.y_max, rec.y_max + grid_size_y };
    
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        
        // we already did middle one
        if (i == 1 && j == 1)
          continue;
      
        rec.x_min = per_x_min[i];
        rec.x_max = per_x_max[i];
        rec.y_min = per_y_min[j];
        rec.y_max = per_y_max[j];
      
        // do same thing as above
      
        if (x_ > rec.x_min && x_ < rec.x_max && y_ > rec.y_min && y_ < rec.y_max) {

           delta_x = std::min( abs(rec.x_min - x_), abs(rec.x_max - x_) );
           delta_y = std::min( abs(rec.y_min - y_), abs(rec.y_max - y_) );
    
        distance_ = std::min(distance_, -std::min(delta_x, delta_y));   
  
        }
        else {
  
          delta_x = (x_ > rec.x_min && x_ < rec.x_max) ? 0 : std::min( abs(rec.x_min - x_), abs(rec.x_max - x_) );    
          delta_y = (y_ > rec.y_min && y_ < rec.y_max) ? 0 : std::min( abs(rec.y_min - y_), abs(rec.y_max - y_) );
  
         distance_ = std::min(distance_, sqrt( pow(delta_x, 2) + pow(delta_y, 2) ));
        }
      }
    }
  }
}



void GridPoint::calculateDistanceToSphere(Sphere sphere, uint32_t grid_size_x, uint32_t grid_size_y, float spacing, BC bc) {
  
  
  double radius = sqrt( pow(x_ - sphere.center_x, 2) + pow(y_ - sphere.center_y, 2)); 

  distance_ = radius - sphere.radius;

  if (bc == BC::periodic) {
  
    float per_center_x[] = { sphere.center_x - grid_size_x + 1, sphere.center_x, sphere.center_x + grid_size_x - 1 };
    float per_center_y[] = { sphere.center_y - grid_size_y + 1, sphere.center_y, sphere.center_y + grid_size_y - 1 };  

    
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        
        // we already did middle one
        if (i == 1 && j == 1)
          continue;
      
        sphere.center_x = per_center_x[i];
        sphere.center_y = per_center_y[j];
      
        double radius = sqrt( pow(x_ - sphere.center_x, 2) + pow(y_ - sphere.center_y, 2)); 

        distance_ = std::min(distance_, radius - sphere.radius);
      }
    }
  }
  
  distance_ = distance_ * spacing;
}

double GridPoint::getDistance() {
  return distance_;
}
