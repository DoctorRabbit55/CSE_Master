#include "Grid.h"

Grid::Grid(int size_x, int size_y, float spacing, BC bc){

  spacing_ = spacing;
  bc_ = bc;
  size_x_ = (int) size_x/spacing;
  size_y_ = (int) size_y/spacing;
 
  grid_points_ = std::vector<std::vector<GridPoint> >(size_x_); 

  for (int i = 0; i < size_x_; i++) {

    grid_points_[i] = std::vector<GridPoint>(size_y_);  
    
    for (int j = 0; j < size_y_; j++) {
      grid_points_[i][j] = GridPoint(i*spacing, j*spacing);

    }
  }
}

/**
------------------------------------------------------------------------------------------------------------------
**/

void Grid::calculateDistancesToRectangle(Rectangle rec) {

  // convert to grid coordinates
  rec.x_min = (int) rec.x_min / spacing_;
  rec.x_max = (int) rec.x_max / spacing_;
  rec.y_min = (int) rec.y_min / spacing_;
  rec.y_max = (int) rec.y_max / spacing_;

  for (size_t x = 0; x < size_x_; x++) {  
    for (size_t y = 0; y < size_y_; y++) {
      grid_points_[x][y].calculateDistanceToRectangle(rec, size_x_, size_y_, spacing_, bc_);
    }
  }
}

void Grid::calculateDistancesToSphere(Sphere sphere) {

  // convert to grid coordinates
  sphere.center_x = (int) sphere.center_x / spacing_;
  sphere.center_y = (int) sphere.center_y / spacing_;
  sphere.radius /= spacing_;


  for (size_t x = 0; x < size_x_; x++) {  
    for (size_t y = 0; y < size_y_; y++) {
      grid_points_[x][y].calculateDistanceToSphere(sphere, size_x_, size_y_, spacing_, bc_);
    }
  }

}

/**
------------------------------------------------------------------------------------------------------------------
**/

std::vector<std::vector<double> > Grid::getDistances() {

  std::vector<std::vector<double> > distances = std::vector<std::vector<double> >(size_x_);

  for (size_t x = 0; x < size_x_; x++) {  
    
    distances[x] = std::vector<double>(size_y_);
    
    for (size_t y = 0; y < size_y_; y++) {
      distances[x][y] = grid_points_[x][y].getDistance();
    }
  }
  
  return distances;
}

/**
------------------------------------------------------------------------------------------------------------------
**/

double Grid::getCurvature(int x, int y) {


  double D_x = (getNormalVector(x+1, y).x - getNormalVector(x-1, y).x) / (2*spacing_);
  double D_y = (getNormalVector(x, y+1).y - getNormalVector(x, y-1).y) / (2*spacing_);

  return D_x + D_y;

}

/**
------------------------------------------------------------------------------------------------------------------
**/

Vector2d Grid::getNormalVector(int x, int y) {

  Vector2d vec;

  vec.x = getDerivative(x, y, Direction::x, Derivative::central);
  vec.y = getDerivative(x, y, Direction::y, Derivative::central); 

  if (vec.x == 0 && vec.y == 0)
    return vec;
    
  double normalization = sqrt( vec.x*vec.x + vec.y*vec.y );

  vec.x /= normalization;
  vec.y /= normalization;

  return vec;
}

/**
------------------------------------------------------------------------------------------------------------------
**/

double Grid::getDerivative(int x, int y, Direction direction, Derivative type) {

  switch(type) {
  
    case Derivative::forwards:
      
      if (direction == Direction::x)
        return (grid_points_[checkXandUpdate(x+1)][y].getDistance() - grid_points_[x][y].getDistance()) / spacing_;
      
      if (direction == Direction::y)
        return (grid_points_[x][checkYandUpdate(y+1)].getDistance() - grid_points_[x][y].getDistance()) / spacing_;
      
      break;

    case Derivative::backwards:
      
      if (direction == Direction::x)
        return (grid_points_[x][y].getDistance() - grid_points_[checkXandUpdate(x-1)][y].getDistance()) / spacing_;
      
      if (direction == Direction::y)
        return (grid_points_[x][y].getDistance() - grid_points_[x][checkYandUpdate(y-1)].getDistance()) / spacing_;
      
      break;
  
    case Derivative::central:
  
      return (getDerivative(x, y, direction, Derivative::forwards) + getDerivative(x, y, direction, Derivative::backwards)) / 2;    
      break;
  }
  
  return 0;
}

int Grid::checkXandUpdate(int x) {

  if (bc_ == BC::periodic) {
  
    if (x == -1)
      return size_x_-1;
      
    if (x == size_x_)
      return 0;
      
  }

  if (bc_ == BC::reflective) {
  
    if (x == -1)
      return 0;
      
    if (x == size_x_)
      return size_x_-1;
      
  }
  
  return x;
}

int Grid::checkYandUpdate(int y) {

  if (bc_ == BC::periodic) {
  
    if (y == -1)
      return size_y_-1;
      
    if (y == size_y_)
      return 0;
      
  }

  if (bc_ == BC::reflective) {
  
    if (y == -1)
      return 0;
      
    if (y == size_y_)
      return size_y_-1;
      
  }
  
  return y;
}
