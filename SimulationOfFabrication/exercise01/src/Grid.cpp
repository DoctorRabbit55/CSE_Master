#include "Grid.h"

Grid::Grid(int size_x, int size_y, float spacing, BC bc){

  spacing_ = spacing;
  bc_ = bc;
  size_x_ = size_x;
  size_y_ = size_y;
 
  grid_points_ = std::vector<std::vector<GridPoint> >(size_x); 

  for (int i = 0; i < size_x; i++) {

    grid_points_[i] = std::vector<GridPoint>(size_y);  
    
    for (int j = 0; j < size_y; j++) {
      grid_points_[i][j] = GridPoint(i*spacing, j*spacing);

    }
  }
}

/**
------------------------------------------------------------------------------------------------------------------
**/

void Grid::calculateDistancesToRectangle(Rectangle rec) {

  for (size_t x = 0; x < size_x_; x++) {  
    for (size_t y = 0; y < size_y_; y++) {
      grid_points_[x][y].calculateDistanceToRectangle(rec, size_x_, size_y_, spacing_, bc_);
    }
  }
}

void Grid::calculateDistancesToSphere(Sphere sphere) {

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

  std::cout << "D_x: " << D_x << std::endl;
  std::cout << "D_y: " << D_y << std::endl;

  return D_x + D_y;

}

/**
------------------------------------------------------------------------------------------------------------------
**/

Vector2d Grid::getNormalVector(int x, int y) {

  Vector2d vec;

  vec.x = getDerivative(x, y, Direction::x, Derivative::central);
  vec.y = getDerivative(x, y, Direction::y, Derivative::central); 

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
      
    if (x == size_y_)
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
