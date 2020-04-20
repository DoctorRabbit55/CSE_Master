#include "Grid.h"

Grid::Grid(int size_x, int size_y, float spacing, bool BC){

  spacing_ = spacing;
  bc_ = BC;
 
  grid_points_ = std::vector<std::vector<GridPoint> >(size_x+2); 

  for (int i = 0; i < size_x+2; i++) {

    grid_points_[i] = std::vector<GridPoint>(size_y+2);  
    
    for (int j = 1; j < size_y+1; j++) {
      grid_points_[i][j] = GridPoint(i*spacing, j*spacing);

    }
  }
  
  // fill in BC
  if (BC == PERIODIC_BC) {
  
    for (int i = 1; i < size_y+1; i++) {
      std::cout << i << std::endl;
      grid_points_[0][i] = grid_points_[size_x][i];
      std::cout << i << std::endl;
      grid_points_[size_x+1][i] = grid_points_[1][i];
    }
  
    for (int i = 1; i < size_x+1; i++) {
      grid_points_[i][0] = grid_points_[i][size_y];
      grid_points_[i][size_y+1] = grid_points_[i][size_y-1];
    }
  }

  else if (BC == REFLECTIVE_BC) {
  
    for (int i = 1; i < size_y+1; i++) {
      std::cout << i << std::endl;
      grid_points_[0][i] = grid_points_[size_x][i];
      std::cout << i << std::endl;
      grid_points_[size_x+1][i] = grid_points_[1][i];
    }
  
    for (int i = 1; i < size_x+1; i++) {
      grid_points_[i][0] = grid_points_[i][size_y];
      grid_points_[i][size_y+1] = grid_points_[i][size_y-1];
    }
  }

}

void Grid::calculateDistancesToRectangle(Rectangle rec) {

  for (size_t x = 1; x < grid_points_.size()-1; x++) {  
    for (size_t y = 1; y < grid_points_[0].size()-1; y++) {
      grid_points_[x][y].calculateDistanceToRectangle(rec, grid_points_.size()-2, grid_points_[0].size()-2, bc_);
    }
  }
}

void Grid::calculateDistancesToSphere(Sphere sphere) {

  for (size_t x = 1; x < grid_points_.size()-1; x++) {  
    for (size_t y = 1; y < grid_points_[0].size()-1; y++) {
      grid_points_[x][y].calculateDistanceToSphere(sphere, grid_points_.size()-2, grid_points_[0].size()-2, bc_);
    }
  }

}

std::vector<std::vector<double> > Grid::getDistances() {

  std::vector<std::vector<double> > distances = std::vector<std::vector<double> >(grid_points_.size());

  for (size_t x = 1; x < grid_points_.size()-1; x++) {  
    
    distances[x] = std::vector<double>(grid_points_[x].size());
    
    for (size_t y = 1; y < grid_points_[0].size()-1; y++) {
      distances[x][y] = grid_points_[x][y].getDistance();
    }
  }
  
  return distances;

}

double Grid::getDerivative(int x, int y, bool direction, uint8_t kind) {

  switch(kind) {
  
    case FORWARD_DER:
      
      if (direction == X_DIR)
        return (grid_points_[x+1][y].getDistance() - grid_points_[x][y].getDistance()) / spacing_;
      
      if (direction == Y_DIR)
        return (grid_points_[x][y+1].getDistance() - grid_points_[x][y].getDistance()) / spacing_;
      
      break;

    case BACKWARD_DER:
      
      if (direction == X_DIR)
        return (grid_points_[x][y].getDistance() - grid_points_[x-1][y].getDistance()) / spacing_;
      
      if (direction == Y_DIR)
        return (grid_points_[x][y].getDistance() - grid_points_[x][y-1].getDistance()) / spacing_;
      
      break;
  
    case CENTRAL_DER:
  
      return (getDerivative(x, y, direction, FORWARD_DER) + getDerivative(x, y, direction, BACKWARD_DER)) / 2;    
      break;
  }
  
  return 0;
}


