#include "Advancer.h"

Advancer::Advancer(std::shared_ptr<Grid> grid_ptr) {

  grid_ptr_ = grid_ptr;

  temp_grid_ = std::vector<std::vector<GridPoint> >(grid_ptr_->size_x_);
  
  for (int x = 0; x < grid_ptr_->size_x_; x++) {
    temp_grid_[x] = std::vector<GridPoint>(grid_ptr_->size_y_);
  }

}

void Advancer::advanceByConstant(double v, bool do_engquist_osher){

  for (int x = 0; x < grid_ptr_->size_x_; x++) {
    for (int y = 0; y < grid_ptr_->size_y_; y++) {
    
      if (do_engquist_osher)
        temp_grid_[x][y].distance_ = grid_ptr_->grid_points_[x][y].distance_ - v * engquistOsherStep(x, y, v);
      else
        temp_grid_[x][y].distance_ = grid_ptr_->grid_points_[x][y].distance_ - v;
   
    }
  }
  
  grid_ptr_->grid_points_.swap(temp_grid_);
}

void Advancer::advanceByVector(Vector2d vec, bool do_engquist_osher) {

  Vector2d normal_vec;
  double v = 0;

  for (int x = 0; x < grid_ptr_->size_x_; x++) {
    for (int y = 0; y < grid_ptr_->size_y_; y++) {
      
      normal_vec = grid_ptr_->getNormalVector(x, y);

      v = vec.x * normal_vec.x + vec.y * normal_vec.y; 
        //std::cout << normal_vec.x << "   " << normal_vec.y <<  "      " << v <<std::endl;    
      if (do_engquist_osher)
        temp_grid_[x][y].distance_ = grid_ptr_->grid_points_[x][y].distance_ - v * engquistOsherStep(x, y, v);
      else
        temp_grid_[x][y].distance_ = grid_ptr_->grid_points_[x][y].distance_ - v;
     
    }
  }   
  
  grid_ptr_->grid_points_.swap(temp_grid_);
  
}

void Advancer::advanceByCurvature(bool do_engquist_osher) {

  double v;

  for (int x = 0; x < grid_ptr_->size_x_; x++) {
    for (int y = 0; y < grid_ptr_->size_y_; y++) {
    
      v = grid_ptr_->getCurvature(x, y);
    
      std::cout << v << std::endl;
    
      if (do_engquist_osher)
        temp_grid_[x][y].distance_ = grid_ptr_->grid_points_[x][y].distance_ - v * engquistOsherStep(x, y, v);
      else
        temp_grid_[x][y].distance_ = grid_ptr_->grid_points_[x][y].distance_ - v;
   
    }
  }
  
  grid_ptr_->grid_points_.swap(temp_grid_);

}


double Advancer::engquistOsherStep(int x, int y, double v) {

  double result = 0;
    
  double D_x_b = grid_ptr_->getDerivative(x,y,Direction::x,Derivative::backwards);
  double D_x_f = grid_ptr_->getDerivative(x,y,Direction::x,Derivative::forwards);
  double D_y_b = grid_ptr_->getDerivative(x,y,Direction::y,Derivative::backwards);
  double D_y_f = grid_ptr_->getDerivative(x,y,Direction::y,Derivative::forwards);
    
  //std::cout << "D_x_b: " << D_x_b << " D_x_f: " << D_x_f << " D_y_b: " << D_y_b << " D_y_f: " << D_y_f << std::endl;
  
  if (v < 0) {
    
    result = sqrt( pow(std::max(-D_x_b, 0.), 2) + 
                   pow(std::max(-D_y_b, 0.), 2) +
                   pow(std::min(-D_x_f, 0.), 2) +
                   pow(std::min(-D_y_f, 0.), 2)  );
    
  }
  else {
  
    result = sqrt( pow(std::max(D_x_b, 0.), 2) + 
                   pow(std::max(D_y_b, 0.), 2) +
                   pow(std::min(D_x_f, 0.), 2) +
                   pow(std::min(D_y_f, 0.), 2)  );
  
  }
  
  //std::cout << result << std::endl;
  //return result;
  return result > 0 ? result : 1;

}
