#include "Advancer.h"

Advancer::Advancer(std::shared_ptr<Grid> grid_ptr) {

  grid_ptr_ = grid_ptr;

  temp_grid_ = std::vector<std::vector<GridPoint> >(grid_ptr_->size_x_);
  
  for (int x = 0; x < grid_ptr_->size_x_; x++) {
    temp_grid_[x] = std::vector<GridPoint>(grid_ptr_->size_y_);
  }

}

void Advancer::advanceByVector(Vector2d vec, bool do_engquist_osher) {

  Vector2d normal_vec;
  double v = 0;

  // sphere
  if (grid_ptr_->has_rectangle_surface_ == false) {
  //if (1) {

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
  }
  
  // rectangle
  else {
  
    double x_min = grid_ptr_->rec_.x_min;
    double x_max = grid_ptr_->rec_.x_max;
    double y_min = grid_ptr_->rec_.y_min;
    double y_max = grid_ptr_->rec_.y_max;
    
    double center_x = x_min + (x_max-x_min)/2.;
    double center_y = y_min + (y_max-y_min)/2.;
    
    double offset_x = 0;
    double offset_y = 0;
    double norm = 0;
    
    Vector2d normal_vec;
  
    for (int x = 0; x < grid_ptr_->size_x_; x++) {
      for (int y = 0; y < grid_ptr_->size_y_; y++) {
          //std::cout << x << "    " << y << std::endl;      
        if (grid_ptr_->grid_points_[x][y].distance_ != 0) {

          // left 
          if ( x <= x_min) {
            //top
            if (y > y_max) {
            
              offset_x = x - x_min;
              offset_y = y - y_max;
              
              norm = sqrt(offset_x*offset_x + offset_y*offset_y);
              
              normal_vec.x = offset_x / norm;
              normal_vec.y = offset_y / norm;
              
            }
            //middle
            else if (y >= y_min) {
              normal_vec = n1;
            }
            //bottom
            else {
            
              offset_x = x - x_min;
              offset_y = y_min - y;
              
              norm = sqrt(offset_x*offset_x + offset_y*offset_y);
              
              normal_vec.x = offset_x / norm;
              normal_vec.y = offset_y / norm;
            }
          }
          // middle
           else if ( x <= x_max) {
             //top
             if (y > y_max) 
               normal_vec = n8;
             //middle
             else if (y >= y_min) {
               if ( abs(x - center_x) > abs(y - center_y) ) {

                 if (x <= center_x)
                   normal_vec = n1;
                 else
                   normal_vec = n2;
               }
               else if ( abs(x - center_x) < abs(y - center_y) ){
                 if (y <= center_y)
                   normal_vec = n7;
                 else
                   normal_vec = n8;
               }
               else {
               
                 if (x <= center_x)
                   normal_vec.x = -1/sqrt(2);
                 else
                   normal_vec.x = 1/sqrt(2);
                 
                 if (y <= center_y)
                   normal_vec.y = -1/sqrt(2);
                 else
                   normal_vec.y = 1/sqrt(2);
               
                 if (x == center_x && y == center_y) {
                   normal_vec.x = -1/2;
                   normal_vec.y = -1/2;
                 }
               }
             //bottom
             }
             else
               normal_vec = n7; 
           }
           // right
           else {
             //top
             if (y > y_max) { 
               
               offset_x  = x - x_max;
               offset_y = y - y_max;
              
               norm = sqrt(offset_x*offset_x + offset_y*offset_y);
              
               normal_vec.x = offset_x / norm;
               normal_vec.y = offset_y / norm;
             }
             //middle
             else if (y >= y_min)
               normal_vec = n2;
             //bottom
             else {
               offset_x  = x - x_max;
               offset_y = y_min - y;
              
               norm = sqrt(offset_x*offset_x + offset_y*offset_y);
              
               normal_vec.x = offset_x / norm;
               normal_vec.y = offset_y / norm;
               
             }
           }
        }
        else {
          normal_vec = grid_ptr_->getNormalVector(x, y);
        }

        v = vec.x * normal_vec.x + vec.y * normal_vec.y; 

        //std::cout << normal_vec.x << "   " << normal_vec.y <<  "      " << v <<std::endl;    
        if (do_engquist_osher)
          temp_grid_[x][y].distance_ = grid_ptr_->grid_points_[x][y].distance_ - v * engquistOsherStep(x, y, v);
        else
          temp_grid_[x][y].distance_ = grid_ptr_->grid_points_[x][y].distance_ - v;
     
      }
    }   
  }
  
  grid_ptr_->grid_points_.swap(temp_grid_);
  
  Rectangle rec_new;
  rec_new.x_min = grid_ptr_->rec_.x_min + vec.x;
  rec_new.y_min = grid_ptr_->rec_.y_min + vec.y;
  rec_new.x_max = grid_ptr_->rec_.x_max + vec.x;
  rec_new.y_max = grid_ptr_->rec_.y_max + vec.x;
  
  grid_ptr_->rec_ = rec_new;
  
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
