#include <stdexcept>
#include <string>

#include "Grid.h"

#include "vtkOutput.hpp"

using namespace std;

int main(int argc, char** argv) {

  int arg_counter = 0;

  // parse input parameters
  if (argc < 4) {
    cerr << "Not enough input parameters!" << endl;
    return -1;
  }

  int x_size;
  int y_size;
  float spacing;

  try {
    x_size = stoi(argv[++arg_counter]);
    y_size = stoi(argv[++arg_counter]);
    spacing = stof(argv[++arg_counter]);
  }
  catch (invalid_argument& e){
    
    cerr << "An exception occured while parsing arguments: " << e.what() << endl;
    return -1;
  }

  Grid* grid = new Grid(x_size, y_size, spacing, BC::periodic);

  if (argc > 4) {
    
    string surface_class = string(argv[++arg_counter]);
      
    if (surface_class.compare("Rectangle") == 0) {
        
      float x_min = stof(argv[++arg_counter]);
      float x_max = stof(argv[++arg_counter]);
      float y_min = stof(argv[++arg_counter]);
      float y_max = stof(argv[++arg_counter]);
      
      Rectangle rec;
      rec.x_min = x_min; rec.x_max = x_max; rec.y_min = y_min; rec.y_max = y_max;
          
      cout << "Adding Rectangle surface: calculating distance function" << endl; 
      grid->calculateDistancesToRectangle(rec);
      cout << "Finished" << endl;
        
    }
    else if (surface_class.compare("Sphere") == 0) {        
        
      float center_x = stof(argv[++arg_counter]);
      float center_y = stof(argv[++arg_counter]);
      float radius = stof(argv[++arg_counter]);
      
      Sphere sphere;
      sphere.center_x = center_x; sphere.center_y = center_y; sphere.radius = radius;
          
      cout << "Adding Rectangle surface: calculating distance function" << endl; 
      grid->calculateDistancesToSphere(sphere);
      cout << "Finished" << endl;
    }
    else {
      cerr << "Encountered invalid surface class: " << surface_class << endl;
    }
  }
  
  std::vector<std::vector<double> > distances = grid->getDistances();
  
  vtkOutput output("test", spacing, distances);
  output.write();

  if (argc - arg_counter == 3) {
  
    int x = std::stoi(argv[++arg_counter]);
    int y = std::stoi(argv[++arg_counter]);
  
    double deriv = grid->getDerivative(x, y, Direction::y, Derivative::central);
    std::cout << "derivative: " << deriv << std::endl;
    
    Vector2d vec = grid->getNormalVector(x, y);
    std::cout << "normal vector: " << vec.x << " " << vec.y << std::endl;
  
    double curv = grid->getCurvature(x, y);
    std::cout << "curvature: " << curv << std::endl;
  }


  free(grid);

  return 1;

}
