#include <stdexcept>
#include <string>

#include "Grid.h"

#include "vtkOutput.hpp"

using namespace std;

int main(int argc, char** argv) {


  // parse input parameters
  if (argc < 4) {
    cerr << "Not enough input parameters!" << endl;
    return -1;
  }

  int x_size;
  int y_size;
  float spacing;

  try {
    x_size = stoi(argv[1]);
    y_size = stoi(argv[2]);
    spacing = stof(argv[3]); 
  }
  catch (invalid_argument& e){
    
    cerr << "An exception occured while parsing arguments: " << e.what() << endl;
    return -1;
  }

  Grid grid = Grid(x_size, y_size, spacing, BC::periodic);

  if (argc > 4) {
    
    string surface_class = string(argv[4]);
      
    if (surface_class.compare("Rectangle") == 0) {
        
      if (argc != 9) {
        cerr << "Invalid number of input parameters!" << endl;
        return -1;
      }        
        
      float x_min = stof(argv[5]);
      float x_max = stof(argv[6]);
      float y_min = stof(argv[7]);
      float y_max = stof(argv[8]);
      
      Rectangle rec;
      rec.x_min = x_min; rec.x_max = x_max; rec.y_min = y_min; rec.y_max = y_max;
          
      cout << "Adding Rectangle surface: calculating distance function" << endl; 
      grid.calculateDistancesToRectangle(rec);
      cout << "Finished" << endl;
        
    }
    else if (surface_class.compare("Sphere") == 0) {
              
      if (argc != 8) {
        cerr << "Invalid number of input parameters!" << endl;
        return -1;
      }        
        
      float center_x = stof(argv[5]);
      float center_y = stof(argv[6]);
      float radius = stof(argv[7]);
      
      Sphere sphere;
      sphere.center_x = center_x; sphere.center_y = center_y; sphere.radius = radius;
          
      cout << "Adding Rectangle surface: calculating distance function" << endl; 
      grid.calculateDistancesToSphere(sphere);
      cout << "Finished" << endl;
    }
    else {
      cerr << "Encountered invalid surface class: " << surface_class << endl;
    }
  }
  
  std::vector<std::vector<double> > distances = grid.getDistances();
  
  vtkOutput output("test", spacing, distances);
  output.write();

  return 1;

}
