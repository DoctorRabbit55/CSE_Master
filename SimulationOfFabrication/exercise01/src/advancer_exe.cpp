#include <stdexcept>
#include <string>

#include "Grid.h"
#include "Advancer.h"

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

  cout << "What boundary conditions should be applied? [ref|per]" << endl;;
  
  string bc_string;
  BC bc;
  
  while (true) {
  
    std::getline(cin, bc_string);
 
    if (bc_string.compare("ref") == 0) {
      bc = BC::reflective;
      break;
    }
    
    if (bc_string.compare("per") == 0) {
      bc = BC::periodic;
      break;
    }
 
    cout << "Invalid choice. Please choose [ref|per]." << endl;
 
  }
 
  std::shared_ptr<Grid> grid_ptr(new Grid(x_size, y_size, spacing, bc));

  if (argc > 4) {
    
    string surface_class = string(argv[++arg_counter]);
      
    if (surface_class.compare("Rectangle") == 0) {
        
      float x_min = stof(argv[++arg_counter]);
      float y_min = stof(argv[++arg_counter]);
      float x_max = stof(argv[++arg_counter]);
      float y_max = stof(argv[++arg_counter]);
      
      Rectangle rec;
      rec.x_min = x_min; rec.x_max = x_max; rec.y_min = y_min; rec.y_max = y_max;
          
      cout << "Adding Rectangle surface: calculating distance function" << endl; 
      grid_ptr->calculateDistancesToRectangle(rec);
      cout << "Finished" << endl;
        
    }
    else if (surface_class.compare("Sphere") == 0) {        
        
      float center_x = stof(argv[++arg_counter]);
      float center_y = stof(argv[++arg_counter]);
      float radius = stof(argv[++arg_counter]);
      
      Sphere sphere;
      sphere.center_x = center_x; sphere.center_y = center_y; sphere.radius = radius;
          
      cout << "Adding Sphere surface: calculating distance function" << endl; 
      grid_ptr->calculateDistancesToSphere(sphere);
      cout << "Finished" << endl;
    }
    else {
      cerr << "Encountered invalid surface class: " << surface_class << endl;
    }
  }
  

  if (argc - arg_counter == 3) {
  
    int x = std::stoi(argv[++arg_counter]);
    int y = std::stoi(argv[++arg_counter]);
    
    Vector2d vec = grid_ptr->getNormalVector(x, y);
    
    std::cout << "normal vector: " << vec.x << " " << vec.y << std::endl;
  
    double curv = grid_ptr->getCurvature(x, y);
    std::cout << "curvature: " << curv << std::endl;
  }
  
  std::vector<std::vector<double> > distances = grid_ptr->getDistances();
  
  vtkOutput output("grid.vtk", spacing, distances);
  output.write();
  cout << "File saved as \"grid.vtk\"" << endl;



  std::shared_ptr<Advancer> adv_ptr(new Advancer(grid_ptr));

  cout << "Choose an advancer method. [constant|vector|curvature]" << endl;;
  
  string advance_method;

  while (true) {
  
    getline(cin, advance_method);
 
    if (advance_method.compare("constant") == 0) {
      
      cout << "Choose velocity. [double]" << endl;
      string s;
      getline(cin, s);
      double v = stod(s);
      
      cout << "Choose time. [int]" << endl;
      getline(cin, s);
      double time = stod(s);
      
      cout << "Use enqguist-ohser-sheme? [1|0]" << endl;
      getline(cin, s);
      bool eos = stoi(s);
      
      double maxV = abs(v);
      double delta_t = spacing / maxV * 0.1;
      
      int iter = time / delta_t;
      
      for (int i = 0; i < iter; i++) {
        adv_ptr->advanceByConstant(v, delta_t, eos);
      }
      
      break;
    }
    
    if (advance_method.compare("vector") == 0) {
      
      Vector2d vec;
      
      cout << "Choose x velocity. [double]" << endl;
      string s;
      getline(cin, s);
      vec.x = stod(s);
      
      cout << "Choose y velocity. [double]" << endl;     
      getline(cin, s);
      vec.y = stod(s);
      
      cout << "Choose time. [int]" << endl;
      getline(cin, s);
      double time = stod(s);
      
      cout << "Use enqguist-ohser-sheme? [1|0]" << endl;
      getline(cin, s);
      bool eos = stoi(s);
      
      double maxV = abs(max(vec.x, vec.y));
      double delta_t = spacing / maxV * 0.1;
      
      int iter = time / delta_t;
      
      for (int i = 0; i < iter; i++) {
        adv_ptr->advanceByVector(vec, delta_t, eos);
      }
      
      break;
    }
    
    if (advance_method.compare("curvature") == 0) {
            
      cout << "Choose time. [int]" << endl;
      string s;
      getline(cin, s);
      double time = stod(s);
      
      cout << "Use enqguist-ohser-sheme? [1|0]" << endl;
      getline(cin, s);
      bool eos = stoi(s);
      
      double delta_t = 0.05 * spacing;
      int iter = time / delta_t;
      
      for (int i = 0; i < iter; i++) {
       adv_ptr->advanceByCurvature(delta_t, eos);
      }
      
      break;
    }
 
    cout << "Invalid choice. Please choose [constant|vector|curvature]." << endl;
 
  }
  
  distances = grid_ptr->getDistances();
  
  vtkOutput output2("advanced.vtk", spacing, distances);
  output2.write();
  cout << "File saved as \"advanced.vtk\"" << endl;


  return 1;

}
