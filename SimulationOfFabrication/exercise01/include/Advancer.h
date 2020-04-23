#pragma once

#include <memory>
#include <math.h>

#include "Grid.h"

class Advancer {

  public:
    Advancer(std::shared_ptr<Grid> grid_ptr);

    void advanceByConstant(double v, bool do_engquist_osher);
    void advanceByVector(Vector2d vec, bool do_engquist_osher);

  private:
    std::shared_ptr<Grid> grid_ptr_;
    std::vector<std::vector<GridPoint> > temp_grid_;

    double engquistOsherStep(int x, int y, double v);

};
