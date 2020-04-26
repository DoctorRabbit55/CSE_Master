#pragma once

#include <memory>
#include <math.h>

#include "Grid.h"

class Advancer {

  public:
    Advancer(std::shared_ptr<Grid> grid_ptr);

    void advanceByConstant(double v, bool do_engquist_osher);
    void advanceByVector(Vector2d vec, bool do_engquist_osher);
    void advanceByCurvature(bool do_engquist_osher);

  private:
    std::shared_ptr<Grid> grid_ptr_;
    std::vector<std::vector<GridPoint> > temp_grid_;

    double engquistOsherStep(int x, int y, double v);

    const Vector2d n1 = Vector2d(-1, 0);
    const Vector2d n2 = Vector2d(1, 0);
    const Vector2d n3 = Vector2d(-1/sqrt(2), 1/sqrt(2));
    const Vector2d n4 = Vector2d(1/sqrt(2), -1/sqrt(2));
    const Vector2d n5 = Vector2d(1/sqrt(2), 1/sqrt(2));
    const Vector2d n6 = Vector2d(-1/sqrt(2), -1/sqrt(2));
    const Vector2d n7 = Vector2d(0, -1);
    const Vector2d n8 = Vector2d(0, 1);
   
};
