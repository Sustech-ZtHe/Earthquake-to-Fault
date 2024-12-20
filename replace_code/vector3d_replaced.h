//
// vector3d.h
//     Class providing common math operations for 3D points
//
// Author:  Tilman Schramke, Christoph Dalitz
// Date:    2017-03-16
// License: see LICENSE
//

#ifndef VECTOR3D_H_
#define VECTOR3D_H_

#include <iostream>

class Vector3d {
public:
  double x;
  double y;
  double z;
  double mag;
  double year;
  double month;
  double day;
  double hour;
  double minute;
  double second;

  Vector3d();
  Vector3d(double a, double b, double c);
  //Vector3d(double a, double b, double c, double mag);
  //修改2
  Vector3d(double a, double b, double c, double mag, double year, double month, double day, double hour, double minute, double second);
  bool operator==(const Vector3d &rhs) const;
  Vector3d& operator=(const Vector3d& other);
  // nicely formatted output
  friend std::ostream& operator<<(std::ostream& os, const Vector3d& vec);
  // Euclidean norm
  double norm() const;
};



// mathematical vector operations

// vector addition
Vector3d operator+(Vector3d x, Vector3d y);
// vector subtraction
Vector3d operator-(Vector3d x, Vector3d y);
// scalar product
double operator*(Vector3d x, Vector3d y);
// scalar multiplication
Vector3d operator*(Vector3d x, double c);
Vector3d operator*(double c, Vector3d x);
Vector3d operator/(Vector3d x, double c);


#endif /* VECTOR3D_H_ */
