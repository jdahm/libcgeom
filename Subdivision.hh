#ifndef SUBDIVISION_HH
#define SUBDIVISION_HH

#include <string>

class Subdivision {
public:
  virtual void Write(const std::string& fileName) = 0;
  virtual ~Subdivision() { };
};

#endif
