#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "field.h"

void apply_bc(Field &f, const Domain &d);
void apply_bc_periodic(Field &f, const Domain &d);
#endif  // _BOUNDARY_H_
