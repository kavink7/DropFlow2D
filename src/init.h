#ifndef _INIT_H_
#define _INIT_H_

#include "field.h"

void setparams(Domain &d, SimulationParameters &sp, PhysicalParameters &pp);

void allocate_fields(Domain &d, Time &t);

void initialize_data(Domain &d, SimulationParameters &sp);

#endif  // _INIT_H_
