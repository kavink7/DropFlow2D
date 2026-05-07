#ifndef _TIMEADV_H_
#define _TIMEADV_H_

#include "field.h"

void rk4_adv(Domain &d, Time &t, SimulationParameters &sp, PhysicalParameters &pp);
void rk4_adv_lv(Domain &d, Time &t, SimulationParameters &sp, PhysicalParameters &pp);

#endif  // _TIMEADV_H_
