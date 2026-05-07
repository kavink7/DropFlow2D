#ifndef _VELOCITY_H_
#define _VELOCITY_H_

#include "field.h"

void set_u_staggered(Domain &d, SimulationParameters &sp, Field *u);
void set_v_staggered(Domain &d, SimulationParameters &sp, Field *v);
void set_sbr_u(Domain &d, Field *u, Field *v);
void set_sbr_v(Domain &d, Field *u, Field *v);
void set_leveque_u(Domain &d, Field *u, double t, double T);
void set_leveque_v(Domain &d, Field *v, double t, double T);
#endif  // _VELOCITY_H_
