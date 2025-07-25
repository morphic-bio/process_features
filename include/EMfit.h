#ifndef EMFIT_H
#define EMFIT_H

#include "common.h"

void determine_signal_cutoff_from_fit(NBSignalCut *fit, int len, double gposterior,
                                double em_cumulative_limit);

NBSignalCut em_nb_signal_cut(const uint32_t *hist, int len, double gposterior,
                 int max_iter, double tol, double em_cumulative_limit);

#endif // EMFIT_H