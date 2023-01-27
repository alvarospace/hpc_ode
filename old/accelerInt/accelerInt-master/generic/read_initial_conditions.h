/**
 * \file
 * \brief definition of the generic initial condition reader
 *
 * \author Nicholas Curtis
 * \date 03/10/2015
 *
 */

#ifndef READ_INITIAL_CONDITIONS_H
#define READ_INITIAL_CONDITIONS_H
void read_initial_conditions(const char* filename, int NUM, double** y_host, double** variable_host);
#endif