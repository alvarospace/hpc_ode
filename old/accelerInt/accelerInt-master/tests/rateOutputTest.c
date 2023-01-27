/* \file rateOutputTest.c
 *
 * \author Nicholas J. Curtis
 * \date 02/04/2015
This is a simple program that outputs the rates and analytical jacobian of the current mechanism file (mechanism.c)
in order to easily compare between different versions of the rates and jacobian generators
*/

#include <stdio.h>
#include "header.h"

int main (int argc, char *argv[]) {
	write_jacobian_and_rates_output(16384);
}