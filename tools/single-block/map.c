#include <stdlib.h>
#include <math.h>

#include "map.h"


const char *map_name[NUM_MAPS] = {"Identity", "Stretched", "Second",
                                  "Divergence Free", "Curl Free", "Rotation", "Annulus","Half-Annulus","Distored"};
const char *map_key[NUM_MAPS] = {"identity", "stretched", "second", "curl", "div", "rot", "annulus","hannulus","distort"};


double map_identity_x(double r, double s, double t)
{
	return r;
}


double map_identity_y(double r, double s, double t)
{
	return s;
}


double map_identity_z(double r, double s, double t)
{
	return t;
}


double map_polar_x(double r, double s, double t)
{
	return (r+.5) * cos(2*M_PI*s);
}


double map_polar_y(double r, double s, double t)
{
	return (r+.5) * sin(2*M_PI*s);
}


double map_polar_z(double r, double s, double t)
{
	return t;
}


double map_second_x(double r, double s, double t)
{
	//return xi;//*xi;
	return sin(.5 * M_PI * r);
	//return 3*xi - 1;
}


double map_second_y(double r, double s, double t)
{
	//return sqrt(eta);
	//return eta;//*eta;
	//return sin(.5 * M_PI * eta);
	return sin(.5 * M_PI * s);
	//return 3*eta-1;
}


double map_second_z(double r, double s, double t)
{
	return sin(.5*M_PI*t);
}


double map_curl_x(double r, double s, double t)
{
	return -s + 1;
}


double map_curl_y(double r, double s, double t)
{
	return r;
}


double map_curl_y_3d(double r, double s, double t)
{
	return -r + 1;
}


double map_curl_z_3d(double r, double s, double t)
{
	return -t + 1;
}


double map_div_x(double r, double s, double t)
{
	return 3*r - 1;
}


double map_div_y(double r, double s, double t)
{
	return 3*s - 1;
}


double map_rot_x(double r, double s, double t)
{
	return r - s;
}


double map_rot_y(double r, double s, double t)
{
	return r + s - 1;
}


double map_hannulus_x(double r, double s, double t)
{
	double r1 = .1;
	double r2 = 1;
	return (r1 + (r2 - r1)*r)*sin(M_PI*s);
}



double map_hannulus_y(double r, double s, double t)
{
	double r1 = .1;
	double r2 = 1;
	return (r1 + (r2 - r1)*r)*cos(M_PI*s);
}


double map_distort_x(double r, double s, double t)
{
	double eps = .05;
	return r + eps*sin(2*M_PI*r)*sin(2*M_PI*s);
}


double map_distort_y(double r, double s, double t)
{
	double eps = .05;
	return s + eps*sin(2*M_PI*r)*sin(2*M_PI*s);
}


double map_distort_x_3d(double r, double s, double t)
{
	double eps = .01;
	return r + eps*sin(2*M_PI*r)*sin(2*M_PI*s)*sin(2*M_PI*t);
}


double map_distort_y_3d(double r, double s, double t)
{
	double eps = .01;
	return s + eps*sin(2*M_PI*r)*sin(2*M_PI*s)*sin(2*M_PI*t);
}


double map_distort_z_3d(double r, double s, double t)
{
	double eps = .01;
	return t + eps*sin(2*M_PI*r)*sin(2*M_PI*s)*sin(2*M_PI*t);
}



mapping *map_create(map_id mid, int nd)
{
	mapping *ret = (mapping*) malloc(sizeof(mapping));

	ret->id = mid;
	switch (mid) {
	case(MAP_IDENTITY):
		ret->x = &map_identity_x;
		ret->y = &map_identity_y;
		ret->z = &map_identity_z;
		break;
	case(MAP_POLAR):
		ret->x = &map_polar_x;
		ret->y = &map_polar_y;
		ret->z = &map_polar_z;
		break;
	case(MAP_SECOND):
		ret->x = &map_second_x;
		ret->y = &map_second_y;
		ret->z = &map_second_z;
		break;
	case(MAP_CURL):
		ret->x = &map_curl_x;
		ret->y = &map_curl_y;
		break;
	case(MAP_DIV):
		ret->x = &map_div_x;
		ret->y = &map_div_y;
		break;
	case(MAP_ROT):
		ret->x = &map_rot_x;
		ret->y = &map_rot_y;
		break;
	case (MAP_HANNULUS):
		ret->x = &map_hannulus_x;
		ret->y = &map_hannulus_y;
		break;
	case(MAP_DISTORT):
		if (nd == 2) {
			ret->x = &map_distort_x;
			ret->y = &map_distort_y;
		} else {
			ret->x = &map_distort_x_3d;
			ret->y = &map_distort_y_3d;
			ret->z = &map_distort_z_3d;
		}
		break;
	default:
		ret->x = &map_identity_x;
		ret->y = &map_identity_y;
		break;
	}

	return ret;
}


void map_destroy(mapping *mp)
{
}
