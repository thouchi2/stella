#ifndef MAP_H
#define MAP_H

#include "base.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NUM_MAPS 9

typedef enum {
	MAP_IDENTITY=0, MAP_STRETCH=1, MAP_SECOND=2,
	MAP_CURL=3, MAP_DIV=4, MAP_ROT=5,
	MAP_POLAR=6, MAP_HANNULUS=7,MAP_DISTORT=8
} map_id;

extern const char *map_name[NUM_MAPS];
extern const char *map_key[NUM_MAPS];

typedef struct {
	func x;
	func y;
	func z;
	int id;
} mapping;


mapping *map_create(map_id, int nd);


void map_destroy(mapping*);

#ifdef __cplusplus
}
#endif

#endif
