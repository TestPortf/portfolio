#include "init.h"

#define GREEN_ON() GPIOD->ODR=0x1000
#define ORANGE_ON() GPIOD->ODR=0x2000
#define RED_ON() GPIOD->ODR=0x4000
#define BLUE_ON() GPIOD->ODR=0x8000
#define LEDS_OFF() GPIOD->ODR=0

#ifndef MAIN_H
#define MAIN_H
//
#endif
