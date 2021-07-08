#include "init.h"

GPIO_InitTypeDef Init_LEDs;

void LEDs_ini()
{
	RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOD, ENABLE);
	Init_LEDs.GPIO_Pin = GPIO_Pin_12 | GPIO_Pin_13 | GPIO_Pin_14 | GPIO_Pin_15;
	Init_LEDs.GPIO_Mode = GPIO_Mode_OUT;
	Init_LEDs.GPIO_Speed = GPIO_Speed_2MHz;
	Init_LEDs.GPIO_OType = GPIO_OType_PP;
	Init_LEDs.GPIO_PuPd = GPIO_PuPd_NOPULL;
	GPIO_Init(GPIOD, &Init_LEDs);
}

