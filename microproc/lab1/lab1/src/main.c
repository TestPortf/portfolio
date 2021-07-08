#include "main.h"

void Delay(volatile uint32_t nCount)
{
	while(nCount--) {}
}

int main(void){
	LEDs_ini();
	while(1){
		GREEN_ON();
		Delay(800000);
		LEDS_OFF();
		Delay(800000);
		ORANGE_ON();
		Delay(800000);
		LEDS_OFF();
		Delay(800000);
		RED_ON();
		Delay(800000);
		LEDS_OFF();
		Delay(800000);
		BLUE_ON();
		Delay(800000);
		LEDS_OFF();
		Delay(800000);
	}
}
