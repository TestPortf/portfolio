#include "main.h"

void Delay(volatile uint32_t nCount)
{
	while(nCount--) {}
}

int main(void){
	LEDs_ini();
	while(1){
		CHANGE_GREEN();
		Delay(800000);
		CHANGE_GREEN();
		Delay(800000);
		CHANGE_ORANGE();
		Delay(800000);
		CHANGE_ORANGE();
		Delay(800000);
		CHANGE_RED();
		Delay(800000);
		CHANGE_RED();
		Delay(800000);
		CHANGE_BLUE();
		Delay(800000);
		CHANGE_BLUE();
		Delay(800000);
	}
}
