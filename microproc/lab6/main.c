#include "main.h"

int duty_cycle_counter;
int duty_cycles_array [NUMBER_OF_DUTY_CYCLES] = {0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80, 0xFF};


void Delay(volatile uint32_t nCount)
{
	while(nCount--) {}
}

int main(void){
	LEDs_ini();
	TIM4_ini();
	
  while(1)
  {
		for (duty_cycle_counter = 0; duty_cycle_counter < NUMBER_OF_DUTY_CYCLES; duty_cycle_counter++){
			WRITE_REG(TIM4->CCR1, duty_cycles_array[duty_cycle_counter]);
			WRITE_REG(TIM4->CCR2, duty_cycles_array[duty_cycle_counter]);
			WRITE_REG(TIM4->CCR3, duty_cycles_array[duty_cycle_counter]);
			WRITE_REG(TIM4->CCR4, duty_cycles_array[duty_cycle_counter]);
			Delay(0xF0000);
		}
  }
}
