#include "main.h"

int duty_cycles_array [NUMBER_OF_DUTY_CYCLES] = {0x2,0x4,0x8,0x10,0x20,0x40,0x80,0xFF};

void Generate_pwm(int number_of_periods, int pwm_duty_cycle){
	int pwm_period_number;
	int impulse_time;
	int pause_time;
	
	for(pwm_period_number = 0; pwm_period_number < number_of_periods; pwm_period_number++){ 
		for (impulse_time = 0; impulse_time <= pwm_duty_cycle; impulse_time++){
			SET_BIT(GPIOD->ODR, GPIO_ODR_ODR_12);
			SET_BIT(GPIOD->ODR, GPIO_ODR_ODR_13);
			SET_BIT(GPIOD->ODR, GPIO_ODR_ODR_14);
			SET_BIT(GPIOD->ODR, GPIO_ODR_ODR_15);
		}
		for (pause_time = 0; pause_time <= duty_cycles_array[NUMBER_OF_DUTY_CYCLES - 1] - pwm_duty_cycle; pause_time++) CLEAR_REG(GPIOD->ODR);
	}
}

int main(void){
	int duty_cycle_counter;
	int current_duty_cycle = 0;
	
	LEDs_ini();
	BUTTON_ini();
	while(1){
		if (READ_BIT(GPIOA->IDR, GPIO_IDR_IDR_0) == 1){
			duty_cycle_counter++;
			current_duty_cycle = duty_cycles_array[duty_cycle_counter % 8]; 
			Generate_pwm(4000, current_duty_cycle);
		}
		Generate_pwm(1, current_duty_cycle);
	}
}		

