#include "main.h"

int led_switch_count = 0;

void TIM7_IRQHandler(void){ // Прерывание от 7-ого таймера
	CLEAR_BIT(TIM7->SR, TIM_SR_UIF); // Очищаем флаг прерывания 7-ого таймера
	led_switch_count++; 
	switch(led_switch_count){ //Каждое последующее срабатывание прерывания зажигает другой диод
		case 1:
			CLEAR_REG(GPIOD->ODR); //очищаем регистр выходных значений диодов
			SET_BIT(GPIOD->ODR, GPIO_ODR_ODR_12);  //Зажигаем оранжевый диод
			break;
		case 2:
			CLEAR_REG(GPIOD->ODR);
			SET_BIT(GPIOD->ODR, GPIO_ODR_ODR_13); // Зажигаем красный диод
			break;
		case 3:
			CLEAR_REG(GPIOD->ODR);
			SET_BIT(GPIOD->ODR, GPIO_ODR_ODR_14);  // Зажигаем синий диод
			break;
		case 4:
			CLEAR_REG(GPIOD->ODR);
			SET_BIT(GPIOD->ODR, GPIO_ODR_ODR_15); // Зажигаем зелёный диод
			led_switch_count = 0;
			break;
	}
}

int main(void){
	LEDs_ini(); // Инициализация диодов
	TIM7_ini(); // Инициализация 7-ого таймера
	while(1);
}
