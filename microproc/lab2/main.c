#include "main.h"

void Delay(volatile uint32_t nCount)
{
	while(nCount--) {}
}

int main(void){
	int led_switch_count = 0; //Переменная считающая количество нажатий кнопки
	
	LEDs_ini(); //Инициализация диодов
	BUTTON_ini(); //Инициализация кнопки
	
	while(1){
		// Фиксируем смену состояния кнопки
		if (READ_BIT(GPIOA->IDR, GPIO_IDR_IDR_0) != switch_state){
// Запоминаем его
switch_state = READ_BIT(GPIOA->IDR, GPIO_IDR_IDR_0);
		if (READ_BIT(GPIOA->IDR, GPIO_IDR_IDR_0) == 1){
			led_switch_count++; 
			switch(led_switch_count){ //На каждое нажатие зажигается свой диод
				case 1:
					CLEAR_REG(GPIOD→ODR); //очищаем регистр выходных значений диодов
					SET_BIT(GPIOD->ODR, GPIO_ODR_ODR_12); //Зажигаем оранжевый диод
					break;
				case 2:
					CLEAR_REG(GPIOD->ODR);
					SET_BIT(GPIOD->ODR, GPIO_ODR_ODR_13); // Зажигаем красный диод
					break;
				case 3:
					CLEAR_REG(GPIOD->ODR);
					SET_BIT(GPIOD->ODR, GPIO_ODR_ODR_14); // Зажигаем синий диод
					break;
				case 4:
					CLEAR_REG(GPIOD->ODR);
					SET_BIT(GPIOD->ODR, GPIO_ODR_ODR_15); // Зажигаем зелёный диод
					led_switch_count = 0; 
					break;
				}
		  	}
		Delay(90000); // Задержка нужна для того, чтобы не происходило повторного срабатывания кнопки
		}
	}
}
