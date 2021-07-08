#include "init.h"

void LEDs_ini(void){
	SET_BIT(RCC->AHB1ENR, RCC_AHB1ENR_GPIODEN); 
	SET_BIT(GPIOD->MODER, GPIO_MODER_MODER12_1 | GPIO_MODER_MODER13_1 | GPIO_MODER_MODER14_1 | GPIO_MODER_MODER15_1); // Установка портов в альтернативный режим работы
	WRITE_REG(GPIOD->AFR[1], 0x22220000); // Подключение вывода 4-х каналов таймера 4 к 4-ом диодам
}

void TIM4_ini(void){
	SET_BIT(RCC->APB1ENR, RCC_APB1ENR_TIM4EN);
	SET_BIT(TIM4->CCER, TIM_CCER_CC1E | TIM_CCER_CC2E | TIM_CCER_CC3E | TIM_CCER_CC4E); // Разрешение использовать порты ввода-вывода к которым подключены каналы таймера для ШИМа
	WRITE_REG(TIM4->ARR, 0xFF); // Установка предела счёта таймера
	// Выставляем на всех каналах таймера 4 режим  работы – ШИМ
	SET_BIT(TIM4->CCMR1, TIM_CCMR1_OC1M_1 | TIM_CCMR1_OC1M_2);
	SET_BIT(TIM4->CCMR1, TIM_CCMR1_OC2M_1 | TIM_CCMR1_OC2M_2);
	SET_BIT(TIM4->CCMR2, TIM_CCMR2_OC3M_1 | TIM_CCMR2_OC3M_2);
	SET_BIT(TIM4->CCMR2, TIM_CCMR2_OC4M_1 | TIM_CCMR2_OC4M_2); 
	SET_BIT(TIM4->CR1, TIM_CR1_CEN | TIM_CR1_ARPE ); // Включаем таймер и заставляем его запоминать значение предела счёта таймера
}

void SPI1_ini(void){
	SET_BIT(RCC->AHB1ENR, RCC_AHB1ENR_GPIOAEN);
	SET_BIT(RCC->AHB1ENR, RCC_AHB1ENR_GPIOEEN);

	SET_BIT(GPIOA->MODER, GPIO_MODER_MODER5_1 | GPIO_MODER_MODER6_1 | GPIO_MODER_MODER7_1); // Установка портов в альтернативный режим работы
	SET_BIT(GPIOA->AFR[0], 0x55500000); // Подключение выводов SPI модуля к 3-м портам
	
	// Настройка порта отвечающего за сигнал о начале и окончании передачи данных
	SET_BIT(GPIOE->MODER, GPIO_MODER_MODER3_0);
	SET_BIT(GPIOE->PUPDR, GPIO_PUPDR_PUPDR3_0);
	// Сигнализируем, что передача данных не ведётся
	SET_BIT(GPIOE->BSRRL, GPIO_ODR_ODR_3);

	SET_BIT(RCC->APB2ENR, RCC_APB2ENR_SPI1EN);

	// Настройка SPI для работы с акселерометром
	WRITE_REG(SPI1->CR1, SPI_CR1_SSM | SPI_CR1_SSI | SPI_CR1_MSTR | SPI_CR1_CPHA | SPI_CR1_CPOL | SPI_CR1_BR_1 | SPI_CR1_BR_2);
	// Включение SPI
	SET_BIT(SPI1->CR1, SPI_CR1_SPE);
}
