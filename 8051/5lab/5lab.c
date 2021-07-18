#include <stdio.h>
#include "C8051F120.h"

sbit LED = P1^6;
cnt = 0x10;
int FLAG = 0x01;

// Периферийные специфические функции инициализации,
// Вызывается из функции Init_Device()
void Reset_Sources_Init(){
    WDTCN     = 0xDE;// Отключаем стороживой таймер
    WDTCN     = 0xAD;
}

void Timer_Init(){
    SFRPAGE   = TIMER01_PAGE;
    TMOD      = 0x01; // Устанавливает таймер в 16-разрядный режим
}

void Port_IO_Init(){
    SFRPAGE   = CONFIG_PAGE;
    P1MDOUT   = 0x40;
    XBR2      = 0x40;
}

void Interrupts_Init(){
    IE        = 0x82; // Разрешение прерывания от таймера Т0
}

// Функция инициализации для устройства,
// Вызов Init_Device () из вашей основной программы
void Init_Device(void){
    Reset_Sources_Init();
    Port_IO_Init();
    Interrupts_Init();
	Timer_Init();
}

Timer_ISR(void) interrupt 1 {
    TF0 = 0;    //чистка флага прерывания
	cnt --; // декрементируем переменную задержки
	if (cnt == 0x00){ // если переменная задержки равна 0, то
	    cnt = 0x10; // Перезаписываем в неё число 16
		FLAG=0x01; // Устанавливаем флаг
	}
}

void main(void){
	Init_Device();
	TR0 = 1; // Включаем таймер 0
	LED = 0; // Выключаем диод
	while (1){
	    if(FLAG==0x01){ // Проверяем не установлен ли флаг, если установлен, то
	        FLAG=0x00; // Обнуляем флагг
	        LED=~LED; // Инвертируем состояние диода
	    }
	}
}
