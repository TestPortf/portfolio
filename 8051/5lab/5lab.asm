$NOMOD51; директива подавляет предварительное определение имен 8051 SFR
$include (C8051F120.inc); подключение библиотеки микропроцессорного семейства

"C8051F120"

	LED BIT P1.6; указание на LED, подключенный к выводу P1.6
	Flag BIT 00h
	
	CSEG AT 0; определяет начало сегмента программы (в нашем случае 0)
	LJMP MAIN; передача управление по адресу MAIN
	
	ORG 000Bh; указывает метку, на которую должна прыгнуть программа при прерывании от нулевого таймера
	AJMP T0_ISR

; Периферийные специфические функции инициализации,
; Вызываются из метки Init_Device
Reset_Sources_Init:
    mov  WDTCN,     #0DEh; отключение стороживого таймера
    mov  WDTCN,     #0ADh
    ret
	
Timer_Init:
    mov  SFRPAGE,   #TIMER01_PAGE
    mov  TMOD,      #01h; Устанавливает таймер в 16-разрядный режим
    ret

Port_IO_Init:
    mov  SFRPAGE,   #CONFIG_PAGE
    mov  P1MDOUT,   #040h
    mov  XBR2,      #040h
    ret

Interrupts_Init:
mov  IE,        #082h; Разрешение прерывания от таймера Т0
    ret

; Функция инициализации для устройства,
; Вызывает Init_Device из основной программы
Init_Device:
    lcall Reset_Sources_Init
    lcall Port_IO_Init
	lcall Timer_Init
    lcall Interrupts_Init
    ret
	
T0_ISR:     ;Обработчик прерываний
	CLR TF0;     Снимаем запрос на прерывание от таймера Т0
	DJNZ R0, OUT; декрементируем переменную задержки и если она равна 0, то
	MOV R0, #10h; загружаем число в счётчик прерываний
	SETB Flag; Устанавливаем флаг
OUT: RETI;	Возвращаемся обратно в MAIN, там где мы были перед прерыванием
		
MAIN:
	CALL Init_Device; инициализация устройства
	CLR LED; выключение светодиода
	MOV R0, #10h; R0 - счётчик числа прерываний
	CLR FLAG
	SETB  TR0; Включаем таймер 0
LOOP:
    JBC Flag, SKIP; Если флаг установлен, то мы его сбрасываем и переходим к метке SKIP 
    SJMP LOOP

SKIP:
    CPL LED;  Инвертируем состояние диода
    SJMP LOOP; Начинаем сначала
	
END
