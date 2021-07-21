$NOMOD51; директива подавляет предварительное определение имен 8051 SFR
$include (C8051F120.inc); подключение библиотеки микропроцессорного семейства

"C8051F120"

	LED BIT P1.6; указание на LED, подключенный к выводу P1.6
	
	CSEG AT 0; определяет начало сегмента программы (в нашем случае 0)
	LJMP MAIN; передача управления по адресу MAIN
	
	ORG 000Bh; указывает метку, на которую должна прыгнуть программа при прерывании от нулевого таймера
	AJMP T0_ISR

Reset_Sources_Init:
    mov  WDTCN,     #0DEh; 
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
mov  IE,        #082h; Разрешение прерывания от таймера T0
    ret

Init_Device:
    lcall Reset_Sources_Init
    lcall Port_IO_Init
	lcall Timer_Init
    lcall Interrupts_Init
    ret
	
T0_ISR:     ; Обработчик прерываний на языке Assembler
	CLR TF0;     Снимаем запрос на прерывания от таймера T0
	RETI;	Возвращаемся обратно туда, где мы были перед прерыванием
		
MAIN:
	CALL Init_Device; инициализация устройства
	MOV R0, #10h; записываем в счётчик числа прерываний число 16
	setb TR0; запуск нулевого таймера
	
IDLE:
	ORL PCON, #01h; Включаем режим ожидания
	DJNZ R0, IDLE; после каждого прерывания из r0 вычитается 1
	MOV R0, #10h
	CPL LED; Меняем состояние диода
	JMP IDLE; 
	
END
