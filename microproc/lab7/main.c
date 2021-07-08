#include "main.h"

void Delay(volatile uint32_t nCount)
{
	while(nCount--) {}
}

int main(void){
	int8_t Accel_x, Accel_y;
	
	// Инициализируемуем требуемую перефирию
	LEDs_ini();
	TIM4_ini();
	SPI1_ini();
	
	// Включаем отслеживание по осям у и х, а так же включаем сам акселерометр 
	SPI_Tx(CTRL_REG1, CTRL_REG1_PD | CTRL_REG1_XEN | CTRL_REG1_YEN);
	
	while(1){
		// Считываем измеренные акселерометром параметры
		Accel_x = SPI_Rx(OUTX);
		Accel_y = SPI_Rx(OUTY);
		
		// В зависимости от полярности параметров меняем скважность ШИМ
		if (Accel_x < -NOISE_LIMIT || Accel_x > NOISE_LIMIT){
			if (Accel_x < 0){ 
				WRITE_REG(TIM4->CCR1, -Accel_x);
				WRITE_REG(TIM4->CCR3, 0);
			}
			else{
				WRITE_REG(TIM4->CCR3, Accel_x);
				WRITE_REG(TIM4->CCR1, 0);
			}
		}
		if (Accel_y < -NOISE_LIMIT || Accel_y > NOISE_LIMIT){
			if (Accel_y < 0){
				WRITE_REG(TIM4->CCR4, -Accel_y);
				WRITE_REG(TIM4->CCR2, 0);
			}
			else{
				WRITE_REG(TIM4->CCR2, Accel_y);
				WRITE_REG(TIM4->CCR4, 0);
			}
		}
	}
}
