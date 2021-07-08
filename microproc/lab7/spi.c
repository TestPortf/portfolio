void SPI_Tx(uint8_t address, uint8_t data){
	// Сигнализируем о начале передачи
	SET_BIT(GPIOE->BSRRH, GPIO_ODR_ODR_3);
	// Ждём готовности передачи
	while (READ_BIT(SPI1->SR, SPI_SR_TXE) == 0);
	// Передаём адрес регистра в который будем записывать данные
	WRITE_REG(SPI1->DR, address);
	// Ждём ответа устройства
	while (READ_BIT(SPI1->SR, SPI_SR_RXNE) == 0);
	// Считываем мусор, чтобы была возможность послать следующее сообщение
	READ_REG(SPI1->DR);
	while (READ_BIT(SPI1->SR, SPI_SR_TXE) == 0);
	// Передаём данные, которые нужно записать
	WRITE_REG(SPI1->DR, data);
	while (READ_BIT(SPI1->SR, SPI_SR_RXNE) == 0);
	// Считываем мусор
	READ_REG(SPI1->DR);
	// Ждём окончания работы SPI
	while (READ_BIT(SPI1->SR, SPI_SR_BSY) == SPI_SR_BSY);
	// Сигнализируем об окончании передачи
	SET_BIT(GPIOE->BSRRL, GPIO_ODR_ODR_3);
}

uint8_t SPI_Rx(uint8_t address){
	uint8_t recived_byte;
	
	// Сигнализируем о начале передачи
	SET_BIT(GPIOE->BSRRH, GPIO_ODR_ODR_3);
	// Ждём готовности передачи
	while (READ_BIT(SPI1->SR, SPI_SR_TXE) == 0);
	// Передаём бит чтения данных и адрес регистра, который мы собирамся прочесть
	WRITE_REG(SPI1->DR, address | TRANSMIT_RW);
	// Ждём ответа устройства
	while (READ_BIT(SPI1->SR, SPI_SR_RXNE) == 0);
	// Считываем мусор, чтобы была возможность принять новые данные
	recived_byte = READ_REG(SPI1->DR);
	while (READ_BIT(SPI1->SR, SPI_SR_TXE) == 0);
	// Посылаем мусор, чтобы устройство послало нам запрашиваемые данные
	WRITE_REG(SPI1->DR, 0x00);
	while (READ_BIT(SPI1->SR, SPI_SR_RXNE) == 0);
	// Наконец считаываем нужное нам
	recived_byte = READ_REG(SPI1->DR);
	// Ждём окончания работы SPI
	while (READ_BIT(SPI1->SR, SPI_SR_BSY) == SPI_SR_BSY);
	// Сигнализируем об окончании передачи
	SET_BIT(GPIOE->BSRRL, GPIO_ODR_ODR_3);
	
	return recived_byte;
}
