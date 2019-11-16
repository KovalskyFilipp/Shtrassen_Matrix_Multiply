#define _CRTDBG_MAP_ALLOC 
#include <stdlib.h>
#include <crtdbg.h> 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFF_SIZE 1024

/* Вспомогательная функция для определения необходимого размера массивов */
int log2(int x)
{
	int result = 0;
	while (x != 0) {
		x = x >> 1;
		result++;
	}
	return result;
}

/* Вспомогательная функция для проверки размеров массивов */
int is_power_of_2(size_t n)
{
	if (n != 0 && n != 1 && (n & (n - 1)) == 0)
		return 1;
	return 0;
}

/* Функция динамического выделения памяти под  массивы*/
float** dynamic_array_alloc(size_t N, size_t M)
{
	float** A = (float**)calloc(N, sizeof(float*));
	for (int i = 0; i < N; i++) {
		A[i] = (float*)calloc(M, sizeof(float));
	}
	return A;
}

/* Функция для освобождения памяти */
void dynamic_array_free(float** A, size_t N)
{
	for (int i = 0; i < N; i++) {
		free(A[i]);
	}
	free(A);
}

/* Функция для вывода двумерных массивов */
void dynamic_array_print(float** A, size_t N, size_t M)
{
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			printf("%.3f ", A[i][j]);
		}
		printf("\n");
	}
}

/* Вспомогательная функция умножения матриц */
float** Multiply(float** res, float** a, float** b)
{
	float m1 = (a[0][0] + a[1][1]) * (b[0][0] + b[1][1]);
	float m2 = (a[1][0] + a[1][1]) * b[0][0];
	float m3 = (b[0][1] - b[1][1]) * a[0][0];
	float m4 = (b[1][0] - b[0][0]) * a[1][1];
	float m5 = (a[0][0] + a[0][1]) * b[1][1];
	float m6 = (a[1][0] - a[0][0]) * (b[0][0] + b[0][1]);
	float m7 = (a[0][1] - a[1][1]) * (b[1][0] + b[1][1]);

	res[0][0] = m1 + m4 - m5 + m7;
	res[0][1] = m3 + m5;
	res[1][0] = m2 + m4;
	res[1][1] = m1 - m2 + m3 + m6;

	return res;
}

/* Вспомогательная функция сложения матриц */
float** Addition(float** res, float** a, float** b, int n, int m)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			res[i][j] = a[i][j] + b[i][j];

	return res;
}

/* Вспомогательная функция вычитания матриц */
float** Subtract(float** res, float** a, float** b, int n, int m)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			res[i][j] = a[i][j] - b[i][j];

	return res;
}

/* Функция умножения матриц по алгоритму Штрассена */
float** MultiShtrassen(float** dest, float** srcA, float** srcB, int length)
{
	if (length == 2)
		return Multiply(dest, srcA, srcB);

	int len = length / 2;

	float** a11 = dynamic_array_alloc(len, len);
	float** a12 = dynamic_array_alloc(len, len);
	float** a21 = dynamic_array_alloc(len, len);
	float** a22 = dynamic_array_alloc(len, len);

	float** b11 = dynamic_array_alloc(len, len);
	float** b12 = dynamic_array_alloc(len, len);
	float** b21 = dynamic_array_alloc(len, len);
	float** b22 = dynamic_array_alloc(len, len);

	float** c11 = dynamic_array_alloc(len, len);
	float** c12 = dynamic_array_alloc(len, len);
	float** c21 = dynamic_array_alloc(len, len);
	float** c22 = dynamic_array_alloc(len, len);

	float** m1 = dynamic_array_alloc(len, len);
	float** m2 = dynamic_array_alloc(len, len);
	float** m3 = dynamic_array_alloc(len, len);
	float** m4 = dynamic_array_alloc(len, len);
	float** m5 = dynamic_array_alloc(len, len);
	float** m6 = dynamic_array_alloc(len, len);
	float** m7 = dynamic_array_alloc(len, len);

	float** temp1 = dynamic_array_alloc(len, len);
	float** temp2 = dynamic_array_alloc(len, len);


	for (int i = 0; i < len; i++) {
		for (int j = 0; j < len; j++) {
			a11[i][j] = srcA[i][j];
			a12[i][j] = srcA[i][j + len];
			a21[i][j] = srcA[i + len][j];
			a22[i][j] = srcA[i + len][j + len];

			b11[i][j] = srcB[i][j];
			b12[i][j] = srcB[i][j + len];
			b21[i][j] = srcB[i + len][j];
			b22[i][j] = srcB[i + len][j + len];
		}
	}


	MultiShtrassen(m1, Addition(temp1, a11, a22, len, len), Addition(temp2, b11, b22, len, len), len);
	MultiShtrassen(m2, Addition(temp1, a21, a22, len, len), b11, len);
	MultiShtrassen(m3, a11, Subtract(temp1, b12, b22, len, len), len);
	MultiShtrassen(m4, a22, Subtract(temp1, b21, b11, len, len), len);
	MultiShtrassen(m5, Addition(temp1, a11, a12, len, len), b22, len);
	MultiShtrassen(m6, Subtract(temp1, a21, a11, len, len), Addition(temp2, b11, b12, len, len), len);
	MultiShtrassen(m7, Subtract(temp1, a12, a22, len, len), Addition(temp2, b21, b22, len, len), len);


	Subtract(c11, Addition(temp1, m1, m4, len, len), Subtract(temp2, m5, m7, len, len), len, len);
	Addition(c12, m3, m5, len, len);
	Addition(c21, m2, m4, len, len);
	Addition(c22, Subtract(temp1, m1, m2, len, len), Addition(temp2, m3, m6, len, len), len, len);


	for (int i = 0; i < len; i++) {
		for (int j = 0; j < len; j++) {
			dest[i][j] = c11[i][j];
			dest[i][j + len] = c12[i][j];
			dest[i + len][j] = c21[i][j];
			dest[i + len][j + len] = c22[i][j];
		}
	}

	dynamic_array_free(a11, len);
	dynamic_array_free(a12, len);
	dynamic_array_free(a21, len);
	dynamic_array_free(a22, len);

	dynamic_array_free(b11, len);
	dynamic_array_free(b12, len);
	dynamic_array_free(b21, len);
	dynamic_array_free(b22, len);

	dynamic_array_free(c11, len);
	dynamic_array_free(c12, len);
	dynamic_array_free(c21, len);
	dynamic_array_free(c22, len);

	dynamic_array_free(m1, len);
	dynamic_array_free(m2, len);
	dynamic_array_free(m3, len);
	dynamic_array_free(m4, len);
	dynamic_array_free(m5, len);
	dynamic_array_free(m6, len);
	dynamic_array_free(m7, len);

	dynamic_array_free(temp1, len);
	dynamic_array_free(temp2, len);

	return dest;
}

/* Функция определения необходимых размеров матриц, для выполнения умножения по алгоритму Штрассена  */
int getNewSize(int rows, int columns)
{
	int m[2] = { rows, columns };
	int max = 0;
	for (int i = 0; i < 2; ++i)
	{
		if (m[i] > max)
		{
			max = m[i];
		}
	}
	if (is_power_of_2(max) == 1)
		return max;
	int result = 1 << log2(max);
	return result;
}

/* Функция считывания данных из файла */
void getData(char buff[], int* m_columns, int* m_rows, float** temp1, float** temp2, int* columns_1, int* rows_1, int* columns_2, int* rows_2, int* counter)
{

	char* end_buff;
	char* token = strtok_s(buff, ";\n", &end_buff);
	int columns = 0;
	int rows = 0;
	while (token != NULL)
	{
		char* end_token;
		char* row = strtok_s(token, "{}", &end_token);
		rows = 0;
		*counter = *counter +1;
		while (row != NULL)
		{
			char* end_row;
			char* column = strtok_s(row, ",", &end_row);
			rows++;
			columns = 0;
			while (column != NULL)
			{
				columns++;
				float buf = atof(column, &end_row);
				if(*counter == 1)
					temp1[rows-1][columns-1] = buf;
				else
					temp2[rows - 1][columns - 1] = buf;
				column = strtok_s(NULL, ",", &end_row);
			}
			row = strtok_s(NULL, "{}", &end_token);
		}
		token = strtok_s(NULL, ";\n", &end_buff);
		if (*counter == 1) {
			*columns_1 = columns;
			*rows_1 = rows;
		}
		else
			*columns_2 = columns;
		*rows_2 = rows;
	}

	if (*m_rows <= rows)
		*m_rows = rows;
	if (*m_columns <= columns)
		*m_columns = columns;
}

int main(void)
{
	float** temp1 = dynamic_array_alloc(BUFF_SIZE, BUFF_SIZE);
	float** temp2 = dynamic_array_alloc(BUFF_SIZE, BUFF_SIZE);
	int size = 0;
	int counter_val = 0;
	int m_columns_val = 0;
	int m_rows_val = 0;
	int columns_1_val = 0;
	int rows_1_val = 0;
	int columns_2_val = 0;
	int rows_2_val = 0;
	int* counter = &counter_val;
	int* m_columns = &m_columns_val;
	int* m_rows = &m_rows_val;
	int* columns_1 = &columns_1_val;
	int* rows_1 = &rows_1_val;
	int* columns_2 = &columns_2_val;
	int* rows_2 = &rows_2_val;
	int c = 0;
	FILE* stream;
	fopen_s(&stream,"Data.csv", "r");
	int count = 0;
	do
	{
		char buff[BUFF_SIZE];
		fgets(buff, BUFF_SIZE, stream);
		getData(buff, m_columns, m_rows, temp1, temp2, columns_1, rows_1, columns_2, rows_2, counter);
		count++;
	} while ((getc(stream)) != EOF);
	int fclose(FILE * stream);
	size = getNewSize(*m_rows, *m_columns);
	float** A = dynamic_array_alloc(size+1,size+1);
	float** B = dynamic_array_alloc(size+1, size+1);
	float** result = dynamic_array_alloc(size+1, size+1);
	for (int i = 0; i <= *columns_1; i++) {
		for (int j = 0; j <= *rows_1; j++) {
			A[i][j] = temp1[i][j];
		}
	}
	for (int i = 0; i <= *columns_2; i++) {
		for (int j = 0; j <= *rows_2; j++) {
			B[i][j] = temp2[i][j];
		}
	}
	dynamic_array_free(temp1, BUFF_SIZE);
	dynamic_array_free(temp2, BUFF_SIZE);
	result = MultiShtrassen(result, A,B, size);
	dynamic_array_print(result, size, size);
	dynamic_array_free(A, size+1);
	dynamic_array_free(B, size+1);
	dynamic_array_free(result, size+1);
	printf("Press any key to close this window");
	_CrtDumpMemoryLeaks();
	char ex = getch();
	return 0;

}