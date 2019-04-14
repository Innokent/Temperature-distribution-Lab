/*
8. Разработать, используя средства многопотокового программирования, параллельную программу решения двумерной нестационарной краевой задачи методом конечных разностей с использованием явной вычислительной схемы. Объект моделирования - прямоугольная пластина постоянной толщины. Подробности постановки подобной задачи даны выше. Возможны граничные условия первого и второго рода в различных узлах расчетной сетки. Количество потоков, временной интервал моделирования и количество (кратное 8) узлов расчетной сетки - параметры программы. См. также замечание. Программа должна демонстрировать ускорение по сравнению с последовательным вариантом. Предусмотреть визуализацию результатов посредством утилиты gnuplot. 
*/
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>

#define k 1 // Коэффицент теплопроводности
#define Gt 0 // Приведенная скорость взаимного превращения тепловой энергии в другие виды энергии
#define dx 2 // Шаг по Х
#define dy 2 // Шаг по Y
#define L 100 // Длинна стержня
#define S 10 // Площадь поперечного сечения стержня
#define boundaryConditionFirstTop 10
#define boundaryConditionSecondTop 0.1
#define boundaryConditionFirstBot 20
#define boundaryConditionSecondBot 0.2
#define boundaryConditionFirstLeft 30
#define boundaryConditionSecondLeft 0.3
#define boundaryConditionFirstRight 40
#define boundaryConditionSecondRight 0.4
#define maxTime 100
#define currMode 1

typedef struct 
{
	pthread_t tid;
	double *mas; // Ссылка на массив начальных значений
	int firstIndexStart;
	int firstIndexEnd;
	int secondIndexStart;
	int secondIndexEnd;
	double dt;
	int n;
	int m;
} Thread_param;

pthread_mutex_t mutx;
Thread_param *threads;
double *array;
int control[1] = {0};
int block[2] = {0, 0};

double calcNext(double T01, double T11, double T21, double T10, double T12, double dt)
{
	double numerator1 = T01 - 2 * T11 + T21;
	double numerator2 = T10 - 2 * T11 + T12;
	double fraction1 = numerator1 / (dx * dx);
	double fraction2 = numerator2 / (dy * dy);
	return (dt * (k * (fraction1 + fraction2) + Gt) + T11);
}

double calcBoundary(int i, int j, int n, int m, double dt, double T, int mode)
{
	if(mode == 1) // Граничные условия первого рода
	{
		
		if(i == 0)
		{
			return boundaryConditionFirstTop;
		}
		if(i == n - 1)
		{
			return boundaryConditionFirstBot;
		}
		if(j == 0)
		{
			return boundaryConditionFirstLeft;
		}
		if(j == m - 1)
		{
			return boundaryConditionFirstRight;
		}
	}
	else // Граничные условия второго рода
	{
		if(i == 0)
		{
			return (boundaryConditionSecondTop * dt / dx + 1) * T;
		}
		if(i == n - 1)
		{
			return (boundaryConditionSecondBot * dt / dx + 1) * T;
		}
		if(j == 0)
		{
			return (boundaryConditionSecondLeft * dt / dy + 1) * T;
		}
		if(j == m - 1)
		{
			return (boundaryConditionSecondRight * dt / dy + 1) * T;
		}
	}
}

void *solver(void *arg_p)
{
	Thread_param *params = (Thread_param *) arg_p;
	double *prevLayer = (double *) calloc (params->n * params->m, sizeof(double));
	double *currLayer = (double *) calloc (params->n * params->m, sizeof(double));
	double T01, T11, T21, T10, T12;
	int K = 0;
	
	for(double t = 0.0 + params->dt; t <= maxTime; t += params->dt)
	{
		for(int i = params->firstIndexStart; i <= params->firstIndexEnd; i++ )
		{
			for(int j = params->secondIndexStart; j <= params->secondIndexEnd; j++)
			{
				if((i != 0) && (i != (params->n - 1)) && (j != 0) && (j != (params->m - 1)))
				{
					T01 = array[params->n * params->m * K + params->n * j + i - 1];
					T11 = array[params->n * params->m * K + params->n * j + i];
					T21 = array[params->n * params->m * K + params->n * j + i + 1];
					T10 = array[params->n * params->m * K + params->n * (j - 1) + i];
					T12 = array[params->n * params->m * K + params->n * (j + 1) + i];
					array[params->n * params->m * (K + 1) + params->n * j + i] = calcNext(T01, T11, T21, T10, T12, params->dt);
				}
				else
				{
					array[params->n * params->m * (K + 1) + params->n * j + i] = calcBoundary(i, j, params->n, params->m, params->dt, array[params->n * params->m * K + params->n * j + i], currMode);
				}
			}
		}
		pthread_mutex_lock(&mutx);
		block[0]++;
		pthread_mutex_unlock(&mutx);
		K++;
		while(block[0] < K * block[1])
		{
		}
		fflush(stdin);
		fflush(stdout);
	}
	control[0]++;	
}

int main(int argc, char* argv[])
{
	pthread_attr_t pattr;
	
	if(argc != 5) // Проверяем количество параметров
	{
		puts("Ошибка: Неверное количесвто параметров");
		exit(1);
	}
	
	int c_threads = atoi(argv[1]);
	double dt = atof(argv[2]);
	int N = atoi(argv[3]);
	int M = atoi(argv[4]);
	int K = (int)(maxTime / dt) + 1;
	block[1] = c_threads;
	
	if((N * M) % 8) // Проверям кратность числа узлов восьми
	{
		puts("Ошибка: Количество узлов сетки должно быть кратно 8");
		exit(3);
	}
	
	if(N % c_threads) // Проверяем условие делимости нацело числа узлов и числа потоков
	{
		puts("Ошибка: Количество узов сетки по одной координте не делится нацело на количество потоков");
		exit(4);
	}

	N += 2; // Учитываем граничные условия для верха и низа пластины
	M += 2; // Учитываем граничные условия для левого и правого края пластины

	array = (double *) calloc (K * N * M, sizeof(double)); // Выделяем память под массив

	if(!array) // Проверяем удалось ли выделить память под массив
	{
		puts("Ошибка: Не удалось выделить память");
		exit(2);
	}
	
	for(int i = 0; i < N; i++) // Заполняем узлы балки значениями в начальный момент времни (t = 0)
	{
		for(int j = 0; j < M; j++)
		{
			if((i != 0) && (i != N - 1) && (j != 0) && (j != M - 1))
			{
				array[N * j + i] = 0.0;
			}
			else
			{
				array[N * j + i] = calcBoundary(i, j, N, M, dt, 0, 1);
			}
		}
	}

	pthread_attr_init(&pattr);
	pthread_attr_setscope(&pattr, PTHREAD_SCOPE_SYSTEM);
	pthread_attr_setdetachstate(&pattr, PTHREAD_CREATE_JOINABLE);
	pthread_mutex_init(&mutx, NULL);

	threads = (Thread_param *) calloc(c_threads, sizeof(Thread_param)); // Выделяем память под потоки

	for(int i = 0; i < c_threads; i++) // Инициализируем атрибуты потоков
	{
		threads[i].mas = array;
		threads[i].n = N;
		threads[i].m = M;
		threads[i].dt = dt;
		threads[i].firstIndexStart = 1 + i * ((N - 2) / c_threads);
		if(i == 0)
		{
			threads[i].firstIndexStart--;
		}
		threads[i].firstIndexEnd = (i + 1) * ((N - 2) / c_threads);
		if(i == c_threads - 1)
		{
			threads[i].firstIndexEnd++;
		}
		threads[i].secondIndexStart = 0;
		threads[i].secondIndexEnd = M - 1;
		
	}
	
	for(int i = 0; i < c_threads; i++)
	{
		pthread_create(&threads[i].tid, &pattr, solver, (void *) &(threads[i])); //Создаем потоки
	}

	while(control[0] != c_threads) // Ждем завершения работы всех потоков
	{
	}

	
	FILE *output;
	output = fopen("output.txt", "w");
	
	for(int i = 0; i < M; i++)
	{
		for(int j = 0; j < N; j++)
		{
			fprintf(output, "%d %d %lf\n", i, j, array[N * M * (K - 2) + N * j + i]);
		}
	}
	
	fclose(output);

	pthread_mutex_destroy(&mutx); // Освобождаем мьютекс

	free(array); // Освобождаем массив данных
	free(threads); // Освобождаем массив потоков

	return 1;
}
