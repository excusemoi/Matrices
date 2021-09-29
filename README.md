Демонстрация работы - 1.txt_tex.pdf

Формулировка задания: 
Разработайте класс квадратной матрицы, наследующий класс
TeX_convertible. Класс должен содержать указатель на место в
памяти, где располагаются элементы матрицы (вещественного типа
данных) и её размерность (целое неотрицательное число). Для
класса реализуйте необходимое число конструкторов (при этом хотя
бы один из них должен принимать аргументы по умолчанию);
перегрузите конструктор копирования; деструктор; оператор
присваивания; арифметические операторы для сложения матриц,
вычитания матриц, умножения матриц, умножения матрицы на
число и числа на матрицу; операторы сравнения матриц на предмет
полного равенства (epsilon принять равным 1e-6); индексатор для
взятия значения из матрицы по индексам строки/столбца (отсчёт с
0); операторы вставки в поток и выгрузки из потока. Реализуйте
дружественные методы вычисления определителя, нахождения
обратной матрицы, нахождения транспонированной матрицы,
вычисление следа матрицы и матричной экспоненты. При невозможности
выполнения операции, должна быть сгенерирована исключительная
ситуация (для каждого типа ошибки - свой тип исключительной
ситуации), которая должна быть перехвачена и обработана в
вызывающем коде.
Продемонстрируйте работу с вашим классом: на вход программе
подаётся файл, содержащий выражения с матрицами (каждое из
выражений содержит одну из операций: сложение / вычитание /
умножение матриц, умножение матрицы на число, умножение числа
на матрицу, сравнение матриц (==, !=); нахождение определителя
матрицы; нахождение обратной матрицы; нахождение
транспонированной матрицы; нахождение следа матрицы;
нахождение матричной экспоненты; формат представления данных в
файле определите самостоятельно). Необходимо вычислить
значения всех выражений и сгенерировать TeX-файл (выражения,
при вычислении которых была сгенерирована исключительная
ситуация, в выходной файл выписываться не должны), где каждое
выражение будет иметь вид
<исходное выражение> = <результат вычисления выражения>.
(optional) После генерации TeX-файла необходимо запустить его
компиляцию и получить на выходе pdf-файл.
Замечания. Арифметические операции необходимо
реализовать с помощью соответствующих операции присваивания:
например, операция + должна быть реализована с помощью
операции +=. Необходимо продемонстрировать передачу аргументов
в функции по значению и по ссылке; возврат объекта из функции.

