#include "TeX_convertible.h"
#include <random>
#define EPS pow(10, -6)

int factorial(int n)
{
    int fact = 1;
    for (int i = 1; i <= n; i++)
        fact *= i;
    return fact;
}

string toString(const double d) {
    std::ostringstream strm;
    strm << d;
    return strm.str();
}

class SqMatrix : protected TeX_convertible {
private:
    double** sqMatrix;
    int order = 0;
    int sqMatrixMemoryAllocator(int order);
    int sqMatrixMemoryDeleter();
public:
    SqMatrix(double** SqMatrix, int order);
    SqMatrix(int order, double elemToFill, bool isDiagonal);
    SqMatrix(const SqMatrix& matrix);
    SqMatrix();
    ~SqMatrix();
    double* operator[](const int rowIndex) const;
    bool operator==(SqMatrix secondMatrix);
    bool operator!=(SqMatrix secondMatrix);
    int getOrder();
    string convert() const;

    friend istream& operator>>(istream& stream, SqMatrix& matrix)
    {
        int order = matrix.order;
        matrix.sqMatrixMemoryDeleter();
        stream >> order;
        matrix.sqMatrixMemoryAllocator(order);

        for (int i = 0; i < order; i++)
            for (int j = 0; j < order; j++)
                stream >> matrix[i][j];

        return stream;
    }

    friend ostream& operator<<(ostream& stream, const SqMatrix& matrix)
    {
        int order = matrix.order;
        for (int i = 0; i < order; i++)
            for (int j = 0; j < order; j++) {
                stream << matrix[i][j] << " ";
                if (j == order - 1)
                    stream << "\n";
            }
        return stream;
    }

    friend SqMatrix operator*(double realNumLeft, const SqMatrix& matrix)
    {
        return SqMatrix(matrix) *= realNumLeft;
    }

    friend SqMatrix& transp(SqMatrix& matrix)
    {
        SqMatrix copy(matrix);
        for (int i = 0; i < matrix.order; i++)
            for (int j = 0; j < matrix.order; j++)
                matrix[i][j] = copy[j][i];
        return matrix;
    }

    friend double det(const SqMatrix& matr)
    {
        double det = 1;
        int length = matr.order;
        SqMatrix matrix(matr);
        for (int i = 0; i < length; ++i) {
            int k = i;
            for (int j = i + 1; j < length; ++j)
                if (abs(matrix[j][i]) > abs(matrix[k][i]))
                    k = j;
            if (abs(matrix[k][i]) < EPS) {
                det = 0;
                return det;
            }
            swap(matrix.sqMatrix[i], matrix.sqMatrix[k]);
            if (i != k)
                det = -det;
            det *= matrix[i][i];
            for (int j = i + 1; j < length; ++j)
                matrix[i][j] /= matrix[i][i];
            for (int j = 0; j < length; ++j)
                if (j != i && abs(matrix[j][i]) > EPS)
                    for (int k = i + 1; k < length; ++k)
                        matrix[j][k] -= matrix[i][k] * matrix[j][i];
        }
        return det;
    }

    friend double algebraicComplement(const SqMatrix& matrix, int rowIndex, int columnIndex)
    {
        int length,
            iCopy,
            jCopy;
        length = matrix.order;
        SqMatrix copy(length - 1, 0, false);
        for (int i = 0, iCopy = 0; i < length; i++) {
            for (int j = 0, jCopy = 0; j < length; j++)
                if (i != rowIndex && j != columnIndex)
                    copy[iCopy][jCopy++] = matrix[i][j];
            jCopy = 0;
            if (i != rowIndex)
                iCopy++;
        }
        return pow(-1, rowIndex + columnIndex) * det(copy);
    }

    friend SqMatrix inverseMatrix(const SqMatrix& matrix)
    {
        int length = matrix.order;
        double determ = det(matrix);
        if (abs(determ) < EPS)
            throw "determinant = 0";
        SqMatrix copy(matrix);
        for (int i = 0; i < length; i++)
            for (int j = 0; j < length; j++)
                copy[j][i] = algebraicComplement(matrix, i, j);
        return copy * (1.0 / determ);
    }

    friend double traceMatrix(const SqMatrix& matrix)
    {
        double trace = 0;
        for (int i = 0; i < matrix.order; i++)
            trace += matrix[i][i];
        return trace;
    }

    friend SqMatrix expMatrix(const SqMatrix& matrix, int n)
    {
        SqMatrix exp(matrix.order, 1, true);
        for (int i = 1; i <= n; i++)
            exp += degreeMatrix(matrix, i) * (1.0 / factorial(i));
        return exp;
    }

    friend SqMatrix degreeMatrix(const SqMatrix& matrix, int degree)
    {
        SqMatrix copy(matrix);
        for (int i = 0; i < degree - 1; i++)
            copy *= matrix;
        return copy;
    }

    friend void swap(double* strMatr, double* strMatrSwap, int length)
    {
        double tmp;
        for (int i = 0; i < length; i++) {
            tmp = strMatr[i];
            strMatr[i] = strMatrSwap[i];
            strMatrSwap[i] = tmp;
        }
    }

    SqMatrix operator*(double realNumRight) const
    {
        return SqMatrix(*this) *= realNumRight;
    }

    SqMatrix& operator*=(double realNumRight)
    {
        for (int i = 0; i < order; i++)
            for (int j = 0; j < order; j++)
                this->sqMatrix[i][j] *= realNumRight;
        return *this;
    }

    SqMatrix operator+(const SqMatrix& secondMatrix)
    {
        return SqMatrix(*this) += secondMatrix;
    }

    SqMatrix& operator+=(const SqMatrix& secondMatrix)
    {
        int length = secondMatrix.order;
        if (this->order != secondMatrix.order)
            throw length_error("orders of matrices are different");
        for (int i = 0; i < length; i++)
            for (int j = 0; j < length; j++)
                this->sqMatrix[i][j] += secondMatrix.sqMatrix[i][j];
        return *this;
    }

    SqMatrix operator-(const SqMatrix& secondMatrix)
    {
        return SqMatrix(*this) -= secondMatrix;
    }

    SqMatrix& operator-=(const SqMatrix& secondMatrix)
    {
        int length = secondMatrix.order;
        if (this->order != secondMatrix.order)
            throw length_error("orders of matrices are different");
        for (int i = 0; i < length; i++)
            for (int j = 0; j < length; j++)
                this->sqMatrix[i][j] -= secondMatrix.sqMatrix[i][j];
        return *this;
    }


    SqMatrix& operator*=(const SqMatrix& secondMatrix)
    {
        if (secondMatrix.order != this->order)
            throw length_error("orders of matrices are different");
        int length = secondMatrix.order;
        const SqMatrix objCopy(*this);
        const SqMatrix secondCopy(secondMatrix);
        *this = SqMatrix(secondMatrix.order, 0, false);
        for (int i = 0; i < length; i++)
            for (int j = 0; j < length; j++) {
                for (int k = 0; k < length; k++)
                    (*this)[i][j] += objCopy[i][k] * secondCopy[k][j];
            }
        return *this;
    }

    SqMatrix operator*(const SqMatrix& secondMatrix) const
    {
        return SqMatrix(*this) *= secondMatrix;
    }

    SqMatrix& operator= (const SqMatrix& secondMatrix)
    {
        int length = secondMatrix.order;
        this->sqMatrixMemoryDeleter();
        this->sqMatrixMemoryAllocator(length);
        this->order = length;
        for (int i = 0; i < length; i++)
            for (int j = 0; j < length; j++)
                this->sqMatrix[i][j] = secondMatrix.sqMatrix[i][j];
        return *this;
    }
};

SqMatrix::~SqMatrix()
{
    sqMatrixMemoryDeleter();
}

SqMatrix::SqMatrix(int order, double elemToFill, bool isDiagonal)
{
    sqMatrixMemoryAllocator(order);
    for (int i = 0; i < order; i++)
        for (int j = 0; j < order; j++)
            if (isDiagonal == true)
                sqMatrix[i][j] = (i == j) ? (elemToFill) : (0);
            else
                sqMatrix[i][j] = elemToFill;
}

SqMatrix::SqMatrix(const SqMatrix& matrix)
{
    int length = matrix.order;
    sqMatrixMemoryAllocator(length);
    for (int i = 0; i < length; i++)
        for (int j = 0; j < length; j++)
            this->sqMatrix[i][j] = matrix.sqMatrix[i][j];
}

SqMatrix::SqMatrix(double** matrix, int order) :order(order)
{
    sqMatrixMemoryAllocator(order);
    for (int i = 0; i < order; i++)
        for (int j = 0; j < order; j++)
            sqMatrix[i][j] = matrix[i][j];
}

SqMatrix::SqMatrix()
{
    sqMatrixMemoryAllocator(2);
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            sqMatrix[i][j] = 0;
}

int SqMatrix::sqMatrixMemoryAllocator(int order)
{
    try {
        sqMatrix = new double* [order];
        this->order = order;
        for (int i = 0; i < order; i++)
            sqMatrix[i] = new double[order];
        return 1;
    }
    catch (bad_alloc& ex) {
        sqMatrixMemoryDeleter();
        cerr << "Can't allocate memory for SqMatrix object\n";
        return 0;
    }
}

int SqMatrix::sqMatrixMemoryDeleter()
{
    try {
        for (int i = 0; i < this->order; i++) {
            delete[] sqMatrix[i];
        }
        delete[] sqMatrix;
    }
    catch (const exception& ex) {
        cerr << "Can't free memory\n";
        return 0;
    }
    return 1;
}

double* SqMatrix::operator[](const int rowIndex) const
{
    if (rowIndex < 0 || rowIndex >= this->order)
        throw out_of_range("out of range");
    return this->sqMatrix[rowIndex];
}

bool SqMatrix::operator==(SqMatrix secondMatrix)
{
    if (order != secondMatrix.order)
        return false;
    int length = order;
    for (int i = 0; i < length; i++)
        for (int j = 0; j < length; j++)
            if (fabs(sqMatrix[i][j] - secondMatrix[i][j]) > EPS)
                return false;
    return true;
}

bool SqMatrix::operator!=(SqMatrix secondMatrix)
{
    return !((*this) == secondMatrix);
}


int SqMatrix::getOrder()
{
    return this->order;
}

string SqMatrix::convert() const
{
    string matrix("\n$\\begin{pmatrix}\n");
    if(this->order < 2)
        throw invalid_argument("Incorrect matrix order");
    for (int i = 0; i < this->order; i++)
        for (int j = 0; j < this->order; j++) {
            matrix += toString((*this)[i][j]);
            if (j != this->order - 1)
                matrix += " & ";
            else
                if (i == this->order - 1)
                    matrix += "\n\\end{pmatrix}$\n%";
                else
                    matrix += "\\\\\n";
        }
    return matrix;
}
